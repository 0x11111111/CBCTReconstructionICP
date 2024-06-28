/**
 * @file CBCTJawRegistration.cpp
 * @brief Implementation of the CBCTJawRegistration class for registering CBCT images of jaws.
 *
 * This file contains the implementation of the CBCTJawRegistration class, which provides methods to register
 * cone-beam computed tomography (CBCT) images of upper and lower jaws with corresponding 3D models of teeth.
 *
 * This project is derived from https://github.com/0x11111111/VisualLab/tree/base_meshlib_lab and is dedicated to generating
 * the `cbct_jaw_registration.pyd` module, which helps in the registration of oral scan data to CBCT images.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

#include <filesystem>
#include <ranges>

#include <MRMesh/MRPython.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRAffineXf.h>
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRICP.h>
#include <MRMesh/MRMeshOrPoints.h>
#include <MRMesh/MRBox.h>
#include <MRMesh/MRObjectVoxels.h>
#include <MRMesh/MRVoxelsLoad.h>
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRMeshSave.h>
#include <MRMesh/MRMeshIntersect.h>
#include <MRMesh/MRLine3.h>
#include <MRMesh/MRMeshComponents.h>
#include <MRMesh/MRMeshPart.h>
#include <MRMesh/MRRingIterator.h>
#include <MRMesh/MRRegionBoundary.h>
#include <MRMesh/MRMeshFillHole.h>
#include <MRMesh/MRDenseBox.h>

#include "CBCTJawRegistration.h"

namespace fs = std::filesystem;

/**
 * @brief Constructor for CBCTJawRegistration class.
 *
 * Initializes the CBCTJawRegistration object by loading the upper and lower teeth meshes,
 * the upper jaw mesh, and the CBCT voxel data from the specified files and directories.
 *
 * @param upper_teeth_file Path to the upper teeth mesh file.
 * @param lower_teeth_file Path to the lower teeth mesh file.
 * @param tiff_folder Path to the directory containing the CBCT TIFF images.
 * @param upper_teeth_jaw Path to the upper jaw mesh file.
 * @param lower_x Lower bound of the x-dimension for the region of interest.
 * @param upper_x Upper bound of the x-dimension for the region of interest.
 * @param lower_y Lower bound of the y-dimension for the region of interest.
 * @param upper_y Upper bound of the y-dimension for the region of interest.
 * @param lower_z Lower bound of the z-dimension for the region of interest.
 * @param upper_z Upper bound of the z-dimension for the region of interest.
 * @param voxel_space Voxel spacing used for loading the CBCT data.
 * @param temp Path to the directory for saving temporary files.
 */
CBCTJawRegistration::CBCTJawRegistration(
    const std::wstring& upper_teeth_file,
    const std::wstring& lower_teeth_file,
    const std::wstring& tiff_folder,
    const std::wstring& upper_teeth_jaw,
    const int lower_x,
    const int upper_x,
    const int lower_y,
    const int upper_y,
    const int lower_z,
    const int upper_z,
    const float voxel_space,
    const std::wstring& temp)
{
    m_upper_teeth_mr = LoadMeshFromFile(upper_teeth_file);
    m_lower_teeth_mr = LoadMeshFromFile(lower_teeth_file);
    m_upper_jaw_mr = LoadMeshFromFile(upper_teeth_jaw);
    m_temp_folder = temp;

    if (!m_temp_folder.empty())
    {
        if (!fs::exists(m_temp_folder))
        {
            fs::create_directories(m_temp_folder);
        }
    }

    auto voxel_load = MR::VoxelsLoad::loadTiffDir(
        {
            fs::path(tiff_folder),
            { voxel_space, voxel_space, voxel_space },
            MR::VoxelsLoad::GridType::DenseGrid,
            {}
        }
    );
    auto vdb = voxel_load.value();

    auto dims = vdb.dims;

    auto dim_x = dims.x;
    auto dim_y = dims.y;
    auto dim_z = dims.z;

    auto min_x = std::max(0, lower_x - 10);
    auto max_x = std::min(dim_x, upper_x + 10);
    auto min_y = std::max(0, lower_y - 10);
    auto max_y = std::min(dim_y, upper_y + 10);
    auto min_z = lower_z;
    auto max_z = std::min(dim_z, upper_z + 10);

    auto upper_active_bounds = MR::Box3i({ min_x, min_y, static_cast<int>((max_z - min_z) * 0.4) + min_z }, { max_x, max_y, max_z });
    auto lower_active_bounds = MR::Box3i({ min_x, min_y, min_z }, { max_x, max_y, max_z });

    MR::ObjectVoxels object_voxels;
    object_voxels.construct(vdb);
    auto& bins = object_voxels.histogram().getBins();

    auto sum = std::accumulate(bins.begin(), bins.end(), 0);
    auto rit = bins.rbegin();
    auto accumulation = 0;

    auto it = bins.begin();
    for (; it != bins.end(); ++it)
    {
        accumulation += *it;
        if (accumulation > sum * 0.99)
        {
            break;
        }
    }

    auto bin_id = std::distance(bins.begin(), it);

    object_voxels.setActiveBounds(upper_active_bounds);
    object_voxels.setIsoValue(bin_id);
    m_cbct_upper_mr = MR::Mesh(*object_voxels.mesh());
    if (!m_temp_folder.empty())
    {
        MR::MeshSave::toAnySupportedFormat(m_cbct_upper_mr, m_temp_folder / "upper_recon.ply");
    }

    object_voxels.setActiveBounds(lower_active_bounds);
    object_voxels.setIsoValue(bin_id);
    m_cbct_all_mr = MR::Mesh(*object_voxels.mesh());

    if (!m_temp_folder.empty())
    {
        MR::MeshSave::toAnySupportedFormat(m_cbct_all_mr, m_temp_folder / "all_recon.ply");
    }
}

/**
 * @brief Get the number of valid vertices in the upper teeth mesh.
 *
 * @return int The number of valid vertices in the upper teeth mesh.
 */
int CBCTJawRegistration::GetUpperTeethVertNum()
{
    return m_upper_teeth_mr.topology.numValidVerts();
}

/**
 * @brief Get the transformation matrix for the upper teeth.
 *
 * @return pybind11::tuple A tuple containing the transformation matrix and translation vector.
 */
pybind11::tuple CBCTJawRegistration::GetUpperTransformation()
{
    return AffineXf3fToPyTuple(m_upper_teeth_transformation);
}

/**
 * @brief Get the transformation matrix for the lower teeth.
 *
 * @return pybind11::tuple A tuple containing the transformation matrix and translation vector.
 */
pybind11::tuple CBCTJawRegistration::GetLowerTransformation()
{
    return AffineXf3fToPyTuple(m_lower_teeth_transformation);
}

/**
 * @brief Perform the registration of the CBCT images with the teeth meshes.
 *
 * This method aligns the upper and lower teeth meshes with the corresponding CBCT images using an iterative
 * closest point (ICP) algorithm. It also saves the transformed meshes to the temporary folder if specified.
 */
void CBCTJawRegistration::Registration()
{
    using namespace MR;
    m_upper_teeth_transformation = ICPHalf(m_cbct_upper_mr, m_upper_teeth_mr, true);
    m_upper_teeth_mr.transform(m_upper_teeth_transformation);
    if (!m_temp_folder.empty())
    {
		MR::MeshSave::toAnySupportedFormat(m_upper_teeth_mr, m_temp_folder / "upper_teeth_icp_transformed.ply");
	}
    
    m_upper_jaw_mr.transform(m_upper_teeth_transformation);
    auto expansion = 1.f;
    for (auto& v : m_upper_jaw_mr.topology.getValidVerts())
    {
        Vector3f average_normal(0.f, 0.f, 0.f);
        int count = 0;
        for (auto& e : orgRing(m_upper_jaw_mr.topology, v))
        {
            auto& dst_v = m_upper_jaw_mr.topology.dest(e);
            average_normal += m_upper_jaw_mr.normal(dst_v);
            count++;
        }
        average_normal /= static_cast<float>(count);
        m_upper_jaw_mr.points[v] += expansion * average_normal;
    }

    auto right_boundary = findRightBoundary(m_upper_jaw_mr.topology);

    FaceBitSet filled_new_faces;

    for (auto& boundary : right_boundary)
    {
        fillHole(m_upper_jaw_mr, boundary[0],
            {
                {},
                &filled_new_faces,
                FillHoleParams::MultipleEdgesResolveMode::None,
                true,
                20
            }
        );
    }

    auto outside_vbs = m_cbct_all_mr.topology.getValidVerts();

    for (auto& v : m_cbct_all_mr.topology.getValidVerts())
    {
        std::vector<MeshIntersectionResult> intersection_results;
        rayMeshIntersectAll(
            m_upper_jaw_mr,
            Line3f(m_cbct_all_mr.points[v], m_cbct_all_mr.normal(v)),
            [&intersection_results](const MeshIntersectionResult& result) -> bool
            {
                intersection_results.push_back(result);
                return true;
            },
            -100.f,
            100.f
        );
        // std::cout << intersection_results.size() << std::endl;

        if (intersection_results.size() > 0)
        {
            auto first_distance = intersection_results[0].distanceAlongLine;
            for (auto& intersection_res : intersection_results | std::ranges::views::drop(1))
            {
                auto distance = intersection_res.distanceAlongLine;
                if (std::signbit(first_distance) != std::signbit(distance))
                {
                    outside_vbs.reset(v);
                    break;
                }
            }
        }
    }


    auto outside_fbs = getInnerFaces(m_cbct_all_mr.topology, outside_vbs);
    auto outside_largest_components_fbs = MeshComponents::getLargestComponent({ m_cbct_all_mr, &outside_fbs });

    m_cbct_lower_mr = ExtractSubMeshByFbs(m_cbct_all_mr, outside_largest_components_fbs);

    m_lower_teeth_transformation = ICPHalf(m_cbct_lower_mr, m_lower_teeth_mr, false);
    m_lower_teeth_mr.transform(m_lower_teeth_transformation);

    // Lower teeth requires second alignment towards the whole CBCT
    MR::MeshOrPoints lower_teeth_second_mop(m_lower_teeth_mr);

    MR::ICPProperties lower_second_icp_properties;
    lower_second_icp_properties.cosTreshold = 2.0f / 3;
    lower_second_icp_properties.method = MR::ICPMethod::PointToPlane;
    lower_second_icp_properties.icpMode = MR::ICPMode::AnyRigidXf;
    lower_second_icp_properties.badIterStopCount = 100;
    lower_second_icp_properties.iterLimit = 300;
    lower_second_icp_properties.farDistFactor = 1.5f;
    lower_second_icp_properties.distThresholdSq = 2500.f;
    lower_second_icp_properties.mutualClosest = false;

    MR::ICP lower_second_icp(lower_teeth_second_mop, m_cbct_all_mr, MR::AffineXf3f(), MR::AffineXf3f(), 0.3);
    lower_second_icp.setParams(lower_second_icp_properties);
    lower_second_icp.updatePointPairs();

    m_lower_teeth_transformation = m_lower_teeth_transformation * lower_second_icp.calculateTransformation();

    if (!m_temp_folder.empty())
    {
        MR::MeshSave::toAnySupportedFormat(m_cbct_lower_mr, m_temp_folder / "lower_recon.ply");
        MR::MeshSave::toAnySupportedFormat(m_lower_teeth_mr, m_temp_folder / "lower_teeth_icp_transformed.ply");
    }
}

/**
 * @brief Load a mesh from a file.
 *
 * @param file_path Path to the mesh file.
 * @return MR::Mesh The loaded mesh.
 * @throws std::runtime_error if the mesh cannot be loaded.
 */
MR::Mesh CBCTJawRegistration::LoadMeshFromFile(const std::wstring& file_path)
{
    auto res = MR::MeshLoad::fromAnySupportedFormat(fs::path(file_path));
    if (!res.has_value())
    {
        throw std::runtime_error(res.error());
    }
    return res.value();
}

/**
 * @brief Perform a half ICP alignment for the given CBCT and teeth meshes.
 *
 * @param cbct_mr The CBCT mesh.
 * @param teeth_mr The teeth mesh.
 * @param upper Boolean indicating if the alignment is for the upper teeth.
 * @return MR::AffineXf3f The transformation matrix resulting from the ICP alignment.
 */
MR::AffineXf3f CBCTJawRegistration::ICPHalf(MR::Mesh cbct_mr, MR::Mesh teeth_mr, bool upper)
{
    using namespace MR;
    auto cbct_bound = cbct_mr.computeBoundingBox();
    auto teeth_bound = teeth_mr.computeBoundingBox();

    auto& cbct_min = cbct_bound.min;
    auto& cbct_max = cbct_bound.max;

    auto& teeth_min = teeth_bound.min;
    auto& teeth_max = teeth_bound.max;

    auto cbct_xz_panel = MR::Box3f({ cbct_min[0], cbct_min[1], cbct_min[2] }, { cbct_max[0], cbct_min[1], cbct_max[2] });
    auto teeth_xz_panel = MR::Box3f({ teeth_min[0], teeth_min[1], teeth_min[2] }, { teeth_max[0], teeth_min[1], teeth_max[2] });

    auto cbct_xz_panel_center = cbct_xz_panel.center();
    auto teeth_xz_panel_center = teeth_xz_panel.center();

    // Align the upper of the teeth to the upper of the CBCT


    auto y_translation_dist = 6.f;
    auto z_translation_dist = 6.f;

    auto cbct_center = cbct_mr.findCenterFromPoints();
    auto teeth_center = teeth_mr.findCenterFromPoints();

    auto first_translation = MR::AffineXf3f::translation(MR::Vector3f(cbct_xz_panel_center - teeth_xz_panel_center) + MR::Vector3f{ 0, y_translation_dist, z_translation_dist });

    //DenseBox cbct_dense_box(cbct_mr);
    //DenseBox teeth_dense_box(teeth_mr);

    //auto cbct_corner_1 = cbct_dense_box.corner({ false, true, false });
    //auto cbct_corner_2 = cbct_dense_box.corner({ true, true, true });
    //auto cbct_center = (cbct_corner_1 + cbct_corner_2) / 2.f;

    //auto teeth_corner_1 = teeth_dense_box.corner({ false, true, false });
    //auto teeth_corner_2 = teeth_dense_box.corner({ true, true, true });
    //auto teeth_center = (teeth_corner_1 + teeth_corner_2) / 2.f;

    //auto first_translation = MR::AffineXf3f::translation(cbct_center - teeth_center);

    auto teeth_first_translation_mr(teeth_mr);
    teeth_first_translation_mr.transform(first_translation);

    std::string first_translated_name = upper ? "upper_teeth_initial_translated.ply" : "lower_teeth_initial_translated.ply";

    MR::MeshSave::toAnySupportedFormat(teeth_first_translation_mr, m_temp_folder / first_translated_name);

    MR::MeshOrPoints cbct_mop(cbct_mr);
    MR::MeshOrPoints teeth_first_mop(teeth_first_translation_mr);

    MR::ICPProperties first_icp_properties;
    first_icp_properties.cosTreshold = -1.f;
    first_icp_properties.method = MR::ICPMethod::PointToPlane;
    first_icp_properties.icpMode = MR::ICPMode::AnyRigidXf;
    first_icp_properties.badIterStopCount = 100;
    first_icp_properties.iterLimit = 300;
    first_icp_properties.farDistFactor = 3.5f;
    first_icp_properties.distThresholdSq = 100000.f;
    first_icp_properties.mutualClosest = false;

    MR::ICP first_icp(teeth_first_mop, cbct_mop, MR::AffineXf3f(), MR::AffineXf3f(), 0.3);
    first_icp.setParams(first_icp_properties);
    first_icp.updatePointPairs();

    auto teeth_first_transformation = first_icp.calculateTransformation();

    // Lower teeth only need first rough alignment
    if (!upper)
    {
        return teeth_first_transformation * first_translation;
    }

    teeth_first_translation_mr.transform(teeth_first_transformation);

    //std::string first_transformed_name = upper ? "upper_teeth_first_transformed.ply" : "lower_teeth_first_transformed.ply";
    //MeshSave::toAnySupportedFormat(teeth_first_translation_mr, m_temp_folder / first_transformed_name);

    auto second_translation = MR::AffineXf3f::translation({ 0.f, 0.f, -5.f});
    teeth_first_translation_mr.transform(second_translation);

    MR::MeshOrPoints teeth_second_mop(teeth_first_translation_mr);

    MR::ICPProperties second_icp_properties;
    second_icp_properties.cosTreshold = 2.0f / 3;
    second_icp_properties.method = MR::ICPMethod::PointToPlane;
    second_icp_properties.icpMode = MR::ICPMode::AnyRigidXf;
    second_icp_properties.badIterStopCount = 100;
    second_icp_properties.iterLimit = 300;
    second_icp_properties.farDistFactor = 1.5f;
    second_icp_properties.distThresholdSq = 15000.f;
    second_icp_properties.mutualClosest = false;

    MR::ICP second_icp(teeth_second_mop, cbct_mop, MR::AffineXf3f(), MR::AffineXf3f(), 0.3);
    second_icp.setParams(second_icp_properties);
    second_icp.updatePointPairs();

    auto teeth_second_transformation = second_icp.calculateTransformation();

    return teeth_second_transformation * second_translation * teeth_first_transformation * first_translation;
}

/**
 * @brief Convert an MR::AffineXf3f transformation to a Python tuple.
 *
 * @param aff The transformation matrix to convert.
 * @return pybind11::tuple A tuple containing the matrix and translation vector.
 */
pybind11::tuple CBCTJawRegistration::AffineXf3fToPyTuple(MR::AffineXf3f aff)
{
    // Convert MR::Matrix3f to NumPy array
    pybind11::array_t<float> A({ 3, 3 });
    auto A_buf = A.request();
    float* A_ptr = static_cast<float*>(A_buf.ptr);

    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            A_ptr[i * 3 + j] = aff.A[i][j];
        }
    }

    // Convert MR::Vector3f to NumPy array
    pybind11::array_t<float> b({ 3 });
    auto b_buf = b.request();
    float* b_ptr = static_cast<float*>(b_buf.ptr);

    for (size_t i = 0; i < 3; ++i)
    {
        b_ptr[i] = aff.b[i];
    }

    return pybind11::make_tuple(A, b);
}

/**
 * @brief Extract a sub-mesh from the given mesh using a face bit set.
 *
 * @param mesh The original mesh.
 * @param fbs The face bit set indicating which faces to include in the sub - mesh.
 * @return MR::Mesh The extracted sub - mesh.
 */
MR::Mesh CBCTJawRegistration::ExtractSubMeshByFbs(const MR::Mesh& mesh, const MR::FaceBitSet& fbs)
{
    std::vector<MR::Triangle3f> pos_triples;
    for (auto& f : mesh.topology.getValidFaces()) {
        if (fbs.test(f)) {
            MR::ThreeVertIds tri_verts;
            mesh.topology.getTriVerts(f, tri_verts);

            MR::Triangle3f triangle3f;
            triangle3f[0] = mesh.points[tri_verts[0]];
            triangle3f[1] = mesh.points[tri_verts[1]];
            triangle3f[2] = mesh.points[tri_verts[2]];

            pos_triples.push_back(triangle3f);
        }
    }
    return MR::Mesh::fromPointTriples(pos_triples, false);
}

/**
 * @brief Expose the CBCTJawRegistration class to Python using pybind11.
 *
 * This macro defines a Python module called `cbct_jaw_registration` and exposes
 * the CBCTJawRegistration class and its methods to Python.
 */
PYBIND11_MODULE(cbct_jaw_registration, m) {
    pybind11::class_<CBCTJawRegistration>(m, "CBCTJawRegistration")
        .def(pybind11::init<>())
        .def(pybind11::init<const std::wstring&, const std::wstring&, const std::wstring&, const std::wstring&, int, int, int, int, int, int, float, const std::wstring&>(),
            pybind11::arg("upper_teeth_file"),
            pybind11::arg("lower_teeth_file"),
            pybind11::arg("tiff_folder"),
            pybind11::arg("upper_jaw_file"),
            pybind11::arg("lower_x"),
            pybind11::arg("upper_x"),
            pybind11::arg("lower_y"),
            pybind11::arg("upper_y"),
            pybind11::arg("lower_z"),
            pybind11::arg("upper_z"),
            pybind11::arg("voxel_space") = 0.3f,
            pybind11::arg("temp") = "")
        .def("Registration", &CBCTJawRegistration::Registration)
        .def("GetUpperTeethVertNum", &CBCTJawRegistration::GetUpperTeethVertNum)
        .def("GetUpperTransformation", &CBCTJawRegistration::GetUpperTransformation)
        .def("GetLowerTransformation", &CBCTJawRegistration::GetLowerTransformation);
}
