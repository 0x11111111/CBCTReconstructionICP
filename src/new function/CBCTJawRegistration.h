#include <MRMesh/MRPython.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRAffineXf.h>
#include <filesystem>
#include <MRMesh/MRObjectVoxels.h>


namespace fs = std::filesystem;

class CBCTJawRegistration
{
public:
    MR::Mesh m_upper_teeth_mr;
    MR::Mesh m_lower_teeth_mr;
    MR::Mesh m_upper_jaw_mr;
    MR::Mesh m_cbct_upper_mr;
    MR::Mesh m_cbct_lower_mr;
    MR::Mesh m_cbct_all_mr;
    MR::AffineXf3f m_upper_teeth_transformation;
    MR::AffineXf3f m_lower_teeth_transformation;
    MR::ObjectVoxels m_object_voxels;
    fs::path m_temp_folder;

    CBCTJawRegistration() = default;
    CBCTJawRegistration(
        const std::wstring& upper_teeth_file,
        const std::wstring& lower_teeth_file,
        const std::wstring& tiff_folder,
        const std::wstring& upper_jaw_file,
        const int lower_x,
        const int upper_x,
        const int lower_y,
        const int upper_y,
        const int lower_z,
        const int upper_z,
        const float voxel_space = 0.3f,
        const std::wstring& temp = std::wstring());


    int GetUpperTeethVertNum();
    void Registration();
    pybind11::tuple GetUpperTransformation();
    pybind11::tuple GetLowerTransformation();

private:
    MR::Mesh LoadMeshFromFile(const std::wstring& file_path);
    MR::AffineXf3f ICPHalf(MR::Mesh cbct_mr, MR::Mesh teeth_mr, bool upper = true);
    pybind11::tuple AffineXf3fToPyTuple(MR::AffineXf3f aff);
    MR::Mesh ExtractSubMeshByFbs(const MR::Mesh& mesh, const MR::FaceBitSet& fbs);
};