#include "stdafx.h"
//#include <pybind11/pybind11.h>
//#include "vtkRenderPipeline.h"
//#include "meshTransform.h"
//#include "simpleRender.h"
//#include "IOManip.hpp"
//#include <MRMesh/MRMeshLoad.h>
//#include <MRMesh/MRId.h>
//#include <MRMesh/MRMesh.h>
//#include <MRMesh/MRBitSetParallelFor.h>
//#include <MRMesh/MRMeshTopology.h>
//#include <MRMesh/MRExpected.h>
//
//#include <MRMesh/MRVoxelsLoad.h>
//#include <MRMesh/MRObjectVoxels.h>
//#include <MRMesh/MRBox.h>
//#include <MRMesh/MRAffineXf.h>
//#include <MRMesh/MRMeshComponents.h>
//#include <MRMesh/MRICP.h>
//#include <MRMesh/MRMeshOrPoints.h>
//#include <MRMesh/MRMeshBoolean.h>
//#include <MRMesh/MREdgePaths.h>
//#include <MRMesh/MRRegionBoundary.h>
//#include <MRMesh/MRFixSelfIntersections.h>
//#include <MRMesh/MRMeshIntersect.h>
//#include <MRMesh/MRLine.h>
//#include <MRMesh/MRMeshDecimate.h>
//
//#include <CGAL/mesh_segmentation.h>
//
//
//#include <filesystem>
//#include <algorithm>
//
//vtkRenderPipeline* pipeline;
//
///**
// * @brief Generate a 4-digit number string with leading zeros.
// *
// * This function takes an integer and converts it into a string representation
// * with a fixed width of 4 characters. If the number has less than 4 digits,
// * leading zeros are added to pad the string to the desired width.
// *
// * @param number The input integer to be converted.
// * @return A string representation of the input number with leading zeros.
// *
// * @note The function uses std::ostringstream, std::setw(), and std::setfill()
// *       to format the output string.
// *
// * @example
// *   int num = 42;
// *   std::string num_str = generate_leading_zero_number_str(num);
// *   // num_str will be "0042"
// */
//std::string generate_leading_zero_number_str(int number)
//{
//	std::ostringstream stream;
//	stream << std::setw(4) << std::setfill('0') << number;
//	return stream.str();
//}
//
//void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//	//std::cout<<"Right Press" << endl;
//}
//void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//	//std::cout<<"Right Released" << endl;
//}
//
//void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//	//std::cout<<"Left Press" << endl;
//}
//
//void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//	
//}
//
//void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
//{
//}
//	
//namespace fs = std::filesystem;
//
//auto ExtractSubMeshByFbs = [&](const Mesh& mesh_copy, const FaceBitSet& fbs) {
//	std::vector<Triangle3f> pos_triples;
//	for (auto& f : mesh_copy.topology.getValidFaces()) {
//		if (fbs.test(f)) {
//			ThreeVertIds tri_verts;
//			mesh_copy.topology.getTriVerts(f, tri_verts);
//
//			Triangle3f triangle3f;
//			triangle3f[0] = mesh_copy.points[tri_verts[0]];
//			triangle3f[1] = mesh_copy.points[tri_verts[1]];
//			triangle3f[2] = mesh_copy.points[tri_verts[2]];
//
//			pos_triples.push_back(triangle3f);
//		}
//	}
//	return Mesh::fromPointTriples(pos_triples, false);
//	};
//
//Mesh SdfFilterJaws(Mesh jaws, double sdf_min_threshold)
//{
//	auto mesh = MRMeshToSurfaceMesh(jaws);
//	SurfaceMesh mesh_copy = mesh;
//
//	// 计算SDF值
//	typedef SurfaceMesh::Property_map<face_descriptor, double> Facet_double_map;
//	Facet_double_map sdf_property_map = mesh.add_property_map<face_descriptor, double>("f:sdf").first;
//	CGAL::sdf_values(mesh, sdf_property_map);
//
//	// 获取SDF值的最小值和最大值，用于归一化
//	double min_sdf = std::numeric_limits<double>::max();
//	double max_sdf = std::numeric_limits<double>::lowest();
//	for (face_descriptor fd : faces(mesh)) {
//		double sdf = sdf_property_map[fd];
//		if (sdf < min_sdf) min_sdf = sdf;
//		if (sdf > max_sdf) max_sdf = sdf;
//	}
//
//	// 为面创建颜色属性映射
//	typedef SurfaceMesh::Property_map<face_descriptor, CGAL::Color> Facet_color_map;
//	Facet_color_map color_property_map = mesh.add_property_map<face_descriptor, CGAL::Color>("f:color").first;
//
//	SurfaceMesh filtered_mesh;
//	std::map<SurfaceMesh::Vertex_index, SurfaceMesh::Vertex_index> vertex_map;
//
//	// 将normalized_sdf低于0.3的面移除
//	for (face_descriptor fd : faces(mesh)) {
//		double sdf = sdf_property_map[fd];
//		double normalized_sdf = (sdf - min_sdf) / (max_sdf - min_sdf);
//		if (normalized_sdf >= 0.3) {
//			std::vector<SurfaceMesh::Vertex_index> vertices;
//			for (auto vd : vertices_around_face(mesh.halfedge(fd), mesh)) {
//				auto it = vertex_map.find(vd);
//				if (it == vertex_map.end()) {
//					SurfaceMesh::Vertex_index new_vd = filtered_mesh.add_vertex(mesh.point(vd));
//					vertex_map[vd] = new_vd;
//					vertices.push_back(new_vd);
//				}
//				else {
//					vertices.push_back(it->second);
//				}
//			}
//			filtered_mesh.add_face(vertices);
//		}
//	}
//
//	// 根据SDF值设置面颜色
//	for (face_descriptor fd : faces(mesh)) {
//		double sdf = sdf_property_map[fd];
//		double normalized_sdf = (sdf - min_sdf) / (max_sdf - min_sdf);
//		unsigned char red = static_cast<unsigned char>(255 * (1.0 - normalized_sdf));
//		unsigned char blue = static_cast<unsigned char>(255 * normalized_sdf);
//		color_property_map[fd] = CGAL::Color(red, 0, blue);
//	}
//
//	//std::string output_file = "data/sdf_reconstructed_surface_cut_down_largest_component.ply";
//	//std::ofstream output(output_file, std::ios::binary);
//	//CGAL::IO::write_PLY(output, mesh, CGAL::parameters::stream_precision(17));
//
//	std::string output_filtered_file = "data/sdf_reconstructed_surface_cut_down_largest_component_filtered.ply";
//	auto filtered_mesh_mr = SurfaceMeshToMRMesh(filtered_mesh);
//
//	auto all_components = MeshComponents::getAllComponents(filtered_mesh_mr);
//	std::sort(all_components.begin(), all_components.end(), [](const FaceBitSet& a, const FaceBitSet& b) { return (a.count() > b.count()); });
//	for (auto& fbs : all_components) {
//		std::cout << fbs.count() << std::endl;
//	}
//	//auto filtered_mesh_mr_fbs = MeshComponents::getLargestComponent(filtered_mesh_mr);
//	auto largest_two_components = all_components[0] | all_components[1];
//	auto filtered_mesh_mr_largest_component = ExtractSubMeshByFbs(filtered_mesh_mr, largest_two_components);
//	
//	return filtered_mesh_mr_largest_component;
//}
//
//AffineXf3f ICPHalf(Mesh cbct_mr, Mesh teeth_mr, bool upper = true)
//{
//	auto cbct_bound = cbct_mr.computeBoundingBox();
//	auto teeth_bound = teeth_mr.computeBoundingBox();
//
//	auto& cbct_min = cbct_bound.min;
//	auto& cbct_max = cbct_bound.max;
//
//	auto& teeth_min = teeth_bound.min;
//	auto& teeth_max = teeth_bound.max;
//
//	auto cbct_xz_panel = Box3f({ cbct_min[0], cbct_min[1], cbct_min[2] }, { cbct_max[0], cbct_min[1], cbct_max[2] });
//	auto teeth_xz_panel = Box3f({ teeth_min[0], teeth_min[1], teeth_min[2] }, { teeth_max[0], teeth_min[1], teeth_max[2] });
//
//	auto cbct_xz_panel_center = cbct_xz_panel.center();
//	auto teeth_xz_panel_center = teeth_xz_panel.center();
//
//	// Align the upper of the teeth to the upper of the CBCT
//
//
//	auto z_cbct_dim = cbct_max[2] - cbct_min[2];
//	auto z_teeth_dim = teeth_max[2] - teeth_min[2];
//	auto z_translation_dist = (z_cbct_dim - z_teeth_dim) / 2 * (upper ? 1 : -1);
//
//	auto y_translation_dist = -1.f;
//
//	auto initial_translation = cbct_xz_panel_center - teeth_xz_panel_center;
//	initial_translation += { 0, y_translation_dist, z_translation_dist };
//	auto initial_transformation = AffineXf3f::translation(initial_translation);
//
//	auto teeth_init_translation_mr(teeth_mr);
//	teeth_init_translation_mr.transform(initial_transformation);
//
//	MeshOrPoints cbct_mop(cbct_mr);
//	MeshOrPoints teeth_mop(teeth_init_translation_mr);
//
//	ICPProperties icp_properties;
//	icp_properties.cosTreshold = -1.f;
//	icp_properties.method = ICPMethod::Combined;
//	icp_properties.icpMode = ICPMode::AnyRigidXf;
//	icp_properties.badIterStopCount = 100;
//	icp_properties.iterLimit = 300;
//	icp_properties.farDistFactor = 1.3f;
//	icp_properties.distThresholdSq = 1000.f;
//	icp_properties.mutualClosest = false;
//
//	ICP icp(teeth_mop, cbct_mop, AffineXf3f(), AffineXf3f(), 0.3);
//	icp.setParams(icp_properties);
//	icp.updatePointPairs();
//
//	auto teeth_transformation = icp.calculateTransformation();
//
//	teeth_init_translation_mr.transform(teeth_transformation);
//
//	return teeth_transformation * initial_transformation;
//}
//
//int main()
//{
//	using namespace MR;
//	pipeline = new vtkRenderPipeline();
//
//	fs::path base_folder_path = "data/04";
//	fs::path cbct_folder_path = base_folder_path / "ct.nii";
//	fs::path intraoral_folder_path = base_folder_path / "口扫";
//
//	fs::path cbct_file_path = cbct_folder_path / "ct.nii";
//	fs::path dicom_file_path = base_folder_path / "dicom";
//	fs::path tiff_file_path = base_folder_path / "tiff";
//	fs::path upper_jaw_file_path = intraoral_folder_path / "UpperJaw.stl";
//	fs::path lower_jaw_file_path = intraoral_folder_path / "LowerJaw.stl";
//
//	auto voxel_space = 0.3f;
//	auto voxel_size = Vector3f(voxel_space, voxel_space, voxel_space);
//	//auto voxel_load = VoxelsLoad::loadDCMFolder(dicom_file_path);
//	auto voxel_load = VoxelsLoad::loadTiffDir(
//		{
//			tiff_file_path,
//			voxel_size,
//			VoxelsLoad::GridType::DenseGrid,
//			{}
//		}
//	);
//	//std::cout << voxel_load.error() << std::endl;
//	auto vdb = voxel_load.value();
//	auto dims = vdb.dims;
//
//	auto dim_x = dims.x;
//	auto dim_y = dims.y;
//	auto dim_z = dims.z;
//
//	auto min_x = 0;
//	auto max_x = dim_x;
//	auto min_y = 0;
//	auto max_y = dim_y / 2;
//	auto min_z = 83;
//	auto max_z = 312;
//	auto middle_z = (min_z + max_z) / 2;
//
//	auto upper_active_bounds = Box3i({ min_x, min_y, static_cast<int>((max_z - min_z) * 0.4) + min_z }, { max_x, max_y, max_z });
//	auto lower_active_bounds = Box3i({ min_x, min_y, min_z }, { max_x, max_y, static_cast<int>((max_z - min_z) * 0.6) + min_z });
//
//
//	ObjectVoxels obj_voxels;
//	obj_voxels.construct(vdb);
//	auto& bins = obj_voxels.histogram().getBins();
//
//	auto sum = std::accumulate(bins.begin(), bins.end(), 0);
//	std::cout << sum << std::endl;
//	
//	std::cout << bins.size() << std::endl;
//	auto rit = bins.rbegin();
//	auto accumulation = 0;
//
//	auto it = bins.begin();
//	for (; it != bins.end(); ++it)
//	{
//		accumulation += *it;
//		if (accumulation > sum * 0.99)
//		{
//			break;
//		}
//	}
//
//	auto bin_id = std::distance(bins.begin(), it);
//	std::cout << bin_id << std::endl;
//	std::cout << std::endl;
//	auto& min_max = obj_voxels.histogram().getBinMinMax(bins.size() / 3);
//	std::cout << min_max.first << " " << min_max.second << std::endl;
//	//obj_voxels.applyScale(1.0f / vdb.voxelSize[0] / 3.5f);
//	//obj_voxels.setDualMarchingCubes(true);
//	obj_voxels.setIsoValue(bin_id);
//	Mesh cbct_recon_full_mr(*obj_voxels.mesh());
//	MeshSave::toAnySupportedFormat(cbct_recon_full_mr, cbct_folder_path / "cbct_recon_full.ply");
//
//	obj_voxels.setActiveBounds(upper_active_bounds);
//	obj_voxels.setIsoValue(bin_id);
//	Mesh cbct_upper_mr(*obj_voxels.mesh());
//	MeshSave::toAnySupportedFormat(cbct_upper_mr, cbct_folder_path / "cbct_upper.ply");
//
//	obj_voxels.setActiveBounds(lower_active_bounds);
//	obj_voxels.setIsoValue(bin_id);
//	Mesh cbct_lower_mr(*obj_voxels.mesh());
//	MeshSave::toAnySupportedFormat(cbct_lower_mr, cbct_folder_path / "cbct_lower.ply");
//
//
//	auto upper_teeth = MeshLoad::fromAnySupportedFormat(intraoral_folder_path / "upper_teeth.ply").value();
//	auto upper_transformation = ICPHalf(cbct_upper_mr, upper_teeth, true);
//	upper_teeth.transform(upper_transformation);
//	MeshSave::toAnySupportedFormat(upper_teeth, cbct_folder_path / "upper_teeth_icp_transformed.ply");
//
//	auto lower_teeth = MeshLoad::fromAnySupportedFormat(intraoral_folder_path / "lower_teeth.ply").value();
//	auto lower_transformation = ICPHalf(cbct_lower_mr, lower_teeth, false);
//	lower_teeth.transform(lower_transformation);
//	MeshSave::toAnySupportedFormat(lower_teeth, cbct_folder_path / "lower_teeth_icp_transformed.ply");
//
//
//	//upper_teeth_icp_transformation.b.z += 2.f;
//	//auto pair_num = upper_teeth_icp.getNumActivePairs();
//	//std::cout << "Number of active pairs: " << pair_num << std::endl;
//
//	//auto rms_p2p_dist = upper_teeth_icp.getMeanSqDistToPoint();
//	//auto last_iter = upper_teeth_icp.getLastICPInfo();
//	//std::cout << "RMS point-to-point distance: " << rms_p2p_dist << std::endl;
//	//std::cout << "Last iteration: " << last_iter << std::endl;
//
//
//	//Mesh upper_teeth_icp_transformed_mr(upper_teeth_initial_translated_mr);
//	//upper_teeth_icp_transformed_mr.transform(upper_teeth_icp_transformation);
//
//	//MeshSave::toAnySupportedFormat(upper_teeth_icp_transformed_mr, cbct_folder_path / "upper_teeth_icp_fisrt_transformed.ply");
//	//MeshOrPoints upper_second_teeth_mop(upper_teeth_icp_transformed_mr);
//
//	//ICPProperties upper_teeth_second_icp_properties;
//	//upper_teeth_second_icp_properties.cosTreshold = -1.f;
//	//upper_teeth_second_icp_properties.method = ICPMethod::Combined;
//	//upper_teeth_second_icp_properties.icpMode = ICPMode::AnyRigidXf;
//	//upper_teeth_second_icp_properties.badIterStopCount = 100;
//	//upper_teeth_second_icp_properties.iterLimit = 300;
//	//upper_teeth_second_icp_properties.farDistFactor = 1.1f;
//	//upper_teeth_second_icp_properties.distThresholdSq = 100.f;
//	//upper_teeth_second_icp_properties.mutualClosest = false;
//
//	//ICP upper_teeth_second_icp(upper_second_teeth_mop, cbct_recon_mop, AffineXf3f(), AffineXf3f(), 0.3);
//	//auto upper_teeth_second_icp_transformation = upper_teeth_second_icp.calculateTransformation();
//	//upper_teeth_icp_transformed_mr.transform(upper_teeth_second_icp_transformation);
//
//	//MeshSave::toAnySupportedFormat(cbct_recon_largest_component_mr, cbct_folder_path / "cbct_recon_mr.ply");
//	//MeshSave::toAnySupportedFormat(upper_teeth_icp_transformed_mr, cbct_folder_path / "upper_teeth_icp_transformed.ply");
//
//	//auto upper_jaw_mr = MeshLoad::fromAnySupportedFormat(intraoral_folder_path / "UpperJaw.stl").value();
//	//upper_jaw_mr.transform(initial_transformation);
//	//upper_jaw_mr.transform(upper_teeth_icp_transformation);
//	//upper_jaw_mr.transform(upper_teeth_second_icp_transformation);
//	//MeshSave::toAnySupportedFormat(upper_jaw_mr, cbct_folder_path / "upper_jaw_icp_transformed.ply");
//
//
//	//auto lower_teeth_mr = MeshLoad::fromAnySupportedFormat(intraoral_folder_path / "lower_teeth.ply").value();
//	//lower_teeth_mr.transform(initial_transformation);
//	//lower_teeth_mr.transform(upper_teeth_icp_transformation);
//	//lower_teeth_mr.transform(upper_teeth_second_icp_transformation);
//
//	//auto lower_teeth_max_z = cbct_min[2] + cbct_z_dim * 0.6;
//
//	//VertBitSet lower_teeth_range_vbs(cbct_recon_largest_component_mr.topology.numValidVerts());
//	//for (auto& v : cbct_recon_largest_component_mr.topology.getValidVerts())
//	//{
//	//	if (cbct_recon_largest_component_mr.points[v][2] < lower_teeth_max_z)
//	//	{
//	//		lower_teeth_range_vbs.set(v);
//	//	}
//	//}
//	//auto lower_teeth_range_fbs = getInnerFaces(cbct_recon_largest_component_mr.topology, lower_teeth_range_vbs);
//	//auto recon_lower_range_mr = ExtractSubMeshByFbs(cbct_recon_largest_component_mr, lower_teeth_range_fbs);
//
//	//MeshSave::toAnySupportedFormat(recon_lower_range_mr, cbct_folder_path / "recon_lower_range_mr.ply");
//
//
//	//MeshOrPoints lower_cbct_recon_mop(recon_lower_range_mr);
//	//MeshOrPoints lower_teeth_mop(lower_teeth_mr);
//
//	////ICPProperties lower_teeth_icp_properties;
//	////lower_teeth_icp_properties.cosTreshold = -1.f;
//	////lower_teeth_icp_properties.method = ICPMethod::PointToPlane;
//	////lower_teeth_icp_properties.icpMode = ICPMode::AnyRigidXf;
//	////lower_teeth_icp_properties.badIterStopCount = 100;
//	////lower_teeth_icp_properties.iterLimit = 200;
//	////lower_teeth_icp_properties.farDistFactor = 1.5f;
//	////lower_teeth_icp_properties.distThresholdSq = 88.f;
//	////lower_teeth_icp_properties.mutualClosest = false;
//
//	//ICP lower_teeth_icp(lower_teeth_mop, lower_cbct_recon_mop, AffineXf3f(), AffineXf3f(), 0.3);
//	//lower_teeth_icp.setParams(upper_teeth_second_icp_properties);
//
//	//auto lower_teeth_icp_transformation = lower_teeth_icp.calculateTransformation();
//
//	//lower_teeth_mr.transform(lower_teeth_icp_transformation);
//	//MeshSave::toAnySupportedFormat(lower_teeth_mr, cbct_folder_path / "lower_teeth_icp_transformed.ply");
//
//
//	//auto lower_jaw_mr = MeshLoad::fromAnySupportedFormat(intraoral_folder_path / "LowerJaw.stl").value();
//	//lower_jaw_mr.transform(initial_transformation);
//	//lower_jaw_mr.transform(upper_teeth_icp_transformation);
//	//lower_jaw_mr.transform(upper_teeth_second_icp_transformation);
//	//lower_jaw_mr.transform(lower_teeth_icp_transformation);
//	//MeshSave::toAnySupportedFormat(lower_jaw_mr, cbct_folder_path / "lower_jaw_icp_transformed.ply");
//
//	//SelfIntersections::fix(recon_lower_range_mr, SelfIntersections::Settings(
//	//	{
//	//		SelfIntersections::Settings::Method::CutAndFill,
//	//		2,
//	//		3,
//	//		0.0f
//	//	}
//	//));
//	//SelfIntersections::fix(upper_teeth_icp_transformed_mr, SelfIntersections::Settings(
//	//	{
//	//		SelfIntersections::Settings::Method::CutAndFill,
//	//		2,
//	//		3,
//	//		0.0f
//	//	}
//	//));
//
//	//AffineXf3f identity = AffineXf3f();
//	//auto boolean_res = MR::boolean(recon_lower_range_mr, upper_teeth_icp_transformed_mr, BooleanOperation::DifferenceBA, &identity);
//	//std::cout << boolean_res.errorString << std::endl;
//	//auto recon_lower_boolean_mr = boolean_res.mesh;
//
//	////MeshSave::toAnySupportedFormat(recon_lower_boolean_mr, cbct_folder_path / "recon_lower_boolean_mr.ply");
//
//	//auto lower_teeth_remaining_vbs(recon_lower_range_mr.topology.getValidVerts());
//	//auto negative_self_detection_length = 20.f;
//	//auto positive_detection_length = 0.0f;
//	//for (auto& v : recon_lower_range_mr.topology.getValidVerts())
//	//{
//	//	auto norm = recon_lower_range_mr.normal(v);
//
//	//	auto self_intersection_segment = Line3f{ recon_lower_range_mr.points[v], recon_lower_range_mr.points[v] - norm * negative_self_detection_length };
//	//	auto self_intersection_res = rayMeshIntersect(recon_lower_range_mr, self_intersection_segment);
//	//	if (self_intersection_res.has_value())
//	//	{
//	//		auto start_point = recon_lower_range_mr.points[v] + positive_detection_length * norm;
//	//		auto end_point = self_intersection_res.value().proj.point;
//	//		auto intersection_res = rayMeshIntersect(upper_teeth_icp_transformed_mr, { start_point, end_point });
//	//		if (intersection_res.has_value())
//	//		{
//	//			lower_teeth_remaining_vbs.reset(v);
//	//		}
//	//	}
//	//}
//
//	//auto lower_teeth_remaining_fbs = getInnerFaces(recon_lower_range_mr.topology, lower_teeth_remaining_vbs);
//	//auto lower_teeth_remaining_mr = ExtractSubMeshByFbs(recon_lower_range_mr, lower_teeth_remaining_fbs);
//
//	//MeshSave::toAnySupportedFormat(lower_teeth_remaining_mr, cbct_folder_path / "lower_teeth_remaining_mr.ply");
//
//	//FaceBitSet lower_teeth_remaining_fbs(recon_lower_range_mr.topology.getValidFaces());
//	//for (auto& v : lower_jaw_mr.topology.getValidVerts())
//	//{
//	//	auto norm = lower_jaw_mr.normal(v);
//	//	auto positive_detection_length = 0.1f;
//	//	auto negative_detection_length = 5.f;
//	//	auto segment3f = Line3f{ lower_jaw_mr.points[v] + norm * positive_detection_length, lower_jaw_mr.points[v] - norm * negative_detection_length };
//	//	auto intersection_res = rayMeshIntersect(recon_lower_range_mr, segment3f);
//	//	if (intersection_res.has_value())
//	//	{
//	//		lower_teeth_remaining_fbs.reset(intersection_res.value().proj.face);
//	//	}
//	//}
//
//	//auto lower_teeth_remaining_mr = ExtractSubMeshByFbs(recon_lower_range_mr, lower_teeth_remaining_fbs);
//	//MeshSave::toAnySupportedFormat(lower_teeth_remaining_mr, cbct_folder_path / "lower_teeth_remaining_mr.ply");
//
//	// Set up the camera and interactor.
//	//pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
//	//pipeline->Renderer->ResetCamera();
//	//// Set up the callback functions for mouse events.
//	//pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
//	//pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
//	//pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
//	//pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
//	//pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);
//
//	//pipeline->RenderWindowInteractor->Start();
//
//	// Clean up the pipeline after each run.
//	delete pipeline;
//}
