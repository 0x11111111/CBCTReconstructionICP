#include "stdafx.h"
#include "vtkRenderPipeline.h"
#include "meshTransform.h"
#include "simpleRender.h"
#include "IOManip.hpp"
#include <MRMesh/MRMeshLoad.h>
#include <MRMesh/MRId.h>
#include <MRMesh/MRMesh.h>
#include <MRMesh/MRBitSetParallelFor.h>
#include <MRMesh/MRMeshTopology.h>
#include <MRMesh/MRExpected.h>

#include <MRMesh/MRVoxelsLoad.h>
#include <MRMesh/MRObjectVoxels.h>
#include <MRMesh/MRBox.h>
#include <MRMesh/MRAffineXf.h>

#include <filesystem>

vtkRenderPipeline* pipeline;

/**
 * @brief Generate a 4-digit number string with leading zeros.
 *
 * This function takes an integer and converts it into a string representation
 * with a fixed width of 4 characters. If the number has less than 4 digits,
 * leading zeros are added to pad the string to the desired width.
 *
 * @param number The input integer to be converted.
 * @return A string representation of the input number with leading zeros.
 *
 * @note The function uses std::ostringstream, std::setw(), and std::setfill()
 *       to format the output string.
 *
 * @example
 *   int num = 42;
 *   std::string num_str = generate_leading_zero_number_str(num);
 *   // num_str will be "0042"
 */
std::string generate_leading_zero_number_str(int number)
{
	std::ostringstream stream;
	stream << std::setw(4) << std::setfill('0') << number;
	return stream.str();
}

void RightPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Press" << endl;
}
void RightRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Right Released" << endl;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	//std::cout<<"Left Press" << endl;
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}
	
namespace fs = std::filesystem;

int main()
{
	using namespace MR;
	pipeline = new vtkRenderPipeline();

	fs::path base_folder_path = "data/02";
	fs::path cbct_folder_path = base_folder_path / "ct.nii";
	fs::path intraoral_folder_path = base_folder_path / "口扫";

	fs::path cbct_file_path = cbct_folder_path / "ct.nii";
	fs::path dicom_file_path = base_folder_path / "dicom";
	fs::path upper_jaw_file_path = intraoral_folder_path / "UpperJaw.stl";
	fs::path lower_jaw_file_path = intraoral_folder_path / "LowerJaw.stl";


	auto voxel_load = VoxelsLoad::loadDCMFolder(dicom_file_path);
	//std::cout << voxel_load.error() << std::endl;
	auto vdb = voxel_load.value().vdbVolume;
	std::cout << vdb.voxelSize[0] << vdb.voxelSize[1] << vdb.voxelSize[2] << std::endl;

	auto dims = vdb.dims;

	auto dim_x = dims.x;
	auto dim_y = dims.y;
	auto dim_z = dims.z;

	auto max_x = dim_x / 2;
	auto max_y = dim_y;
	auto max_z = dim_z / 16 * 5;
	auto min_z = dim_z / 16 * 3;

	auto active_bounds = Box3i(Vector3i(0, 0, min_z), Vector3i(max_x, max_y, max_z));

	ObjectVoxels obj_voxels;
	obj_voxels.construct(vdb);
	obj_voxels.applyScale(1.0f / vdb.voxelSize[0] / 3.5f);
	obj_voxels.setActiveBounds(active_bounds);
	obj_voxels.setDualMarchingCubes(true);
	obj_voxels.setIsoValue(.444f);
	Mesh cbct_recon_mr(*obj_voxels.mesh());

	auto rot_matrix = Matrix3f::rotation(Vector3f{ 0, 0, 1. }, PI / 2);
	auto rot_xf = AffineXf3f::linear(rot_matrix);

	cbct_recon_mr.transform(rot_xf);
	MeshSave::toAnySupportedFormat(cbct_recon_mr, cbct_folder_path / "cbct_recon_mr.ply");



	// Set up the camera and interactor.
	//pipeline->Renderer->GetActiveCamera()->SetParallelProjection(1);
	//pipeline->Renderer->ResetCamera();
	//// Set up the callback functions for mouse events.
	//pipeline->addObserver(vtkCommand::LeftButtonPressEvent, LeftPress);
	//pipeline->addObserver(vtkCommand::MouseMoveEvent, MouseMove);
	//pipeline->addObserver(vtkCommand::LeftButtonReleaseEvent, LeftRelease);
	//pipeline->addObserver(vtkCommand::RightButtonPressEvent, RightPress);
	//pipeline->addObserver(vtkCommand::RightButtonReleaseEvent, RightRelease);

	//pipeline->RenderWindowInteractor->Start();

	// Clean up the pipeline after each run.
	delete pipeline;
}