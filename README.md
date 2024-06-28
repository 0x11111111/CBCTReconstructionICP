# CBCTReconstructionICP: Reconstruction of CBCT and ICP Alignment of Jaws by MeshLib

## This is a project dedicated to compiling pyd module for [this python project](https://github.com/0x11111111/MeshSegNet-DentalArchAligner)

### This project is derived from [this branch](https://github.com/0x11111111/VisualLab/tree/base_meshlib_lab)

### Instructions for deployment (Only Windows and MSVC supported)

1. git clone [This Project](../../)
2. Launching project via "VTK_CGAL_Libigl.sln".
3. Switch the "Solution Configurations" to Release.
4. RUN (the project itself).

### Core

In [this Class CBCTJawRegistration](https://github.com/0x11111111/CBCTReconstructionICP/blob/main/src/new%20function/CBCTJawRegistration.cpp). This class exposes the following to Python:

1. Constructor

```C++
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
        const std::wstring& temp = std::wstring())
```

2. `void Registration()`

3. `pybind11::tuple GetUpperTransformation()`

4. `pybind11::tuple GetLowerTransformation()`

### Acknowledgement

1. This Project is designed and compiled by [Schwarz Solomon](https://github.com/SchwarzSolomon).
2. Tidied up and uploaded by [0x11111111](https://github.com/0x11111111).
