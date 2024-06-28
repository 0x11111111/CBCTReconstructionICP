import os
from teeth_segmentation import extract_teeth_from_jaw

def registration(upper_jaw_mesh: str, lower_jaw_mesh: str, cbct_path: str, temp_path: str):
    if not os.path.exists(temp_path):
        os.makedirs(temp_path)
    
    print('upper_jaw_mesh:', upper_jaw_mesh)

    output_upper = os.path.join(temp_path, 'upper_teeth.ply')
    output_lower = os.path.join(temp_path, 'lower_teeth.ply')
    extract_teeth_from_jaw(upper_jaw=upper_jaw_mesh, lower_jaw=lower_jaw_mesh, output_upper_teeth=output_upper, output_lower_teeth=output_lower)


if __name__ == '__main__':
    data_rel_path = 'data/01/'
    data_abs_path = os.path.abspath(data_rel_path)
    print(data_abs_path)

    temp_path = os.path.join(data_abs_path, 'temp')
    upper_jaw_path = os.path.join(data_abs_path, '口扫/UpperJaw.stl')
    lower_jaw_path = os.path.join(data_abs_path, '口扫/LowerJaw.stl')
    cbct_path = os.path.join(data_abs_path, 'ct.nii/ct.nii')

    print(temp_path)

    registration(upper_jaw_path, lower_jaw_path, cbct_path, temp_path)
