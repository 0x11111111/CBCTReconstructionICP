# Test cbct_jaw_registration module
from meshlib import mrmeshpy as mm

import cbct_jaw_registration
import os
os.environ['YOLO_VERBOSE'] = str(False)
from ultralytics import YOLOv10
import cv2
import numpy as np
import shutil
from pathlib import Path
from tqdm import tqdm
import glob

data_rel_path = '../../data/01'.replace('/', os.path.sep)
data_abs_path = os.path.abspath(data_rel_path)
intraoral_scan_path = os.path.join(data_abs_path, '口扫')
cbct_path = os.path.join(data_abs_path, 'ct.nii')
tiff_foler = os.path.join(data_abs_path, 'tiff')
adaptive_thres_dir = os.path.join(data_abs_path, 'adaptive_threshold_jpg')

print(data_abs_path)

def apply_adaptive_threshold_and_filter(input_dir, output_dir):

    # 检查目录是否存在
    if os.path.exists(output_dir):
    # 删除目录及其内容
        shutil.rmtree(output_dir)
        print(f"目录 '{output_dir}' 已删除")
    else:
        print(f"目录 '{output_dir}' 不存在")
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)

    # 获取输入目录中的所有图像文件
    image_files = [f for f in os.listdir(input_dir) if os.path.join(input_dir, f).endswith('.tiff')]
    
    for img_name in image_files:
        input_path = os.path.join(input_dir, img_name)
        output_path = os.path.join(output_dir, img_name.replace('.tiff', '.jpg'))
        
        # 读取图像
        image = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)  # 读取为灰度图像
        if image is None:
            print(f"无法读取图像: {input_path}")
            continue
        
        # 应用自适应阈值处理
        adaptive_thresh = cv2.adaptiveThreshold(
            image, 255, cv2.ADAPTIVE_THRESH_MEAN_C, 
            cv2.THRESH_BINARY_INV, 17, 2
        )
        
        # 定义形态学操作的内核
        kernel = np.ones((3, 3), np.uint8)
        
        # 进行形态学开运算去除噪点
        morph_open = cv2.morphologyEx(adaptive_thresh, cv2.MORPH_OPEN, kernel)
        
        # 进行形态学闭运算填充线段间的空隙
        morph_close = cv2.morphologyEx(morph_open, cv2.MORPH_CLOSE, kernel)

        # 连通组件分析去除小的联通组件
        num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(morph_close, connectivity=8)
        min_area = 1100  # 设定最小区域面积阈值
        
        filtered_image = np.zeros_like(morph_close)
        
        for i in range(1, num_labels):  # 从1开始，跳过背景
            if stats[i, cv2.CC_STAT_AREA] >= min_area:
                filtered_image[labels == i] = 255
        
        # 保存处理后的图像
        cv2.imwrite(output_path, filtered_image)
        print(f"已处理并保存图像: {output_path}")


# 应用自适应阈值处理和过滤
apply_adaptive_threshold_and_filter(tiff_foler, adaptive_thres_dir)


os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

# Train the model using the 'coco8.yaml' dataset for 3 epochs
model = YOLOv10('D:/Code/CBCT_Reconstruction_ICP/runs/detect/train14/weights/best.pt')

def find_valid_indices(positive_indices):
    n = len(positive_indices)
    if n < 4:
        return [], None, None  # 如果列表长度小于4，无法找到四个连续的数值

    # 找到第一个包含四个连续递增数值的子序列
    start = None
    for i in range(n - 3):
        if (positive_indices[i + 1] == positive_indices[i] + 1 and
                positive_indices[i + 2] == positive_indices[i] + 2 and
                positive_indices[i + 3] == positive_indices[i] + 3):
            start = i
            break

    # 找到第一个包含四个连续递减数值的子序列
    end = None
    for i in range(n - 1, 2, -1):
        if (positive_indices[i] == positive_indices[i - 1] + 1 and
                positive_indices[i - 1] == positive_indices[i - 2] + 1 and
                positive_indices[i - 2] == positive_indices[i - 3] + 1):
            end = i
            break

    if start is not None and end is not None and start < end:
        complete_subsequence = positive_indices[start:end + 1]
        return complete_subsequence, min(complete_subsequence), max(complete_subsequence)
    else:
        return [], None, None

def plot_predictions(img, preds, output_path):
    for pred in preds:
        x1, y1, x2, y2 = map(int, pred['bbox'])
        label = f"{pred['label']} {pred['confidence']:.2f}"

        # 绘制矩形框
        cv2.rectangle(img, (x1, y1), (x2, y2), color=(0, 255, 0), thickness=2)

        # 绘制标签
        cv2.putText(img, label, (x1, y1 - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 2)

    # 保存绘制好预测结果的图像
    cv2.imwrite(output_path, cv2.cvtColor(img, cv2.COLOR_RGB2BGR))

def extract_number_from_filename(filename):
    # 假设文件名格式为 'prefix_0000.extension'
    # 提取文件名中四位前补0的数字部分
    return int(filename.split('_')[-1][:4])

# 批量预测函数
def predict_folder(input_folder, output_folder, batch_size=16):
    # 获取输入文件夹中的所有图像文件路径
    image_paths = glob.glob(os.path.join(input_folder, '*.*'))
    print(image_paths)
    # Path(output_folder).mkdir(parents=True, exist_ok=True)  # 创建输出文件夹

    positive_indices = []

    for i in tqdm(range(0, len(image_paths), batch_size), desc="Processing batches"):
        batch_paths = image_paths[i:i + batch_size]
        imgs = [cv2.imread(img_path) for img_path in batch_paths]

        # 将 BGR 图像（OpenCV 读取的格式）转换为 RGB
        imgs_rgb = [cv2.cvtColor(img, cv2.COLOR_BGR2RGB) for img in imgs]

        # 进行预测
        results = model(imgs_rgb)  # 你可以调整 size 参数

        for img_path, img, result in zip(batch_paths, imgs, results):
            preds = []
            has_box = False
            for box in result.boxes:
                confidence = box.conf.item()  # 获取置信度
                if confidence <= 0.40:
                    continue
                has_box = True
                bbox = box.xyxy[0].tolist()  # 获取 bbox 的 [x1, y1, x2, y2] 坐标
                label = model.names[int(box.cls)]  # 获取标签名称
                pred = {
                    'bbox': bbox,
                    'label': label,
                    'confidence': confidence
                }
                preds.append(pred)

            if has_box:
                # output_img_path = os.path.join(output_folder, os.path.basename(img_path))
                # plot_predictions(img, preds, output_img_path)  # 保存带框的图片
                positive_indices.append(extract_number_from_filename(os.path.basename(img_path)))

    positive_indices.sort()
    return positive_indices

predict_output_folder = os.path.join(data_abs_path, 'predict_output')
positive_indices = predict_folder(adaptive_thres_dir, predict_output_folder, batch_size=16)

print(positive_indices)
ls, min_index, max_index = find_valid_indices(positive_indices)

lower_teeth_path = os.path.join(intraoral_scan_path, 'lower_teeth.ply')
upper_teeth_path = os.path.join(intraoral_scan_path, 'upper_teeth.ply')
upper_jaw_path = os.path.join(intraoral_scan_path, 'UpperJaw.stl')
lower_jaw_path = os.path.join(intraoral_scan_path, 'LowerJaw.stl')
print(lower_teeth_path)
print(upper_teeth_path)
registration = cbct_jaw_registration.CBCTJawRegistration(upper_teeth_path, lower_teeth_path, tiff_foler, max_index, min_index, 0.3, 'temp/')
registration.Registration()
upper_a, upper_b = registration.GetUpperTransformation()
lower_a, lower_b = registration.GetLowerTransformation()
print(upper_a)
print(lower_b)

affinex3f_upper = mm.AffineXf3f()
upper_a_matrix = mm.Matrix3f()
upper_a_matrix.x = mm.Vector3f(*upper_a[0])
upper_a_matrix.y = mm.Vector3f(*upper_a[1])
upper_a_matrix.z = mm.Vector3f(*upper_a[2])

affinex3f_upper.A = upper_a_matrix
affinex3f_upper.b = mm.Vector3f(*upper_b)


affinex3f_lower = mm.AffineXf3f()
lower_a_matrix = mm.Matrix3f()
lower_a_matrix.x = mm.Vector3f(*lower_a[0])
lower_a_matrix.y = mm.Vector3f(*lower_a[1])
lower_a_matrix.z = mm.Vector3f(*lower_a[2])

affinex3f_lower.A = lower_a_matrix
affinex3f_lower.b = mm.Vector3f(*lower_b)


upper_jaw = mm.loadMesh(upper_jaw_path)
lower_jaw = mm.loadMesh(lower_jaw_path)

upper_jaw.transform(affinex3f_upper)
lower_jaw.transform(affinex3f_lower)


mm.saveMesh(upper_jaw, os.path.join(cbct_path, 'upper_jaw_transformed.ply'))
mm.saveMesh(lower_jaw, os.path.join(cbct_path, 'lower_jaw_transformed.ply'))
