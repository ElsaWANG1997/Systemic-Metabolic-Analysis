import os
import numpy as np
import nibabel as nib
import pandas as pd
import pydicom
from pathlib import Path
from skimage.filters.rank import entropy as rank_entropy
from skimage.morphology import ball
from skimage.segmentation import slic
from skimage.measure import regionprops
import warnings
warnings.filterwarnings('ignore')


def load_mapping_and_info():
    """
    加载患者映射表和临床信息
    """
    # 加载映射表
    mapping_path = r" "
    df_mapping = pd.read_excel(mapping_path)
    
    # 创建映射字典: 024 -> S3
    pat_to_sn = {}
    for idx, row in df_mapping.iterrows():
        pat_id_raw = str(row.iloc[1]).strip()  # 第二列：024
        # 格式化为3位有效数字
        try:
            pat_id = f"{int(pat_id_raw):03d}"
        except:
            pat_id = pat_id_raw
        sn_id = str(row.iloc[2]).strip()   # 第三列：S3
        pat_to_sn[pat_id] = sn_id
    
    # 加载临床信息
    info_path = r" "
    df_info = pd.read_excel(info_path)
    
    # 创建信息字典: 024 -> {sex, BMI}
    pat_info = {}
    for idx, row in df_info.iterrows():
        pat_id_raw = str(row.iloc[1]).strip()  # 第二列：患者名
        # 格式化为3位有效数字
        try:
            pat_id = f"{int(pat_id_raw):03d}"
        except:
            pat_id = pat_id_raw
        sex = str(row.iloc[2]).strip()      # 第三列：性别
        bmi = float(row.iloc[6])            # 第7列：BMI
        pat_info[pat_id] = {'sex': sex, 'BMI': bmi}
    
    return pat_to_sn, pat_info


def calculate_SUL(pet_data, dicom_pet, sex, bmi):
    """
    计算PET SUL
    
    参数:
        pet_data: PET图像数据 (numpy array)
        dicom_pet: PET的DICOM头文件
        sex: 性别 ('M' or 'F')
        bmi: BMI值
    返回:
        PET_SUL: SUL标准化后的PET图像
    """
    # 提取时间信息
    acquisition_time = dicom_pet.AcquisitionTime
    radio_info = dicom_pet.RadiopharmaceuticalInformationSequence[0]
    start_time = radio_info.RadiopharmaceuticalStartTime
    half_life = radio_info.RadionuclideHalfLife
    
    # 计算时间（秒）
    def time_to_seconds(time_str):
        if len(time_str) >= 6:
            hours = int(time_str[0:2])
            minutes = int(time_str[2:4])
            seconds = int(time_str[4:6])
            return hours * 3600 + minutes * 60 + seconds
        else:
            raise ValueError(f"Invalid time format: {time_str}")
    
    acq_time_sec = time_to_seconds(acquisition_time)
    start_time_sec = time_to_seconds(start_time)
    decay_time = acq_time_sec - start_time_sec
    
    # 计算衰减剂量
    total_dose = radio_info.RadionuclideTotalDose
    decay_dose = total_dose * np.power(2, -decay_time / half_life) / 1000.0
    
    # 计算瘦体重 (LBM)
    patient_weight = float(dicom_pet.PatientWeight)
    
    if sex == 'F':
        LBM = (9.27 * 1000 * patient_weight) / (8.78 * 1000 + 244 * bmi)
    elif sex == 'M':
        LBM = (9.27 * 1000 * patient_weight) / (6.68 * 1000 + 216 * bmi)
    else:
        raise ValueError(f"Unknown sex: {sex}")
    
    # 计算SUL
    rescale_slope = float(dicom_pet.RescaleSlope)
    rescale_intercept = float(dicom_pet.RescaleIntercept)
    
    PET_SUL = (pet_data * rescale_slope + rescale_intercept) / (decay_dose / LBM)
    
    return PET_SUL


def lung_window_normalize(hu, WL=-600.0, WW=1500.0, out_range=(0.0, 1.0)):
    """
    肺窗归一化
    
    参数:
        hu: HU值数组
        WL: 窗位
        WW: 窗宽
        out_range: 输出范围
    返回:
        归一化后的数组
    """
    minHU = WL - WW / 2.0
    maxHU = WL + WW / 2.0
    hu = np.asarray(hu, dtype=np.float32)
    
    # 裁剪
    hu = np.clip(hu, minHU, maxHU)
    
    # [0,1] 归一化
    x = (hu - minHU) / (maxHU - minHU)
    
    if out_range == (0.0, 1.0):
        return x
    elif out_range == (-1.0, 1.0):
        return 2.0 * x - 1.0
    else:
        lo, hi = out_range
        return lo + x * (hi - lo)


def local_entropy_map(label_img, img_data, radius=2, is_quantized=False):
    """
    计算局部熵图
    
    参数:
        label_img: 二值肿瘤mask
        img_data: 输入图像数据
        radius: 半径 (2 -> 5x5x5)
        is_quantized: 是否已经量化
    返回:
        熵图
    """
    # 如果未量化，先量化到uint8
    if not is_quantized:
        # 归一化到[0, 255]
        img_min = np.nanmin(img_data[label_img > 0])
        img_max = np.nanmax(img_data[label_img > 0])
        if img_max > img_min:
            img_q = ((img_data - img_min) / (img_max - img_min) * 255).astype(np.uint8)
        else:
            img_q = np.zeros_like(img_data, dtype=np.uint8)
    else:
        img_q = img_data.astype(np.uint8)
    
    # 创建3D结构元素
    selem_3d = ball(radius)
    selem_2d = selem_3d[:, :, selem_3d.shape[2] // 2]  # 取中间切片作为2D结构元素
    
    # 逐层计算熵
    ent = np.zeros_like(img_q, dtype=np.float32)
    for z in range(img_q.shape[2]):  # 注意：这里假设z是第三个维度
        slice_data = img_q[:, :, z]
        ent[:, :, z] = rank_entropy(slice_data, selem=selem_2d)
    
    # Mask外设为NaN
    ent[label_img == 0] = np.nan
    
    return ent


def align_volumes(volume, mask):
    """
    对齐volume和mask，使用尾部对齐策略
    
    参数:
        volume: 图像volume
        mask: ROI mask
    返回:
        对齐后的volume
    """
    vol_shape = volume.shape
    mask_shape = mask.shape
    
    if vol_shape == mask_shape:
        return volume
    
    # 检查x, y维度
    if vol_shape[0] != mask_shape[0] or vol_shape[1] != mask_shape[1]:
        raise ValueError(f"X或Y维度不匹配: volume {vol_shape}, mask {mask_shape}")
    
    # Z轴尾部对齐
    vol_z = vol_shape[2]
    mask_z = mask_shape[2]
    
    if vol_z < mask_z:
        raise ValueError(f"Volume的z轴({vol_z})小于mask的z轴({mask_z})")
    
    # 使用尾部切片
    start_z = vol_z - mask_z
    aligned_volume = volume[:, :, start_z:vol_z]
    
    print(f"  对齐: 原始volume shape {vol_shape}, 使用z切片[{start_z+1}:{vol_z}], 对齐后shape {aligned_volume.shape}")
    
    return aligned_volume


def extract_tumor_features(patient_folder, pat_to_sn, pat_info):
    """
    提取单个患者的肿瘤特征
    
    参数:
        patient_folder: 患者文件夹路径
        pat_to_sn: 患者ID到SN的映射
        pat_info: 患者临床信息
    返回:
        features_dict: 包含4个特征的字典
    """
    patient_name = patient_folder.name
    print(f"\n{'='*80}")
    print(f"处理患者: {patient_name}")
    print(f"{'='*80}")
    
    # 加载ROI mask
    roi_path = patient_folder / "PET_ROI.nii.gz"
    if not roi_path.exists():
        print(f"  错误: 未找到 {roi_path}")
        return None
    
    roi_img = nib.load(roi_path)
    roi_mask = roi_img.get_fdata().astype(bool)
    print(f"  ROI mask shape: {roi_mask.shape}, 体素数: {np.sum(roi_mask)}")
    
    # 加载PET数据
    pet_path = patient_folder / "PET.nii"
    if not pet_path.exists():
        print(f"  错误: 未找到 {pet_path}")
        return None
    
    pet_img = nib.load(pet_path)
    pet_data = pet_img.get_fdata()
    print(f"  PET shape: {pet_data.shape}")
    
    # 加载CT数据
    ct_path = patient_folder / "resample_CT.nii"
    if not ct_path.exists():
        print(f"  错误: 未找到 {ct_path}")
        return None
    
    ct_img = nib.load(ct_path)
    ct_data = ct_img.get_fdata()
    print(f"  CT shape: {ct_data.shape}")
    
    # 对齐数据
    pet_data = align_volumes(pet_data, roi_mask)
    ct_data = align_volumes(ct_data, roi_mask)
    
    # 获取映射的SN ID
    if patient_name not in pat_to_sn:
        print(f"  警告: 患者 {patient_name} 在映射表中未找到")
        return None
    
    sn_id = pat_to_sn[patient_name]
    print(f"  映射到SN: {sn_id}")
    
    # 获取临床信息
    if patient_name not in pat_info:
        print(f"  警告: 患者 {patient_name} 的临床信息未找到")
        return None
    
    sex = pat_info[patient_name]['sex']
    bmi = pat_info[patient_name]['BMI']
    print(f"  性别: {sex}, BMI: {bmi}")
    
    # 加载PET DICOM头文件
    dicom_pet_dir = Path(r" ") / sn_id / "AC-PT"
    dicom_pet_path = dicom_pet_dir / "00000001.dcm"
    
    if not dicom_pet_path.exists():
        print(f"  错误: 未找到PET DICOM {dicom_pet_path}")
        return None
    
    dicom_pet = pydicom.dcmread(dicom_pet_path)
    
    # 加载CT DICOM头文件
    dicom_ct_dir = Path(r" ") / sn_id / "AC-CT"
    dicom_ct_path = dicom_ct_dir / "00000001.dcm"
    
    if not dicom_ct_path.exists():
        print(f"  错误: 未找到CT DICOM {dicom_ct_path}")
        return None
    
    dicom_ct = pydicom.dcmread(dicom_ct_path)
    
    # 1. 计算PET SUL
    print(f"\n  计算PET SUL...")
    pet_sul = calculate_SUL(pet_data, dicom_pet, sex, bmi)
    pet_sul_tumor = pet_sul[roi_mask].flatten()
    print(f"    SUL范围: [{pet_sul_tumor.min():.2f}, {pet_sul_tumor.max():.2f}], mean={pet_sul_tumor.mean():.2f}")
    
    # 2. 计算PET熵
    print(f"  计算PET局部熵 (5x5x5)...")
    pet_entropy = local_entropy_map(roi_mask, pet_sul, radius=2)
    pet_entropy_tumor = pet_entropy[roi_mask].flatten()
    pet_entropy_tumor = pet_entropy_tumor[~np.isnan(pet_entropy_tumor)]
    print(f"    熵范围: [{pet_entropy_tumor.min():.2f}, {pet_entropy_tumor.max():.2f}], mean={pet_entropy_tumor.mean():.2f}")
    
    # 3. 计算CT HU
    print(f"  计算CT HU...")
    rescale_slope_ct = float(dicom_ct.RescaleSlope)
    rescale_intercept_ct = float(dicom_ct.RescaleIntercept)
    ct_hu = ct_data * rescale_slope_ct + rescale_intercept_ct
    
    # 肺窗归一化
    ct_normalized = lung_window_normalize(ct_hu, WL=-600.0, WW=1500.0, out_range=(0.0, 1.0))
    ct_normalized_tumor = ct_normalized[roi_mask].flatten()
    print(f"    归一化CT范围: [{ct_normalized_tumor.min():.2f}, {ct_normalized_tumor.max():.2f}], mean={ct_normalized_tumor.mean():.2f}")
    
    # 4. 计算CT熵
    print(f"  计算CT局部熵 (5x5x5)...")
    ct_entropy = local_entropy_map(roi_mask, ct_normalized, radius=2)
    ct_entropy_tumor = ct_entropy[roi_mask].flatten()
    ct_entropy_tumor = ct_entropy_tumor[~np.isnan(ct_entropy_tumor)]
    print(f"    熵范围: [{ct_entropy_tumor.min():.2f}, {ct_entropy_tumor.max():.2f}], mean={ct_entropy_tumor.mean():.2f}")
    
    # 确保所有特征长度一致
    n_voxels = np.sum(roi_mask)
    features_dict = {
        'patient': patient_name,
        'PET_SUL': pet_sul,
        'PET_Entropy': pet_entropy,
        'CT_HU': ct_normalized,
        'CT_Entropy': ct_entropy,
        'ROI_mask': roi_mask,
        'n_voxels': n_voxels
    }
    
    print(f"  特征提取完成，共 {n_voxels} 个肿瘤体素")
    
    return features_dict


def perform_slic_segmentation(features_dict, n_segments=100, compactness=10.0):
    """
    对提取的特征进行SLIC超像素分割
    
    参数:
        features_dict: 特征字典
        n_segments: 超像素数量
        compactness: 紧凑度参数
    返回:
        slic_labels: SLIC分割标签
        superpixel_features: 超像素级别的特征
    """
    patient_name = features_dict['patient']
    print(f"\n  执行SLIC超像素分割...")
    print(f"    参数: n_segments={n_segments}, compactness={compactness}")
    
    # 构建4D特征volume
    pet_sul = features_dict['PET_SUL']
    pet_entropy = features_dict['PET_Entropy']
    ct_hu = features_dict['CT_HU']
    ct_entropy = features_dict['CT_Entropy']
    roi_mask = features_dict['ROI_mask']
    
    # 归一化每个特征到[0, 1]
    def normalize_feature(feat, mask):
        feat_masked = feat[mask]
        feat_min = np.nanmin(feat_masked)
        feat_max = np.nanmax(feat_masked)
        if feat_max > feat_min:
            feat_norm = (feat - feat_min) / (feat_max - feat_min)
        else:
            feat_norm = np.zeros_like(feat)
        feat_norm[~mask] = 0
        return feat_norm
    
    pet_sul_norm = normalize_feature(pet_sul, roi_mask)
    pet_entropy_norm = normalize_feature(pet_entropy, roi_mask)
    ct_hu_norm = ct_hu  # 已经归一化
    ct_entropy_norm = normalize_feature(ct_entropy, roi_mask)
    
    # 堆叠成4D volume (x, y, z, 4)
    feature_volume = np.stack([
        pet_sul_norm,
        pet_entropy_norm,
        ct_hu_norm,
        ct_entropy_norm
    ], axis=-1)
    
    # 只在ROI区域进行SLIC
    # 创建一个包含ROI的边界框
    coords = np.where(roi_mask)
    x_min, x_max = coords[0].min(), coords[0].max()
    y_min, y_max = coords[1].min(), coords[1].max()
    z_min, z_max = coords[2].min(), coords[2].max()
    
    # 裁剪到ROI边界框
    feature_volume_cropped = feature_volume[x_min:x_max+1, y_min:y_max+1, z_min:z_max+1, :]
    roi_mask_cropped = roi_mask[x_min:x_max+1, y_min:y_max+1, z_min:z_max+1]
    
    print(f"    裁剪后volume shape: {feature_volume_cropped.shape}")
    
    # 执行SLIC (3D版本)
    # 注意：skimage的slic支持多通道3D图像
    slic_labels_cropped = slic(
        feature_volume_cropped,
        n_segments=n_segments,
        compactness=compactness,
        mask=roi_mask_cropped,
        start_label=1
    )
    
    # 恢复到原始尺寸
    slic_labels = np.zeros_like(roi_mask, dtype=int)
    slic_labels[x_min:x_max+1, y_min:y_max+1, z_min:z_max+1] = slic_labels_cropped
    
    # 统计超像素
    unique_labels = np.unique(slic_labels[slic_labels > 0])
    n_superpixels = len(unique_labels)
    
    print(f"    生成 {n_superpixels} 个超像素")
    
    # 计算每个超像素的平均特征
    superpixel_features = []
    
    for label in unique_labels:
        sp_mask = (slic_labels == label)
        sp_size = np.sum(sp_mask)
        
        sp_feat = {
            'patient': patient_name,
            'superpixel_id': label,
            'n_voxels': sp_size,
            'PET_SUL_mean': np.mean(pet_sul[sp_mask]),
            'PET_SUL_std': np.std(pet_sul[sp_mask]),
            'PET_Entropy_mean': np.nanmean(pet_entropy[sp_mask]),
            'PET_Entropy_std': np.nanstd(pet_entropy[sp_mask]),
            'CT_HU_mean': np.mean(ct_hu[sp_mask]),
            'CT_HU_std': np.std(ct_hu[sp_mask]),
            'CT_Entropy_mean': np.nanmean(ct_entropy[sp_mask]),
            'CT_Entropy_std': np.nanstd(ct_entropy[sp_mask])
        }
        
        superpixel_features.append(sp_feat)
    
    return slic_labels, superpixel_features


def save_results(all_superpixel_features, output_dir):
    """
    保存结果
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 保存超像素特征
    df_features = pd.DataFrame(all_superpixel_features)
    feature_path = output_dir / "SLIC_Superpixel_Features.xlsx"
    df_features.to_excel(feature_path, index=False)
    print(f"\n超像素特征已保存至: {feature_path}")
    
    return feature_path


if __name__ == "__main__":
    # 设置路径
    tumor_seg_dir = Path(r" ")
    output_dir = Path(r" ")
    
    print("=" * 80)
    print("SLIC 超像素分割 - 肿瘤特征提取")
    print("=" * 80)
    print(f"输入目录: {tumor_seg_dir}")
    print(f"输出目录: {output_dir}")
    print("=" * 80)
    
    # 加载映射和信息
    print("\n加载患者映射和临床信息...")
    pat_to_sn, pat_info = load_mapping_and_info()
    print(f"  映射表加载: {len(pat_to_sn)} 个患者")
    print(f"  临床信息加载: {len(pat_info)} 个患者")
    
    # 获取所有患者文件夹
    patient_folders = [f for f in tumor_seg_dir.iterdir() if f.is_dir()]
    print(f"\n找到 {len(patient_folders)} 个患者文件夹")
    
    # 处理每个患者
    all_superpixel_features = []
    
    for patient_folder in sorted(patient_folders):
        try:
            # 提取特征
            features_dict = extract_tumor_features(patient_folder, pat_to_sn, pat_info)
            
            if features_dict is None:
                continue
            
            # SLIC分割
            slic_labels, superpixel_features = perform_slic_segmentation(
                features_dict,
                n_segments=100,
                compactness=10.0
            )
            
            all_superpixel_features.extend(superpixel_features)
            
            # 保存SLIC超像素图像到指定目录
            superpixels_dir = Path(r" ")
            superpixels_dir.mkdir(parents=True, exist_ok=True)
            
            slic_nii = nib.Nifti1Image(slic_labels.astype(np.int16), affine=nib.load(patient_folder / "PET_ROI.nii.gz").affine)
            slic_path = superpixels_dir / f"{patient_folder.name}_tumor_superpixels.nii"
            nib.save(slic_nii, slic_path)
            print(f"  超像素图像已保存: {slic_path}")
            
        except Exception as e:
            print(f"\n错误: 处理患者 {patient_folder.name} 时出错")
            print(f"  {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # 保存结果
    print("\n" + "=" * 80)
    print("保存结果...")
    print("=" * 80)
    
    feature_path = save_results(all_superpixel_features, output_dir)
    
    print("\n" + "=" * 80)
    print("处理完成!")
    print("=" * 80)
    print(f"\n生成的文件:")
    print(f"  1. {feature_path} - 超像素特征表")
    superpixels_dir = Path(r" ")
    print(f"  2. {superpixels_dir / '[患者名]_tumor_superpixels.nii'} - 超像素分割图像")
    print(f"\n共处理 {len(all_superpixel_features)} 个超像素")
    print("=" * 80)

