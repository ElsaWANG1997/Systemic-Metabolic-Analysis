% Systemic metabolic network analysis: Comparing network connectivity patterns under different smoking states
% Validation using Bootstrap and Combat correction of data


%% 读取Combat校正后的数据
disp('=== 读取Combat校正和GLM调整后的数据 ===');

% 读取SUL校正数据
SUL_harmonized = readtable(' ');
disp(['SUL harmonized data: ' num2str(size(SUL_harmonized,1)) ' patients × ' num2str(size(SUL_harmonized,2)) ' columns']);

% 提取基本信息
patient_names = SUL_harmonized{:, 1};  % 患者姓名
center_info = SUL_harmonized{:, 2};    % 中心信息 (1=RJ, 2=HN, 3=QFS)
disease_info = SUL_harmonized{:, 9};   % 疾病信息 (0=HC, 1=Benign, 2=ADC, 3=SCC, 4=SCLC)
smoking_info = SUL_harmonized{:, 8};   % 吸烟史信息 (0=不吸烟, 1=吸烟)

% 提取ROI数据（第10列开始，因为第8、9列是smoking和disease）
roi_cols = 10:size(SUL_harmonized, 2);
SUL_data_all = table2array(SUL_harmonized(:, roi_cols)); %原始的SUL值

disp(['ROI数量: ' num2str(length(roi_cols))]);

% 从Excel表格的列名中提取ROI名称
roi_column_names = SUL_harmonized.Properties.VariableNames(roi_cols);
ROInames = roi_column_names;

%% 计算代谢率（相对于肝脏的SUL值）
disp('=== 计算代谢率（相对于肝脏） ===');

% 肝脏是第28个ROI（在原始38个ROI中）
liver_roi_index = 28;
fprintf('肝脏ROI索引: %d\n', liver_roi_index);
fprintf('肝脏ROI名称: %s\n', ROInames{liver_roi_index});

% 提取肝脏的SUL值
liver_SUL = SUL_data_all(:, liver_roi_index);

% 检查肝脏SUL值
fprintf('肝脏SUL统计:\n');
fprintf('  均值: %.4f\n', mean(liver_SUL, 'omitnan'));
fprintf('  标准差: %.4f\n', std(liver_SUL, 'omitnan'));
fprintf('  范围: [%.4f, %.4f]\n', min(liver_SUL), max(liver_SUL));
fprintf('  零值数量: %d\n', sum(liver_SUL == 0));
fprintf('  NaN数量: %d\n', sum(isnan(liver_SUL)));

% 处理肝脏SUL为零或NaN的情况
liver_SUL_safe = liver_SUL;
liver_SUL_safe(liver_SUL_safe == 0 | isnan(liver_SUL_safe)) = NaN;

if sum(isnan(liver_SUL_safe)) > 0
    fprintf('警告: %d个患者的肝脏SUL值为零或NaN，将在代谢率计算中排除\n', sum(isnan(liver_SUL_safe)));
end

% 计算代谢率：每个ROI的SUL / 肝脏SUL
fprintf('计算代谢率...\n');
metabolic_rates = zeros(size(SUL_data_all));

for roi = 1:size(SUL_data_all, 2)
    metabolic_rates(:, roi) = SUL_data_all(:, roi) ./ liver_SUL_safe;
end

% 检查代谢率计算结果
fprintf('代谢率统计:\n');
fprintf('  总体均值: %.4f\n', mean(metabolic_rates(:), 'omitnan'));
fprintf('  总体标准差: %.4f\n', std(metabolic_rates(:), 'omitnan'));
fprintf('  NaN比例: %.2f%%\n', sum(isnan(metabolic_rates(:)))/numel(metabolic_rates)*100);

% 显示各ROI的代谢率统计
fprintf('\n各ROI代谢率统计（前10个）:\n');
for roi = 1:min(10, size(metabolic_rates, 2))
    roi_rates = metabolic_rates(:, roi);
    fprintf('  %s: 均值=%.3f, 标准差=%.3f\n', ROInames{roi}, mean(roi_rates, 'omitnan'), std(roi_rates, 'omitnan'));
end

% 更新数据：使用代谢率替代原始SUL值
SUL_data_all = metabolic_rates;
fprintf('\n✓ 已将SUL数据替换为相对于肝脏的代谢率数据\n');

%% 排除Aorta（第21号器官）
disp('=== 排除Aorta器官 ===');

% Aorta是第21个ROI，在全身器官分析中需要排除
aorta_roi_index = 21;
fprintf('排除器官: 第%d号ROI - %s\n', aorta_roi_index, ROInames{aorta_roi_index});

% 创建排除Aorta的索引
roi_indices_to_keep = 1:size(SUL_data_all, 2);
roi_indices_to_keep(aorta_roi_index) = [];  % 移除第21号ROI

% 更新数据和ROI名称
SUL_data_all = SUL_data_all(:, roi_indices_to_keep);
ROInames = ROInames(roi_indices_to_keep);

fprintf('排除Aorta后的ROI数量: %d\n', length(ROInames));
fprintf('✓ 已从分析中排除Aorta器官\n');

%% 按吸烟史和疾病分组数据
disp('=== 按吸烟史和疾病分组数据 ===');

% 检查数据类型并转换
fprintf('检查数据类型...\n');
fprintf('disease_info类型: %s\n', class(disease_info));
fprintf('smoking_info类型: %s\n', class(smoking_info));

% 转换为数值类型（如果是cell或其他类型）
if iscell(disease_info)
    % 尝试不同的转换方法
    try
        disease_info = cell2mat(disease_info);
    catch
        % 如果cell2mat失败，尝试逐个转换
        fprintf('使用逐个转换方法处理disease_info...\n');
        disease_numeric = zeros(size(disease_info));
        for i = 1:length(disease_info)
            if isnumeric(disease_info{i})
                disease_numeric(i) = disease_info{i};
            elseif ischar(disease_info{i}) || isstring(disease_info{i})
                disease_numeric(i) = str2double(disease_info{i});
            else
                disease_numeric(i) = NaN;
            end
        end
        disease_info = disease_numeric;
    end
end

if iscell(smoking_info)
    % 尝试不同的转换方法
    try
        smoking_info = cell2mat(smoking_info);
    catch
        % 如果cell2mat失败，尝试逐个转换
        fprintf('使用逐个转换方法处理smoking_info...\n');
        smoking_numeric = zeros(size(smoking_info));
        for i = 1:length(smoking_info)
            if isnumeric(smoking_info{i})
                smoking_numeric(i) = smoking_info{i};
            elseif ischar(smoking_info{i}) || isstring(smoking_info{i})
                smoking_numeric(i) = str2double(smoking_info{i});
            else
                smoking_numeric(i) = NaN;
            end
        end
        smoking_info = smoking_numeric;
    end
end

% 如果是字符串，尝试转换为数值
if isstring(disease_info) || ischar(disease_info)
    disease_info = str2double(disease_info);
end
if isstring(smoking_info) || ischar(smoking_info)
    smoking_info = str2double(smoking_info);
end

fprintf('转换后数据类型:\n');
fprintf('disease_info类型: %s, 唯一值: %s\n', class(disease_info), mat2str(unique(disease_info)));
fprintf('smoking_info类型: %s, 唯一值: %s\n', class(smoking_info), mat2str(unique(smoking_info)));

% 定义分组索引
% HC组（疾病状态=0）
hc_never_smoking_idx = (disease_info == 0) & (smoking_info == 0);  % HC不吸烟
hc_smoking_idx = (disease_info == 0) & (smoking_info == 1);        % HC吸烟

% LC内部组（疾病状态1-4，中心2）
lc_internal_never_smoking_idx = (disease_info >= 1) & (disease_info <= 4) & (center_info == 2) & (smoking_info == 0);  % LC内部不吸烟
lc_internal_smoking_idx = (disease_info >= 1) & (disease_info <= 4) & (center_info == 2) & (smoking_info == 1);        % LC内部吸烟

% LC外部组（疾病状态1-4，中心3）
lc_external_never_smoking_idx = (disease_info >= 1) & (disease_info <= 4) & (center_info == 3) & (smoking_info == 0);  % LC外部不吸烟
lc_external_smoking_idx = (disease_info >= 1) & (disease_info <= 4) & (center_info == 3) & (smoking_info == 1);        % LC外部吸烟

% 提取分组数据
HC_never_smoking_data = SUL_data_all(hc_never_smoking_idx, :);
HC_smoking_data = SUL_data_all(hc_smoking_idx, :);
LC_internal_never_smoking_data = SUL_data_all(lc_internal_never_smoking_idx, :);
LC_internal_smoking_data = SUL_data_all(lc_internal_smoking_idx, :);
LC_external_never_smoking_data = SUL_data_all(lc_external_never_smoking_idx, :);
LC_external_smoking_data = SUL_data_all(lc_external_smoking_idx, :);

% 获取患者名称
HC_never_smoking_names = patient_names(hc_never_smoking_idx);
HC_smoking_names = patient_names(hc_smoking_idx);
LC_internal_never_smoking_names = patient_names(lc_internal_never_smoking_idx);
LC_internal_smoking_names = patient_names(lc_internal_smoking_idx);
LC_external_never_smoking_names = patient_names(lc_external_never_smoking_idx);
LC_external_smoking_names = patient_names(lc_external_smoking_idx);

% 显示分组统计
fprintf('\n=== 分组统计 ===\n');
fprintf('HC组:\n');
fprintf('  不吸烟: %d patients\n', size(HC_never_smoking_data, 1));
fprintf('  吸烟: %d patients\n', size(HC_smoking_data, 1));

fprintf('LC内部组 (Center=2):\n');
fprintf('  不吸烟: %d patients\n', size(LC_internal_never_smoking_data, 1));
fprintf('  吸烟: %d patients\n', size(LC_internal_smoking_data, 1));

fprintf('LC外部组 (Center=3):\n');
fprintf('  不吸烟: %d patients\n', size(LC_external_never_smoking_data, 1));
fprintf('  吸烟: %d patients\n', size(LC_external_smoking_data, 1));

% 检查各组样本数量
groups_to_check = {
    {'HC不吸烟', HC_never_smoking_data};
    {'HC吸烟', HC_smoking_data};
    {'LC内部不吸烟', LC_internal_never_smoking_data};
    {'LC内部吸烟', LC_internal_smoking_data};
    {'LC外部不吸烟', LC_external_never_smoking_data};
    {'LC外部吸烟', LC_external_smoking_data}
};

fprintf('\n=== 样本充足性检查 ===\n');
for i = 1:length(groups_to_check)
    group_name = groups_to_check{i}{1};
    group_data = groups_to_check{i}{2};
    n_samples = size(group_data, 1);
    
    if n_samples >= 5
        fprintf('✓ %s: %d patients (充足)\n', group_name, n_samples);
    else
        fprintf('⚠ %s: %d patients (样本不足，建议≥5)\n', group_name, n_samples);
    end
end

%% Bootstrap参数设置
disp('=== Bootstrap参数设置 ===');
n_boot = 1000; % Bootstrap重采样次数
alpha = 0.05; % 显著性水平
consistency_thresh = 0.95; % 一致性阈值（95%的Bootstrap样本中显著）
r_thresh = 0.5; % 相关性强度阈值

fprintf('Bootstrap参数:\n');
fprintf('  n_boot = %d\n', n_boot);
fprintf('  alpha = %.3f\n', alpha);
fprintf('  consistency_thresh = %.3f\n', consistency_thresh);
fprintf('  r_thresh = %.3f (相关性强度阈值)\n', r_thresh);

%% 计算各组协方差网络
disp('===============================================');
disp('=== 计算各组协方差网络 ===');
disp('===============================================');

% 为每个组分别建立协方差网络矩阵
% 每个网络都通过Bootstrap重采样和FDR校正来确保稳定性和统计显著性
networks_results = struct();

% 定义要分析的组
groups_to_analyze = {
    {'HC_never_smoking', HC_never_smoking_data, 'HC不吸烟组'};
    {'HC_smoking', HC_smoking_data, 'HC吸烟组'};
    {'LC_internal_never_smoking', LC_internal_never_smoking_data, 'LC内部不吸烟组'};
    {'LC_internal_smoking', LC_internal_smoking_data, 'LC内部吸烟组'};
    {'LC_external_never_smoking', LC_external_never_smoking_data, 'LC外部不吸烟组'};
    {'LC_external_smoking', LC_external_smoking_data, 'LC外部吸烟组'}
};

% 为每组建立独立的协方差网络
for i = 1:length(groups_to_analyze)
    group_name = groups_to_analyze{i}{1};
    group_data = groups_to_analyze{i}{2};
    group_desc = groups_to_analyze{i}{3};
    
    if size(group_data, 1) >= 5
        fprintf('\n建立%s协方差网络...\n', group_desc);
        fprintf('  样本数: %d, ROI数: %d\n', size(group_data, 1), size(group_data, 2));
        
        % 使用Bootstrap方法建立稳定的协方差网络
        network_matrix = calculate_group_network_bootstrap(group_data, n_boot, alpha, consistency_thresh, r_thresh);
        networks_results.(group_name) = network_matrix;
        
        % 显示网络统计
        nonzero_connections = sum(abs(network_matrix(:)) > 0) / 2;
        max_strength = max(abs(network_matrix(:)));
        fprintf('  %s协方差网络: %d个显著连接, 最大相关强度: %.4f\n', group_desc, nonzero_connections, max_strength);
    else
        fprintf('\n跳过%s（样本不足: %d < 5）\n', group_desc, size(group_data, 1));
        networks_results.(group_name) = [];
    end
end

%% 网络统计分析
disp('===============================================');
disp('=== 网络统计分析 ===');
disp('===============================================');

% 输出各组网络的统计信息
fprintf('\n=== 各组网络统计 ===\n');

for i = 1:length(groups_to_analyze)
    group_name = groups_to_analyze{i}{1};
    group_desc = groups_to_analyze{i}{3};
    
    if ~isempty(networks_results.(group_name))
        network_matrix = networks_results.(group_name);
        
        % 计算网络统计指标
        nonzero_connections = sum(abs(network_matrix(:)) > 0) / 2;
        total_possible_connections = (size(network_matrix,1) * (size(network_matrix,1) - 1)) / 2;
        network_density = nonzero_connections / total_possible_connections;
        max_strength = max(abs(network_matrix(:)));
        mean_strength = mean(abs(network_matrix(network_matrix ~= 0)));
        
        fprintf('%s:\n', group_desc);
        fprintf('  显著连接数: %d / %d (%.2f%%)\n', nonzero_connections, total_possible_connections, network_density*100);
        fprintf('  最大连接强度: %.4f\n', max_strength);
        fprintf('  平均连接强度: %.4f\n', mean_strength);
        fprintf('  网络密度: %.4f\n', network_density);
        fprintf('\n');
    else
        fprintf('%s: 未计算（样本不足）\n\n', group_desc);
    end
end

%% 保存结果
disp('===============================================');
disp('=== 保存结果 ===');
disp('===============================================');

% 创建保存结构
results_struct = struct();
results_struct.ROInames = ROInames;
results_struct.patient_info = table(patient_names, center_info, disease_info, smoking_info, ...
    'VariableNames', {'PatientName', 'Center', 'Disease', 'Smoking'});

% 保存网络数据
field_names = fieldnames(networks_results);
for i = 1:length(field_names)
    field_name = field_names{i};
    if ~isempty(networks_results.(field_name))
        results_struct.(field_name) = networks_results.(field_name);
    end
end

% 保存MAT文件
save('Internal_External_LC_Networks_Analysis_Results.mat', '-struct', 'results_struct');
fprintf('结果已保存到: Internal_External_LC_Networks_Analysis_Results.mat\n');

% %% 网络差异分析和可视化
% disp('===============================================');
% disp('=== 网络差异分析和可视化 ===');
% disp('===============================================');
% 
% % 可视化参数设置
% vis_threshold = 0.3; % 可视化阈值，过滤掉很小的连接
% top_n_connections = 15; % 显示前N个最强差异连接
% 
% % 定义网络差异比较
% difference_analyses = {
%     % 1. Never smoking internal LC - Never smoking HC
%     {'LC_internal_never_smoking', 'HC_never_smoking', 'Never smoking internal LC - Never smoking HC', 'Never_smoking_internal_LC_vs_Never_smoking_HC'};
%     % 2. Smoking internal LC - Never smoking internal LC
%     {'LC_internal_smoking', 'LC_internal_never_smoking', 'Smoking internal LC - Never smoking internal LC', 'Smoking_internal_LC_vs_Never_smoking_internal_LC'};
%     % 3. Smoking internal LC - Smoking HC
%     {'LC_internal_smoking', 'HC_smoking', 'Smoking internal LC - Smoking HC', 'Smoking_internal_LC_vs_Smoking_HC'};
%     % 4. Never smoking external LC - Never smoking HC
%     {'LC_external_never_smoking', 'HC_never_smoking', 'Never smoking external LC - Never smoking HC', 'Never_smoking_external_LC_vs_Never_smoking_HC'};
%     % 5. Smoking external LC - Never smoking external LC
%     {'LC_external_smoking', 'LC_external_never_smoking', 'Smoking external LC - Never smoking external LC', 'Smoking_external_LC_vs_Never_smoking_external_LC'};
%     % 6. Smoking external LC - Smoking HC
%     {'LC_external_smoking', 'HC_smoking', 'Smoking external LC - Smoking HC', 'Smoking_external_LC_vs_Smoking_HC'}
% };
% 
% % 存储差异网络结果
% difference_networks = struct();

for i = 1:length(difference_analyses)
    network1_name = difference_analyses{i}{1};  % 被减数网络
    network2_name = difference_analyses{i}{2};  % 减数网络
    comparison_title = difference_analyses{i}{3};
    save_name = difference_analyses{i}{4};
    
    fprintf('\n=== %s ===\n', comparison_title);
    
    % 检查两个网络是否都存在
    if isfield(networks_results, network1_name) && isfield(networks_results, network2_name) && ...
       ~isempty(networks_results.(network1_name)) && ~isempty(networks_results.(network2_name))
        
        network1 = networks_results.(network1_name);
        network2 = networks_results.(network2_name);
        
        % 计算网络差异
        diff_network = network1 - network2;
        difference_networks.(save_name) = diff_network;
        
        % 统计差异网络信息
        total_connections = sum(abs(diff_network(:)) > 0) / 2;
        max_positive_diff = max(diff_network(:));
        max_negative_diff = min(diff_network(:));
        mean_abs_diff = mean(abs(diff_network(:)), 'omitnan');
        
        fprintf('差异网络统计:\n');
        fprintf('  非零差异连接数: %d\n', total_connections);
        fprintf('  最大正差异: %.4f\n', max_positive_diff);
        fprintf('  最大负差异: %.4f\n', max_negative_diff);
        fprintf('  平均绝对差异: %.4f\n', mean_abs_diff);
        
        % 找出差异最大的前N个连接
        fprintf('\n提取前%d个最强差异连接...\n', top_n_connections);
        [top_connections, top_values] = find_top_network_differences(diff_network, ROInames, top_n_connections);
        
        % 显示前N个最强差异连接
        fprintf('前%d个最强差异连接:\n', min(top_n_connections, length(top_values)));
        for j = 1:min(top_n_connections, length(top_values))
            fprintf('  %d. %s <-> %s: %.4f\n', j, top_connections{j,1}, top_connections{j,2}, top_values(j));
        end
        
        % 创建只包含前N个最强连接的网络用于可视化
        top_diff_network = zeros(size(diff_network));
        for j = 1:min(top_n_connections, length(top_values))
            roi1_idx = strcmp(ROInames, top_connections{j,1});
            roi2_idx = strcmp(ROInames, top_connections{j,2});
            roi1_idx = find(roi1_idx);
            roi2_idx = find(roi2_idx);
            if ~isempty(roi1_idx) && ~isempty(roi2_idx)
                top_diff_network(roi1_idx, roi2_idx) = diff_network(roi1_idx, roi2_idx);
                top_diff_network(roi2_idx, roi1_idx) = diff_network(roi2_idx, roi1_idx);
    end
end

        % 可视化差异网络（前N个最强连接）
        fprintf('创建%s差异网络可视化...\n', comparison_title);
        try
            figure('Name', sprintf('%s (前%d个最强差异)', comparison_title, top_n_connections));
            A = top_diff_network;
            symmetricMatrix = A + A.' - diag(diag(A));
            circularGraph(symmetricMatrix, symmetricMatrix, 'Label', ROInames);
            title(sprintf('%s\n(前%d个最强差异连接)', comparison_title, top_n_connections));
        catch ME
            fprintf('  警告: %s差异网络可视化失败: %s\n', comparison_title, ME.message);
            % 创建简单的热力图
            figure('Name', sprintf('%s差异热力图', comparison_title));
            imagesc(top_diff_network);
            colorbar;
            title(sprintf('%s差异网络', comparison_title));
            axis equal tight;
        end
        
    else
        fprintf('跳过%s（缺少必要的网络数据）\n', comparison_title);
        if ~isfield(networks_results, network1_name) || isempty(networks_results.(network1_name))
            fprintf('  缺少网络: %s\n', network1_name);
        end
        if ~isfield(networks_results, network2_name) || isempty(networks_results.(network2_name))
            fprintf('  缺少网络: %s\n', network2_name);
        end
    end
end

% %% 原始网络可视化（可选）
% fprintf('\n=== 原始网络可视化 ===\n');
% 
% % 为每个有效网络创建热力图
% visualization_groups = {
%     {'HC_never_smoking', 'HC不吸烟组网络'};
%     {'HC_smoking', 'HC吸烟组网络'};
%     {'Disease_never_smoking', '疾病不吸烟组网络'};
%     {'Disease_smoking', '疾病吸烟组网络'}
% };
% 
% for i = 1:length(visualization_groups)
%     group_name = visualization_groups{i}{1};
%     group_title = visualization_groups{i}{2};
%     
%     if ~isempty(networks_results.(group_name))
%         fprintf('创建%s热力图...\n', group_title);
%         
%         % 过滤小的连接
%         filtered_matrix = networks_results.(group_name);
%         filtered_matrix(abs(filtered_matrix) < vis_threshold) = 0;
%         
%         try
%             figure('Name', group_title);
%             A = filtered_matrix;
%             symmetricMatrix = A + A.' - diag(diag(A));
%             circularGraph(symmetricMatrix, symmetricMatrix, 'Label', ROInames);
%             title(group_title);
%         catch ME
%             fprintf('  警告: %s可视化失败: %s\n', group_title, ME.message);
%         end
%     else
%         fprintf('跳过%s（数据不足）\n', group_title);
%     end
% end

% %% 按中心分组的协方差网络分析
% disp('===============================================');
% disp('=== 按中心分组的协方差网络分析 ===');
% disp('===============================================');
% 
% % 健康对照组（disease=0）按中心分组
% hc_center1_idx = (disease_info == 0) & (center_info == 1);  % HC中心1
% hc_center2_idx = (disease_info == 0) & (center_info == 2);  % HC中心2
% 
% % 疾病组（disease≠0）按中心分组
% disease_center2_idx = (disease_info ~= 0) & (center_info == 2);  % 疾病中心2
% disease_center3_idx = (disease_info ~= 0) & (center_info == 3);  % 疾病中心3
% 
% % 提取分组数据
% HC_center1_data = SUL_data_all(hc_center1_idx, :);
% HC_center2_data = SUL_data_all(hc_center2_idx, :);
% Disease_center2_data = SUL_data_all(disease_center2_idx, :);
% Disease_center3_data = SUL_data_all(disease_center3_idx, :);
% 
% % 显示分组统计
% fprintf('\n=== 按中心分组统计 ===\n');
% fprintf('健康对照组:\n');
% fprintf('  中心1: %d patients\n', size(HC_center1_data, 1));
% fprintf('  中心2: %d patients\n', size(HC_center2_data, 1));
% 
% fprintf('疾病组:\n');
% fprintf('  中心2: %d patients\n', size(Disease_center2_data, 1));
% fprintf('  中心3: %d patients\n', size(Disease_center3_data, 1));
% 
% % 定义要分析的中心组
% center_groups_to_analyze = {
%     {'HC_center1', HC_center1_data, '健康对照组中心1'};
%     {'HC_center2', HC_center2_data, '健康对照组中心2'};
%     {'Disease_center2', Disease_center2_data, '疾病组中心2'};
%     {'Disease_center3', Disease_center3_data, '疾病组中心3'}
% };
% 
% % 存储中心网络结果
% center_networks_results = struct();
% 
% % 为每个中心组建立协方差网络
% for i = 1:length(center_groups_to_analyze)
%     group_name = center_groups_to_analyze{i}{1};
%     group_data = center_groups_to_analyze{i}{2};
%     group_desc = center_groups_to_analyze{i}{3};
%     
%     if size(group_data, 1) >= 5
%         fprintf('\n建立%s协方差网络...\n', group_desc);
%         fprintf('  样本数: %d, ROI数: %d\n', size(group_data, 1), size(group_data, 2));
%         
%         % 使用Bootstrap方法建立稳定的协方差网络
%         network_matrix = calculate_group_network_bootstrap(group_data, n_boot, alpha, consistency_thresh, r_thresh);
%         center_networks_results.(group_name) = network_matrix;
%         
%         % 显示网络统计
%         nonzero_connections = sum(abs(network_matrix(:)) > 0) / 2;
%         max_strength = max(abs(network_matrix(:)));
%         fprintf('  %s协方差网络: %d个显著连接, 最大相关强度: %.4f\n', group_desc, nonzero_connections, max_strength);
%     else
%         fprintf('\n跳过%s（样本不足: %d < 5）\n', group_desc, size(group_data, 1));
%         center_networks_results.(group_name) = [];
%         end
%     end
%     
% % 保存中心网络结果
% fprintf('\n=== 保存按中心分组的网络结果 ===\n');
% 
% % 创建保存结构
% center_results_struct = struct();
% center_results_struct.ROInames = ROInames;
% center_results_struct.description = '按中心分组的全身器官代谢协方差网络矩阵';
% center_results_struct.groups = {
%     'HC_center1: 健康对照组中心1 (RJ)';
%     'HC_center2: 健康对照组中心2 (HN)';
%     'Disease_center2: 疾病组中心2 (HN)';
%     'Disease_center3: 疾病组中心3 (QFS)'
% };
% 
% % 保存网络数据
% field_names = fieldnames(center_networks_results);
% for i = 1:length(field_names)
%     field_name = field_names{i};
%     if ~isempty(center_networks_results.(field_name))
%         center_results_struct.(field_name) = center_networks_results.(field_name);
%         fprintf('✓ 保存%s网络矩阵\n', field_name);
%     else
%         fprintf('⚠ 跳过%s（数据不足）\n', field_name);
%         end
%     end
%     
% % 保存MAT文件
% save('Center_Based_Metabolic_Networks.mat', '-struct', 'center_results_struct');
% fprintf('\n结果已保存到: Center_Based_Metabolic_Networks.mat\n');

% % 显示保存的网络统计
% fprintf('\n=== 保存的网络统计 ===\n');
% for i = 1:length(field_names)
%     field_name = field_names{i};
%     if ~isempty(center_results_struct.(field_name))
%         network_matrix = center_results_struct.(field_name);
%         nonzero_connections = sum(abs(network_matrix(:)) > 0) / 2;
%         total_possible_connections = (size(network_matrix,1) * (size(network_matrix,1) - 1)) / 2;
%         network_density = nonzero_connections / total_possible_connections;
%         
%         fprintf('%s:\n', field_name);
%         fprintf('  网络大小: %dx%d\n', size(network_matrix));
%         fprintf('  显著连接数: %d / %d (%.2f%%)\n', nonzero_connections, total_possible_connections, network_density*100);
%         fprintf('  网络密度: %.4f\n', network_density);
%         fprintf('\n');
%     end
% end

%% 内外部网络拓扑相似性分析
disp('===============================================');
disp('=== 内外部网络拓扑相似性分析 ===');
disp('===============================================');

% 定义内外部对应比较
topology_comparisons = {
    {'Never_smoking_internal_LC_vs_Never_smoking_HC', 'Never_smoking_external_LC_vs_Never_smoking_HC', 'Never smoking LC vs HC (Internal vs External)'};
    {'Smoking_internal_LC_vs_Never_smoking_internal_LC', 'Smoking_external_LC_vs_Never_smoking_external_LC', 'Smoking effect in LC (Internal vs External)'};
    {'Smoking_internal_LC_vs_Smoking_HC', 'Smoking_external_LC_vs_Smoking_HC', 'Smoking LC vs HC (Internal vs External)'}
};

% 存储拓扑分析结果
topology_results = struct();

for i = 1:length(topology_comparisons)
    internal_network_name = topology_comparisons{i}{1};
    external_network_name = topology_comparisons{i}{2};
    comparison_desc = topology_comparisons{i}{3};
    
    fprintf('\n=== %s ===\n', comparison_desc);
    
    % 检查网络是否存在
    if isfield(difference_networks, internal_network_name) && isfield(difference_networks, external_network_name)
        internal_network = difference_networks.(internal_network_name);
        external_network = difference_networks.(external_network_name);
        
        % 计算Hub分析和Rich-club分析
        fprintf('计算Hub分析和Rich-club分析...\n');
        
        % 内部网络Hub和Rich-club分析
        internal_analysis = calculate_hub_and_richclub_analysis(internal_network, ROInames, sprintf('Internal_%d', i));
        external_analysis = calculate_hub_and_richclub_analysis(external_network, ROInames, sprintf('External_%d', i));
        
        % 存储结果
        comparison_key = sprintf('comparison_%d', i);
        topology_results.(comparison_key) = struct();
        topology_results.(comparison_key).description = comparison_desc;
        topology_results.(comparison_key).internal_network_name = internal_network_name;
        topology_results.(comparison_key).external_network_name = external_network_name;
        topology_results.(comparison_key).internal_analysis = internal_analysis;
        topology_results.(comparison_key).external_analysis = external_analysis;
        
        % 显示比较结果
        fprintf('\n--- Hub分析比较 ---\n');
        fprintf('Hub器官数量:\n');
        fprintf('  Internal: %d hubs\n', length(internal_analysis.hub_organs.combined_hubs));
        fprintf('  External: %d hubs\n', length(external_analysis.hub_organs.combined_hubs));
        
        fprintf('Hub器官重叠:\n');
        common_hubs = intersect(internal_analysis.hub_organs.combined_hubs, external_analysis.hub_organs.combined_hubs);
        fprintf('  共同Hub器官: %d (%s)\n', length(common_hubs), strjoin(common_hubs, ', '));
        
        fprintf('\n--- Rich-club分析比较 ---\n');
        fprintf('Rich-club强度 (平均):\n');
        fprintf('  Internal: %.4f\n', internal_analysis.rich_club.mean_coefficient);
        fprintf('  External: %.4f\n', external_analysis.rich_club.mean_coefficient);
        fprintf('  Difference: %.4f\n', abs(internal_analysis.rich_club.mean_coefficient - external_analysis.rich_club.mean_coefficient));
        
        fprintf('Rich-club成员器官:\n');
        fprintf('  Internal: %s\n', strjoin(internal_analysis.rich_club.member_organs, ', '));
        fprintf('  External: %s\n', strjoin(external_analysis.rich_club.member_organs, ', '));
        
        % 计算网络相似性
        network_similarity = calculate_network_similarity(internal_network, external_network);
        topology_results.(comparison_key).network_similarity = network_similarity;
        
        fprintf('Network Similarity (Pearson r): %.4f (p = %.4f)\n', network_similarity.pearson_r, network_similarity.pearson_p);
        
    else
        fprintf('跳过比较（缺少必要的网络数据）\n');
        if ~isfield(difference_networks, internal_network_name)
            fprintf('  缺少内部网络: %s\n', internal_network_name);
        end
        if ~isfield(difference_networks, external_network_name)
            fprintf('  缺少外部网络: %s\n', external_network_name);
        end
    end
end

% 保存拓扑分析结果
fprintf('\n=== 保存拓扑分析结果 ===\n');
topology_results.ROInames = ROInames;
topology_results.analysis_description = '内外部LC网络Hub分析和Rich-club分析';
topology_results.analysis_methods = {
    'Hub Analysis: 通过Degree、Betweenness、Strength识别关键器官节点';
    'Rich-club Analysis: 分析高连接度器官间的连接倾向'
};

save('Internal_External_Topology_Analysis_Results.mat', 'topology_results');
fprintf('拓扑分析结果已保存到: Internal_External_Topology_Analysis_Results.mat\n');

% 生成Hub和Rich-club分析总结
fprintf('\n=== Hub和Rich-club分析总结 ===\n');
for i = 1:length(topology_comparisons)
    comparison_key = sprintf('comparison_%d', i);
    if isfield(topology_results, comparison_key)
        result = topology_results.(comparison_key);
        fprintf('\n%d. %s\n', i, result.description);
        
        % Hub分析总结
        internal_hubs = result.internal_analysis.hub_organs.combined_hubs;
        external_hubs = result.external_analysis.hub_organs.combined_hubs;
        common_hubs = intersect(internal_hubs, external_hubs);
        
        fprintf('   Hub器官分析:\n');
        fprintf('     内部Hub数量: %d\n', length(internal_hubs));
        fprintf('     外部Hub数量: %d\n', length(external_hubs));
        fprintf('     共同Hub数量: %d (一致性: %.2f%%)\n', length(common_hubs), ...
            length(common_hubs)/max(length(internal_hubs), length(external_hubs))*100);
        
        % Rich-club分析总结
        fprintf('   Rich-club分析:\n');
        fprintf('     内部Rich-club强度: %.4f\n', result.internal_analysis.rich_club.mean_coefficient);
        fprintf('     外部Rich-club强度: %.4f\n', result.external_analysis.rich_club.mean_coefficient);
        fprintf('     Rich-club强度差异: %.4f\n', abs(result.internal_analysis.rich_club.mean_coefficient - result.external_analysis.rich_club.mean_coefficient));
        
        
        % 网络相似性
        if isfield(result, 'network_similarity')
            fprintf('   网络相似性: r=%.4f, p=%.4f\n', result.network_similarity.pearson_r, result.network_similarity.pearson_p);
        end
    end
end

disp('===============================================');
disp('=== 内外部LC全身代谢网络分析完成 ===');
disp('===============================================');

%% 核心函数定义

% 群体协方差网络计算函数（Bootstrap版本）
function network_matrix = calculate_group_network_bootstrap(group_data, n_boot, alpha, consistency_thresh, r_thresh)
    % 计算群体的协方差网络矩阵，使用Bootstrap确保稳定性
    % 每次Bootstrap采样后计算Pearson相关性矩阵，然后进行FDR校正
    % 输入: group_data - 患者数据矩阵 (patients × ROIs)
    %      n_boot - Bootstrap重采样次数
    %      alpha - 显著性水平
    %      consistency_thresh - 一致性阈值（多少比例的Bootstrap样本中显著）
    %      r_thresh - 相关性强度阈值
    % 输出: network_matrix - 协方差网络矩阵 (ROIs × ROIs)
    
    fprintf('  开始Bootstrap协方差网络分析...\n');
    fprintf('    数据: %dx%d (患者x ROI)\n', size(group_data));
    
    % 数据质量检查
    nan_ratio = sum(isnan(group_data(:))) / numel(group_data);
    if nan_ratio > 0.1
        warning('数据中包含较多缺失值: %.1f%%', nan_ratio*100);
    end
    
    n_rois = size(group_data, 2);
    z_boot = @(r) 0.5 * log((1 + r) ./ (1 - r)); % Fisher z变换
    inv_z = @(z) (exp(2*z) - 1) ./ (exp(2*z) + 1); % Fisher z逆变换
    
    % 存储相关性矩阵
    r_stack = zeros(n_rois, n_rois, n_boot);
    sig_count = zeros(n_rois);
    valid_boot_count = 0;
    
    for b = 1:n_boot
        if mod(b, 200) == 0
            fprintf('    Bootstrap进度: %d/%d\n', b, n_boot);
        end
        
        % Bootstrap采样
        sample_idx = randi(size(group_data, 1), size(group_data, 1), 1);
        sample_data = group_data(sample_idx, :);
        
        % 检查样本质量
        if sum(isnan(sample_data(:))) > numel(sample_data) * 0.5
            continue;
        end
        
        % 计算协方差网络（Pearson相关性矩阵）
        try
            [R, P] = corr(sample_data, 'Rows', 'pairwise', 'Type', 'Pearson');
            
            % 处理NaN值
            R(isnan(R)) = 0; 
            P(isnan(P)) = 1;
            
            % 提取上三角部分进行显著性检验（避免重复计算对称连接）
            upper_idx = find(triu(true(n_rois), 1));
            
            % FDR校正：控制假发现率
            [h, ~, ~, ~] = fdr_bh(P(upper_idx), alpha);
            
            % 双重筛选：既要统计显著，又要达到强度阈值
            strong_corr = abs(R(upper_idx)) > r_thresh;
            
            % 构建显著性掩码矩阵
            sig = false(n_rois);
            sig(upper_idx) = h & strong_corr;  % 同时满足显著性和强度要求
            sig = sig | sig';  % 确保矩阵对称性
            
            % 累计显著性计数（用于后续一致性检验）
            sig_count = sig_count + sig;
            
            % Fisher z变换存储相关性值（提高数值稳定性）
            z_r = z_boot(R);
            z_r(isnan(z_r) | isinf(z_r)) = 0;
            
            r_stack(:, :, b) = z_r;
            valid_boot_count = valid_boot_count + 1;
            
        catch ME
            warning('Bootstrap %d 失败: %s', b, ME.message);
            continue;
        end
    end
    
    fprintf('    完成 %d/%d 有效Bootstrap\n', valid_boot_count, n_boot);
    
    if valid_boot_count == 0
        warning('没有有效的Bootstrap样本');
        network_matrix = zeros(n_rois);
        return;
    end
    
    % 计算Bootstrap平均协方差网络
    r_z_mean = mean(r_stack(:,:,1:valid_boot_count), 3, 'omitnan');  % 平均Fisher z值
    r_mean = inv_z(r_z_mean);  % 转换回相关性系数
    r_mean(isnan(r_mean) | isinf(r_mean)) = 0;  % 清理无效值
    
    % 一致性筛选：只保留在足够比例的Bootstrap样本中都显著的连接
    consistency_mask = sig_count >= round(valid_boot_count * consistency_thresh);
    
    % 构建最终协方差网络：平均相关性 × 一致性掩码
    network_matrix = r_mean .* consistency_mask;
    
    % 确保网络矩阵的数学性质
    network_matrix = (network_matrix + network_matrix') / 2;  % 对称性
    network_matrix(logical(eye(size(network_matrix)))) = 0;   % 对角线为0
    
    % 输出网络统计信息
    total_sig_edges = sum(consistency_mask(:)) / 2;  % 一致显著连接数
    final_edges = sum(abs(network_matrix(:)) > 0) / 2;  % 最终网络连接数
    
    fprintf('    结果: 一致性连接%d, 最终显示%d\n', total_sig_edges, final_edges);
    fprintf('    相关性强度范围: [%.4f, %.4f]\n', min(network_matrix(:)), max(network_matrix(:)));
end

% FDR校正函数（Benjamini-Hochberg方法）
function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)
    if nargin < 2, q = 0.05; end
    if nargin < 3, method = 'pdep'; end
    if nargin < 4, report = 'no'; end
    
    pvals = pvals(:);
    s = length(pvals);
    
    if s == 0
        h = []; crit_p = []; adj_ci_cvrg = []; adj_p = [];
        return;
    end
    
    [pvals_sorted, sort_ids] = sort(pvals);
    
    if strcmpi(method, 'pdep')
        thresh = (1:s)' * q / s;
    else
        thresh = (1:s)' * q / (s * sum(1./(1:s)));
    end
    
    wtd_p = pvals_sorted ./ thresh;
    
    if all(wtd_p > 1)
        h = zeros(s, 1);
        crit_p = 0;
        adj_ci_cvrg = 1 - q;
    else
        rej = find(pvals_sorted <= thresh);
        max_id = rej(end);
        
        h = zeros(s, 1);
        h(sort_ids(1:max_id)) = 1;
        crit_p = thresh(max_id);
        
        if strcmpi(method, 'pdep')
            adj_ci_cvrg = 1 - q;
        else
            adj_ci_cvrg = 1 - q * sum(1./(1:s));
        end
    end
    
    adj_p = zeros(s, 1);
    adj_p(sort_ids) = cummin(flipud(cummax(flipud(wtd_p))));
    adj_p = min(adj_p, 1);
end

% 网络相关性计算函数
function [pearson_r, pearson_p, valid_n] = calculate_network_correlation(network1, network2)
    % 计算两个网络之间的相关性
    
    if ~isequal(size(network1), size(network2))
        error('网络尺寸必须相同');
    end
    
    % 提取上三角
    upper_idx = find(triu(true(size(network1)), 1));
    net1_vec = network1(upper_idx);
    net2_vec = network2(upper_idx);
    
    % 找到非零连接（至少一个网络中非零）
    valid_idx = ((net1_vec ~= 0) | (net2_vec ~= 0)) & ~isnan(net1_vec) & ~isnan(net2_vec);
    valid_n = sum(valid_idx);
    
    if valid_n < 3
        fprintf('    警告: 有效连接太少(%d)，无法计算可靠相关性\n', valid_n);
        pearson_r = NaN; 
        pearson_p = NaN;
        return;
    end
    
    valid_net1 = net1_vec(valid_idx);
    valid_net2 = net2_vec(valid_idx);
    
    % Pearson相关性
    try
        [pearson_r, pearson_p] = corr(valid_net1, valid_net2, 'Type', 'Pearson');
    catch
        pearson_r = NaN; 
        pearson_p = NaN;
    end
end

% 查找网络差异中最强的N个连接
function [top_connections, top_values] = find_top_network_differences(diff_network, roi_names, top_n)
    % 找出差异网络中绝对值最大的前N个连接
    % 输入: diff_network - 差异网络矩阵
    %      roi_names - ROI名称列表
    %      top_n - 要提取的连接数量
    % 输出: top_connections - 前N个连接的ROI名称对 {N×2 cell}
    %      top_values - 对应的差异值 [N×1 double]
    
    % 提取上三角部分（避免重复）
    [row_idx, col_idx] = find(triu(true(size(diff_network)), 1));
    
    % 获取对应的差异值
    diff_values = zeros(length(row_idx), 1);
    for i = 1:length(row_idx)
        diff_values(i) = diff_network(row_idx(i), col_idx(i));
    end
    
    % 按绝对值排序
    [~, sort_idx] = sort(abs(diff_values), 'descend');
    
    % 提取前N个
    n_extract = min(top_n, length(sort_idx));
    top_indices = sort_idx(1:n_extract);
    
    % 构建结果
    top_connections = cell(n_extract, 2);
    top_values = zeros(n_extract, 1);
    
    for i = 1:n_extract
        idx = top_indices(i);
        roi1_idx = row_idx(idx);
        roi2_idx = col_idx(idx);
        
        top_connections{i, 1} = roi_names{roi1_idx};
        top_connections{i, 2} = roi_names{roi2_idx};
        top_values(i) = diff_values(idx);
    end
end

% Hub分析和Rich-club分析函数
function analysis_results = calculate_hub_and_richclub_analysis(network_matrix, roi_names, network_name)
    % 计算网络的Hub分析和Rich-club分析
    % 输入: network_matrix - 网络矩阵
    %      roi_names - ROI名称列表
    %      network_name - 网络名称（用于调试）
    % 输出: analysis_results - 包含Hub和Rich-club分析结果的结构体
    
    fprintf('  计算%s的Hub和Rich-club分析...\n', network_name);
    
    % 初始化结果结构
    analysis_results = struct();
    
    % 将网络矩阵转换为绝对值矩阵（保留连接强度）
    abs_network = abs(network_matrix);
    binary_network = abs_network > 0;  % 二进制网络用于度计算
    
    n_nodes = size(abs_network, 1);
    
    % === 1. Hub Analysis ===
    fprintf('    进行Hub分析...\n');
    
    % 计算各种中心性指标
    hub_analysis = struct();
    
    % Degree Centrality (度中心性)
    degrees = sum(binary_network, 2);
    hub_analysis.degrees = degrees;
    
    % Strength Centrality (强度中心性) - 基于连接权重
    strengths = sum(abs_network, 2);
    hub_analysis.strengths = strengths;
    
    % Betweenness Centrality (介数中心性) - 简化计算
    betweenness = calculate_betweenness_centrality(binary_network);
    hub_analysis.betweenness = betweenness;
    
    % 设定Hub阈值（前20%的节点）
    hub_threshold_percentile = 80;  % 前20%
    
    % 基于不同指标识别Hub
    degree_threshold = prctile(degrees, hub_threshold_percentile);
    strength_threshold = prctile(strengths, hub_threshold_percentile);
    betweenness_threshold = prctile(betweenness, hub_threshold_percentile);
    
    degree_hubs = find(degrees >= degree_threshold);
    strength_hubs = find(strengths >= strength_threshold);
    betweenness_hubs = find(betweenness >= betweenness_threshold);
    
    % 综合Hub识别（至少满足两个条件）
    all_hub_candidates = unique([degree_hubs; strength_hubs; betweenness_hubs]);
    combined_hubs = [];
    
    for i = 1:length(all_hub_candidates)
        node = all_hub_candidates(i);
        hub_count = 0;
        if any(degree_hubs == node), hub_count = hub_count + 1; end
        if any(strength_hubs == node), hub_count = hub_count + 1; end
        if any(betweenness_hubs == node), hub_count = hub_count + 1; end
        
        if hub_count >= 2  % 至少满足两个条件
            combined_hubs = [combined_hubs; node];
        end
    end
    
    % 存储Hub结果
    hub_organs = struct();
    hub_organs.degree_hubs = roi_names(degree_hubs);
    hub_organs.strength_hubs = roi_names(strength_hubs);
    hub_organs.betweenness_hubs = roi_names(betweenness_hubs);
    hub_organs.combined_hubs = roi_names(combined_hubs);
    
    analysis_results.hub_organs = hub_organs;
    analysis_results.hub_metrics = hub_analysis;
    
    fprintf('      识别到%d个Hub器官: %s\n', length(combined_hubs), strjoin(roi_names(combined_hubs), ', '));
    
    % === 2. Rich-club Analysis ===
    fprintf('    进行Rich-club分析...\n');
    
    rich_club_analysis = struct();
    
    % 计算Rich-club系数
    [rich_club_coeffs, degree_levels] = calculate_rich_club_coefficient_advanced(abs_network);
    
    % 找到显著的Rich-club成员（高度节点）
    high_degree_threshold = prctile(degrees, 75);  % 前25%高度节点
    rich_club_members = find(degrees >= high_degree_threshold);
    
    % 计算Rich-club内部连接强度
    if length(rich_club_members) >= 2
        rich_club_subnetwork = abs_network(rich_club_members, rich_club_members);
        rich_club_density = (sum(rich_club_subnetwork(:)) - trace(rich_club_subnetwork)) / ...
                           (length(rich_club_members) * (length(rich_club_members) - 1));
        
        rich_club_analysis.members = rich_club_members;
        rich_club_analysis.member_organs = roi_names(rich_club_members);
        rich_club_analysis.internal_density = rich_club_density;
        rich_club_analysis.coefficients = rich_club_coeffs;
        rich_club_analysis.degree_levels = degree_levels;
        rich_club_analysis.mean_coefficient = mean(rich_club_coeffs(~isnan(rich_club_coeffs)));
        
        fprintf('      Rich-club成员器官: %s\n', strjoin(roi_names(rich_club_members), ', '));
        fprintf('      Rich-club内部密度: %.4f\n', rich_club_density);
    else
        rich_club_analysis.members = [];
        rich_club_analysis.member_organs = {};
        rich_club_analysis.internal_density = 0;
        rich_club_analysis.mean_coefficient = 0;
        fprintf('      Rich-club成员不足，无法进行分析\n');
    end
    
    analysis_results.rich_club = rich_club_analysis;
    
    % 分析完成
    fprintf('    Hub和Rich-club分析完成\n');
end

% 模块性计算函数
function Q = calculate_modularity(adj_matrix)
    % 使用Newman-Girvan算法计算模块性
    % 这里使用简化版本的模块性计算
    
    n = size(adj_matrix, 1);
    m = sum(adj_matrix(:)) / 2;  % 总边数
    
    if m == 0
        Q = 0;
        return;
    end
    
    % 计算节点度
    k = sum(adj_matrix, 2);
    
    % 使用贪婪算法进行社区检测
    % 初始化每个节点为一个社区
    communities = 1:n;
    
    % 计算初始模块性
    Q = 0;
    for i = 1:n
        for j = 1:n
            if communities(i) == communities(j)
                Q = Q + (adj_matrix(i,j) - k(i)*k(j)/(2*m));
            end
        end
    end
    Q = Q / (2*m);
    
    % 简化的社区合并过程（只进行几轮迭代）
    max_iterations = min(10, n-1);
    for iter = 1:max_iterations
        best_Q = Q;
        best_merge = [];
        
        % 尝试合并相邻的社区
        unique_communities = unique(communities);
        for c1_idx = 1:length(unique_communities)
            for c2_idx = (c1_idx+1):length(unique_communities)
                c1 = unique_communities(c1_idx);
                c2 = unique_communities(c2_idx);
                
                % 临时合并社区
                temp_communities = communities;
                temp_communities(communities == c2) = c1;
                
                % 计算新的模块性
                temp_Q = 0;
                for i = 1:n
                    for j = 1:n
                        if temp_communities(i) == temp_communities(j)
                            temp_Q = temp_Q + (adj_matrix(i,j) - k(i)*k(j)/(2*m));
                        end
                    end
                end
                temp_Q = temp_Q / (2*m);
                
                if temp_Q > best_Q
                    best_Q = temp_Q;
                    best_merge = [c1, c2];
                end
            end
        end
        
        % 如果找到更好的合并，执行合并
        if ~isempty(best_merge)
            communities(communities == best_merge(2)) = best_merge(1);
            Q = best_Q;
        else
            break;  % 没有改进，停止迭代
        end
    end
end

% 介数中心性计算函数（简化版本）
function betweenness = calculate_betweenness_centrality(adj_matrix)
    % 计算介数中心性（简化版本，基于度和强度的近似）
    % 为了避免复杂的路径计算，使用简化的近似方法
    
    n = size(adj_matrix, 1);
    
    % 计算度中心性
    degrees = sum(adj_matrix > 0, 2);
    
    % 计算邻接矩阵的平方，用于估算2步连接
    adj_squared = adj_matrix * adj_matrix;
    
    % 简化的介数中心性估算：基于节点在局部网络中的重要性
    betweenness = zeros(n, 1);
    
    for i = 1:n
        % 找到节点i的邻居
        neighbors = find(adj_matrix(i, :) > 0);
        
        if length(neighbors) > 1
            % 计算通过节点i连接的邻居对数量
            neighbor_pairs = 0;
            for j = 1:length(neighbors)
                for k = (j+1):length(neighbors)
                    neighbor1 = neighbors(j);
                    neighbor2 = neighbors(k);
                    % 如果两个邻居不直接相连，则节点i在它们之间起到桥梁作用
                    if adj_matrix(neighbor1, neighbor2) == 0
                        neighbor_pairs = neighbor_pairs + 1;
                    end
                end
            end
            
            % 介数中心性与能够桥接的邻居对数量成正比
            betweenness(i) = neighbor_pairs / (degrees(i) + 1);
        end
    end
    
    % 标准化
    if max(betweenness) > 0
        betweenness = betweenness / max(betweenness);
    end
end

% 高级Rich-club系数计算函数
function [rich_club_coeffs, degree_levels] = calculate_rich_club_coefficient_advanced(weight_matrix)
    % 计算基于权重的Rich-club系数
    
    % 计算加权度（strength）
    strengths = sum(abs(weight_matrix), 2);
    unique_strengths = unique(strengths);
    unique_strengths = sort(unique_strengths(unique_strengths > 0), 'descend');
    
    % 只取前10个强度水平以避免计算过载
    degree_levels = unique_strengths(1:min(10, length(unique_strengths)));
    rich_club_coeffs = zeros(size(degree_levels));
    
    for i = 1:length(degree_levels)
        k = degree_levels(i);
        
        % 找到强度≥k的节点
        high_strength_nodes = find(strengths >= k);
        n_high_strength = length(high_strength_nodes);
        
        if n_high_strength < 2
            rich_club_coeffs(i) = NaN;
            continue;
        end
        
        % 计算这些高强度节点之间的连接权重总和
        subgraph = abs(weight_matrix(high_strength_nodes, high_strength_nodes));
        actual_weight_sum = (sum(subgraph(:)) - trace(subgraph)) / 2;  % 排除对角线，除以2避免重复计算
        
        % 计算可能的最大权重总和（使用网络中最强的连接）
        all_weights = abs(weight_matrix(:));
        all_weights = all_weights(all_weights > 0);
        all_weights = sort(all_weights, 'descend');
        
        max_possible_edges = n_high_strength * (n_high_strength - 1) / 2;
        if max_possible_edges <= length(all_weights)
            max_possible_weight_sum = sum(all_weights(1:max_possible_edges));
        else
            max_possible_weight_sum = sum(all_weights);
        end
        
        % Rich-club系数
        if max_possible_weight_sum > 0
            rich_club_coeffs(i) = actual_weight_sum / max_possible_weight_sum;
        else
            rich_club_coeffs(i) = NaN;
        end
    end
end


% 网络相似性计算函数
function similarity_metrics = calculate_network_similarity(network1, network2)
    % 计算两个网络之间的相似性
    
    if ~isequal(size(network1), size(network2))
        error('网络尺寸必须相同');
    end
    
    % 提取上三角部分（避免重复和对角线）
    upper_idx = find(triu(true(size(network1)), 1));
    net1_vec = network1(upper_idx);
    net2_vec = network2(upper_idx);
    
    % 找到非零连接（至少一个网络中非零）
    valid_idx = ((net1_vec ~= 0) | (net2_vec ~= 0)) & ~isnan(net1_vec) & ~isnan(net2_vec);
    valid_n = sum(valid_idx);
    
    similarity_metrics = struct();
    similarity_metrics.valid_connections = valid_n;
    
    if valid_n < 3
        fprintf('    警告: 有效连接太少(%d)，无法计算可靠相似性\n', valid_n);
        similarity_metrics.pearson_r = NaN;
        similarity_metrics.pearson_p = NaN;
        return;
    end
    
    valid_net1 = net1_vec(valid_idx);
    valid_net2 = net2_vec(valid_idx);
    
    % Pearson相关性
    try
        [r, p] = corr(valid_net1, valid_net2, 'Type', 'Pearson');
        similarity_metrics.pearson_r = r;
        similarity_metrics.pearson_p = p;
    catch
        similarity_metrics.pearson_r = NaN;
        similarity_metrics.pearson_p = NaN;
    end
end


