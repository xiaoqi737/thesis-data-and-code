% 主脚本：根据标签选择调用函数并循环增加u_alpha
clc; clear; close all;

% 假设您的数据存储在变量data中（31行6列）
% 前5列：函数参数，最后一列：标签(0或1)
% 请将您的数据加载到data变量中

% 参数初始化
delta = 1;  % 设置delta值，根据实际情况调整
max_p = 0.99;  % 设置u_alpha的最大值，防止无限循环
p_step = 0.01;  % u_alpha的步长

% 存储结果的结构体数组
results = struct();
data_chanchu_juece = readmatrix("D:\内农大\学位论文\原值.xlsx",'Sheet','Sheet5');
data_touru_juece = readmatrix("D:\内农大\学位论文\原值.xlsx",'Sheet','Sheet3');
data_touru_sample_ex=readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet1');%不动
data_touru_sample_sig=readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet2');%不动
data_chanchu_sample_ex=readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet3');%不动
data_chanchu_sample_sig=readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet4');%不动
data = [data_touru_juece(2:32,2:4),data_chanchu_juece(2:32,5:6)];
label_zong = readmatrix("D:\内农大\学位论文\最终效率值.xlsx",'Sheet','Sheet2');
% 循环处理每一行数据
for i = 1:31
    fprintf('处理第%d行数据...\n', i);
    results(i).p_tried = [];   % 记录尝试过的 p 值
    results(i).exitflag = [];   % 记录尝试过的 p 值
    % 提取当前行的数据和标签
    row_data = data(i, 1:5);
    label = label_zong(i,2);
     x_ip = row_data(1:3)';     % 前两列，转换为列向量
     y_rp = row_data(4:5)';     % 第3-4列，转换为列向量
    % 注意：这里假设数据列的对应关系，请根据您的实际数据结构调整
        m = length(x_ip);
        l = length(y_rp);
        
        % 确保 mu_bar 和 sigma2 有正确的行数
       mu_bar = [data_touru_sample_ex(2:32,2:4),data_chanchu_sample_ex(2:32,5:6)]';
       sigma2 = [data_touru_sample_sig(2:32,2:4),data_chanchu_sample_sig(2:32,5:6)]';
        % 第5列可能是其他参数，根据您的实际情况调整
    if label == 0
        % 标签为0，调用solve_theta函数
        fprintf('标签为0，调用solve_theta函数\n');
        
        % 提取函数1需要的参数
        
        
        % 初始化u_alpha
        p = 0.5;
        %exitflag = 1;  % 初始化为有效解
        
        % 循环增加u_alpha直到无解
        while p > 0 && p <= max_p
            [lambda, omega_minus, omega_plus, theta, exitflag, z, fval] = solve_theta_chuantong(mu_bar, sigma2, x_ip, y_rp, p, delta);
             results(i).p_tried(end+1) = p;
            results(i).exitflag(end+1) = exitflag;
            if (0 < theta) && (theta < 1)||(exitflag==1)
                % 有解，存储当前结果
                results(i).label = 0;
                results(i).p = p;
                results(i).lambda = lambda;
                results(i).omega_minus = omega_minus;
                results(i).omega_plus = omega_plus;
                results(i).theta = theta;
                results(i).exitflag = exitflag;
                results(i).z = z;
                results(i).fval = fval;
                
                % 增加u_alpha继续循环
                fprintf('  u_alpha=%.4f 有解，theta=%.4f\n', p, theta);
                p = p + p_step;
            else
                % 无解，退出循环
                break;
                
            end
            %fprintf('  u_alpha=%.4f 无解，退出循环\n', p);
            
        end
        
    else
        % 标签为1，调用solve_rho函数
        fprintf('  标签为1，调用solve_rho函数\n');
        
        % 提取函数2需要的参数
        % 注意：这里假设数据列的对应关系，请根据您的实际数据结构调整
       % mu_bar = row_data(1);
        %sigma2_bar = row_data(2);
       % mu_ip = row_data(3);
        %sigma_ip = row_data(4);
        % 第5列可能是其他参数，根据您的实际情况调整
        
        % 初始化缺失的参数（根据您的实际情况调整）
        %mu_mrp_p = 0;  % 请替换为实际值
        %sigma_mrp_p = 0;  % 请替换为实际值
        %d = 0;  % 请替换为实际值
        
        % 初始化u_alpha
        p = 0.5;
        %exitflag = 2;  % 初始化为有效解
        %rho = 1.1;  % 初始化rho
        
        % 循环增加u_alpha直到theta=1或无解
            while p > 0 && p <= max_p
            [lambda, omega_minus, omega_plus, rho, exitflag, z, fval] = solve_theta(mu_bar, sigma2, x_ip, y_rp, p, delta);
            results(i).p_tried(end+1) = p;
            results(i).exitflag(end+1) = exitflag;
                if rho > 1||(exitflag==1)
                % 有解，存储当前结果
                results(i).label = 1;
                results(i).p = p;
                results(i).lambda = lambda;
                results(i).omega_minus = omega_minus;
                results(i).omega_plus = omega_plus;
                results(i).rho = rho;
                results(i).exitflag = exitflag;
                results(i).z = z;
                results(i).fval = fval;
                
                fprintf('  u_alpha=%.4f 有解，rho=%.4f\n', p, rho);
                p = p - p_step;
                else
                    % rho=1，退出循环
                    break;
                end
            
                % 无解，退出循环
                %fprintf('  u_alpha=%.4f 无解，退出循环\n', p);
            end
            
     end
 end
    
    % 记录循环终止时的u_alpha
    results(i).final_p = p;


% 显示汇总结果
fprintf('\n========== 处理完成 ==========\n');
fprintf('总处理行数: %d\n', size(data, 1));

% 分别统计标签0和标签1的结果
label0_count = 0;
label1_count = 0;

for i = 1:length(results)
    if isfield(results(i), 'label')
        if results(i).label == 0
            label0_count = label0_count + 1;
            fprintf('第%d行（标签0）: 最终u_alpha=%.4f, 最终theta=%.4f\n', ...
                i, results(i).p, results(i).theta);
        else
            label1_count = label1_count + 1;
            fprintf('第%d行（标签1）: 最终u_alpha=%.4f, 最终rho=%.4f\n', ...
                i, results(i).p, results(i).rho);
        end
    end
end

fprintf('\n标签统计: 0标签=%d行, 1标签=%d行\n', label0_count, label1_count);