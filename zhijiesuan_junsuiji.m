% 主脚本：根据标签选择调用函数并循环增加u_alpha
clc; clear; close all;

% 假设您的数据存储在变量data中（31行6列）
% 前5列：函数参数，最后一列：标签(0或1)
% 请将您的数据加载到data变量中

% 参数初始化
% 假设您的数据存储在变量 data 中，31行6列，最后一列为标签
% 假设数据列含义：
%   第1-2列: x_ip 的两个元素（假设 m=2）
%   第3-4列: y_rp 的两个元素（假设 l=2）
%   第5列: 其他参数（如初始 p 值或 delta）
%   第6列: 标签（0或1）

% 预定义的 mu_bar 和 sigma2 矩阵
data_juece_mu = readmatrix("D:\内农大\学位论文\算例数据.xlsx",'Sheet','Sheet1');
data_juece_sig = readmatrix("D:\内农大\学位论文\算例数据.xlsx",'Sheet','Sheet3');
data_sample_mu=readmatrix("D:\内农大\学位论文\算例数据.xlsx",'Sheet','Sheet2');
data_sample_sig=readmatrix("D:\内农大\学位论文\算例数据.xlsx",'Sheet','Sheet4');
%data_chanchu_juece_mu = readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet13');
%data_touru_juece_mu = readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet11');
%data_chanchu_juece_sig = readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet14');
%data_touru_juece_sig = readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet12');
%data_chanchu_sample_mu=readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet3');%不动
%data_touru_sample_mu=readmatrix("D:\内农大\学位论文\期望值.xlsx",'Sheet','Sheet1');%不动
%data_chanchu_sample_sig=readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet4');%不动
%data_touru_sample_sig=readmatrix("D:\内农大\学位论文\方差值.xlsx",'Sheet','Sheet2');%不动

% 参数
delta = 1;
d = 1;
results = struct();
% 处理每一行数据
for i = 1:4
    fprintf('\n=== 处理第%d行数据 ===\n', i);
    %row_data = data(i+1, 2:6);
    mu_ip =  data_juece_mu(i,1:2)';
    sigma_ip = sqrt(data_juece_sig(i,1:2))';
    mu_mrp_p = data_juece_mu(i,3:4)';
    sigma_mrp_p = sqrt(data_juece_sig(i,3:4))';
    % 提取参数
        %mu_ip = row_data(1:3)';     % 前两列，转换为列向量
        %y_rp = row_data(4:5)';     % 第3-4列，转换为列向量
        
        % 初始 p 值可以从数据获取或使用固定值
        p = 0.95;  % 或从数据中获取: row_data(5)
        
        % 检查维度
        m = length(mu_ip);
        l = length(mu_mrp_p);
        
        % 确保 mu_bar 和 sigma2 有正确的行数
        mu_bar = data_sample_mu';
        sigma2_bar = data_sample_sig';
        
        
            
            [lambda, omega_minus, omega_plus, rho, exitflag, z, fval] = solve_rho(...
    mu_bar, sigma2_bar, mu_ip, sigma_ip, mu_mrp_p, sigma_mrp_p, p, d, delta);
            results(i).row_index = i;
            results(i).p = p;
            results(i).theta = rho;
            results(i).exitflag = exitflag;
            results(i).fval = fval;
            results(i).lambda = lambda;
            results(i).omega_minus = omega_minus;
            results(i).omega_plus = omega_plus;
            results(i).z = z;
            
            
end