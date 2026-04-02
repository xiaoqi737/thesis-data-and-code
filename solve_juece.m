% solve_rho.m
% 使用 fmincon 求解图中给定的优化问题（Image 2）
%决策单元随机
% 变量顺序: z = [lambda (nbar x 1); omega_minus (m x 1); omega_plus (l x 1)]
%
% 输入说明:
%   x_bar        : (m + l) x nbar 矩阵， mū_{i,j}（前 m 行对应 i=1..m，后 l 行对应 m+r）
%   y_bar    : (m + l) x nbar 矩阵， σ̄_{i,j}^2
%   mu_ip         : m x 1 向量， μ_{i,p}
%   sigma_ip      : m x 1 向量， σ_{i,p}
%   mu_mrp_p      : l x 1 向量， μ_{m+r,p}
%   sigma_mrp_p   : l x 1 向量， σ_{m+r,p}
%   u_alpha       : 标量 u_α
%   d             : 标量 d（出现在第二组约束）
%   delta         : 标量 δ（若 δ~=0 则等式约束为 sum(lambda)=1）
%
% 输出:
%   lambda, omega_minus, omega_plus, rho (目标值), exitflag, z, fval
%
% 使用示例见文件末尾。

function [lambda, omega_minus, omega_plus, rho, exitflag, z, fval] = solve_juece(...
    x_bar, y_bar, mu_ip, sigma_ip, mu_mrp_p, sigma_mrp_p, p, delta)

if nargin < 9, delta = 1; end

% 尺寸检查
[nRows, nbar] = size(x_bar);
 m = numel(mu_ip);
 l = numel(mu_mrp_p);
% if nRows ~= m + l
%     error('x_bar 的行数应为 m + l');
% end
% if any(size(y_bar) ~= size(x_bar))
%     error('sigma2_bar 应与 mu_bar 尺寸相同');
% end
% if numel(sigma_ip) ~= m || numel(sigma_mrp_p) ~= l
%     error('sigma_ip 或 sigma_mrp_p 维度错误');
% end

% 检查分母正性（在目标的分数里需要）
 denom_num = mu_ip - norminv(p) .* sigma_ip; % 应该为正
 denom_den = mu_mrp_p + norminv(p) .* sigma_mrp_p; % 应该为正
% if any(denom_num <= 0)
%     error('存在 i 使得 mu_ip - u_alpha*sigma_ip <= 0，需调整输入参数');
% end
% if any(denom_den <= 0)
%     error('存在 r 使得 mu_{m+r,p} + u_alpha*sigma_{m+r,p} <= 0，需调整输入参数');
% end

% 变量数
nvar = nbar + m + l;

% 初始点
lambda0_rand = rand(nbar, 1); lambda0 = lambda0_rand / sum(lambda0_rand);
omega_minus0 = 0.01 + 0.1 * rand(m, 1);
omega_plus0 = 0.01 + 0.1 * rand(l, 1);
x0 = [lambda0; omega_minus0; omega_plus0];

% 下界
lb = zeros(nvar,1);
ub = [];

% fmincon 选项
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',1000);

% 目标和约束句柄
obj = @(z) obj_rho(z, m, l, nbar, denom_num, denom_den);
nonlcon = @(z) nlcons_rho(z, x_bar, y_bar, mu_ip, sigma_ip, mu_mrp_p, sigma_mrp_p, p, delta, m, l, nbar);

% 调用 fmincon
[z,fval,exitflag,output] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, options);

% 解析结果
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : end);
rho = obj(z);

fprintf('fmincon exitflag = %d, final objective = %.8g\n', exitflag, fval);
end

%% 目标函数
function r = obj_rho(z, m, l, nbar, denom_num, denom_den)
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : nbar + m + l);

% 分母稳定性保证
eps_den = 1e-9;
num_term = 1 + (1/m) * sum( omega_minus ./ denom_num );
den_term = 1 - (1/l) * sum( omega_plus ./ denom_den );
if den_term <= eps_den
    r = 1e20 + (1/l)*sum(max(0, omega_plus ./ denom_den));
    return;
end
r = num_term / den_term;
end

%% 非线性不等式/等式约束: c(z) <= 0, ceq(z) == 0
function [c, ceq] = nlcons_rho(z, x_bar, y_bar, mu_ip, sigma_ip, mu_mrp_p, sigma_mrp_p, p, delta, m, l, nbar)
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : nbar + m + l);

% 约束数： m + l + 1 (最后一项用于保证目标分母正)
c = zeros(m + l + l + 1, 1);

% i = 1..m: omega_i^- + mu_ip - sum_j mū_{i,j} lambda_j >= u_alpha * sqrt( sum_j sigmā_{i,j}^2 lambda_j^2 + sigma_ip^2 )
for i = 1:m
    left = omega_minus(i) + mu_ip(i) - x_bar(i,:) * lambda;
    rhs = -norminv(p) * sqrt( sigma_ip(i)^2 );
    % left >= rhs  <=> rhs - left <= 0
    c(i) = rhs - left;
end

% r = 1..l: omega_r^+ - mu_{m+r,p} + d * sum_j mū_{m+r,j} lambda_j >= u_alpha * sqrt( sum_j sigmā_{m+r,j}^2 lambda_j^2 + sigma_{m+r,p}^2 )
for r = 1:l
    
    left = omega_plus(r) - mu_mrp_p(r) + ( y_bar(r,:) * lambda );
    rhs = norminv(p) * sqrt( sigma_mrp_p(r)^2 );
    c(m + r) = rhs - left;
end
for rr = 1:l
    denom_den = mu_mrp_p + norminv(p) .* sigma_mrp_p;
    %left = - y_rp(rr) + mu_bar(idx,:) * lambda + omega_plus(rr);
    %rhs = u_alpha * sqrt( sum( sigma2(idx,:) .* (lambda'.^2) ) );
    % 要求 left >= rhs -> rhs - left <= 0
    c(m +l+ rr) = omega_plus(rr) - denom_den(rr);
end
% 额外约束：保证目标分母不为零且为正
denom_den = mu_mrp_p + norminv(p) .* sigma_mrp_p;
eps_den = 1e-8;
c(end) = (1/l) * sum( omega_plus ./ denom_den ) - (1 - eps_den);

% 等式约束： δ * sum(lambda) = δ  -> 当 δ~=0 等价于 sum(lambda)=1
if abs(delta) < eps
    ceq = [];
else
    ceq = sum(lambda) - 1;
end
end