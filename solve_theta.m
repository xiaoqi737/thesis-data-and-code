% 用 fmincon 求解图像中给出的非线性优化问题
% 变量顺序: z = [lambda( nbar x 1 ); omega_minus (m x 1); omega_plus (l x 1)]
%样本单元随机
% 输入:
%   mu_bar   : (m + l) x nbar 矩阵，mū_{i,j}
%   sigma2   : (m + l) x nbar 矩阵，σ̂^2_{i,j}
%   x_ip     : m x 1 向量，x_{i,p}
%   y_rp     : l x 1 向量，y_{r,p}
%   u_alpha  : 标量，u_α
%   delta    : 标量（公式中为 δ），若 δ>0 可等价为 sum(lambda)=1
%
% 输出:
%   lambda, omega_minus, omega_plus, theta 最终目标值, exitflag
%
% 说明:
%   - 约束形式已按 fmincon 的非线性不等式 c(z) <= 0 组织
%   - 所有下界 (lambda, omegas) 均为 0
%   - 为避免目标分母为 0，加入了(可选)非线性不等式：(1/l)*sum(omega_plus./y_rp) < 1 - eps_den
%
function [lambda, omega_minus, omega_plus, rho, exitflag, z, fval] = solve_theta(mu_bar, sigma2, x_ip, y_rp, p, delta)

if nargin < 6
    delta = 1; % 默认使 sum(lambda)=1
end

% 参数尺寸
[nRows, nbar] = size(mu_bar);
m = numel(x_ip);
l = numel(y_rp);
if nRows ~= m + l
    error('mu_bar 的行数应为 m + l');
end
if any(size(sigma2) ~= size(mu_bar))
    error('sigma2 应与 mu_bar 尺寸相同');
end

% 变量数
nvar = nbar + m + l;

% 初始点（可根据问题调整）
lambda0 = ones(nbar,1) / max(1,sum(ones(nbar,1))); % 均匀分布并且和为1
omega_minus0 = max(1e-6, 0.1 * x_ip); % 小正数，避免零带来的奇异
omega_plus0  = max(1e-6, 0.1 * y_rp);
x0 = [lambda0; omega_minus0; omega_plus0];

% 下界（非负）
lb = zeros(nvar,1);
ub = []; % 无上界

% fmincon 选项
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',1000);

% 目标函数句柄
obj = @(z) obj_theta(z, m, l, nbar, x_ip, y_rp);

% 非线性约束句柄
nonlcon = @(z) constraints(z, mu_bar, sigma2, x_ip, y_rp, p, delta, m, l, nbar);

% 调用 fmincon
[z,fval,exitflag,output] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, options);

% 解析结果
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : end);
rho = obj(z);

fprintf('fmincon exitflag = %d, final objective = %.6g\n', exitflag, fval);
end

%% 目标函数实现
function t = obj_theta(z, m, l, nbar, x_ip, y_rp)
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : nbar + m + l);

% 为数值稳定，确保分母不为0或负
den = 1 - (1/l) * sum(omega_plus ./ y_rp);
eps_den = 1e-8;
if den <= eps_den
    % 若分母接近或小于0，返回一个很大的目标值（不可行）
    t = 1e20 + (1/l)*sum(max(0,omega_plus./y_rp)) ;
    return;
end

t = (1 + (1/m) * sum(omega_minus ./ x_ip)) / den;
end

%% 非线性约束实现： c(z) <= 0, ceq(z) == 0
function [c, ceq] = constraints(z, mu_bar, sigma2, x_ip, y_rp, p, delta, m, l, nbar)
lambda = z(1:nbar);
omega_minus = z(nbar+1 : nbar + m);
omega_plus  = z(nbar + m + 1 : nbar + m + l);

% 约束个数： m + l
c = zeros(m + l + l + 1, 1); % 最后一行用于目标分母 < 1 的约束（可选）
% i = 1..m:
for i = 1:m
    left = x_ip(i) - mu_bar(i,:) * lambda + omega_minus(i);
    % 右侧 u_alpha * sqrt( sum_j sigma2(i,j) * lambda_j^2 )
    rhs = norminv(p) * sqrt( sum( sigma2(i,:) .* (lambda'.^2) ) );
    % 要求 left >= rhs  -> rhs - left <= 0
    c(i) = rhs - left;
end

% r = 1..l (行索引在 mu_bar 中为 m + r)
for r = 1:l
    idx = m + r;
    left = - y_rp(r) + mu_bar(idx,:) * lambda + omega_plus(r);
    rhs = norminv(p) * sqrt( sum( sigma2(idx,:) .* (lambda'.^2) ) );
    % 要求 left >= rhs -> rhs - left <= 0
    c(m + r) = rhs - left;
end
for rr = 1:l
    
    %left = - y_rp(rr) + mu_bar(idx,:) * lambda + omega_plus(rr);
    %rhs = u_alpha * sqrt( sum( sigma2(idx,:) .* (lambda'.^2) ) );
    % 要求 left >= rhs -> rhs - left <= 0
    c(m +l+ rr) = omega_plus(rr) - y_rp(rr);
end

% 额外约束： 保证目标分母为正且不接近0
% 1 - (1/l)*sum(omega_plus./y_rp) > eps_den  -> (1/l)*sum(omega_plus./y_rp) - (1 - eps_den) <= 0
eps_den = 1e-6;
c(end) = (1/l) * sum( omega_plus ./ y_rp ) - (1 - eps_den);

% 等式约束： delta * sum(lambda) = delta  -> sum(lambda) - 1 = 0 (当 delta>0)
if abs(delta) < eps
    % 若 delta = 0，则约束退化为 0 = 0，不添加
    ceq = [];
else
    ceq = sum(lambda) - 1;
end
end