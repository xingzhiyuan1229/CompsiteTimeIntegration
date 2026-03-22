# 公式解读文档 —— 第 2 章：基于矩阵指数有理近似的时步格式综述

本文档对论文 *High-order composite implicit time integration schemes based on rational approximations for elastodynamics*（Song & Zhang, 2024）第 2 章中出现的所有公式，逐一给出库中最底层对应计算代码并进行解读，格式为**代码在前、公式在后**。

> **GitHub 仓库根地址**：`https://github.com/xingzhiyuan1229/CompsiteTimeIntegration`

---

## 2.1 时间积分方法——矩阵指数

---

### 公式 (1)：结构动力学运动方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L18-L23)

```matlab
% src/TimeSolverRExpn.m，第 18–23 行
ndof = size(K,1);

K =  dt*dt*sparse(K);   % 对应 Δt²K
C = dt*sparse(C);       % 对应 ΔtC
M = sparse(M);          % 质量矩阵 M
F =  dt*dt*F;           % 对应 Δt²f
```

**代码解读**：入口函数 `TimeSolverRExpn` 接收物理质量矩阵 \(\mathbf{M}\)、阻尼矩阵 \(\mathbf{C}\)、刚度矩阵 \(\mathbf{K}\) 及外力向量 \(\mathbf{F}\)。为适应后续无量纲时间框架，代码首先将 \(\mathbf{K}\) 乘以 \(\Delta t^{2}\)，\(\mathbf{C}\) 乘以 \(\Delta t\)，\(\mathbf{F}\) 乘以 \(\Delta t^{2}\)。这些矩阵正是公式 (5) 中各项所需的比例形式，也是公式 (1) 在算法中的直接载体。

**对应公式 (1)**：

\[
\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{f}(t) \tag{1}
\]

其中 \(\mathbf{M}\)、\(\mathbf{C}\)、\(\mathbf{K}\) 分别为质量矩阵、阻尼矩阵和刚度矩阵，\(\mathbf{f}(t)\) 为外力向量，\(\mathbf{u}\)、\(\dot{\mathbf{u}}\)、\(\ddot{\mathbf{u}}\) 分别为位移、速度和加速度向量。

---

### 公式 (2a)：初始条件

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L25-L53)

```matlab
% src/TimeSolverRExpn.m，第 25 行与第 50–53 行
% Initial velocity --> normalized with dt
v0 = dt*v0;           % u̇°(0) = Δt · u̇(0)

z = zeros(2*ndof,p+1);

% Initial conditions
z(:,1) = [v0; u0];    % z(s=0) = {u̇°₀; u₀}
```

**代码解读**：物理初速度 `v0`（即 \(\dot{\mathbf{u}}(t=0)\)）乘以 \(\Delta t\) 转化为无量纲速度 \(\dot{\mathbf{u}}^{\circ}(0) = \Delta t \cdot \dot{\mathbf{u}}(0)\)（见公式 (4)）。随后构造初始状态向量 \(\mathbf{z}(s=0) = \{\dot{\mathbf{u}}^{\circ}_0;\, \mathbf{u}_0\}\) 并存入 `z(:,1)`，这直接对应公式 (2a) 中的初始条件。

**对应公式 (2a)**：

\[
\mathbf{u}(t=0) = \mathbf{u}_0,\quad \dot{\mathbf{u}}(t=0) = \dot{\mathbf{u}}_0 \tag{2a}
\]

---

### 公式 (3)：无量纲时间变量

**对应代码文件**：[`src/ForceSeries.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/ForceSeries.m#L18-L24)

```matlab
% src/ForceSeries.m，第 18–24 行
function [ft] = ForceSeries(order,fHist,ns,dt)

pf = order + 1;
np = pf + 0;
s = forceSamplingPoints(np);   % 采样点 s ∈ [0,1]
t = dt*reshape(repmat((0:ns-1), np,1)+repmat(s,1,ns),[],1);
                               % t = Δt·( n-1 + s )，n=1,...,ns
f = reshape(fHist(t),np,[]);
T = transMtxPointsToPoly(s, pf);
ft=f'*T;
```

另见 [`src/forceSamplingPoints.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/forceSamplingPoints.m#L18-L22)：

```matlab
% src/forceSamplingPoints.m，第 18–22 行
function [s] = forceSamplingPoints(np)

xi = lglnodes(np-1);
xi = flip(xi);
% Mapping to the interval [0,1]
scl = 1 - 1.0d-12;
s  = 1/2*(scl*xi + 1);   % s ∈ [0,1]
```

**代码解读**：`ForceSeries.m` 中的变量 `s` 正是公式 (3) 中的无量纲时间 \(s\in[0,1]\)。`forceSamplingPoints` 将 Legendre–Gauss–Lobatto（LGL）节点从区间 \([-1,1]\) 线性映射到 \([0,1]\)，得到在当前时步内的采样位置。进而 `t = dt*(n-1+s)` 精确实现了公式 (3)：\(t(s) = t_{n-1} + s\Delta t\)。

**对应公式 (3)**：

\[
t(s) = t_{n-1} + s\Delta t,\quad 0 \le s \le 1 \tag{3}
\]

---

### 公式 (4)：无量纲时间下的速度与加速度

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L25)，[第 82–83 行](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L82-L83)

```matlab
% src/TimeSolverRExpn.m，第 25 行
v0 = dt*v0;   % 入口归一化：u̇°₀ = Δt · u̇₀

% src/TimeSolverRExpn.m，第 82–83 行（循环结束后恢复物理量）
vel = vel/dt;           % u̇ = u̇°/Δt
acc = acc/(dt*dt);      % ü = ü°/Δt²
```

**代码解读**：输入速度 \(\dot{\mathbf{u}}_0\) 在进入算法前被乘以 \(\Delta t\)，使得内部的状态向量始终存储无量纲量 \(\dot{\mathbf{u}}^{\circ} = \Delta t\,\dot{\mathbf{u}}\)（即公式 (4) 中 \(\dot{\mathbf{u}} = \dot{\mathbf{u}}^{\circ}/\Delta t\) 的变形）。计算完成后，`vel` 除以 \(\Delta t\)、`acc` 除以 \(\Delta t^2\) 以恢复物理速度和加速度，这正对应公式 (4) 的逆变换。

**对应公式 (4)**：

\[
\dot{\mathbf{u}} = \frac{1}{\Delta t}\frac{\mathrm{d}\mathbf{u}}{\mathrm{d}s} = \frac{1}{\Delta t}\dot{\mathbf{u}}^{\circ}, \qquad \ddot{\mathbf{u}} = \frac{1}{\Delta t^2}\frac{\mathrm{d}^2\mathbf{u}}{\mathrm{d}s^2} = \frac{1}{\Delta t^2}\ddot{\mathbf{u}}^{\circ} \tag{4}
\]

---

### 公式 (5)：无量纲时间中的运动方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L20-L23)

```matlab
% src/TimeSolverRExpn.m，第 20–23 行
K =  dt*dt*sparse(K);   % → Δt²K（公式(5)中 Δt²K·u 项）
C = dt*sparse(C);       % → ΔtC （公式(5)中 ΔtC·u̇° 项）
M = sparse(M);          % → M   （公式(5)中 M·ü°° 项）
F =  dt*dt*F;           % → Δt²f（公式(5)右端项）
```

**代码解读**：这四行是公式 (5) 在代码层面的直接体现。将物理矩阵乘以适当的 \(\Delta t\) 幂次后，所有后续计算都在无量纲时间 \(s\) 框架下进行。代码中的 `K`、`C`、`M` 分别对应公式 (5) 中的 \(\Delta t^2\mathbf{K}\)、\(\Delta t\mathbf{C}\)、\(\mathbf{M}\)，而 `F` 对应 \(\Delta t^2\mathbf{f}\)。

**对应公式 (5)**：

\[
\mathbf{M}\ddot{\mathbf{u}}^{\circ\circ} + \Delta t\mathbf{C}\dot{\mathbf{u}}^{\circ} + \Delta t^2\mathbf{K}\mathbf{u} = \Delta t^2\mathbf{f} \tag{5}
\]

---

### 公式 (6)：状态空间向量

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L50-L56)

```matlab
% src/TimeSolverRExpn.m，第 50–56 行
z = zeros(2*ndof,p+1);

% Initial conditions
z(:,1) = [v0; u0];          % z = { u̇° ; u }，即公式(6)

for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
end
```

另见 [`src/SolverPadeAx.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/SolverPadeAx.m#L18-L21)：

```matlab
% src/SolverPadeAx.m，第 18–21 行
function y = SolverPadeAx(dM, K, C, x)

n = length(x)/2;
tmp = C*x(1:n)+K*x(n+1:end);  % 上半部分作用在 u̇° 和 u 上
y = [-dM\tmp ; x(1:n)];        % y = A·x，x 分为上下两块
```

**代码解读**：代码中的 `z`（形状 `2*ndof × (p+1)`）的每一列就是当前时步在不同 \(\mathbf{A}_r\) 幂次下的状态向量。`z(:,1)` 存储 \(\mathbf{z} = \{\dot{\mathbf{u}}^{\circ};\, \mathbf{u}\}\)，正是公式 (6) 所定义的状态向量。`SolverPadeAx` 的参数 `x` 也遵循同样的分块结构（上半 \(\dot{\mathbf{u}}^{\circ}\)，下半 \(\mathbf{u}\)）。

**对应公式 (6)**：

\[
\mathbf{z} = \begin{Bmatrix} \dot{\mathbf{u}}^{\circ} \\ \mathbf{u} \end{Bmatrix} \tag{6}
\]

---

### 公式 (7)：一阶 ODE 系统

**对应代码文件**：[`src/SolverPadeAx.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/SolverPadeAx.m#L18-L21)

```matlab
% src/SolverPadeAx.m，第 18–21 行（完整函数体）
function y = SolverPadeAx(dM, K, C, x)

n = length(x)/2;
tmp = C*x(1:n)+K*x(n+1:end);
y = [-dM\tmp ; x(1:n)];

end
```

**代码解读**：`SolverPadeAx` 计算矩阵 \(\mathbf{A}\)（公式 (8)）与向量 \(\mathbf{x}\) 的乘积 \(\mathbf{y} = \mathbf{A}\mathbf{x}\)，是公式 (7) 右端齐次部分 \(\mathbf{A}\mathbf{z}\) 的最底层实现。函数将输入向量 `x` 分为上半部分 \(\dot{\mathbf{u}}^{\circ}\)（`x(1:n)`）和下半部分 \(\mathbf{u}\)（`x(n+1:end)`），分别按照矩阵 \(\mathbf{A}\) 的分块结构进行计算。该函数被 `TimeSolverRExpn` 反复调用以推进子步。

**对应公式 (7)**：

\[
\dot{\mathbf{z}} \equiv \frac{\mathrm{d}\mathbf{z}}{\mathrm{d}s} = \mathbf{A}\mathbf{z} + \mathbf{F} \tag{7}
\]

---

### 公式 (8)：矩阵 A 的定义

**对应代码文件**：[`src/SolverPadeAx.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/SolverPadeAx.m#L18-L21)

```matlab
% src/SolverPadeAx.m，第 18–21 行
function y = SolverPadeAx(dM, K, C, x)

n = length(x)/2;
tmp = C*x(1:n)+K*x(n+1:end);
%  C 此时已是 ΔtC，K 已是 Δt²K（见 TimeSolverRExpn 第 20–21 行）
%  tmp = ΔtC·x₁ + Δt²K·x₂
y = [-dM\tmp ; x(1:n)];
% dM = decomposition(M)，故 dM\tmp = M⁻¹·tmp
% y(1:n)   = -M⁻¹(ΔtC·x₁ + Δt²K·x₂) = (-ΔtM⁻¹C)·x₁ + (-ΔtM⁻¹·ΔtK)·x₂
%           = A[1,1]·x₁ + A[1,2]·x₂
% y(n+1:2n)= x₁ = I·x₁ + 0·x₂ = A[2,1]·x₁ + A[2,2]·x₂

end
```

在 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L35-L36) 中：

```matlab
% src/TimeSolverRExpn.m，第 35–36 行
dM = decomposition(M);   % 对 M 进行分解，供 dM\tmp 调用
Fb = (dM\F);             % Fb = M⁻¹·(Δt²f) = Δt²M⁻¹f
```

**代码解读**：`SolverPadeAx` 以**矩阵-向量乘法**的形式隐式实现了矩阵 \(\mathbf{A}\)，而无需显式存储这个 \(2n\times2n\) 的大矩阵。代码中 `C` 已被预缩放为 \(\Delta t\mathbf{C}\)，`K` 为 \(\Delta t^2\mathbf{K}\)，因此 `tmp = (ΔtC)·x₁ + (Δt²K)·x₂`，再左乘 \(\mathbf{M}^{-1}\) 得到 \(-\Delta t\mathbf{M}^{-1}\mathbf{C}\cdot x_1 - \Delta t\mathbf{M}^{-1}(\Delta t\mathbf{K})\cdot x_2\)，下半部分直接返回 \(x_1 = \mathbf{I}\cdot x_1\)，完整实现了 \(\mathbf{A}\mathbf{x}\)。

**对应公式 (8)**：

\[
\mathbf{A} = \begin{bmatrix} -\Delta t\mathbf{M}^{-1}\mathbf{C} & -\Delta t\mathbf{M}^{-1}\mathbf{K} \\ \mathbf{I} & \mathbf{0} \end{bmatrix} \tag{8}
\]

---

### 公式 (9)：非齐次项 F

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L23-L44)

```matlab
% src/TimeSolverRExpn.m，第 23、35–36、43–44 行
F =  dt*dt*F;           % → Δt²f（标量力向量缩放）

dM = decomposition(M);
Fb = (dM\F);            % Fb = M⁻¹·(Δt²f) = Δt²M⁻¹f（非齐次项上半）

f = zeros(2*ndof,p);
f(:,1)  = [Fb; zeros(ndof,1)];   % F = { Δt²M⁻¹f ; 0 }，即公式(9)
for ii = 2:p
    f(:,ii) = SolverPadeAx(dM,K,C,f(:,ii-1));  % 计算 A^(ii-1)·F
end
Pf = f*cf';             % 力的贡献积分
```

**代码解读**：`Fb = M⁻¹·(Δt²f)` 是公式 (9) 上半块 \(\Delta t^2\mathbf{M}^{-1}\mathbf{f}\) 的直接实现。将其与零向量拼接 `[Fb; zeros(ndof,1)]` 得到完整的 \(\mathbf{F}\) 向量（形状 \(2n\times1\)），与公式 (9) 的分块结构完全一致。后续循环利用 `SolverPadeAx` 计算 \(\mathbf{A}^k\mathbf{F}\)，用于构造力的积分矩阵 `Pf`（对应公式 (12) 中的 \(\sum \mathbf{B}_k\mathbf{F}_{mn}^{(k)}\) 项）。

**对应公式 (9)**：

\[
\mathbf{F} = \begin{Bmatrix} \Delta t^2\mathbf{M}^{-1}\mathbf{f} \\ \mathbf{0} \end{Bmatrix} \tag{9}
\]

---

### 公式 (10)：矩阵指数通解

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L65-L79)

```matlab
% src/TimeSolverRExpn.m，第 65–79 行（时间推进主循环）
for it = 2:ns
    is = it - 1;

    % 公式(18)：z_n = P_r·z_{n-1} + C_k·F_mn^(k)（经 P_r 变换）
    z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';

    % 公式(85)/(88)/(89)：递推子步求解（近似 Q^{-1} 的作用）
    for ip = p:-1:1
        z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
        z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    end

    % 存储响应
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(ndof+pDOF,1);
    fn = fHist(is*dt);
    acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);
end
```

**代码解读**：公式 (10) 给出了系统的精确解析解，但对大型问题无法直接实现（计算矩阵指数代价过高）。整个 `TimeSolverRExpn` 函数正是以有理近似替代精确矩阵指数的工程实现：`prcoe` 编码了有理近似的分子多项式（在 \(\mathbf{A}_r\) 基下），`Pf` 编码了外力积分，子步循环则近似实现了 \(\mathbf{Q}^{-1}\) 的作用（见公式 (85)）。整体结构与公式 (10) 的解析结构一一对应。

**对应公式 (10)**：

\[
\mathbf{z}(s) = e^{\mathbf{A}s}\mathbf{z}_{n-1} + e^{\mathbf{A}s}\int_0^s e^{-\mathbf{A}\tau}\mathbf{F}(\tau)\,\mathrm{d}\tau \tag{10}
\]

---

### 公式 (11)：外力的多项式展开

**对应代码文件**：[`src/ForceSeries.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/ForceSeries.m#L18-L24)，[`src/transMtxPointsToPoly.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/transMtxPointsToPoly.m#L18-L25)

```matlab
% src/ForceSeries.m，第 18–24 行（完整函数体）
function [ft] = ForceSeries(order,fHist,ns,dt)

pf = order + 1;         % 多项式项数 = p_t + 1
np = pf + 0;
s = forceSamplingPoints(np);   % LGL 采样点 s ∈ [0,1]
t = dt*reshape(repmat((0:ns-1), np,1)+repmat(s,1,ns),[],1);
f = reshape(fHist(t),np,[]);   % 在采样时刻求值
T = transMtxPointsToPoly(s, pf);  % 点到多项式系数的变换矩阵
ft=f'*T;                       % ft(is,:) = [F_mn^(0), F_mn^(1), ..., F_mn^(p_t)]

end
```

```matlab
% src/transMtxPointsToPoly.m，第 18–25 行（完整函数体）
function [T] = transMtxPointsToPoly(s, nf)

s  = reshape(s,[],1);
if length(s) >= nf
    a = (s-0.5).^(0:nf-1);   % Vandermonde 矩阵：(s-0.5)^k，k=0,...,nf-1
    T = a/(a'*a);             % T = (Aᵀ·A)⁻¹·Aᵀ（最小二乘拟合矩阵）
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ' ...
        num2str(nf)]);
end
```

**代码解读**：`transMtxPointsToPoly` 构造了以 \((s-0.5)^k\) 为基底的 Vandermonde 矩阵 `a`，然后计算最小二乘拟合矩阵 `T = (aᵀa)⁻¹aᵀ`。当采样点数 ≥ 多项式阶数时，拟合为精确插值。`ForceSeries` 调用该矩阵将每个时步的外力采样值 `f` 转化为多项式展开系数 `ft`，即公式 (11) 中的 \(\mathbf{F}_{mn}^{(k)}\)（\(k=0,1,\ldots,p_t\)）。

**对应公式 (11)**：

\[
\mathbf{F}_n(s) = \sum_{k=0}^{p_t} \mathbf{F}_{mn}^{(k)} (s-0.5)^k \tag{11}
\]

---

### 公式 (12)：时步推进方程（矩阵指数版本）

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L43-L48)，[第 65–79 行](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L65-L79)

```matlab
% src/TimeSolverRExpn.m，第 43–48 行（预计算力的积分贡献）
f = zeros(2*ndof,p);
f(:,1)  = [Fb; zeros(ndof,1)];    % F 向量（公式(9)）
for ii = 2:p
    f(:,ii) = SolverPadeAx(dM,K,C,f(:,ii-1));  % A^(ii-1)·F
end
Pf = f*cf';    % Pf = [B₀F, B₁F, ..., B_{p_t}F] × cf（力积分矩阵）

% src/TimeSolverRExpn.m，第 68 行（每步时间推进核心）
z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
% 第一项：P_r·z_{n-1}（对应 e^A·z_{n-1}）
% 第二项：Σ C_k·F_mn^(k)（对应 Σ B_k·F_mn^(k)，经 Q 缩放）
```

**代码解读**：`Pf*ft(is,:)'` 计算 \(\sum_{k=0}^{p_t} \mathbf{C}_k \mathbf{F}_{mn}^{(k)}\)（公式 (18)/(12) 的力项），`(z(:,1:p+1))*prcoe'` 计算 \(\mathbf{P}_r \mathbf{z}_{n-1}\)（矩阵指数作用于前一步状态，对应 \(e^\mathbf{A}\mathbf{z}_{n-1}\)）。两者之和 `z(:,p+1)` 是经 \(\mathbf{Q}=\mathbf{A}_r^M\) 加权后的右端项，再经子步递推（公式 (85)）得到最终的 \(\mathbf{z}_n\)。

**对应公式 (12)**：

\[
\mathbf{z}_n = e^{\mathbf{A}}\mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{B}_k \mathbf{F}_{mn}^{(k)} \tag{12}
\]

---

### 公式 (13)：B_k 的递推计算

**对应代码文件**：[`src/TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeIntgCoeffForce.m#L18-L26)

```matlab
% src/TimeIntgCoeffForce.m，第 18–26 行（完整函数体）
function [C] = TimeIntgCoeffForce(p,q)

M = length(q) - 1;
tmp = p - q;
C(1,:) = tmp(2:end);  % C_0 = A^{-1}(P-Q)，公式(20)

for k = 1:M
    tmp = ((-1/2)^k)*(p-((-1)^k)*q); % (-0.5)^k·(P-(-1)^k·Q)
    tmp(1:M) = tmp(1:M) + k*C(k,:);  % + k·C_{k-1}
    C(k+1,:) =  tmp(2:end);           % C_k = A^{-1}(·)
end
```

**代码解读**：循环体实现了公式 (19)：`tmp = (-0.5)^k * (P - (-1)^k * Q) + k * C_{k-1}`，然后 `tmp(2:end)` 相当于对多项式左移一位（除以 \(\lambda\) 即乘以 \(\mathbf{A}^{-1}\)）。由于代码中 \(\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k\)，因此 `C(k+1,:)` 存储的是 \(\mathbf{C}_k\)（以多项式系数表示，最低次项已被 \(\mathbf{A}^{-1}\) 消除）。这是公式 (13) 的计算实现，其中用 \(\mathbf{C}_k\) 替换 \(\mathbf{B}_k\) 以避免显式计算矩阵指数。

**对应公式 (13)**：

\[
\mathbf{B}_k = \mathbf{A}^{-1}\!\left(k\mathbf{B}_{k-1} + (-0.5)^k \bigl(e^{\mathbf{A}} - (-1)^k\mathbf{I}\bigr)\right),\quad k=0,1,2,\ldots,p_t \tag{13}
\]

---

### 公式 (14)：B_0 的计算

**对应代码文件**：[`src/TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeIntgCoeffForce.m#L18-L21)

```matlab
% src/TimeIntgCoeffForce.m，第 18–21 行
function [C] = TimeIntgCoeffForce(p,q)

M = length(q) - 1;
tmp = p - q;                  % 多项式 P - Q 的系数
C(1,:) = tmp(2:end);          % C_0 = A^{-1}(P-Q)，公式(20)/(14)
```

**代码解读**：`tmp = p - q` 计算多项式 \(\mathbf{P} - \mathbf{Q}\) 的系数向量，`tmp(2:end)` 丢弃常数项并将系数整体左移，等价于乘以 \(\mathbf{A}^{-1}\)（因为对多项式 \(f(\mathbf{A}) = a_1\mathbf{A} + a_2\mathbf{A}^2 + \cdots\)，有 \(\mathbf{A}^{-1}f(\mathbf{A}) = a_1\mathbf{I} + a_2\mathbf{A} + \cdots\)，在系数数组中体现为右移一位并丢弃常数项）。结果 `C(1,:)` 即为 \(\mathbf{C}_0 = \mathbf{Q}\mathbf{B}_0 = \mathbf{A}^{-1}(\mathbf{P}-\mathbf{Q})\)，是公式 (14) 替换矩阵指数 \(e^\mathbf{A}\) 为有理近似 \(\mathbf{R}=\mathbf{P}/\mathbf{Q}\) 后的结果（见公式 (20)）。

**对应公式 (14)**：

\[
\mathbf{B}_0 = e^{\mathbf{A}}\int_0^1 e^{-\mathbf{A}\tau}\,\mathrm{d}\tau = \mathbf{A}^{-1}(e^{\mathbf{A}} - \mathbf{I}) \tag{14}
\]

---

## 2.2 矩阵指数的有理近似

---

### 公式 (15)：有理近似 R = P/Q

**对应代码文件**：[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L18-L35)

```matlab
% src/InitSchemeRExpn.m，第 18–35 行（完整函数体）
function [r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho )

if contains(scheme,'MP1')
    r =  MP1SchemeRoot(p);    % (M+1)-格式：求最优根 r
else
    r = MschemeRoot(p, rho);  % M-格式：由 ρ_∞ 确定根 r
end

pcoe = pCoefficients(p, r);  % 计算分子多项式 P 的系数 p_i

disp(['r = ', num2str(r)]);
disp(['Selected rhoInfty = ', num2str(abs(pcoe(end)), '%.5f')]);

qrcoe = [zeros(1,p) 1];      % Q 在 A_r 基下：A_r^M（仅最高次项为1）
qcoe = shiftPolycoe(qrcoe,r); % 将 Q 从 A_r 基转换回 A 基

cf = TimeIntgCoeffForce(pcoe,qcoe);  % 计算力积分系数 C_k（公式(19)(20)）

prcoe = shiftPolycoe(pcoe,r);  % 将 P 从 A 基转换到 A_r 基

end
```

**代码解读**：`InitSchemeRExpn` 是整个有理近似框架的初始化入口，对应公式 (15)。它完成三件事：(1) 由 \(\rho_\infty\) 或最优精度条件确定分母多项式的单重根 \(r\)（公式 (34)/(35)）；(2) 计算分子多项式 \(\mathbf{P}\) 的系数 `pcoe`；(3) 将两个多项式从不同基下互相转换，并计算外力积分系数。其输出 `r`、`prcoe`、`cf` 完整编码了有理近似 \(\mathbf{R} = \mathbf{P}/\mathbf{Q}\)。

**对应公式 (15)**：

\[
e^{\mathbf{A}} \approx \mathbf{R} = \frac{\mathbf{P}}{\mathbf{Q}} \tag{15}
\]

---

### 公式 (16a/b)：多项式 P 与 Q 的定义

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24)，[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L29-L30)

```matlab
% src/pCoefficients.m，第 18–24 行（完整函数体，计算分子多项式 P 系数）
function pcoe = pCoefficients(M, r)

pcoe = zeros(1,M+1);
for ii = 0:M
    j = 0:ii;
    % 公式(40)：p_i(r) = Σⱼ (-1)^j · M!/(M-j)!/j!/(i-j)! · r^(M-j)
    p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
    pcoe(ii+1) = p*(r.^(M-j))';   % p_i = p(r)，按公式(40)
end
end
```

```matlab
% src/InitSchemeRExpn.m，第 29–30 行（分母多项式 Q 在 A_r 基下）
qrcoe = [zeros(1,p) 1];      % Q = A_r^M = (rI-A)^M，系数 [0,...,0,1]
qcoe = shiftPolycoe(qrcoe,r); % 转换至 A 的幂次基表示
```

另见 [`src/shiftPolycoe.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/shiftPolycoe.m#L18-L24) 用于多项式基变换：

```matlab
% src/shiftPolycoe.m，第 18–24 行（完整函数体）
function [prcoe] = shiftPolycoe(pcoe,r)

p = size(pcoe,2) - 1;
zc = zeros(size(pcoe,1),1);
prcoe = pcoe;
for ii = p:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
                  - [zc prcoe(:,ii+1:end)];
end

end
```

**代码解读**：`pCoefficients.m` 直接实现公式 (40)，利用二项式展开和 Taylor 级数系数对公式逐项求和，输出 \(p_i\)（\(i=0,1,\ldots,M\)）即公式 (16a) 中的分子多项式 \(\mathbf{P}\) 的系数。分母多项式 \(\mathbf{Q} = (r\mathbf{I}-\mathbf{A})^M\) 在 \(\mathbf{A}_r\) 基下为 \(\mathbf{A}_r^M\)（系数向量 `[0,...,0,1]`），由 `shiftPolycoe` 转换回 \(\mathbf{A}\) 的幂次基（公式 (16b)）。

**对应公式 (16a)**：

\[
\mathbf{P} = \sum_{i=0}^{N_p} p_i \mathbf{A}^i = p_0\mathbf{I} + p_1\mathbf{A} + \cdots + p_{N_p}\mathbf{A}^{N_p} \tag{16a}
\]

**对应公式 (16b)**：

\[
\mathbf{Q} = \sum_{i=0}^{N_q} q_i \mathbf{A}^i = q_0\mathbf{I} + q_1\mathbf{A} + \cdots + q_{N_q}\mathbf{A}^{N_q} \tag{16b}
\]

---

### 公式 (17)：稳定性条件 N_q ≥ N_p

**对应代码文件**：[`src/MschemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MschemeRoot.m#L18-L24)

```matlab
% src/MschemeRoot.m，第 18–24 行（完整函数体）
function [r] = MschemeRoot(M, rhoInfty)

RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;  % 各阶次 M 对应的右端项符号
ir  = [1,  2,  2,  2,  3,  3];           % 对应从小到大排序后的根序号

j = 0:M;
% 公式(42)的多项式系数：p_M(r) = Σⱼ (-1)^j·M!/j!/(M-j)!² · r^(M-j) = ±ρ_∞
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);  % p_M(r) - (±ρ_∞) = 0

rs = sort( roots(pMcoe) );    % 求所有根并升序排列
r  = rs(ir(M));               % 选取满足 |R(λ)| ≤ 1 且相期误差最小的根
end
```

**代码解读**：公式 (17) 是稳定隐式格式的充要条件（\(N_q \ge N_p\)）。在本代码框架中，对 M 子步格式，分子次数 \(N_p = M-1\)（M-格式）或 \(N_p = M\)（M+1-格式），分母次数 \(N_q = M\)，故始终满足 \(N_q \ge N_p\)。`MschemeRoot` 通过求解多项式方程（公式 (42)）来确定满足所需 \(\rho_\infty\) 和稳定性条件的根 \(r\)，其选根策略（`RHS` 符号和 `ir` 序号）正是公式 (17) 约束在数值选优中的体现。

**对应公式 (17)**：

\[
N_q \ge N_p \tag{17}
\]

---

### 公式 (18)：有理近似替代矩阵指数后的时步方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L65-L72)

```matlab
% src/TimeSolverRExpn.m，第 65–72 行（时间推进核心）
for it = 2:ns
    is = it - 1;

    % 右端项计算：P_r·z_{n-1} + Σ C_k·F_mn^(k)（公式(82)/(18)右端）
    z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';

    % 递推子步求解：A_r^M · z_n = z^(0)（公式(85)/(86)）
    for ip = p:-1:1
        z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
        z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    end
    ...
end
```

在 [`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L32-L34) 中预计算系数：

```matlab
% src/InitSchemeRExpn.m，第 32–34 行
cf = TimeIntgCoeffForce(pcoe,qcoe);  % 计算 C_k（公式(19)/(20)）
prcoe = shiftPolycoe(pcoe,r);        % P 在 A_r 基下的系数（P_r）
```

**代码解读**：第 68 行 `z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)'` 是公式 (18)（即 \(\mathbf{Q}\mathbf{z}_n = \mathbf{P}\mathbf{z}_{n-1} + \sum\mathbf{C}_k\mathbf{F}_{mn}^{(k)}\)）的直接计算：
- `(z(:,1:p+1))*prcoe'` 计算 \(\mathbf{P}_r\mathbf{z}_{n-1}\)（P 以 A_r 基表示，与 z 的各阶 A_r 幂次分量点积求和）；
- `Pf*ft(is,:)'` 计算 \(\sum_{k=0}^{p_t}\mathbf{C}_k\mathbf{F}_{mn}^{(k)}\)；
- 结果存入 `z(:,p+1)` 是右端项，后续子步循环对其作用 \(\mathbf{A}_r^{-M}\) 得到 \(\mathbf{z}_n\)。

**对应公式 (18)**：

\[
\mathbf{Q}\mathbf{z}_n = \mathbf{P}\mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{C}_k \mathbf{F}_{mn}^{(k)} \tag{18}
\]

---

### 公式 (19)：C_k 的递推计算

**对应代码文件**：[`src/TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeIntgCoeffForce.m#L22-L26)

```matlab
% src/TimeIntgCoeffForce.m，第 22–26 行（递推循环）
for k = 1:M
    tmp = ((-1/2)^k)*(p-((-1)^k)*q);
    % 计算 (-0.5)^k·(P - (-1)^k·Q)，对应公式(19)中的 (-0.5)^k(P-(-1)^k Q)
    tmp(1:M) = tmp(1:M) + k*C(k,:);
    % 加上 k·C_{k-1}（其中 C(k,:) 即 C_{k-1} 的系数）
    C(k+1,:) =  tmp(2:end);
    % A^{-1} 作用：去掉常数项并左移（因 C_k = A^{-1}(...)）
end
```

**代码解读**：此循环精确实现了公式 (19)。代码中 `p` 和 `q` 为多项式 \(\mathbf{P}\) 和 \(\mathbf{Q}\) 的系数向量，`C(k,:)` 存储 \(\mathbf{C}_{k-1}\) 的多项式系数。`tmp(2:end)` 截取掉零次项后等价于 \(\mathbf{A}^{-1}\) 作用（在多项式系数空间中，乘以 \(\mathbf{A}^{-1}\) 等效于将系数数组右移一位，即 `tmp(2:end)` 使得系数从 \(A^1\) 开始）。

**对应公式 (19)**：

\[
\mathbf{C}_k = \mathbf{Q}\mathbf{B}_k = \mathbf{A}^{-1}\!\left(k\mathbf{C}_{k-1} + (-0.5)^k \bigl(\mathbf{P} - (-1)^k\mathbf{Q}\bigr)\right),\quad k=1,2,\ldots,p_t \tag{19}
\]

---

### 公式 (20)：C_0 的计算

**对应代码文件**：[`src/TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeIntgCoeffForce.m#L18-L21)

```matlab
% src/TimeIntgCoeffForce.m，第 18–21 行
function [C] = TimeIntgCoeffForce(p,q)

M = length(q) - 1;
tmp = p - q;          % P - Q（多项式系数相减）
C(1,:) = tmp(2:end);  % C_0 = A^{-1}(P-Q)
%  tmp(1)（常数项，对应 A^0）= p_0 - q_0 = 0（因 p_0 = q_0 保证 e^0 = I）
%  所以 tmp(2:end) 截取后得到 A^{-1}(P-Q) 的系数，无需担心 A^{-1} 奇异问题
```

**代码解读**：公式 (20) 要求 \(\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P}-\mathbf{Q})\)。代码中 `tmp = p - q` 是分子分母多项式之差的系数，其常数项（对应 \((\mathbf{P}-\mathbf{Q})|_{\mathbf{A}=0} = p_0 - q_0 = 0\)，因为归一化条件 \(p_0 = q_0\)）恰好为零，故 `tmp(2:end)` 相当于将多项式除以 \(\mathbf{A}\)（即乘以 \(\mathbf{A}^{-1}\)），结果正是 \(\mathbf{C}_0\) 的多项式系数。

**对应公式 (20)**：

\[
\mathbf{C}_0 = \mathbf{Q}\mathbf{B}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q}) \tag{20}
\]

---

*（第 2.3 节：数值耗散与频散性质中的公式 (21)–(33) 见 [公式解读_Part2.md](./公式解读_Part2.md)）*
