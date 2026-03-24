# 代码详解文档（第一部分）

## 按主函数运行顺序逐行解析

本文档以 `ExampleThreeDOFs.m` 主脚本为入口，按照实际运行顺序对库中所有代码逐行进行详细解释。涉及计算的代码均给出对应数学公式，代码在前，公式在后，可分段对应。

---

## 一、主脚本 `ExampleThreeDOFs.m`

### 1.1 初始化与路径设置

```matlab
% 第 19 行
clear; close all;

% 第 20 行
dbstop if error;

% 第 22 行
addpath(".\src\");
```

- **第 19 行**：`clear` 清除工作区所有变量，`close all` 关闭所有图形窗口，确保每次运行从干净状态开始。
- **第 20 行**：`dbstop if error` 是调试指令，当程序发生运行时错误时自动进入调试断点，便于排查问题，不参与实际计算。
- **第 22 行**：将 `src` 子目录加入 MATLAB 搜索路径，使得后续调用 `TimeSolverRExpn`、`InitSchemeRExpn` 等函数时 MATLAB 能正确找到对应文件。

---

### 1.2 复合时间积分参数设置

```matlab
% 第 25–27 行
scheme   = 'M';   % 格式类型：'M' 或 'MP1'
nSubStep = 4;     % 子步数 M，M格式支持1~6，(M+1)格式支持1,2,3,5
rhoInfty = 0;     % 高频极限谱半径 ρ_∞，仅 M 格式使用
```

- `scheme = 'M'`：选择 M 阶格式（M sub-steps，M-th order accuracy）。若为 `'MP1'` 则选 (M+1) 阶格式。
- `nSubStep = 4`：本例选 \(M = 4\) 个子步，即四阶精度格式。
- `rhoInfty = 0`：用户指定高频极限谱半径

\[
\rho_\infty = \lim_{\omega \to \infty} |R(\lambda)| = 0
\]

\(\rho_\infty = 0\) 表示高频分量将被完全耗散（最大数值阻尼），\(\rho_\infty = 1\) 则无数值阻尼。

---

### 1.3 三自由度系统参数

```matlab
% 第 30–38 行
tmax = 5010;
dt   = 0.14;
k1 = 10^7;
k2 = 1;
m1 = 0;
m2 = 1;
m3 = 1;
```

本例模拟一个由弹簧-质量组成的串联系统：地基（固定）→ 弹簧 \(k_1\) → 质量 \(m_2\) → 弹簧 \(k_2\) → 质量 \(m_3\)。\(m_1 = 0\) 表示地基节点无质量。两个频率尺度相差极大（\(\sqrt{k_1/m_2} \approx 3162\,\text{rad/s}\)，\(\sqrt{k_2/m_3} = 1\,\text{rad/s}\)），用于检验格式处理刚性（stiff）问题的能力。

- `tmax = 5010`：总仿真时长（秒）。
- `dt = 0.14`：时间步长 \(\Delta t = 0.14\,\text{s}\)。相对于高频模态（\(T_1 \approx 0.002\,\text{s}\)），\(\Delta t / T_1 \approx 70\)，远大于 1，需要格式具有强数值耗散来压制该高频成分。

---

### 1.4 有限元矩阵组装

```matlab
% 第 40–47 行
np = 2;
M  = [m2 0; 0, m3];            % 质量矩阵
C  = zeros(2);                  % 阻尼矩阵（无阻尼）
K  = [k1+k2, -k2; -k2, k2];   % 刚度矩阵
F0 = [1; 0]*k1;                % 力幅向量
u0 = [0; 0];                    % 初始位移
v0 = u0;                        % 初始速度
BC_Accl = [];
```

这四个矩阵对应如下运动方程（两自由度系统）：

\[
\mathbf{M}\ddot{\mathbf{u}}(t) + \mathbf{C}\dot{\mathbf{u}}(t) + \mathbf{K}\mathbf{u}(t) = \mathbf{F}_0 f(t)
\]

其中

\[
\mathbf{M} = \begin{bmatrix} m_2 & 0 \\ 0 & m_3 \end{bmatrix} = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}, \quad
\mathbf{C} = \mathbf{0}, \quad
\mathbf{K} = \begin{bmatrix} k_1+k_2 & -k_2 \\ -k_2 & k_2 \end{bmatrix} = \begin{bmatrix} 10^7+1 & -1 \\ -1 & 1 \end{bmatrix}
\]

\[
\mathbf{F}_0 = \begin{bmatrix} k_1 \\ 0 \end{bmatrix} = \begin{bmatrix} 10^7 \\ 0 \end{bmatrix}, \quad
\mathbf{u}(0) = \dot{\mathbf{u}}(0) = \mathbf{0}
\]

`BC_Accl = []` 声明加速度边界条件为空（本例不使用），无实际计算意义。

---

### 1.5 输出自由度与激励函数

```matlab
% 第 50–52 行
pDOF  = [1;2];
fHist = @(t) sin(1.2*t);
```

- `pDOF = [1;2]`：指定需要保存位移/速度/加速度历史的自由度编号，这里保存全部两个自由度。
- 激励信号为正弦函数，激励频率 \(\Omega = 1.2\,\text{rad/s}\)，通过 \(\mathbf{F}(t) = \mathbf{F}_0 \cdot f(t)\) 作用于系统：

\[
f(t) = \sin(1.2\, t)
\]

---

### 1.6 时间序列与时步数

```matlab
% 第 54–55 行
ns = floor(tmax/dt)+1;
tp = (0:0.02:tmax);
```

- `ns`：数值求解所需的总时步数，等于 \(\lfloor t_{\max}/\Delta t \rfloor + 1\)（包含初始时刻 \(t=0\)）。
- `tp`：用于计算参考解的细密时间轴（时间间隔 0.02 s），用于生成高精度参考曲线以供对比。

---

### 1.7 参考解与反力计算

```matlab
% 第 57–58 行
[uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp);
R1Ref = k1*(fHist(tp)-uRef(1,:));
```

- 调用局部函数 `ThreeDOFsRefSln`（见第二节）计算精确解析参考解。
- `R1Ref`：弹簧 \(k_1\) 的反力（弹力）。由胡克定律，弹簧两端的相对位移为 \(f(t_p) - u_1(t_p)\)（地基位移 = \(f(t)\)，质量 2 位移 = \(u_1\)），故弹力为：

\[
R_1(t) = k_1 \bigl(f(t) - u_1(t)\bigr)
\]

---

### 1.8 数值求解与绘图

```matlab
% 第 60–63 行
[dsp,vel,acc] = TimeSolverRExpn(scheme,nSubStep, rhoInfty, ...
                             ns,dt,fHist, K,M,C,F0,u0,v0,pDOF);
tn = (0:ns-1)*dt;
R1 = k1*((fHist(tn))'-dsp(:,1));
```

- 调用主求解器 `TimeSolverRExpn`（见第三节）得到数值位移 `dsp`、速度 `vel`、加速度 `acc`（均为 `ns × length(pDOF)` 矩阵）。
- `tn`：数值解对应的时间序列 \(t_n = (n-1)\Delta t\)，\(n = 1, \ldots, \text{ns}\)。
- `R1`：数值解对应的弹簧反力 \(R_1(t_n) = k_1(f(t_n) - u_1(t_n))\)。

```matlab
% 第 65–133 行：7 个 figure 绘图块
```

这七组绘图分别对比：质量 2 位移、质量 3 位移、质量 2 速度、质量 3 速度、质量 2 加速度、质量 3 加速度和弹簧反力的参考解（红色实线）与数值解（蓝色虚线）。每个图的结构相同，此处不逐一展开。

---

## 二、局部函数 `ThreeDOFsRefSln`（第 136–155 行）

```matlab
function [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp)
```

该函数利用模态叠加法求解无阻尼受迫振动的精确解析解。

---

### 2.1 模态分析

```matlab
% 第 137–140 行
[Vec,D] = eig(full(K));
Mg = Vec'*M*Vec;
Kg = Vec'*K*Vec;
Fg = Vec'*F0;
```

- **第 137 行**：`eig(full(K))` 对刚度矩阵 \(\mathbf{K}\) 进行特征值分解，返回特征向量矩阵 \(\mathbf{V}\) 和特征值对角矩阵 \(\mathbf{D}\)，满足：

\[
\mathbf{K}\,\mathbf{V} = \mathbf{M}\,\mathbf{V}\,\mathbf{\Lambda}, \quad \mathbf{\Lambda} = \operatorname{diag}(\omega_1^2, \omega_2^2)
\]

- **第 138–140 行**：计算广义（模态）坐标下的质量矩阵 \(\mathbf{M}_g\)、刚度矩阵 \(\mathbf{K}_g\)、广义力向量 \(\mathbf{F}_g\)：

\[
\mathbf{M}_g = \mathbf{V}^T\mathbf{M}\mathbf{V}, \quad
\mathbf{K}_g = \mathbf{V}^T\mathbf{K}\mathbf{V}, \quad
\mathbf{F}_g = \mathbf{V}^T\mathbf{F}_0
\]

对于质量归一化的特征向量，\(\mathbf{M}_g = \mathbf{I}\)，\(\mathbf{K}_g = \mathbf{\Lambda}\)（对角）。由于代码使用 `eig(K)` 而非广义特征值分解，\(\mathbf{M}_g\) 不一定是单位矩阵，但由于 \(\mathbf{M} = \mathbf{I}\)（本例 \(m_2=m_3=1\)），实际上 \(\mathbf{M}_g = \mathbf{I}\)，\(\mathbf{K}_g = \mathbf{D}\)。

---

### 2.2 固有频率

```matlab
% 第 141–142 行
o1 = sqrt(D(1,1));
o2 = sqrt(D(2,2));
```

提取两个固有频率：

\[
\omega_1 = \sqrt{\lambda_1}, \quad \omega_2 = \sqrt{\lambda_2}
\]

本例 \(\omega_1 \approx 1\,\text{rad/s}\)（低频，由 \(k_2\) 主导），\(\omega_2 \approx 3162\,\text{rad/s}\)（高频，由 \(k_1\) 主导）。

---

### 2.3 精确位移响应函数

```matlab
% 第 143–148 行
Uex  = @(o,t)  1/(o*o-1.2^2)*(sin(1.2*t) - 1.2/o*sin(o*t));
Uref = @(o,t)  1/(o*o-1.2^2)*(sin(1.2*t)                  );
Vex  = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t) - cos(o*t));
Vref = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t)            );
Aex  = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t) + o*sin(o*t));
Aref = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t)              );
```

对单自由度无阻尼系统 \(\ddot{q} + \omega^2 q = F_g \sin(\Omega t)\)（初始条件为零），精确位移解为：

\[
q(t) = \frac{F_g}{\omega^2 - \Omega^2}\!\left(\sin\Omega t - \frac{\Omega}{\omega}\sin\omega t\right)
\]

其中 \(\Omega = 1.2\,\text{rad/s}\)。  

- **`Uex`**：包含自由振动项 \(-\frac{\Omega}{\omega}\sin\omega t\)，用于第一模态（\(\omega_1 \approx 1\)），自由振动不可忽略。  
- **`Uref`**：仅保留受迫振动项 \(\sin\Omega t\)，用于第二模态（\(\omega_2 \gg \Omega\)），由于 \(\omega_2\Delta t \gg 1\)，高频自由振动已被数值耗散压制，参考解中的高频项 \(-(\Omega/\omega_2)\sin\omega_2 t\) 幅值极小，可忽略不计。

速度是位移对时间的导数，加速度是速度对时间的导数：

\[
\dot{q}(t) = \frac{\Omega}{\omega^2-\Omega^2}\!\left(\cos\Omega t - \cos\omega t\right), \quad
\ddot{q}(t) = \frac{\Omega}{\omega^2-\Omega^2}\!\left(-\Omega\sin\Omega t + \omega\sin\omega t\right)
\]

---

### 2.4 模态叠加还原物理坐标

```matlab
% 第 149–154 行
U    = [Fg(1)*Uex(o1,tp);  Fg(2)*Uref(o2,tp)];
uRef = Vec*U;
V    = [Fg(1)*Vex(o1,tp);  Fg(2)*Vref(o2,tp)];
vRef = Vec*V;
A    = [Fg(1)*Aex(o1,tp);  Fg(2)*Aref(o2,tp)];
aRef = Vec*A;
```

先计算广义坐标响应，再通过特征向量矩阵 \(\mathbf{V}\) 还原到物理坐标：

\[
\mathbf{u}(t) = \mathbf{V}\,\mathbf{q}(t), \quad
\mathbf{q}(t) = \begin{bmatrix} F_{g,1}\cdot q_1(t;\,\omega_1) \\ F_{g,2}\cdot q_2(t;\,\omega_2) \end{bmatrix}
\]

---

## 三、主求解器 `TimeSolverRExpn.m`

```matlab
function [dsp,vel,acc] = TimeSolverRExpn(scheme,p,rho,ns,dt,fHist,K,M,C,F,u0,v0,pDOF)
```

**输入参数说明：**

| 参数 | 含义 |
|---|---|
| `scheme` | 格式类型字符串（`'M'` 或 `'MP1'`） |
| `p` | 子步数 \(M\) |
| `rho` | 用户指定的 \(\rho_\infty\)（仅 M 格式） |
| `ns` | 总时步数 |
| `dt` | 时间步长 \(\Delta t\) |
| `fHist` | 激励函数句柄 \(f(t)\) |
| `K`, `M`, `C` | 刚度、质量、阻尼矩阵 |
| `F` | 力幅向量 \(\mathbf{F}_0\) |
| `u0`, `v0` | 初始位移和速度 |
| `pDOF` | 需要输出的自由度编号 |

---

### 3.1 矩阵归一化（无量纲时间变换）

```matlab
% 第 18–25 行
ndof = size(K,1);

K  = dt*dt*sparse(K);
C  = dt*sparse(C);
M  = sparse(M);
F  = dt*dt*F;
% Initial velocity --> normalized with dt
v0 = dt*v0;
```

**关键设计**：引入无量纲时间变量 \(s\)，使得时步内 \(s \in [0,1]\)，对应物理时间 \(t = t_{n-1} + s\Delta t\)。用 \(\circ\) 表示对 \(s\) 的导数：

\[
\dot{\mathbf{u}} = \frac{1}{\Delta t}\mathbf{u}^\circ, \quad
\ddot{\mathbf{u}} = \frac{1}{\Delta t^2}\mathbf{u}^{\circ\circ}
\]

代入物理运动方程，得无量纲运动方程：

\[
\mathbf{M}\,\mathbf{u}^{\circ\circ} + \Delta t\,\mathbf{C}\,\mathbf{u}^\circ + \Delta t^2\mathbf{K}\,\mathbf{u} = \Delta t^2\,\mathbf{F}_0 f(s)
\]

因此代码中将矩阵预乘 \(\Delta t\) 的相应幂次：
- `K ← Δt²K`（刚度矩阵），
- `C ← ΔtC`（阻尼矩阵），
- `F ← Δt²F₀`（力幅），
- `v0 ← Δt·v0`（初速度归一化为 \(\dot{\mathbf{u}}^\circ(0) = \Delta t \dot{\mathbf{u}}(0)\)）。

`ndof = 2`（本例），`sparse(·)` 将矩阵转为稀疏格式以节省内存并加速后续线性求解。

---

### 3.2 外力多项式展开预计算

```matlab
% 第 27 行
ft = ForceSeries(p, fHist, ns-1, dt);
```

调用 `ForceSeries`（见第四节）将激励函数 \(f(t)\) 在每个时步上展开为 \(p+1\) 项的多项式（以时步中点为展开中心）：

\[
f_n(s) \approx \sum_{k=0}^{p} c_{n,k}\,(s - 0.5)^k, \quad s \in [0,1]
\]

返回矩阵 `ft`，形状为 \((n_s-1) \times (p+1)\)，其中 `ft(is,k+1)` 是第 `is` 个时步的第 \(k\) 阶多项式系数 \(c_{n,k}\)。

---

### 3.3 输出数组初始化

```matlab
% 第 31–33 行
dsp = zeros(ns, length(pDOF));
vel = dsp;
acc = dsp;
```

预分配 `ns × 2`（本例）的零矩阵，分别存储所有时步的位移、速度、加速度响应历史。

---

### 3.4 质量矩阵分解与 M⁻¹F 预计算

```matlab
% 第 35–36 行
dM = decomposition(M);
Fb = (dM\F);
```

- **第 35 行**：`decomposition(M)` 对质量矩阵 \(\mathbf{M}\) 做 LU 分解（或 Cholesky 分解），返回分解对象 `dM`，可高效重复求解线性方程组。
- **第 36 行**：`Fb = M⁻¹F = M⁻¹(Δt²F₀)`，即归一化后的广义力向量。对应无量纲运动方程右端：

\[
\mathbf{F}_b = \Delta t^2 \mathbf{M}^{-1}\mathbf{F}_0
\]

后续加速度的计算以及状态向量中的力贡献都需要用到 `Fb`。

---

### 3.5 格式初始化（调用 InitSchemeRExpn）

```matlab
% 第 38 行
[r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho);
```

调用 `InitSchemeRExpn`（见第五节）返回三个关键量：

| 输出 | 含义 |
|---|---|
| `r` | 有理近似分母多项式的单重根 \(r\)（实数标量） |
| `prcoe` | 分子多项式 \(\mathbf{P}\) 在 \(\mathbf{A}_r\) 基下的系数，形状 \(1 \times (M+1)\) |
| `cf` | 外力积分系数矩阵，形状 \((M+1) \times M\)，对应 \(\mathbf{C}_k\) 的多项式表示 |

---

### 3.6 有效刚度矩阵组装与分解

```matlab
% 第 40–41 行
Kd  = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);
```

**核心矩阵**：有效刚度矩阵（每个子步需求解的线性方程组系数矩阵）：

\[
\mathbf{K}_d = r^2\mathbf{M} + r\Delta t\,\mathbf{C} + \Delta t^2\mathbf{K}
\]

（代码中 `M`、`C`、`K` 已含 \(\Delta t\) 因子，故实际代码为 `r²M + rC + K`，等价于上式）

\(\mathbf{K}_d\) 是所有子步共享的同一有效刚度矩阵——这正是本方法相对于一般 Padé 展开（每个子步刚度矩阵不同）的关键优势：只需分解一次，所有子步重复利用。`decomposition(Kd)` 完成 LU 分解，后续每步仅需三角回代。

---

### 3.7 力贡献矩阵预计算

```matlab
% 第 43–48 行
f        = zeros(2*ndof, p);
f(:,1)   = [Fb; zeros(ndof,1)];
for ii = 2:p
    f(:,ii) = SolverPadeAx(dM,K,C,f(:,ii-1));
end
Pf = f*cf';
```

**目标**：预计算矩阵 \(\mathbf{P}_f\)，使得每个时步的力贡献可以简单地用 `Pf * ft(is,:)'` 计算。

- 状态向量的结构为 \(\mathbf{z} = \{\mathbf{u}^\circ;\, \mathbf{u}\}\)（上半为无量纲速度，下半为位移）。
- 力项的状态向量形式为：

\[
\tilde{\mathbf{F}} = \begin{bmatrix} \Delta t^2\mathbf{M}^{-1}\mathbf{F}_0 \\ \mathbf{0} \end{bmatrix} = \begin{bmatrix} \mathbf{F}_b \\ \mathbf{0} \end{bmatrix}
\]

- **第 44 行**：`f(:,1) = [Fb; zeros(ndof,1)]` 初始化 \(\tilde{\mathbf{F}}\)。
- **第 45–47 行**：递推计算 \(\mathbf{A}^{k-1}\tilde{\mathbf{F}}\)，其中 `SolverPadeAx(dM,K,C,x)` 计算 \(\mathbf{A}\cdot\mathbf{x}\)：

\[
\mathbf{f}_{:,k} = \mathbf{A}^{k-1}\tilde{\mathbf{F}}, \quad k = 1, 2, \ldots, M
\]

- **第 48 行**：`Pf = f * cf'`，矩阵乘法将力状态向量与积分系数矩阵 `cf` 组合，得到 \((2n) \times M\) 的矩阵 \(\mathbf{P}_f\)，满足：

\[
\mathbf{P}_f \cdot \tilde{\mathbf{c}}_n = \sum_{k=0}^{M-1} \mathbf{C}_k\,\mathbf{F}_{mn}^{(k)}\cdot f_0
\]

（其中 \(\tilde{\mathbf{c}}_n\) 是当前时步的力多项式系数向量 `ft(is,:)'`）

---

### 3.8 状态向量初始化（A_r 幂次展开）

```matlab
% 第 50–56 行
z = zeros(2*ndof, p+1);

% Initial conditions
z(:,1) = [v0; u0];
for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
end
```

引入矩阵 \(\mathbf{A}_r = r\mathbf{I} - \mathbf{A}\)（其中 \(\mathbf{A}\) 为状态矩阵），则：

\[
\mathbf{A}_r\,\mathbf{z} = (r\mathbf{I} - \mathbf{A})\,\mathbf{z} = r\,\mathbf{z} - \mathbf{A}\,\mathbf{z}
\]

代码第 55 行正是在计算 \(\mathbf{A}_r\,\mathbf{z}_{:,ii} = r\,\mathbf{z}_{:,ii} - \mathbf{A}\,\mathbf{z}_{:,ii}\)，
其中 `SolverPadeAx(dM,K,C,z(:,ii))` 返回 \(\mathbf{A}\,\mathbf{z}_{:,ii}\)。

因此初始化循环完成：

\[
\mathbf{z}_{:,1} = \mathbf{z}_{n-1}, \quad
\mathbf{z}_{:,k+1} = \mathbf{A}_r^k\,\mathbf{z}_{n-1}, \quad k = 1, \ldots, M
\]

这 \(M+1\) 列共同编码了初始状态在 \(\mathbf{A}_r\) 幂次基下的分量，供后续时步推进使用。

---

### 3.9 第一时步的响应存储

```matlab
% 第 59–63 行
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
fn = fHist(0);
acc(it,:) = r*z(pDOF,1) - z(pDOF,2) + fn*Fb(pDOF);
```

- 第 60–61 行存储 \(t=0\) 的初始位移和速度（已归一化速度 `v0` 此时存入 `vel`，后续第 82 行再除以 `dt` 还原为物理速度）。
- 初始加速度由运动方程和格式的递推关系得到。`r*z(pDOF,1) - z(pDOF,2)` 利用 `z(:,1) = z₀` 和 `z(:,2) = A_r·z₀` 的关系，实现：

\[
\mathbf{r}\,\mathbf{z}_{:,1} - \mathbf{z}_{:,2} = r\,\mathbf{z}_0 - \mathbf{A}_r\,\mathbf{z}_0 = \mathbf{A}\,\mathbf{z}_0
\]

取上半部分 \(\mathbf{A}\mathbf{z}_0\) 的分量（即 \(\dot{\mathbf{u}}^{\circ\circ}(0)\)），加上力贡献 \(f(0)\cdot\mathbf{F}_b\)，即可得到初始加速度（无量纲形式，第 83 行再除以 `dt²` 还原）：

\[
\mathbf{u}^{\circ\circ}(0) = [\mathbf{A}\,\mathbf{z}_0]_{\text{上半}} + f(0)\,\mathbf{F}_b
\]

---

### 3.10 时间推进主循环

```matlab
% 第 65–80 行
for it = 2:ns
    is = it - 1;

    % （a）计算右端项（公式82的右端）
    z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';

    % （b）递推子步求解（公式85）
    for ip = p:-1:1
        z(1:ndof,ip)       = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
        z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    end

    % （c）存储响应
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(ndof+pDOF,1);
    fn = fHist(is*dt);
    acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);
end
```

这是整个求解器最核心的部分，详细解析如下：

#### （a）右端项计算（第 68 行）

时步方程（将有理近似代入矩阵指数通解后乘以 \(\mathbf{Q}\)）为：

\[
\mathbf{Q}\,\mathbf{z}_n = \mathbf{P}\,\mathbf{z}_{n-1} + \sum_{k=0}^{M} \mathbf{C}_k\,\mathbf{F}_{mn}^{(k)}
\]

将 \(\mathbf{P}\) 和 \(\mathbf{Q} = \mathbf{A}_r^M\) 转换到 \(\mathbf{A}_r\) 基后，方程变为：

\[
\mathbf{A}_r^M\,\mathbf{z}_n = \mathbf{P}_r\,\mathbf{z}_{n-1} + \sum_{k=0}^{M-1} \mathbf{C}_k\,\mathbf{F}_{mn}^{(k)}
\]

定义右端项为 \(\mathbf{z}^{(0)}\)，存于 `z(:,p+1)`：

\[
\mathbf{z}^{(0)} = \underbrace{(z(:,1:p+1) \cdot \text{prcoe}')}_{\mathbf{P}_r\,\mathbf{z}_{n-1}} + \underbrace{\mathbf{P}_f \cdot \tilde{\mathbf{c}}_n}_{\text{力项}}
\]

其中 `prcoe` 的每列是 \(\mathbf{P}_r\) 在 \(\mathbf{A}_r\) 幂次基下的系数，`z(:,1:p+1)` 的各列是 \(\mathbf{A}_r^0, \ldots, \mathbf{A}_r^M\) 作用在 \(\mathbf{z}_{n-1}\) 上的结果（在第 54–56 行和每步递推后更新）。

#### （b）递推子步求解（第 69–72 行）

将 \(\mathbf{A}_r^M\,\mathbf{z}_n = \mathbf{z}^{(0)}\) 分解为 \(M\) 个递推子步，每步求解：

\[
\mathbf{A}_r\,\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)}, \quad k = 1, 2, \ldots, M
\]

展开 \(\mathbf{A}_r = r\mathbf{I} - \mathbf{A}\) 并代入 \(\mathbf{A}\) 的分块结构（见公式 (8)），每个子步分两行完成：

**（b1）第 70 行——求解上半部分（无量纲速度）：**

\[
(r\mathbf{I} - \mathbf{A})\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)}
\]

展开矩阵方程的第一行（\(\text{ndof}\) 个方程）：

\[
(r^2\mathbf{M} + r\Delta t\mathbf{C} + \Delta t^2\mathbf{K})\,\mathbf{z}_1^{(k)} = r\mathbf{M}\,\mathbf{z}_1^{(k-1)} - \Delta t^2\mathbf{K}\,\mathbf{z}_2^{(k-1)}
\]

即 \(\mathbf{K}_d\,\mathbf{z}_1^{(k)} = r\mathbf{M}\,\mathbf{z}_1^{(k-1)} - \mathbf{K}\,\mathbf{z}_2^{(k-1)}\)，其中 \(\mathbf{K}_d = r^2\mathbf{M}+r\mathbf{C}+\mathbf{K}\)。

**（b2）第 71 行——由代数关系更新下半部分（位移）：**

展开矩阵方程第二行（由 \(\mathbf{A}\) 的结构给出 \(\dot{\mathbf{z}}_2 = \mathbf{z}_1\)）：

\[
\mathbf{z}_2^{(k)} = \frac{1}{r}\!\left(\mathbf{z}_1^{(k)} + \mathbf{z}_2^{(k-1)}\right)
\]

循环从 `ip = p` 到 `ip = 1`（逆序），即从 \(k=M\) 到 \(k=1\) 递推，最终 `z(:,1)` 存储 \(\mathbf{z}_n\)（当前时步结束状态）。

#### （c）响应提取与加速度（第 75–78 行）

- `vel(it,:) = z(pDOF,1)`：提取上半部分（无量纲速度 \(\dot{\mathbf{u}}^\circ\)）在输出自由度上的值。
- `dsp(it,:) = z(ndof+pDOF,1)`：提取下半部分（位移 \(\mathbf{u}\)）。
- 加速度由与初始时刻相同的方法计算（第 78 行）。

注意这里存储的 `vel` 是 \(\Delta t \cdot \dot{\mathbf{u}}\)（无量纲速度），循环结束后统一处理。

---

### 3.11 量纲还原

```matlab
% 第 82–83 行
vel = vel/dt;
acc = acc/(dt*dt);
```

将无量纲速度和加速度还原为物理量：

\[
\dot{\mathbf{u}} = \frac{\dot{\mathbf{u}}^\circ}{\Delta t}, \quad
\ddot{\mathbf{u}} = \frac{\ddot{\mathbf{u}}^{\circ\circ}}{\Delta t^2}
\]

位移 `dsp` 已是物理位移，无需换算。

---

## 四、外力展开链

### 4.1 `ForceSeries.m`（第 1–26 行）

```matlab
function [ft] = ForceSeries(order, fHist, ns, dt)
```

**作用**：将激励函数 \(f(t)\) 在每个时步上做多项式拟合，返回每个时步的多项式系数矩阵。

```matlab
% 第 18–24 行
pf = order + 1;
np = pf + 0;
s  = forceSamplingPoints(np);
t  = dt*reshape(repmat((0:ns-1), np,1)+repmat(s,1,ns),[],1);
f  = reshape(fHist(t),np,[]);
T  = transMtxPointsToPoly(s, pf);
ft = f'*T;
```

- **第 18 行**：`pf = order + 1`：多项式项数（包含常数项）。`order = p`（子步数），所以展开为 \(p+1\) 项，截断阶次与分子多项式 \(\mathbf{P}\) 的次数一致。
- **第 19 行**：`np = pf`，采样点个数等于多项式项数（恰好插值，无超定）。
- **第 20 行**：调用 `forceSamplingPoints(np)` 获取 \(np\) 个在 \([0,1]\) 内的 LGL 采样点 \(s_j\)（见第 4.2 节）。
- **第 21 行**：将采样点映射到物理时间。对每个时步 `is = 0, 1, ..., ns-1`，在每个采样点 \(s_j\) 处的物理时刻为：

\[
t_{n,j} = \Delta t\,(\text{is} + s_j), \quad j = 1, \ldots, n_p
\]

`repmat((0:ns-1), np,1)` 构造 \(n_p \times n_s\) 矩阵（每列重复时步序号），`repmat(s,1,ns)` 构造同形的采样偏移矩阵，相加后乘以 `dt`，再 `reshape` 展平为长向量。

- **第 22 行**：在所有时刻一次性求值 `fHist(t)`，然后 `reshape(·, np, [])` 重排为 \(n_p \times n_s\) 矩阵，`f(:,is+1)` 是第 `is` 个时步在各采样点的函数值。
- **第 23 行**：调用 `transMtxPointsToPoly(s, pf)` 计算从采样点值到多项式系数的变换矩阵 \(\mathbf{T}\)（见第 4.4 节）。
- **第 24 行**：`ft = f' * T`。`f'` 是 \(n_s \times n_p\)，`T` 是 \(n_p \times (p+1)\)，结果 `ft` 是 \(n_s \times (p+1)\)，`ft(is,:)` 即第 `is` 时步激励函数的多项式展开系数向量：

\[
\bigl[c_{n,0},\, c_{n,1},\, \ldots,\, c_{n,p}\bigr]
\]

满足近似：\(f_n(s) \approx \sum_{k=0}^{p} c_{n,k}(s-0.5)^k\)

---

### 4.2 `forceSamplingPoints.m`（第 1–24 行）

```matlab
function [s] = forceSamplingPoints(np)
```

**作用**：生成 \(n_p\) 个分布在 \([0,1]\) 的 Legendre-Gauss-Lobatto（LGL）采样点。

```matlab
% 第 18–22 行
xi  = lglnodes(np-1);
xi  = flip(xi);
scl = 1 - 1.0d-12;
s   = 1/2*(scl*xi + 1);
```

- **第 18 行**：`lglnodes(np-1)` 返回 \(n_p\) 个标准区间 \([-1,1]\) 上的 LGL 节点 \(\xi_j\)（从大到小排列）。
- **第 19 行**：`flip(xi)` 将排列颠倒为从小到大（升序）。
- **第 21 行**：`scl = 1 - 10^{-12}` 是微小缩放因子，将端点从 \(\pm1\) 缩进一个机器精度量级，避免恰好在时步边界取值时出现数值奇异。
- **第 22 行**：线性映射 \([-1,1] \to [0,1]\)：

\[
s_j = \frac{1}{2}(\text{scl}\cdot\xi_j + 1) \approx \frac{\xi_j + 1}{2}
\]

LGL 节点是用于 Gauss-Lobatto 积分的最优节点，它包含端点 \(\pm1\)（对应 \(s=0\) 和 \(s=1\)），在高阶多项式插值时具有最小化 Runge 震荡的性质，是高阶谱方法的标准选择。

---

### 4.3 `lglnodes.m`（第 1–50 行）

```matlab
function [x,w,P] = lglnodes(N)
```

**作用**：用 Newton-Raphson 迭代计算 \(N+1\) 个 LGL 节点、积分权重和 Legendre Vandermonde 矩阵。

LGL 节点是下述方程的零点：

\[
(1 - x^2) P'_N(x) = 0
\]

即 \(x = \pm1\) 加上 \(N-1\) 个 \(P'_N(x)\) 的内部零点（共 \(N+1\) 个节点）。

```matlab
% 第 20 行
N1 = N+1;
```

`N1 = N+1` 为节点总数。

```matlab
% 第 23 行
x = cos(pi*(0:N)/N)';
```

以 Chebyshev-Gauss-Lobatto（CGL）节点作为 Newton 迭代的初始猜测：

\[
x_j^{(0)} = \cos\!\left(\frac{j\pi}{N}\right), \quad j = 0, 1, \ldots, N
\]

CGL 节点分布与 LGL 节点相近，收敛快。

```matlab
% 第 26 行
P = zeros(N1,N1);
```

\((N+1)\times(N+1)\) 矩阵，`P(:,k+1)` 存储所有节点处 \(P_k(x)\) 的值（Legendre Vandermonde 矩阵）。

```matlab
% 第 32 行
xold = 2;
```

令初始旧值为 2（超出 \([-1,1]\) 范围），确保 `while` 条件 `max(abs(x-xold)) > eps` 在第一次迭代时成立。

```matlab
% 第 34–46 行
while max(abs(x-xold)) > eps
    xold = x;
    P(:,1) = 1;    P(:,2) = x;
    for k = 2:N
        P(:,k+1) = ((2*k-1)*x.*P(:,k) - (k-1)*P(:,k-1))/k;
    end
    x = xold - (x.*P(:,N1) - P(:,N)) ./ (N1*P(:,N1));
end
```

**Bonnet 递推公式**（第 41 行），用于高效计算各阶 Legendre 多项式：

\[
P_{k+1}(x) = \frac{(2k-1)\,x\,P_k(x) - (k-1)\,P_{k-1}(x)}{k}
\]

**Newton-Raphson 更新**（第 44 行），利用 \((1-x^2)P'_N(x) = N(P_{N-1}(x) - x P_N(x))\) 的恒等式：

\[
x \leftarrow x - \frac{x\,P_N(x) - P_{N-1}(x)}{(N+1)\,P_N(x)}
\]

此公式等价于对方程 \(x\,P_N(x) - P_{N-1}(x) = 0\)（乘以 \(x\) 重写的形式）做 Newton 迭代，迭代至机器精度 `eps` 收敛。

```matlab
% 第 48 行
w = 2./(N*N1*P(:,N1).^2);
```

LGL 积分权重公式：

\[
w_j = \frac{2}{N(N+1)\,[P_N(x_j)]^2}
\]

此权重使得 Gauss-Lobatto 积分公式 \(\int_{-1}^1 f(x)\,\mathrm{d}x \approx \sum_{j=0}^N w_j f(x_j)\) 对所有次数不超过 \(2N-1\) 的多项式精确成立。

---

### 4.4 `transMtxPointsToPoly.m`（第 1–26 行）

```matlab
function [T] = transMtxPointsToPoly(s, nf)
```

**作用**：构造从函数采样值 \(f(s_j)\) 到多项式系数 \([c_0, c_1, \ldots, c_{n_f-1}]\) 的变换矩阵 \(\mathbf{T}\)，使得：

\[
f(s) \approx \sum_{k=0}^{n_f-1} c_k\,(s-0.5)^k
\]

```matlab
% 第 18 行
s = reshape(s,[],1);
```

确保 `s` 为列向量。

```matlab
% 第 19–21 行
if length(s) >= nf
    a = (s-0.5).^(0:nf-1);
    T = a/(a'*a);
```

- **第 20 行**：构造 Vandermonde 矩阵 \(\mathbf{a}\)，形状 \(n_p \times n_f\)：

\[
\mathbf{a}_{j,k} = (s_j - 0.5)^{k-1}, \quad j=1,\ldots,n_p;\; k=1,\ldots,n_f
\]

- **第 21 行**：`T = a / (a'*a)`，即计算：

\[
\mathbf{T} = \mathbf{a}\,(\mathbf{a}^T\mathbf{a})^{-1}
\]

当 \(n_p = n_f\)（恰好插值）时，\(\mathbf{a}\) 为方阵，\(\mathbf{T} = (\mathbf{a}^T)^{-1} = (\mathbf{a}^{-1})^T\)，此时变换矩阵精确。  
当 \(n_p > n_f\)（超定）时，\(\mathbf{T}\) 是最小二乘意义下的右伪逆，满足 \(\hat{\mathbf{c}} = \mathbf{T}^T \mathbf{f}\)（系数 = 伪逆 × 采样值）。  
最终 `ft(is,:) = f(:,is+1)' * T`，即将 \(n_p\) 个采样值映射为 \(n_f\) 个多项式系数（行向量）。

```matlab
% 第 22–26 行
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ', num2str(nf)]);
end
```

当采样点数不足时给出错误提示（欠定情况无法稳健拟合）。

---

*（续见 [代码详解_Part2.md](./代码详解_Part2.md)，包含 `InitSchemeRExpn` 调用链及时步主循环的子函数详解）*
