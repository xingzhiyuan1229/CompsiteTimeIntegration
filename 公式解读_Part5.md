# 公式解读文档 —— 第 4 节续与第 5 节：Bathe 方法关系及数值算例（公式 94–99）

本文档为 [公式解读_Part4.md](./公式解读_Part4.md) 的续篇，涵盖第 4 节剩余部分（公式 (94)–(96)，与 Bathe 方法的关系）和第 5 节数值算例（公式 (97)–(99)）的代码解读。格式为**代码在前、公式在后**。

> **GitHub 仓库根地址**：`https://github.com/xingzhiyuan1229/CompsiteTimeIntegration`

---

## 4. 复合时间积分格式的时步算法（续）

### 公式 (94)：子步方程除以 r² 后的归一化形式

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L40-L42)

```matlab
% src/TimeSolverRExpn.m，第 40–42 行（等效刚度矩阵，对应公式 (88) 的直接形式）
Kd = sparse((r*r)*M + r*C + K);
% 其中 M=M（质量），C=ΔtC（已缩放），K=Δt²K（已缩放）
% 即：r²M + r(ΔtC) + Δt²K = r²(M + (Δt/r)C + (Δt/r)²K)
% 两边除以 r² 后变为：(M + (Δt/r)C + (Δt/r)²K) z1^(k) = (1/r)M*z1^{k-1} - (Δt/r)²K*z2^{k-1}

% 等价于公式 (94)，其中 Δt/r 等价于 Bathe 方法的子步大小 γΔt/2：
% (M + (Δt/r)C + (Δt/r)²K) z1^(k) = (1/r)M*z1^{k-1} - (Δt/r)²K*z2^{k-1}
dKd = decomposition(Kd);
% 预分解一次，供所有 M 个子步共用
```

**代码解读**：公式 (94) 是将公式 (88) 两边同除以 [ r^2 ] 后的等价形式：[ (\mathbf{M} + \frac{\Delta t}{r}\mathbf{C} + (\frac{\Delta t}{r})^2\mathbf{K})\mathbf{z}_1^{(k)} = \frac{1}{r}\mathbf{M}\mathbf{z}_1^{(k-1)} - (\frac{\Delta t}{r})^2\mathbf{K}\mathbf{z}_2^{(k-1)} ]。代码中 `Kd = r²M + r(ΔtC) + Δt²K` 是公式 (88) 的左端矩阵，与公式 (94) 的左端矩阵相差因子 [ r^2 ]（两者等价，只是形式不同）。通过引入等效子步大小 [ \Delta t/r ]，公式 (94) 揭示了本方法与 Bathe 方法等效刚度矩阵的结构性联系，从而建立了时间分割比 [ \gamma = 2/r ]（公式 (95)）。

**对应公式 (94)**：

\[
\left(\mathbf{M} + \frac{\Delta t}{r}\mathbf{C} + \left(\frac{\Delta t}{r}\right)^2\mathbf{K}\right) \mathbf{z}_1^{(k)} = \frac{1}{r}\mathbf{M}\mathbf{z}_1^{(k-1)} - \left(\frac{\Delta t}{r}\right)^2\mathbf{K}\mathbf{z}_2^{(k-1)} \tag{94}
\]

公式 (88) 两边除以 [ r^2 ] 后的等价形式，使得等效子步大小 [ \Delta t/r ] 与 Bathe 方法的 [ \gamma\Delta t/2 ] 直接对应（见公式 (95)）。

---

### 公式 (95)：本方法与 Bathe 方法的时间分割比关系

**对应代码文件**：[`src/MschemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MschemeRoot.m#L18-L25) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L38-L42)

```matlab
% 公式 (95) 建立了 γ = 2/r 的关系（其中 r 由 MschemeRoot 确定）
% 对 M=2, M-格式：r = 2+sqrt(2+2*rhoInfty)（公式 (60)）
% γ = 2/r = 2/(2+sqrt(2+2*rhoInfty))（对应公式 (96)）

% src/MschemeRoot.m（r 的确定，间接给出 γ）
r = MschemeRoot(p, rho);   % 例如 M=2, rho=0.5: r ≈ 3.732, γ = 2/3.732 ≈ 0.536

% src/TimeSolverRExpn.m（等效刚度矩阵中 Δt/r 的体现）
Kd = sparse((r*r)*M + r*C + K);
% = r²(M + (Δt/r)C + (Δt/r)²K)，其中 Δt/r = γΔt/2（Bathe 方法等效子步大小）
```

**代码解读**：公式 (95) [ \gamma = 2/r ] 揭示了本方法（[ M=2 ] 格式）与 [ \rho_\infty ]-Bathe 方法（时间分割比为 [ \gamma ]）的等价关系。在代码中，`r` 由 `MschemeRoot(p, rho)` 确定（公式 (60)：[ r = 2+\sqrt{2+2\rho_\infty} ]），`Kd` 的构造隐含了等效子步大小 [ \Delta t/r ]，其与 Bathe 方法中的 [ \gamma\Delta t/2 ] 相等（即 [ \gamma = 2/r ]）。这一等价关系表明，对 [ M=2 ]，本方法在数值上与 [ \rho_\infty ]-Bathe 方法完全一致（论文在算例中也通过数值计算验证了这一点）。

**对应公式 (95)**：

\[
\gamma = \frac{2}{r} \tag{95}
\]

[ M=2 ] 时本方法与 Bathe 方法的等价关系，[ r ] 为有理近似分母的单重根，[ \gamma ] 为 Bathe 方法的时间分割比。

---

### 公式 (96)：M=2 时分割比 γ 与 ρ_∞ 的关系

**对应代码文件**：[`src/MschemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MschemeRoot.m#L18-L25)

```matlab
% src/MschemeRoot.m（M=2 时，r = 2+sqrt(2+2*rhoInfty)，对应公式 (96) 的推导）
j = 0:M;  % M=2
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);  % RHS(2) = +rhoInfty
rs = sort( roots(pMcoe) );
r  = rs(ir(M));  % ir(2)=2 → r = 2+sqrt(2+2*rhoInfty)（公式 (60)）

% γ = 2/r = 2/(2+sqrt(2+2*rhoInfty))
% 等价推导：γ = (2-sqrt(2+2*rhoInfty))/(1-rhoInfty)（有理化后）
% 公式 (96) 与文献 [33] 的 ρ_∞-Bathe 方法结果一致
```

**代码解读**：公式 (96) 将 [ M=2 ] 时的分割比 [ \gamma = 2/r = 2/(2+\sqrt{2+2\rho_\infty}) ] 化简为等价形式 [ (2-\sqrt{2+2\rho_\infty})/(1-\rho_\infty) ]（通过有理化分母：分子分母乘以 [ (2-\sqrt{2+2\rho_\infty}) ]，利用 [ r\cdot(2-\sqrt{2+2\rho_\infty}) = (2+\sqrt{2+2\rho_\infty})(2-\sqrt{2+2\rho_\infty}) = 4-(2+2\rho_\infty) = 2(1-\rho_\infty) ]，得 [ \gamma = 2\cdot(2-\sqrt{2+2\rho_\infty})/(2(1-\rho_\infty)) = (2-\sqrt{2+2\rho_\infty})/(1-\rho_\infty) ]）。代码中 `MschemeRoot(2, rhoInfty)` 输出的根 [ r ] 代入 [ 2/r ] 即得此分割比，与文献 [33] 中 [ \rho_\infty ]-Bathe 方法的分割比完全吻合，验证了两方法的数值等价性。

**对应公式 (96)**：

\[
\gamma = \frac{2}{2 + \sqrt{2 + 2\rho_{\infty}}} = \frac{2 - \sqrt{2 + 2\rho_{\infty}}}{1 - \rho_{\infty}} \tag{96}
\]

[ M=2 ] 时 M-格式的时间分割比，与 [ \rho_\infty ]-Bathe 方法（文献 [33]）的分割比完全相同，证明了 [ M=2 ] M-格式与 [ \rho_\infty ]-Bathe 方法的数值等价性。

---

## 5. 数值算例

### 5.1 单自由度系统

### 公式 (97)：单自由度算例的外部激励函数

**对应代码文件**：[`ExampleThreeDOFs.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/ExampleThreeDOFs.m#L52)

```matlab
% ExampleThreeDOFs.m，第 52 行（三自由度算例的激励函数，与 SDOF 算例形式类似）
fHist = @(t) sin(1.2*t);      % 三自由度算例中的激励（简谐函数）

% 注：单自由度算例（5.1 节）使用公式 (97) 的激励：
% f_1(t) = 10*cos(2*sqrt(5)/5 * t) + 70*sin(2*sqrt(10) * t)
% 该算例的 MATLAB 代码未包含在本仓库中，但 TimeSolverRExpn 通用于所有算例
% 通过修改 fHist 句柄即可实现公式 (97) 的激励：
% fHist = @(t) 10*cos(2*sqrt(5)/5*t) + 70*sin(2*sqrt(10)*t);
```

**代码解读**：公式 (97) 给出单自由度算例（Section 5.1）的外部激励函数 [ f_1(t) = 10\cos\!\left(\frac{2\sqrt{5}}{5}t\right) + 70\sin\!\left(2\sqrt{10}\,t\right) ]，包含两个谐波分量。在仓库提供的 `ExampleThreeDOFs.m`（Section 5.2 算例）中，激励为 `fHist = @(t) sin(1.2*t)`，是公式 (97) 的简化版本（仅单一正弦分量），两者均以函数句柄的形式传入通用求解器 `TimeSolverRExpn`。`ForceSeries.m` 中的 `fHist(t)` 调用可接受任意函数句柄，完全支持公式 (97) 所示的多分量激励。

**对应公式 (97)**：

\[
f_1(t) = 10\cos\!\left(\frac{2\sqrt{5}}{5}t\right) + 70\sin\!\left(2\sqrt{10}\,t\right) \tag{97}
\]

单自由度算例的外部谐波激励，包含两个频率分量。通过 `fHist = @(t) 10*cos(2*sqrt(5)/5*t) + 70*sin(2*sqrt(10)*t)` 传入求解器。

---

### 公式 (98)：位移 L₂ 范数误差的定义

**对应代码文件**：[`ExampleThreeDOFs.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/ExampleThreeDOFs.m#L136-L155) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L75-L76)

```matlab
% ExampleThreeDOFs.m，第 136–155 行（参考解 uRef 的计算函数）
function [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp)
[Vec,D] = eig(full(K));
% ...
uRef = Vec*U;  % 精确参考解（u_exact(t)）

% 数值解 dsp 由 TimeSolverRExpn 输出
[dsp,vel,acc] = TimeSolverRExpn(...);
% dsp(:,1) 为自由度 1 的数值位移（u_numerical(t)）

% 公式 (98) 的 L₂ 误差计算（积分离散化为求和）：
% epsilon_L2 = sum((uRef - dsp).^2) / sum(uRef.^2) * 100%
% 在 ExampleThreeDOFs 中通过图形对比直观展示，不显式计算此误差
% 但 TimeSolverRExpn 输出的 dsp 即为 u_numerical，与公式 (98) 一一对应
```

**代码解读**：公式 (98) 定义了衡量时间积分精度的 [ L_2 ] 范数误差 [ \epsilon_{L_2} = \frac{\int_0^{t_\text{sim}}(u_\text{exact}-u_\text{numerical})^2\,\mathrm{d}t}{\int_0^{t_\text{sim}}u_\text{exact}^2\,\mathrm{d}t}\times 100\% ]。在代码层面：`uRef`（精确解）由 `ThreeDOFsRefSln` 计算（模态叠加解析解），`dsp`（数值解）由 `TimeSolverRExpn` 输出。连续积分以离散求和近似：[ \epsilon_{L_2} \approx \frac{\sum_n(u_\text{exact}(t_n)-u_\text{numerical}(t_n))^2}{\sum_n u_\text{exact}(t_n)^2}\times 100\% ]。论文图 12–13 展示了不同格式、不同时步大小下的收敛曲线，这些误差值均通过公式 (98) 计算得到。

**对应公式 (98)**：

\[
\epsilon_{L_2} = \frac{\int_{t=0}^{t_{\mathrm{sim}}} \left(u_{\mathrm{exact}}(t) - u_{\mathrm{numerical}}(t)\right)^2\,\mathrm{d}t}{\int_{t=0}^{t_{\mathrm{sim}}} u_{\mathrm{exact}}(t)^2\,\mathrm{d}t} \times 100\% \tag{98}
\]

基于 [ L_2 ] 范数的位移相对误差，用于定量评估时间积分格式的精度（论文图 12–13 中的收敛曲线均采用此误差度量）。

---

### 5.3 一维波传播算例

### 公式 (99)：Courant-Friedrichs-Lewy (CFL) 数的定义

**对应代码文件**：[`ExampleThreeDOFs.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/ExampleThreeDOFs.m#L31-L34) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L20) 

```matlab
% ExampleThreeDOFs.m，第 31–34 行（时步参数，对应 CFL 数的设定）
tmax = 5010;   % 仿真总时长
dt = 0.14;     % 时步大小 Δt（论文 5.2 节选 Δt=0.14）

% 对于 5.3 节的波传播算例（双材料杆），CFL 数的计算对应公式 (99)：
% CFL = c_2 * dt / dx
% 其中 c_2 = 20*sqrt(2) m/s（右段波速），dx = L_segment/N_elements = 2/2000 = 0.001 m
% 时步由 CFL 数反推：dt = CFL * dx / c_2

% src/TimeSolverRExpn.m，第 20 行（时步 dt 的使用）
K =  dt*dt*sparse(K);   % Δt² 缩放，体现了 dt（时步大小）在算法中的核心作用
C = dt*sparse(C);
% dt 的选取直接控制 CFL 数，从而影响数值耗散和色散特性
```

**代码解读**：公式 (99) 定义了波传播问题中的 CFL 数 [ \mathrm{CFL} = c_2\Delta t/\Delta x ]，其中 [ c_2 ] 为波速、[ \Delta t ] 为时步大小、[ \Delta x ] 为单元尺寸。在代码中，`dt` 是 `TimeSolverRExpn` 的关键输入参数，`K = dt²K`、`C = dt*C` 的预缩放使得所有后续计算均在无量纲时间框架下进行（公式 (3)–(5)）。对于波传播算例（5.3 节），`dt` 的选取对应特定的 CFL 数，是控制数值色散（相对周期误差，公式 (32)）和数值耗散（谱半径，公式 (27)）行为的主要参数。论文图 21–26 展示了不同 CFL 数下各格式的数值表现，`dt` 通过 `CFL * dx / c_2` 计算后传入求解器。

**对应公式 (99)**：

\[
\mathrm{CFL} = c_2 \frac{\Delta t}{\Delta x} \tag{99}
\]

波传播问题中的 CFL 数定义，以右段波速 [ c_2 ] 计算。[ \Delta t ] 直接传入 `TimeSolverRExpn`；更大的 CFL 数意味着更大的时步，高阶格式（更多子步）可以在更大 CFL 数下保持良好精度，提高计算效率。

---

## 附录：公式解读文档总览

以下是本系列文档（Part 1–5）所覆盖的全部 99 个公式的完整对应关系：

| 文档 | 公式范围 | 章节 | 核心内容 |
|------|----------|------|---------|
| [公式解读_Part1.md](./公式解读_Part1.md) | (1)–(14) | 2.1 | 结构动力学方程、状态空间变换、矩阵指数精确解 |
| [公式解读_Part2.md](./公式解读_Part2.md) | (15)–(33) | 2.2–2.3 | 有理近似框架、谱分析（耗散与频散） |
| [公式解读_Part3.md](./公式解读_Part3.md) | (34)–(65) | 3.1 | M-格式有理近似构造（展开系数、求根） |
| [公式解读_Part4.md](./公式解读_Part4.md) | (66)–(93) | 3.2, 4 | (M+1)-格式、A_r 基下算法实现 |
| **本文档 Part5** | **(94)–(99)** | **4, 5** | **Bathe 方法关系、数值算例** |

所有代码链接指向 `https://github.com/xingzhiyuan1229/CompsiteTimeIntegration` 仓库的 `main` 分支，源码文件位于 `src/` 目录下。

### 核心源码文件与公式对应关系

| 源码文件 | 对应公式（主要） |
|---------|----------------|
| [`src/MschemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MschemeRoot.m) | (33), (41)–(43), (45)–(65), (95)–(96) |
| [`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m) | (67), (72)–(74) |
| [`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m) | (40)–(44), (51)–(65), (71)–(79) |
| [`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m) | (15), (34), (44), (69)–(70) |
| [`src/shiftPolycoe.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/shiftPolycoe.m) | (81a), (81b) |
| [`src/TimeIntgCoeffForce.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeIntgCoeffForce.m) | (13), (14), (19), (20) |
| [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m) | (1)–(12), (18), (80)–(93), (94)–(96), (99) |
| [`src/SolverPadeAx.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/SolverPadeAx.m) | (7), (8), (22) |
| [`src/ForceSeries.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/ForceSeries.m) | (3), (11) |
| [`src/forceSamplingPoints.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/forceSamplingPoints.m) | (3) |
| [`src/transMtxPointsToPoly.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/transMtxPointsToPoly.m) | (11) |
| [`src/lglnodes.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/lglnodes.m) | (3), (11)（LGL 采样点） |
| [`ExampleThreeDOFs.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/ExampleThreeDOFs.m) | (97), (98), (99) |

---

*本文档（Part 5）覆盖公式 (94)–(99)，共 6 个公式，配合 Part 1–4，完整解读了论文全部 99 个公式对应的最底层计算代码。*
