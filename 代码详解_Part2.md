# 代码详解文档（第二部分）

## InitSchemeRExpn 调用链及核心子函数

本文档为 [代码详解_Part1.md](./代码详解_Part1.md) 的续篇，详解格式初始化调用链（`InitSchemeRExpn` 及其所有子函数）和 `SolverPadeAx`。

---

## 五、格式初始化链

### 5.1 `InitSchemeRExpn.m`（第 1–36 行）

```matlab
function [r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho)
```

**作用**：根据格式类型和参数，完整初始化有理近似格式，返回时间推进所需的全部系数。

有理近似的一般形式为：

\[
e^{\mathbf{A}} \approx \mathbf{R}(\mathbf{A}) = \frac{\mathbf{P}(\mathbf{A})}{\mathbf{Q}(\mathbf{A})} = \frac{p_0\mathbf{I} + p_1\mathbf{A} + \cdots + p_M\mathbf{A}^M}{(r\mathbf{I} - \mathbf{A})^M}
\]

分母 \(\mathbf{Q} = (r\mathbf{I}-\mathbf{A})^M\) 具有单重根 \(r\)，这保证了所有子步的有效刚度矩阵相同。

---

#### 5.1.1 选取根 r（第 18–22 行）

```matlab
if contains(scheme,'MP1')
    r = MP1SchemeRoot(p);
else
    r = MschemeRoot(p, rho);
end
```

- 若格式为 `'MP1'`：调用 `MP1SchemeRoot(p)` 返回预设的根 \(r\)（见第 5.3 节）。
- 否则为 M 格式：调用 `MschemeRoot(p, rho)` 根据用户指定的 \(\rho_\infty\) 数值求解得到 \(r\)（见第 5.2 节）。

---

#### 5.1.2 计算分子多项式系数 pcoe（第 24 行）

```matlab
pcoe = pCoefficients(p, r);
```

调用 `pCoefficients(p, r)` 计算多项式 \(\mathbf{P}(\mathbf{A}) = \sum_{i=0}^{M} p_i\mathbf{A}^i\) 的系数向量 `pcoe`（见第 5.4 节）。

---

#### 5.1.3 输出诊断信息（第 26–27 行）

```matlab
disp(['r = ', num2str(r)]);
disp(['Selected rhoInfty = ', num2str(abs(pcoe(end)), '%.5f')]);
```

- 打印选定的根 \(r\)。
- `pcoe(end) = p_M` 是分子多项式的最高次项系数。由公式 (33)，高频极限谱半径为：

\[
\rho_\infty = \left|\frac{p_M}{q_M}\right| = |p_M|
\]

（因为 \(q_M = 1\)，分母多项式最高次项系数为 1）。`abs(pcoe(end))` 即为实际确认的 \(\rho_\infty\)，打印供用户验证。

---

#### 5.1.4 分母多项式系数（第 29–30 行）

```matlab
qrcoe = [zeros(1,p) 1];
qcoe  = shiftPolycoe(qrcoe, r);
```

- **第 29 行**：分母 \(\mathbf{Q} = \mathbf{A}_r^M = (r\mathbf{I}-\mathbf{A})^M\) 在 \(\mathbf{A}_r\) 基下的系数为 `[0, 0, ..., 0, 1]`（只有最高次 \(\mathbf{A}_r^M\) 项系数为 1），形状 \(1\times(M+1)\)。
- **第 30 行**：调用 `shiftPolycoe(qrcoe, r)` 将 \(\mathbf{Q}\) 从 \(\mathbf{A}_r\) 基转换到 \(\mathbf{A}\) 基，得到 \(\mathbf{Q}(\mathbf{A}) = (r\mathbf{I}-\mathbf{A})^M\) 的系数向量 `qcoe`（见第 5.5 节）。

---

#### 5.1.5 外力积分系数（第 32 行）

```matlab
cf = TimeIntgCoeffForce(pcoe, qcoe);
```

调用 `TimeIntgCoeffForce(pcoe, qcoe)` 计算矩阵 \(\mathbf{C}_k\)（\(k=0,1,\ldots,M\)）的多项式表示，形状 \((M+1)\times M\)（见第 5.6 节）。

---

#### 5.1.6 分子多项式在 A_r 基下的系数（第 34 行）

```matlab
prcoe = shiftPolycoe(pcoe, r);
```

调用 `shiftPolycoe(pcoe, r)` 将分子多项式 \(\mathbf{P}(\mathbf{A})\) 从 \(\mathbf{A}\) 基转换到 \(\mathbf{A}_r\) 基，得到 `prcoe`（见第 5.5 节）。

这样 `prcoe` 使得以下计算成立（主循环第 68 行）：

\[
\mathbf{P}(\mathbf{A})\,\mathbf{z}_{n-1} = \sum_{k=0}^{M} \text{prcoe}(k+1) \cdot \mathbf{A}_r^k\,\mathbf{z}_{n-1}
\]

代码中 `z(:,1:p+1)` 的各列已预先计算了 \(\mathbf{A}_r^0\mathbf{z}_{n-1}, \ldots, \mathbf{A}_r^M\mathbf{z}_{n-1}\)，因此点积 `(z(:,1:p+1)) * prcoe'` 即为 \(\mathbf{P}_r(\mathbf{A}_r)\,\mathbf{z}_{n-1}\)。

---

### 5.2 `MschemeRoot.m`（第 1–26 行）

```matlab
function [r] = MschemeRoot(M, rhoInfty)
```

**作用**：给定子步数 \(M\) 和用户指定的高频极限谱半径 \(\rho_\infty\)，数值求解分母多项式的根 \(r\)。

---

#### 5.2.1 右端项符号表（第 18–19 行）

```matlab
RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
ir  = [1,  2,  2,  2,  3,  3];
```

由理论分析，对 M-格式（M 个子步，M 阶精度），方程 \(p_M(r) = \pm\rho_\infty\) 的右端符号和对应根的排序序号如下查表确定：

| M | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| 右端 | \(+\rho_\infty\) | \(+\rho_\infty\) | \(-\rho_\infty\) | \(+\rho_\infty\) | \(-\rho_\infty\) | \(-\rho_\infty\) |
| 选根序号 | 1 | 2 | 2 | 2 | 3 | 3 |

这些值是通过数值实验（检验 \(|R(\lambda)|\le1\) 且相期误差最小）预先确定的（论文表 1）。

---

#### 5.2.2 构造 p_M(r) 的多项式系数（第 20–22 行）

```matlab
j      = 0:M;
pMcoe  = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);
```

**关键计算**：`pMcoe` 是多项式 \(p_M(r)\) 关于变量 \(r\) 的系数向量（降幂排列），第 \(j+1\) 个元素对应 \(r^{M-j}\) 的系数：

\[
p_M(r) = \sum_{j=0}^{M} (-1)^j \frac{M!}{j!\,(M-j)!^2}\,r^{M-j}
\]

这是由论文公式 (40) 在 \(i = M\) 时代入得到的：

\[
p_M(r) = \sum_{j=0}^{M} (-1)^j \frac{M!}{j!(M-j)!} \cdot \frac{1}{(M-j)!}\,r^{M-j}
= \sum_{j=0}^{M} (-1)^j \frac{M!}{j!\,(M-j)!^2}\,r^{M-j}
\]

- **第 22 行**：将方程右端项从系数向量的最后一项（常数项，\(r^0\) 的系数，对应 \(j=M\)）中减去 \(\text{RHS}(M)\)，等价于求解方程：

\[
p_M(r) - (\pm\rho_\infty) = 0
\]

---

#### 5.2.3 求解并选根（第 23–24 行）

```matlab
rs = sort(roots(pMcoe));
r  = rs(ir(M));
```

- **第 23 行**：`roots(pMcoe)` 求多项式方程 \(p_M(r) = \pm\rho_\infty\) 的所有根（复数向量），`sort` 按升序排列。
- **第 24 行**：`rs(ir(M))` 按查表选取第 `ir(M)` 个根，确保所选根 \(r\) 满足：
  1. 谱半径条件：\(|R(i\omega\Delta t)| \le 1\)（无条件稳定）
  2. 相期误差最小（精度最优）

对本例 \(M=4\)，`ir(4)=2`，选取升序排列后的第 2 个根。

---

### 5.3 `MP1SchemeRoot.m`（第 1–27 行）

```matlab
function [r] = MP1SchemeRoot(M)
```

**作用**：返回 \((M+1)\) 阶精度格式（M 个子步）的预设根 \(r\)（硬编码值）。

```matlab
switch M
    case 2
        r = 3 - sqrt(3);         % r ≈ 1.267949192431123
    case 3
        r = 0.9358222275240879;
    case 5
        r = 2.1129659585785242;
end
```

\((M+1)\)-格式的根通过令 \(p_{M+1}(r) = 0\)（使截断误差阶次从 \(M+1\) 提升到 \(M+2\)）求得：

\[
p_{M+1}(r) = \sum_{j=0}^{M} (-1)^j\frac{M!}{j!(M-j)!}\cdot\frac{1}{(M+1-j)!}\,r^{M-j} = 0
\]

只有 \(M = 2, 3, 5\) 时，方程的实数根满足 \(|R(\lambda)|\le1\)（无条件稳定），其他 \(M\) 值不适用。根值由 Mathematica 预先计算并硬编码于此。

---

### 5.4 `pCoefficients.m`（第 1–24 行）

```matlab
function pcoe = pCoefficients(M, r)
```

**作用**：根据已确定的根 \(r\)，计算分子多项式 \(P(\lambda)\) 的系数 \(p_0, p_1, \ldots, p_M\)。

```matlab
% 第 18–23 行
pcoe = zeros(1,M+1);
for ii = 0:M
    j = 0:ii;
    p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
    pcoe(ii+1) = p*(r.^(M-j))';
end
```

**关键公式**（论文公式 (40)）：

\[
p_i(r) = \sum_{j=0}^{\min(i,M)} (-1)^j \frac{M!}{j!(M-j)!} \cdot \frac{1}{(i-j)!}\,r^{M-j}, \quad i = 0, 1, \ldots, M
\]

- **第 18 行**：`pcoe` 初始化为零向量，长度 \(M+1\)，`pcoe(i+1)` 存储 \(p_i\)。
- **第 19–23 行**：外循环对 \(i = 0, 1, \ldots, M\) 逐一计算 \(p_i(r)\)：
  - `j = 0:ii`：求和指标（因为当 \(j > i\) 时 \(1/(i-j)! = 0\)，无贡献）。
  - 第 21 行按公式逐元素计算求和项系数：`(-1)^j * M!/（(M-j)!*j!*(i-j)!）`，其中 `factorial(ii-j)` 对应 \(1/(i-j)!\)。
  - 第 22 行点乘 \(r^{M-j}\) 并求和得到 \(p_i\)。

**物理含义**：\(p_i\) 是 Taylor 展开 \((r-\lambda)^M e^\lambda\) 在 \(\lambda = 0\) 处 \(\lambda^i\) 项的系数。当选定根 \(r\) 使 \(p_M = \pm\rho_\infty\) 时，格式具有 M 阶精度和指定的数值耗散特性。

---

### 5.5 `shiftPolycoe.m`（第 1–26 行）

```matlab
function [prcoe] = shiftPolycoe(pcoe, r)
```

**作用**：将多项式从 \(\mathbf{A}\) 的幂次基 \(\{1, \mathbf{A}, \mathbf{A}^2, \ldots\}\) 转换到 \(\mathbf{A}_r\) 的幂次基 \(\{1, \mathbf{A}_r, \mathbf{A}_r^2, \ldots\}\)，其中 \(\mathbf{A}_r = r\mathbf{I} - \mathbf{A}\)，即变量替换 \(\mathbf{A} = r\mathbf{I} - \mathbf{A}_r\)。

若输入多项式为

\[
P(\mathbf{A}) = \sum_{i=0}^{M} \text{pcoe}(i+1)\cdot\mathbf{A}^i
\]

则输出多项式为

\[
P_r(\mathbf{A}_r) = P(r\mathbf{I} - \mathbf{A}_r) = \sum_{i=0}^{M} \text{prcoe}(i+1)\cdot\mathbf{A}_r^i
\]

```matlab
% 第 18–24 行
p     = size(pcoe,2) - 1;
zc    = zeros(size(pcoe,1),1);
prcoe = pcoe;
for ii = p:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1)  r*prcoe(:,ii+2:end)  zc] ...
                    - [zc  prcoe(:,ii+1:end)];
end
```

- **第 18 行**：`p = M`（多项式次数）。
- **第 19 行**：`zc = zeros(行数, 1)` 是一个零列向量，用于在数组拼接时填充维度。
- **第 20 行**：`prcoe = pcoe`，将输出初始化为输入，随后就地更新。
- **第 21–24 行（核心循环）**：从最高次 `ii = p` 到 `ii = 1` 逐次执行多项式基变换。

**算法本质**：这是一种迭代的 Horner 型算法，每一步执行变量替换 \(A \to r - A_r\)，相当于对当前多项式做一步"移位"。每步迭代：

拼接右侧：`[prcoe(:,ii)+r*prcoe(:,ii+1),  r*prcoe(:,ii+2:end),  zc]`  
减去左侧：`[zc,  prcoe(:,ii+1:end)]`

实际效果等价于对从第 `ii` 列到末尾的部分多项式执行：

\[
\tilde{p}_{ii} \leftarrow p_{ii} + r\,p_{ii+1}, \quad
\tilde{p}_{k} \leftarrow r\,p_{k+1} - p_{k+1} = (r-1)\,p_{k+1}, \quad k > ii
\]

经过完整循环后，系数向量完成变量替换 \(A^i \to (r - A_r)^i\) 的展开。

**用法一**（第 30 行）：`qcoe = shiftPolycoe([0,...,0,1], r)`  
将 \(Q_r(\mathbf{A}_r) = \mathbf{A}_r^M\) 转为 \(Q(\mathbf{A}) = (r\mathbf{I}-\mathbf{A})^M\) 在 \(\mathbf{A}\) 基下的展开系数。

**用法二**（第 34 行）：`prcoe = shiftPolycoe(pcoe, r)`  
将 \(P(\mathbf{A})\) 转为 \(P_r(\mathbf{A}_r)\) 的展开系数，供主循环中计算 \(\mathbf{P}_r\mathbf{z}_{n-1}\) 使用。

---

### 5.6 `TimeIntgCoeffForce.m`（第 1–26 行）

```matlab
function [C] = TimeIntgCoeffForce(p, q)
```

**作用**：给定分子多项式系数 `p`（对应 \(\mathbf{P}\)）和分母多项式系数 `q`（对应 \(\mathbf{Q}\)），计算外力积分系数矩阵 \(\mathbf{C}\)，使得：

\[
\mathbf{C}_k\,\tilde{\mathbf{F}} = \mathbf{C}(k+1,:)\cdot[\text{状态向量各 }\mathbf{A}\text{ 幂次分量}]
\]

其中 \(\mathbf{C}_k\)（\(k=0,1,\ldots,M\)）由下述递推公式（论文公式 (19)/(20)）定义：

\[
\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P} - \mathbf{Q})
\]

\[
\mathbf{C}_k = \mathbf{A}^{-1}\!\left(k\,\mathbf{C}_{k-1} + (-0.5)^k\bigl(\mathbf{P} - (-1)^k\mathbf{Q}\bigr)\right), \quad k = 1, 2, \ldots, M
\]

---

```matlab
% 第 18–21 行
M   = length(q) - 1;
tmp = p - q;
C   = zeros(M+1, M);
C(1,:) = tmp(2:end);
```

- **第 18 行**：`M = length(q)-1` 是多项式次数（子步数）。
- **第 19 行**：`tmp = p - q` 计算多项式 \(\mathbf{P} - \mathbf{Q}\) 的系数向量。
- **第 20 行**：`C` 为 \((M+1)\times M\) 矩阵，`C(k+1,:)` 存储 \(\mathbf{C}_k\) 的多项式系数（长度 \(M\)，对应 \(\mathbf{A}^1, \mathbf{A}^2, \ldots, \mathbf{A}^M\) 各次项，常数项已被 \(\mathbf{A}^{-1}\) 消除）。
- **第 21 行**：计算 \(\mathbf{C}_0 = \mathbf{A}^{-1}(\mathbf{P}-\mathbf{Q})\)。

**关键解释**：由于 \(p_0 = q_0\)（归一化条件 \(e^0 = \mathbf{I}\) 要求 \(\mathbf{P}(0) = \mathbf{Q}(0)\)），故 `tmp(1) = p(1) - q(1) = 0`。乘以 \(\mathbf{A}^{-1}\) 等价于在多项式系数数组中丢弃常数项：

\[
\mathbf{A}^{-1}(c_1\mathbf{A} + c_2\mathbf{A}^2 + \cdots) = c_1\mathbf{I} + c_2\mathbf{A} + \cdots
\]

在系数数组中体现为 `tmp(2:end)`（去掉对应 \(\mathbf{A}^0\) 的第一个元素）。因此 `C(1,:) = tmp(2:end)` 存储的是次数降低一阶的多项式系数（从 \(\mathbf{A}^0\) 开始，长度 \(M\)）。

---

```matlab
% 第 22–25 行
for k = 1:M
    tmp = ((-1/2)^k)*(p - ((-1)^k)*q);
    tmp(1:M) = tmp(1:M) + k*C(k,:);
    C(k+1,:) = tmp(2:end);
end
```

- **第 23 行**：计算 \((-0.5)^k(\mathbf{P} - (-1)^k\mathbf{Q})\) 的系数向量：

\[
\text{tmp} = (-0.5)^k(p - (-1)^k q)
\]

- **第 24 行**：加上 \(k\,\mathbf{C}_{k-1}\) 的贡献（`k*C(k,:)` 即 \(k\) 倍的 \(\mathbf{C}_{k-1}\) 系数向量，注意 `C(k,:)` 对应 \(\mathbf{C}_{k-1}\)）。但 `C(k,:)` 的长度为 \(M\)，对应 \(\mathbf{A}^0, \ldots, \mathbf{A}^{M-1}\)，加到 `tmp(1:M)`（对应 \(\mathbf{A}^0, \ldots, \mathbf{A}^{M-1}\)）：

\[
\text{tmp}(1:M) \leftarrow \text{tmp}(1:M) + k\cdot\mathbf{C}_{k-1}\text{的系数}
\]

这实现了公式中的 \(k\,\mathbf{C}_{k-1}\) 加法。

- **第 25 行**：`C(k+1,:) = tmp(2:end)`，同样通过去掉常数项实现 \(\mathbf{A}^{-1}\) 运算，得到 \(\mathbf{C}_k\) 的系数。

**整体效果**：`C` 矩阵的 \((k+1)\) 行存储了 \(\mathbf{C}_k\)（以多项式系数形式，乘以状态向量的 \(\mathbf{A}\) 幂次分量即得矩阵作用结果）。主循环第 48 行 `Pf = f * cf'` 利用这一结构高效组装力的贡献。

---

## 六、核心计算子函数

### 6.1 `SolverPadeAx.m`（第 1–24 行）

```matlab
function y = SolverPadeAx(dM, K, C, x)
```

**作用**：计算矩阵-向量乘积 \(\mathbf{y} = \mathbf{A}\,\mathbf{x}\)，其中 \(\mathbf{A}\) 为无量纲状态矩阵（隐式定义，从不显式构造）。

状态矩阵的定义为（论文公式 (8)）：

\[
\mathbf{A} = \begin{bmatrix} -\mathbf{M}^{-1}\mathbf{C} & -\mathbf{M}^{-1}\mathbf{K} \\ \mathbf{I} & \mathbf{0} \end{bmatrix}
\]

（其中 \(\mathbf{C}\) 已含因子 \(\Delta t\)，\(\mathbf{K}\) 已含因子 \(\Delta t^2\)）

状态向量分块：\(\mathbf{x} = [\mathbf{x}_1;\, \mathbf{x}_2]\)，其中 \(\mathbf{x}_1 \in \mathbb{R}^n\)（无量纲速度），\(\mathbf{x}_2 \in \mathbb{R}^n\)（位移）。

```matlab
% 第 18–20 行
n   = length(x)/2;
tmp = C*x(1:n) + K*x(n+1:end);
y   = [-dM\tmp ; x(1:n)];
```

- **第 18 行**：`n = ndof`（自由度数，本例为 2）。
- **第 19 行**：计算上半部分的分子：

\[
\text{tmp} = \mathbf{C}\,\mathbf{x}_1 + \mathbf{K}\,\mathbf{x}_2 = \Delta t\mathbf{C}\,\mathbf{x}_1 + \Delta t^2\mathbf{K}\,\mathbf{x}_2
\]

- **第 20 行**：组装结果向量：

\[
\mathbf{y} = \begin{bmatrix} -\mathbf{M}^{-1}(\mathbf{C}\,\mathbf{x}_1 + \mathbf{K}\,\mathbf{x}_2) \\ \mathbf{x}_1 \end{bmatrix} = \begin{bmatrix} -\mathbf{M}^{-1}\mathbf{C} & -\mathbf{M}^{-1}\mathbf{K} \\ \mathbf{I} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{x}_1 \\ \mathbf{x}_2 \end{bmatrix} = \mathbf{A}\,\mathbf{x}
\]

其中 `dM\tmp` = \(\mathbf{M}^{-1}(\mathbf{C}\mathbf{x}_1 + \mathbf{K}\mathbf{x}_2)\)（利用预计算的 LU 分解 `dM`）。

**关键设计**：`SolverPadeAx` 从不显式构造 \(2n \times 2n\) 的 \(\mathbf{A}\) 矩阵，而是通过分块矩阵乘法实现，避免了大矩阵存储和计算的开销。这使得算法对大规模系统（\(n \gg 1\)）同样高效。

**调用场景**：

1. **主循环第 55 行**（初始化）：`z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii))` 计算 \((r\mathbf{I} - \mathbf{A})\mathbf{z} = \mathbf{A}_r\mathbf{z}\)。
2. **第 46 行**（力预计算）：`f(:,ii) = SolverPadeAx(dM,K,C,f(:,ii-1))` 计算 \(\mathbf{A}^{ii-1}\tilde{\mathbf{F}}\)。

---

## 七、时步主循环深度解析

回到 `TimeSolverRExpn.m` 主循环（第 65–80 行），给出更深入的数学解释。

---

### 7.1 矩阵 z 的物理含义

在进入主循环之前，`z` 的各列含义如下（初始化完成后）：

\[
\mathbf{z}_{:,k+1} = \mathbf{A}_r^k\,\mathbf{z}_{n-1}, \quad k = 0, 1, \ldots, M
\]

在每步循环结束时，`z(:,1)` 被更新为新时步的状态向量 \(\mathbf{z}_n\)，其余列被重新计算为 \(\mathbf{A}_r^k\mathbf{z}_n\)。这 \(M+1\) 列状态向量构成了在 \(\mathbf{A}_r\) 幂次基下的完整表示，是格式递推结构的数学基础。

---

### 7.2 右端项计算（第 68 行）详细推导

```matlab
z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
```

设 \(p = M\)，则：

**左端** `z(:,p+1)` 即 \(\mathbf{z}^{(0)}\)（方程 \(\mathbf{A}_r^M\mathbf{z}_n = \mathbf{z}^{(0)}\) 的右端项）。

**第一项** `(z(:,1:p+1))*prcoe'`：

矩阵 `z(:,1:p+1)` 的列分别为 \(\mathbf{A}_r^0\mathbf{z}_{n-1}, \ldots, \mathbf{A}_r^M\mathbf{z}_{n-1}\)；`prcoe` 是 \(1\times(M+1)\) 行向量（\(\mathbf{P}\) 在 \(\mathbf{A}_r\) 基下的系数）。乘积为：

\[
\sum_{k=0}^{M} \text{prcoe}(k+1)\cdot\mathbf{A}_r^k\,\mathbf{z}_{n-1} = P_r(\mathbf{A}_r)\,\mathbf{z}_{n-1} = \mathbf{P}(\mathbf{A})\,\mathbf{z}_{n-1}
\]

**第二项** `Pf*ft(is,:)'`：

`Pf` 为 \((2n)\times M\) 矩阵，`ft(is,:)'` 为 \(M\times1\) 向量（当前时步力多项式系数 \([c_0, c_1, \ldots, c_{M-1}]^T\)），乘积为：

\[
\mathbf{P}_f\,\tilde{\mathbf{c}}_n = \sum_{k=0}^{M-1} \mathbf{C}_k\,\tilde{\mathbf{F}}\cdot c_{n,k} = \sum_{k=0}^{M-1} \mathbf{C}_k\,\mathbf{F}_{mn}^{(k)}
\]

其中 \(\mathbf{F}_{mn}^{(k)} = c_{n,k}\tilde{\mathbf{F}}\) 是时步中点处展开的 \(k\) 阶力系数。

**合并两项**：

\[
\mathbf{z}^{(0)} = \mathbf{P}(\mathbf{A})\,\mathbf{z}_{n-1} + \sum_{k=0}^{M-1}\mathbf{C}_k\,\mathbf{F}_{mn}^{(k)}
\]

这正是时步方程（论文公式 (82)）的右端：

\[
\mathbf{A}_r^M\,\mathbf{z}_n = \mathbf{P}_r(\mathbf{A}_r)\,\mathbf{z}_{n-1} + \sum_{k=0}^{M} \mathbf{C}_k\,\tilde{\mathbf{F}}_{mn}^{(k)}
\]

---

### 7.3 子步递推（第 69–72 行）详细推导

```matlab
for ip = p:-1:1
    z(1:ndof,ip)        = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
    z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
end
```

求解系统 \(\mathbf{A}_r^M\mathbf{z}_n = \mathbf{z}^{(0)}\) 等价于顺序求解 \(M\) 个线性系统：

\[
\mathbf{A}_r\,\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)}, \quad k=M,M-1,\ldots,1
\]

（循环从 `ip=p` 到 `ip=1`，注意 `z(:,ip+1) = z^{(k-1)}`，`z(:,ip) = z^{(k)}`）

展开 \(\mathbf{A}_r = r\mathbf{I}-\mathbf{A}\)，将方程分块：

**上半块**（速度分量，\(n\) 个方程）：

\[
\bigl(r\mathbf{I} - \mathbf{A}\bigr)\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)}
\]

其上半行为：
\[
r\,\mathbf{z}_1^{(k)} - \bigl(-\mathbf{M}^{-1}\mathbf{C}\,\mathbf{z}_1^{(k)} - \mathbf{M}^{-1}\mathbf{K}\,\mathbf{z}_2^{(k)}\bigr) = \mathbf{z}_1^{(k-1)}
\]
\[
(r\mathbf{I} + \mathbf{M}^{-1}\mathbf{C})\,\mathbf{z}_1^{(k)} + \mathbf{M}^{-1}\mathbf{K}\,\mathbf{z}_2^{(k)} = \mathbf{z}_1^{(k-1)}
\]

乘以 \(\mathbf{M}\)：

\[
(r\mathbf{M} + \mathbf{C})\,\mathbf{z}_1^{(k)} + \mathbf{K}\,\mathbf{z}_2^{(k)} = \mathbf{M}\,\mathbf{z}_1^{(k-1)}
\]

但此式含未知量 \(\mathbf{z}_2^{(k)}\)，需先由下半块关系消去。

**下半块**（位移分量，\(n\) 个方程）：

\[
r\,\mathbf{z}_2^{(k)} - \mathbf{z}_1^{(k)} = \mathbf{z}_2^{(k-1)}
\quad\Rightarrow\quad
\mathbf{z}_2^{(k)} = \frac{\mathbf{z}_1^{(k)} + \mathbf{z}_2^{(k-1)}}{r}
\]

（第 71 行正是此式）

将 \(\mathbf{z}_2^{(k)}\) 代入上半块方程：

\[
(r\mathbf{M} + \mathbf{C})\,\mathbf{z}_1^{(k)} + \mathbf{K}\cdot\frac{\mathbf{z}_1^{(k)} + \mathbf{z}_2^{(k-1)}}{r} = \mathbf{M}\,\mathbf{z}_1^{(k-1)}
\]

乘以 \(r\) 整理：

\[
\bigl(r^2\mathbf{M} + r\mathbf{C} + \mathbf{K}\bigr)\,\mathbf{z}_1^{(k)} = r\mathbf{M}\,\mathbf{z}_1^{(k-1)} - \mathbf{K}\,\mathbf{z}_2^{(k-1)}
\]

即：

\[
\mathbf{K}_d\,\mathbf{z}_1^{(k)} = r\mathbf{M}\,\mathbf{z}_1^{(k-1)} - \mathbf{K}\,\mathbf{z}_2^{(k-1)}
\]

这正是第 70 行所实现的方程，`dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1))` 求解此线性方程组得到 \(\mathbf{z}_1^{(k)}\)（上半部分无量纲速度），再由第 71 行的代数关系得到 \(\mathbf{z}_2^{(k)}\)（下半部分位移）。

**每个子步只需一次三角回代**（利用预分解的 `dKd`），且所有子步共享同一 \(\mathbf{K}_d\)，总计算量为 \(M\) 次三角回代，相比原始矩阵指数方法节省了巨大计算量。

---

### 7.4 加速度计算（第 78 行）

```matlab
acc(it,:) = r*z(pDOF,1) - z(pDOF,2) + fn*Fb(pDOF);
```

从定义 \(\mathbf{z}_{:,2} = \mathbf{A}_r\mathbf{z}_{:,1}\) 出发：

\[
\mathbf{z}_{:,2} = (r\mathbf{I} - \mathbf{A})\,\mathbf{z}_{:,1}
\quad\Rightarrow\quad
[\mathbf{A}\mathbf{z}_{:,1}]_{\text{上半}} = r\,\mathbf{z}_1 - [\mathbf{z}_{:,2}]_{\text{上半}}
\]

由 \(\mathbf{A}\) 的定义，其作用在状态向量上产生的上半部分为 \(-\mathbf{M}^{-1}(\mathbf{C}\mathbf{z}_1 + \mathbf{K}\mathbf{z}_2)\)，即 \(\dot{\mathbf{u}}^{\circ\circ}\)（无量纲加速度）。加上力的贡献 \(f(t_n)\mathbf{F}_b\)（其上半部分即 \(\Delta t^2\mathbf{M}^{-1}\mathbf{F}_0 f(t_n)\)），得无量纲加速度：

\[
\mathbf{u}^{\circ\circ}(t_n) \approx r\,\mathbf{z}_1 - [\mathbf{z}_{:,2}]_{\text{上半}} + f(t_n)\,\mathbf{F}_b
\]

---

## 八、完整调用关系图

```
ExampleThreeDOFs.m（主脚本）
│
├── ThreeDOFsRefSln（局部函数，精确参考解）
│     └── eig, 匿名函数
│
└── TimeSolverRExpn（主求解器）
      │
      ├── ForceSeries（外力多项式展开）
      │     ├── forceSamplingPoints（采样点生成）
      │     │     └── lglnodes（LGL 节点 Newton 迭代）
      │     └── transMtxPointsToPoly（变换矩阵构造）
      │
      ├── InitSchemeRExpn（格式系数初始化）
      │     ├── MschemeRoot（M格式根求解）
      │     │     └── roots（MATLAB 内置多项式根）
      │     ├── MP1SchemeRoot（(M+1)格式预设根）
      │     ├── pCoefficients（分子多项式系数）
      │     ├── shiftPolycoe（多项式基变换，调用两次）
      │     └── TimeIntgCoeffForce（外力积分系数）
      │
      └── SolverPadeAx（A 矩阵向量乘积，多次调用）
```

---

## 九、关键设计哲学总结

| 设计特点 | 实现位置 | 意义 |
|---|---|---|
| 单重根分母多项式 | `MschemeRoot`/`MP1SchemeRoot` | 保证所有子步有效刚度矩阵相同 |
| 无量纲时间变量 | `TimeSolverRExpn` 第 18–25 行 | 使算法与 \(\Delta t\) 解耦，便于分析 |
| LGL 采样点 | `forceSamplingPoints` + `lglnodes` | 高阶多项式插值的最优节点分布 |
| A 矩阵不显式构造 | `SolverPadeAx` | 避免大矩阵存储，适合大规模系统 |
| 有效刚度矩阵仅分解一次 | `TimeSolverRExpn` 第 41 行 | 大幅降低每步计算成本 |
| \(\mathbf{A}_r\) 幂次预存储 | `z` 的多列结构 | 将 \(\mathbf{P}_r\mathbf{z}_{n-1}\) 化为简单向量内积 |
