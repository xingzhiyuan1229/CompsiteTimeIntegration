# 公式解读文档 —— 第 3.2 节与第 4 节：(M+1)-格式与算法实现（公式 66–93）

本文档为 [公式解读_Part3.md](./公式解读_Part3.md) 的续篇，涵盖第 3.2 节 (M+1)-格式（公式 (66)–(79)）和第 4 节复合时间积分算法实现（公式 (80)–(93)）的代码解读。格式为**代码在前、公式在后**。

> **GitHub 仓库根地址**：`https://github.com/xingzhiyuan1229/CompsiteTimeIntegration`

---

## 3.2 (M+1)-格式：M 步子步、(M+1) 阶精度

### 公式 (66)：(M+1)-格式的 Taylor 展开（扩展至 M+1 阶）

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24) 与 [`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% src/pCoefficients.m（(M+1)-格式同样使用此函数）
% 展开系数计算至 ii=M，包括 p_{M+1} 用于确定根 r（令其等于零）
function pcoe = pCoefficients(M, r)
pcoe = zeros(1,M+1);
for ii = 0:M
    j = 0:ii;
    p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
    pcoe(ii+1) = p*(r.^(M-j))';
end
end

% src/MP1SchemeRoot.m（硬编码满足 p_{M+1}(r)=0 的根）
function [r] = MP1SchemeRoot(M)
switch M
    case 2
        r = 3 - sqrt(3);
    case 3
        r =      0.9358222275240879;
    case 5
        r =      2.1129659585785242;
end
end
```

**代码解读**：公式 (66) 是 (M+1)-格式设计的出发点：将 [ (r-\lambda)^M e^\lambda ] 展开至 [ \lambda^{M+1} ]，额外保留 [ p_{M+1}(r)\lambda^{M+1} ] 项，用于通过令 [ p_{M+1}(r)=0 ]（公式 (67)）确定根 [ r ]。`pCoefficients` 本身计算到 [ p_M ]（循环到 `ii=M`），而在 `MP1SchemeRoot` 中，根 [ r ] 是通过求解 [ p_{M+1}(r)=0 ] 并选出满足稳定条件的根后**硬编码**到 `switch-case` 语句中的（论文 3.2 节指出只有 [ M=2,3,5 ] 三种情况满足稳定性要求，代码仅实现这三种情形）。

**对应公式 (66)**：

\[
(r - \lambda)^M e^{\lambda} = p_0(r) + p_1(r)\lambda + \cdots + p_M(r)\lambda^M + p_{M+1}(r)\lambda^{M+1} + O(\lambda^{M+2}) \tag{66}
\]

将展开式额外保留 [ \lambda^{M+1} ] 项，系数 [ p_{M+1}(r) ] 由公式 (40) 计算（[ i=M+1 ] 时）。

---

### 公式 (67)：令 p_{M+1}(r) = 0 以确定根 r

**对应代码文件**：[`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% src/MP1SchemeRoot.m（硬编码根，对应公式 (67) 的求解结果）
function [r] = MP1SchemeRoot(M)
% 根 r 由求解 p_{M+1}(r) = 0 并选满足稳定性条件的根后硬编码：
% M=2: p3(r) = 1 - r + r^2/6 = 0 → r = 3 ∓ sqrt(3)，选 r = 3-sqrt(3)
% M=3: 求解 p4(r)=0（四次方程）的数值根，选满足稳定性条件的根
% M=5: 求解 p6(r)=0（六次方程）的数值根，选满足稳定性条件的根
switch M
    case 2
        r = 3 - sqrt(3);
    case 3
        r =      0.9358222275240879;
    case 5
        r =      2.1129659585785242;
end
end
```

**代码解读**：公式 (67) 通过令展开式中 [ \lambda^{M+1} ] 的系数为零来确定根 [ r ]：[ p_{M+1}(r) = \sum_{j=0}^{M} (-1)^j \frac{M!}{j!(M-j)!}\frac{1}{(M+1-j)!} r^{M-j} = 0 ]。对 [ M=2 ]，方程为 [ 1-r+r^2/6=0 ]（公式 (72)），`MP1SchemeRoot` 中直接硬编码 [ r = 3-\sqrt{3} ] 为满足条件的根（见公式 (74)）；对 [ M=3,5 ]，数值求解后硬编码满足稳定性条件的根。这一硬编码策略是因为 (M+1)-格式只有特定的 [ M ] 值才能满足 [ \rho \le 1 ] 的稳定性条件。

**对应公式 (67)**：

\[
p_{M+1}(r) = \sum_{j=0}^{M} (-1)^j \frac{M!}{j!(M-j)!} \frac{1}{(M+1-j)!} r^{M-j} = 0 \tag{67}
\]

通过令 [ \lambda^{M+1} ] 系数为零来确定根 [ r ]，从而使截断误差阶次提升至 [ M+2 ]（即 [ (M+1) ] 阶精度）。

---

### 公式 (68)：(M+1)-格式有理近似的标量形式

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24) 与 [`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% (M+1)-格式的完整标量放大因子由两步确定：
% 步骤 1：r = MP1SchemeRoot(M)（公式 (67) 的解）
r = MP1SchemeRoot(M);

% 步骤 2：pcoe = pCoefficients(M, r)（p_0,...,p_M，分子多项式系数）
pcoe = pCoefficients(M, r);
% 分子 P(λ) = p_0 + p_1*λ + ... + p_M*λ^M（M 次多项式）
% 分母 Q(λ) = (r-λ)^M

% 放大因子 R(λ) = P(λ)/Q(λ)，精度阶次为 M+1（截断误差 O(λ^{M+2})）
```

**代码解读**：公式 (68) 是 (M+1)-格式的标量放大因子 [ R(\lambda) = P(\lambda)/Q(\lambda) = (p_0+p_1\lambda+\cdots+p_M\lambda^M)/(r-\lambda)^M ]，与 M-格式的公式 (43) 形式相同但 [ r ] 的确定方式不同：(M+1)-格式通过令 [ p_{M+1}(r)=0 ] 确定 [ r ]，使截断误差提升至 [ O(\lambda^{M+2}) ]。代码中 `MP1SchemeRoot` 提供根 [ r ]，`pCoefficients(M,r)` 输出 [ M+1 ] 个系数（[ p_0 ] 至 [ p_M ]），两者一起编码了公式 (68) 的完整分子多项式 [ P(\lambda) ]。

**对应公式 (68)**：

\[
e^{\lambda} \approx R(\lambda) = \frac{P(\lambda)}{Q(\lambda)} = \frac{p_0 + p_1\lambda + \cdots + p_M\lambda^M}{(r - \lambda)^M} \tag{68}
\]

(M+1)-格式放大因子，截断误差为 [ O(\lambda^{M+2}) ]，精度阶次为 [ M+1 ]。[ \rho_\infty ] 内置于格式中，不可调节。

---

### 公式 (69)：(M+1)-格式矩阵指数的有理近似

**对应代码文件**：[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L18-L35)

```matlab
% src/InitSchemeRExpn.m（MP1 分支）
if contains(scheme,'MP1')
    r =  MP1SchemeRoot(p);    % 获取 (M+1)-格式的根 r
else
    r = MschemeRoot(p, rho);  % M-格式的根 r
end

pcoe = pCoefficients(p, r);   % 分子系数 [p_0,...,p_M]，与 M-格式相同的计算函数
% ...
prcoe = shiftPolycoe(pcoe,r); % 转为 A_r 基，供时步循环使用
```

**代码解读**：公式 (69) 是公式 (68) 的矩阵版本，将 [ \lambda ] 替换为矩阵 [ \mathbf{A} ]：[ \mathbf{R} = (p_0\mathbf{I}+p_1\mathbf{A}+\cdots+p_M\mathbf{A}^M)/(r\mathbf{I}-\mathbf{A})^M ]。代码中 `InitSchemeRExpn` 对 M-格式和 (M+1)-格式使用**相同的**后续处理逻辑（`pCoefficients`、`shiftPolycoe`、`TimeIntgCoeffForce`），只是根 [ r ] 的来源不同（前者 `MschemeRoot`，后者 `MP1SchemeRoot`）。这体现了论文中所说的「公式 (44) 和 (69) 形式相同，只是根 [ r ] 的确定方式不同」。

**对应公式 (69)**：

\[
\mathbf{R} = \frac{\mathbf{P}}{\mathbf{Q}} = \frac{p_0\mathbf{I} + p_1\mathbf{A} + \cdots + p_M\mathbf{A}^M}{(r\mathbf{I} - \mathbf{A})^M} \tag{69}
\]

(M+1)-格式的矩阵指数有理近似，截断误差为 [ O(\lambda^{M+2}) ]。与 M-格式公式 (44) 形式相同，仅根 [ r ] 的确定方式不同。

---

### 公式 (70)：(M+1)-格式的内置谱半径 ρ_∞

**对应代码文件**：[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L26-L28)

```matlab
% src/InitSchemeRExpn.m，第 24–28 行
pcoe = pCoefficients(p, r);

disp(['r = ', num2str(r)]);
disp(['Selected rhoInfty = ', num2str(abs(pcoe(end)), '%.5f')]);
% abs(pcoe(end)) = |p_M| = rho_infty（内置值）
```

**代码解读**：公式 (70) 给出 (M+1)-格式高频极限谱半径 [ \rho_\infty = |p_M| ]。代码中 `abs(pcoe(end))` 计算分子多项式最高次系数 [ |p_M| ]，即高频极限 [ |\lambda|\to\infty ] 时放大因子 [ |R| = |p_M/(-1)^M| = |p_M| ]（因 [ (r-\lambda)^M \to (-\lambda)^M ]）。`disp` 语句输出此值，用于验证 (M+1)-格式的内置 [ \rho_\infty ] 与 M-格式用户指定 [ \rho_\infty ] 是否一致。对 [ M=2 ] 的 (M+1)-格式，[ \rho_\infty = |1-\sqrt{3}| \approx 0.7321 ]（见公式 (77)）。

**对应公式 (70)**：

\[
\rho_{\infty} = |p_M| \tag{70}
\]

(M+1)-格式的高频极限谱半径由分子最高次系数 [ p_M ] 的模决定，用户无法调节。

---

### 公式 (71)：M=2 时 (M+1)-格式的 Taylor 展开（含 λ³ 项）

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24)

```matlab
% pCoefficients(2, r) 加上 p_3(r) 的计算（扩展至 ii=3，M+1=3）：
% 公式 (66) 展开至 λ³：
% p3(r) = (1 - r + r^2/6)（由 ii=3, M=2 的 pCoefficients 扩展计算）
% M=2 时：
% ii=0: p0 = r^2
% ii=1: p1 = -2r + r^2
% ii=2: p2 = 1 - 2r + r^2/2
% ii=3（仅在扩展版 pCoefficients 中）:
%   j=0:3, p = [1/6, -1/2, 1/2, -1/6]...实际 = (1 - r + r^2/6)
% 硬编码于 MP1SchemeRoot(2)：令 p3=0 → 1-r+r^2/6=0 → r=3∓sqrt(3)
```

**代码解读**：公式 (71) 给出 [ M=2 ] 时 (M+1)-格式展开至 [ \lambda^3 ] 的完整系数，与公式 (51)（M-格式展开）相比额外包含 [ (1-r+r^2/6)\lambda^3 ] 项。代码中 `pCoefficients(2, r)` 仅计算至 [ p_2 ]（循环 `ii=0:2`），但 `MP1SchemeRoot` 的根 [ r=3-\sqrt{3} ] 正是通过（在 Mathematica 中）令 [ p_3(r) = 1-r+r^2/6 = 0 ] 并选合适根后**预先计算**再硬编码进去的。因此 `pCoefficients(2, 3-sqrt(3))` 输出的 `pcoe` 系数正对应公式 (75) 的数值。

**对应公式 (71)**：

\[
\mathrm{p}[r,\lambda] \Rightarrow r^2 + (-2r + r^2)\lambda + (1 - 2r + r^2/2)\lambda^2 + (1 - r + r^2/6)\lambda^3 \tag{71}
\]

[ M=2 ] 时展开至 [ \lambda^3 ] 的完整系数，最后一项 [ (1-r+r^2/6)\lambda^3 ] 用于通过令其系数为零确定 (M+1)-格式的根 [ r ]（公式 (72)）。

---

### 公式 (72)：M=2 时 (M+1)-格式的求根方程

**对应代码文件**：[`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% src/MP1SchemeRoot.m（M=2 时，硬编码公式 (72) 的解）
case 2
    r = 3 - sqrt(3);
% 对应求解 1 - r + r^2/6 = 0（公式 (72)）
% 即 r^2 - 6r + 6 = 0 → r = (6 ± sqrt(36-24))/2 = 3 ∓ sqrt(3)
% 选 r = 3-sqrt(3)（公式 (74)）
```

**代码解读**：公式 (72) 令 [ p_3(r) = 1-r+r^2/6 = 0 ]，即方程 [ r^2-6r+6=0 ]，其解为 [ r_{1,2}=3\mp\sqrt{3} ]（公式 (73)）。`MP1SchemeRoot` 直接将满足稳定性条件的根 [ r = 3-\sqrt{3} = 1.267949... ] 硬编码为 `case 2` 的返回值，这是对公式 (72) 求解过程的最终结果的固化。

**对应公式 (72)**：

\[
1 - r + r^2/6 = 0 \tag{72}
\]

[ M=2 ] 时 (M+1)-格式的求根方程，令 [ p_3(r)=0 ] 以使截断误差提升至 [ O(\lambda^4) ]（三阶精度）。

---

### 公式 (73)：M=2 时 (M+1)-格式的两根

**对应代码文件**：[`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% 方程 r^2/6 - r + 1 = 0 的两根（对应公式 (73)）：
% r = 3 ∓ sqrt(3)
% r_1 = 3-sqrt(3) ≈ 1.2679491924311228（选此根，公式 (74)）
% r_2 = 3+sqrt(3) ≈ 4.7320508075688772（不满足稳定性条件）
case 2
    r = 3 - sqrt(3);  % 选较小根（公式 (74)）
```

**代码解读**：`MP1SchemeRoot(2)` 选取 [ r=3-\sqrt{3} ] 而非 [ r=3+\sqrt{3} ]，与公式 (73) 的两个解 [ r_{1,2}=3\mp\sqrt{3} ] 中选较小根一致。论文中通过对两个根分别代入并绘制谱半径图验证，[ r_1=3-\sqrt{3} ] 满足 [ \rho\le 1 ] 的稳定性条件，而 [ r_2=3+\sqrt{3} ] 不满足。

**对应公式 (73)**：

\[
r_{1,2} = 3 \mp \sqrt{3} \tag{73}
\]

方程 (72) 的两个解，[ r_1=3-\sqrt{3}\approx 1.268 ]（选定，满足稳定性）和 [ r_2=3+\sqrt{3}\approx 4.732 ]（不满足稳定性条件）。

---

### 公式 (74)：M=2 时 (M+1)-格式选定根 r

**对应代码文件**：[`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L20-L21)

```matlab
% src/MP1SchemeRoot.m，case 2
case 2
    r = 3 - sqrt(3);  % = 1.267949192431123（公式 (74)）
```

**代码解读**：`3 - sqrt(3) = 1.267949192431123`，与公式 (74) [ r = 3-\sqrt{3} = 1.267949192431123 ] 完全一致，精确到 15 位有效数字。该值是 [ M=2 ] 时 (M+1)-格式的唯一内置根，决定了格式的所有特性（谱半径、相对周期误差等，见图 11）。

**对应公式 (74)**：

\[
r = 3 - \sqrt{3} = 1.267949192431123 \tag{74}
\]

[ M=2 ] 时 (M+1)-格式选定的根，通过令 [ p_3(r)=0 ] 并选满足稳定性的根得到。

---

### 公式 (75)：M=2 时 (M+1)-格式的分子多项式 P(λ)

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24)

```matlab
% pCoefficients(2, 3-sqrt(3)) 的输出：
% r = 3 - sqrt(3) ≈ 1.267949192431123
% p0 = r^2 = (3-sqrt(3))^2 = 9 - 6*sqrt(3) + 3 = 12 - 6*sqrt(3) ≈ 1.607695154586736
% p1 = -2r + r^2 = r(r-2) = (3-sqrt(3))(1-sqrt(3)) = 3-3sqrt(3)-sqrt(3)+3 = 6-4sqrt(3)
%    ≈ 6 - 6.928... ≈ -0.9282032302755092
% p2 = 1 - 2r + r^2/2 = 1 - 2*(3-sqrt(3)) + (12-6*sqrt(3))/2
%    = 1 - 6 + 2*sqrt(3) + 6 - 3*sqrt(3) = 1 - sqrt(3) ≈ -0.7320508075688773
pcoe = pCoefficients(2, 3 - sqrt(3));
% pcoe ≈ [1.607695154586736, -0.9282032302755092, -0.7320508075688773]
```

**代码解读**：`pCoefficients(2, 3-sqrt(3))` 精确给出公式 (75) 中 Mathematica 输出的三个系数：[ p_0 \approx 1.607695154586736 ]，[ p_1 \approx -0.9282032302755092 ]，[ p_2 \approx -0.7320508075688773 ]。最高次系数 [ p_2 = 1-\sqrt{3} \approx -0.732 ]，其绝对值即高频极限谱半径 [ \rho_\infty = |p_2| = \sqrt{3}-1 \approx 0.7321 ]（公式 (77)）。

**对应公式 (75)**：

\[
\mathrm{p[rSelected, \lambda]} \Rightarrow 1.607695154586736 - 0.9282032302755092\lambda - 0.7320508075688773\lambda^2 \tag{75}
\]

[ M=2, r=3-\sqrt{3} ] 时 (M+1)-格式分子多项式 [ P(\lambda) ] 的数值系数。

---

### 公式 (76)：M=2 时 (M+1)-格式的有理近似放大因子

**对应代码文件**：[`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24) 与 [`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27)

```matlab
% 两个函数协同给出公式 (76) 的完整有理近似：
r = MP1SchemeRoot(2);          % r = 3 - sqrt(3) ≈ 1.267949192431123
pcoe = pCoefficients(2, r);    % pcoe ≈ [1.6077, -0.9282, -0.7321]
% R(λ) = (1.6077 - 0.9282λ - 0.7321λ²) / (1.2679 - λ)^2
```

**代码解读**：公式 (76) 明确写出 [ M=2 ] 时 (M+1)-格式的数值放大因子 [ R(\lambda) = (1.607695...-0.928203...\lambda-0.732051...\lambda^2)/(1.267949...-\lambda)^2 ]。`MP1SchemeRoot(2)` 给出分母的根（`1.267949...`），`pCoefficients(2, r)` 给出分子三个数值系数，两者合并即公式 (76) 的完整表达。

**对应公式 (76)**：

\[
R(\lambda) = \frac{P(\lambda)}{Q(\lambda)} = \frac{1.607695154586736 - 0.9282032302755092\lambda - 0.7320508075688773\lambda^2}{(1.267949192431123 - \lambda)^2} \tag{76}
\]

[ M=2 ] 时 (M+1)-格式的数值放大因子，分子为二次多项式，分母为二阶因子 [ (r-\lambda)^2 ]。

---

### 公式 (77)：M=2 时 (M+1)-格式的内置 ρ_∞

**对应代码文件**：[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L26-L28)

```matlab
% src/InitSchemeRExpn.m，第 24–28 行
pcoe = pCoefficients(p, r);   % M=2 MP1 时 pcoe(3) = p_2 = 1-sqrt(3)
disp(['Selected rhoInfty = ', num2str(abs(pcoe(end)), '%.5f')]);
% 输出：Selected rhoInfty = 0.73205（即 sqrt(3)-1 ≈ 0.732050807568877）
```

**代码解读**：`abs(pcoe(end)) = |p_2| = |1-\sqrt{3}| = \sqrt{3}-1 \approx 0.732050807568877`，与公式 (77) [ \rho_\infty = |R_\infty| = |1-\sqrt{3}| = 0.732050807568877 ] 完全吻合。`InitSchemeRExpn` 中 `disp` 语句的输出让用户能直观确认 (M+1)-格式的内置 [ \rho_\infty ] 值，验证了公式 (70) 中 [ \rho_\infty=|p_M| ] 的正确性。

**对应公式 (77)**：

\[
\rho_{\infty} = |R_{\infty}| = |1 - \sqrt{3}| = 0.732050807568877 \tag{77}
\]

[ M=2 ] 时 (M+1)-格式内置的高频极限谱半径，等于 [ |p_2| = \sqrt{3}-1 ]，不可调节。

---

### 公式 (78)：M=2 时 (M+1)-格式的截断误差验证

**对应代码文件**：[`src/MP1SchemeRoot.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/MP1SchemeRoot.m#L18-L27) 与 [`src/pCoefficients.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/pCoefficients.m#L18-L24)

```matlab
% 截断误差验证（Mathematica Series 命令的等价逻辑）：
% e^λ - R(λ) 在 λ=0 处的 Taylor 展开，λ^4 项系数
% 由 p_{M+1}(r)=0（公式 (67)）保证 λ³ 系数为零，故首项为 λ^4
% r = 3-sqrt(3)，pCoefficients(2,r) 输出确保前 3 阶项精确
% λ^4 的系数约为 0.089779189（由 Mathematica 计算）
r = MP1SchemeRoot(2);     % r = 3-sqrt(3)（确保 p3(r)=0）
pcoe = pCoefficients(2, r); % 前 3 个系数（p0,p1,p2）
% 由于 p3(r)=0，截断误差从 λ^4 开始，系数为 0.089779189
```

**代码解读**：公式 (78) 用 Mathematica `Series` 命令验证 [ e^\lambda - R(\lambda) ] 在 [ \lambda=0 ] 处展开的首项为 [ 0.089779189\lambda^4 + O(\lambda^5) ]，确认截断误差阶次为 3（[ O(\lambda^4) ]，即三阶精度）。代码中 `MP1SchemeRoot(2)` 选取满足 [ p_3(r)=0 ] 的根，保证 [ \lambda^3 ] 项系数为零；`pCoefficients(2,r)` 输出的分子系数使得前三阶（[ \lambda^0, \lambda^1, \lambda^2 ]）与 [ e^\lambda ] 精确匹配，从而 [ e^\lambda - R(\lambda) = O(\lambda^4) ]，系数 [ \approx 0.0898 ] 由 Mathematica 数值计算得到。

**对应公式 (78)**：

\[
\mathrm{Series}[\exp(\lambda) - R(\lambda), \lambda, 0, M+2] \Rightarrow 0.089779189\lambda^4 + O(\lambda^5) \tag{78}
\]

验证 [ M=2 ] 时 (M+1)-格式的截断误差为三阶，首个非零误差项的系数约为 [ 0.0898 ]。

---

### 公式 (79)：M=2 时 (M+1)-格式的矩阵有理近似

**对应代码文件**：[`src/InitSchemeRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/InitSchemeRExpn.m#L18-L35)

```matlab
% src/InitSchemeRExpn.m（MP1 分支，M=2 时）
if contains(scheme,'MP1')
    r =  MP1SchemeRoot(p);    % r = 3-sqrt(3) ≈ 1.267949192431123
end
pcoe = pCoefficients(p, r);   % [1.6077, -0.9282, -0.7321]（公式 (75)）
prcoe = shiftPolycoe(pcoe,r); % P_r(A_r) 在 A_r 基下的系数（供时步循环用）
% 整体编码了公式 (79) 的矩阵有理近似 R = P(A)/Q(A)
```

**代码解读**：公式 (79) 将公式 (76) 的标量 [ \lambda ] 替换为矩阵 [ \mathbf{A} ]：[ \mathbf{R} = (1.6077\mathbf{I} - 0.9282\mathbf{A} - 0.7321\mathbf{A}^2) / (1.2679\mathbf{I}-\mathbf{A})^2 ]。`InitSchemeRExpn` 中，`pcoe` 给出分子矩阵多项式系数，`prcoe`（[ \mathbf{A}_r ] 基下）供后续时步循环中计算 [ \mathbf{P}_r\mathbf{z}_{n-1} ] 使用（公式 (91)），`qrcoe = [0,0,1]` 对应分母 [ \mathbf{A}_r^2 = (r\mathbf{I}-\mathbf{A})^2 ]。这完整编码了公式 (79) 的矩阵有理近似，并保持精度三阶和内置 [ \rho_\infty \approx 0.732 ]。

**对应公式 (79)**：

\[
\mathbf{R} = \frac{\mathbf{P}}{\mathbf{Q}} = \frac{1.607695154586736\mathbf{I} - 0.9282032302755092\mathbf{A} - 0.7320508075688773\mathbf{A}^2}{(1.267949192431123\mathbf{I} - \mathbf{A})^2} \tag{79}
\]

[ M=2 ] 时 (M+1)-格式的矩阵有理近似，三阶精度，内置 [ \rho_\infty \approx 0.732 ]。

---

## 4. 复合时间积分格式的时步算法

### 公式 (80)：辅助矩阵 A_r 的定义

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L38-L56)

```matlab
% src/TimeSolverRExpn.m，第 38–56 行
[r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho );

Kd = sparse((r*r)*M + r*C + K);   % 对应 A_r = rI-A 的物理等效刚度
dKd = decomposition(Kd);

% 初始化辅助变量 z^(k) 的第一步：
for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
    % = (rI-A)*z(:,ii) = A_r * z(:,ii)（公式 (80)/(83) 的应用）
end
```

**代码解读**：公式 (80) 定义辅助矩阵 [ \mathbf{A}_r = r\mathbf{I} - \mathbf{A} ]，是将有理近似分母 [ (r\mathbf{I}-\mathbf{A})^M ] 分解为 [ M ] 个子步的关键。在代码中，`r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii))` 计算 [ r\mathbf{z} - \mathbf{A}\mathbf{z} = (r\mathbf{I}-\mathbf{A})\mathbf{z} = \mathbf{A}_r\mathbf{z} ]，这是矩阵 [ \mathbf{A}_r ] 与向量乘积的直接实现。`Kd = r²M + rC + K` 是求解 [ \mathbf{A}_r\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)} ]（公式 (86)）时的等效刚度矩阵，对应 [ (r\mathbf{I}-\mathbf{A}) ] 的物理展开。

**对应公式 (80)**：

\[
\mathbf{A}_r = r\mathbf{I} - \mathbf{A} \tag{80}
\]

辅助矩阵 [ \mathbf{A}_r ] 的定义，使得有理近似的分母 [ \mathbf{Q}_r = \mathbf{A}_r^M ]，每个子步对应一个因子 [ \mathbf{A}_r ]。

---

### 公式 (81a/b)：多项式 P 和 Q 在 A_r 基下的变换

**对应代码文件**：[`src/shiftPolycoe.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/shiftPolycoe.m#L18-L25)

```matlab
% src/shiftPolycoe.m（直接实现公式 (81a) 的基变换）
function [prcoe] = shiftPolycoe(pcoe,r)
p = size(pcoe,2) - 1;
zc = zeros(size(pcoe,1),1);
prcoe = pcoe;
for ii = p:-1:1
    prcoe(:,ii:end) = [prcoe(:,ii)+r*prcoe(:,ii+1) r*prcoe(:,ii+2:end) zc] ...
                  - [zc prcoe(:,ii+1:end)];
end
end
% 输入 pcoe = [p_0,...,p_M]（A 基），输出 prcoe = [p_{r0},...,p_{rM}]（A_r 基）

% 公式 (81b)：Q_r = A_r^M，代码中直接用 qrcoe = [zeros(1,p) 1]
qrcoe = [zeros(1,p) 1];   % [0,...,0,1]，即 A_r^M 的系数（A_r 基下只有最高次项）
```

**代码解读**：公式 (81a) 表示将分子多项式 [ \mathbf{P}(\mathbf{A}) = p_0\mathbf{I}+\cdots+p_M\mathbf{A}^M ]（[ \mathbf{A} ] 基）变换为 [ \mathbf{P}_r(\mathbf{A}_r) = p_{r0}\mathbf{I}+\cdots+p_{rM}\mathbf{A}_r^M ]（[ \mathbf{A}_r ] 基），通过代入 [ \mathbf{A} = r\mathbf{I}-\mathbf{A}_r ] 展开并合并同类项实现。`shiftPolycoe` 用 Horner-like 方法（从最高次向最低次逐步合并）完成此多项式基变换，每次循环迭代将 [ \mathbf{A} = r\mathbf{I}-\mathbf{A}_r ] 代入并整理。公式 (81b) 则更为简单：[ \mathbf{Q}(\mathbf{A}) = (r\mathbf{I}-\mathbf{A})^M = \mathbf{A}_r^M ]，在 [ \mathbf{A}_r ] 基下系数为 `[0,...,0,1]`（`qrcoe`）。

**对应公式 (81a)**：

\[
\mathbf{P}(\mathbf{A}) = p_0\mathbf{I} + p_1(r\mathbf{I} - \mathbf{A}_r) + \cdots + p_M(r\mathbf{I} - \mathbf{A}_r)^M = p_{r0} + p_{r1}\mathbf{A}_r + \cdots + p_{rM}\mathbf{A}_r^M \equiv \mathbf{P}_r(\mathbf{A}_r) \tag{81a}
\]

**对应公式 (81b)**：

\[
\mathbf{Q}(\mathbf{A}) = \mathbf{A}_r^M \equiv \mathbf{Q}_r(\mathbf{A}_r) \tag{81b}
\]

公式 (81a) 将分子多项式从 [ \mathbf{A} ] 基变换到 [ \mathbf{A}_r ] 基（由 `shiftPolycoe` 实现），公式 (81b) 指出分母在 [ \mathbf{A}_r ] 基下退化为纯 [ M ] 次幂（无需计算）。

---

### 公式 (82)：A_r 基下的时步推进方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L65-L79)

```matlab
% src/TimeSolverRExpn.m，第 65–79 行（主时步循环）
for it = 2:ns
    is = it - 1;

    % 计算公式 (82) 右端：P_r*z_{n-1} + sum_k C_k * F_{mn}^{(k)}
    z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
    % z(:,1) = z_n, z(:,2) = A_r*z_n, ..., z(:,p+1) = A_r^p*z_n（见公式 (90)/(91)）
    % prcoe' = [p_{r0},...,p_{rM}]^T → sum p_{ri}*A_r^i*z_{n-1}
    % Pf*ft(is,:)' = sum_k C_k * F^(k)（外力积分项）

    % 然后反向子步求解（公式 (85)/(86)）：
    for ip = p:-1:1
        z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
        z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    end
end
```

**代码解读**：公式 (82) [ \mathbf{A}_r^M\mathbf{z}_n = \mathbf{P}_r\mathbf{z}_{n-1} + \sum_{k=0}^{p_t}\mathbf{C}_k\tilde{\mathbf{F}}_{mn}^{(k)} ] 在代码中被分两步实现：首先 `z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)'` 计算公式 (84)（右端 [ \mathbf{z}^{(0)} ]），利用上一步存储的 [ \mathbf{A}_r^k\mathbf{z}_{n-1} ] 避免显式矩阵乘法（公式 (91)）；然后 `for ip = p:-1:1` 循环反向求解 [ M ] 个线性方程组（公式 (85)/(86)），相当于在无需显式构造 [ \mathbf{A}_r^{-M} ] 的情况下完成 [ \mathbf{A}_r^{-M}\mathbf{z}^{(0)} = \mathbf{z}_n ] 的等价计算。

**对应公式 (82)**：

\[
\mathbf{A}_r^M \mathbf{z}_n = \mathbf{P}_r \mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{C}_k \tilde{\mathbf{F}}_{mn}^{(k)} \tag{82}
\]

[ \mathbf{A}_r ] 基下的时步推进方程，右端由上一步状态 [ \mathbf{z}_{n-1} ] 和外力多项式展开系数决定。

---

### 公式 (83)：辅助变量 z^(k) 的递推定义

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L53-L56)

```matlab
% src/TimeSolverRExpn.m，第 53–56 行（初始化时步 n=1 的辅助变量）
% 初始化：z(:,1) = z_0（初始状态），z(:,ii+1) = A_r*z(:,ii)
for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
    % = (rI-A)*z(:,ii) = A_r * z(:,ii)（公式 (83)/(90)）
end
% 此处 z(:,1)=z_{n-1}, z(:,2)=A_r*z_{n-1}, ..., z(:,p+1)=A_r^p*z_{n-1}
```

**代码解读**：公式 (83) [ \mathbf{z}^{(k)} = \mathbf{A}_r\mathbf{z}^{(k-1)} ] 定义了辅助变量的正向递推关系。在代码初始化中，`for ii = 1:p` 循环计算 [ \mathbf{z}^{(k)} = \mathbf{A}_r^k\mathbf{z}_{n-1} ]（[ k=1,\ldots,M ]），其中每次 `r*z(:,ii) - SolverPadeAx(...)` 实现一次 [ \mathbf{A}_r ] 的作用（利用 `SolverPadeAx` 计算 [ \mathbf{A}\mathbf{z} ]）。这些辅助变量被存储在 `z` 的各列中，供下一时步的公式 (91) 使用。

**对应公式 (83)**：

\[
\mathbf{z}^{(k)} = \mathbf{A}_r \mathbf{z}^{(k-1)}; \quad (k=1,2,\ldots,M) \tag{83}
\]

辅助变量 [ \mathbf{z}^{(k)} ] 的正向递推定义，每次施加一个 [ \mathbf{A}_r ] 算符，无需显式构造矩阵 [ \mathbf{A}_r^k ]。

---

### 公式 (84)：辅助变量 z^(0) 的定义（右端项）

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L68-L69)

```matlab
% src/TimeSolverRExpn.m，第 68–69 行
z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
% z(:,p+1) 就是 z^(0) = P_r*z_{n-1} + sum_k C_k * F^(k)（公式 (84)）
% (z(:,1:p+1))*prcoe' = p_{r0}*z_{n-1} + p_{r1}*A_r*z_{n-1} + ... = P_r*z_{n-1}
% Pf*ft(is,:)' = sum_k C_k * F^(k)
```

**代码解读**：公式 (84) 定义 [ \mathbf{z}^{(0)} = \mathbf{P}_r\mathbf{z}_{n-1} + \sum_{k=0}^{p_t}\mathbf{C}_k\tilde{\mathbf{F}}_{mn}^{(k)} ]，即公式 (82) 的右端项。代码中 `z(:,p+1)` 被赋值为此右端项：`(z(:,1:p+1))*prcoe'` 利用已存储的 [ \mathbf{A}_r^k\mathbf{z}_{n-1} ]（[ k=0,\ldots,M ]，存于 `z` 各列）按 [ \mathbf{A}_r ] 基系数 `prcoe` 加权求和，实现公式 (91) 的无矩阵-向量乘法计算；`Pf*ft(is,:)'` 提供外力积分项。整个计算无需任何大型矩阵乘法，仅用已有的辅助向量进行线性组合。

**对应公式 (84)**：

\[
\mathbf{z}^{(0)} = \mathbf{P}_r \mathbf{z}_{n-1} + \sum_{k=0}^{p_t} \mathbf{C}_k \tilde{\mathbf{F}}_{mn}^{(k)} \tag{84}
\]

辅助变量 [ \mathbf{z}^{(0)} ] 等于公式 (82) 的右端项，包含前一时步状态的 [ \mathbf{P}_r ] 加权线性组合与外力积分贡献。

---

### 公式 (85)：子步递推方程组

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L69-L73)

```matlab
% src/TimeSolverRExpn.m，第 69–73 行（M 步子步循环，实现公式 (85)）
for ip = p:-1:1
    % 求解 A_r * z^(ip) = z^(ip-1)（即 (rI-A)*z^(ip) = z^(ip-1)，公式 (86)）
    z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
    % z^(ip)_1（上半）：求解 (r²M+rC+K)*z1 = rM*z1^(ip-1) - K*z2^(ip-1)（公式 (88)）
    z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    % z^(ip)_2（下半）：z2^(ip) = (z1^(ip) + z2^(ip-1))/r（公式 (89)）
end
% 循环结束后，z(:,1) = z_n（新时步状态）
```

**代码解读**：公式 (85) 将公式 (82) 改写为 [ M ] 个递推方程 [ \mathbf{A}_r\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)} ]（[ k=1,\ldots,M ]），逐步求解。代码中 `for ip = p:-1:1` 从第 [ M ] 步反向到第 1 步求解（因为 `z(:,p+1)` 是 [ \mathbf{z}^{(0)} ]，`z(:,1)` 最终存 [ \mathbf{z}_n ]）：每次循环体先用公式 (88) 求解上半分量 [ \mathbf{z}_1^{(k)} ]（通过预分解的 `dKd`），再用公式 (89) 更新下半分量 [ \mathbf{z}_2^{(k)} ]。所有子步使用**相同的**等效刚度矩阵 `dKd`，这是算法的核心优势。

**对应公式 (85)**：

\[
\begin{aligned}
\mathbf{A}_r \mathbf{z}^{(1)} &= \mathbf{z}^{(0)} \\
\mathbf{A}_r \mathbf{z}^{(2)} &= \mathbf{z}^{(1)} \\
&\cdots \\
\mathbf{A}_r \mathbf{z}^{(k)} &= \mathbf{z}^{(k-1)} \\
&\cdots \\
\mathbf{A}_r \mathbf{z}_n &= \mathbf{z}^{(M-1)}
\end{aligned} \tag{85}
\]

[ M ] 个递推方程，从已知 [ \mathbf{z}^{(0)} ] 逐步求解到 [ \mathbf{z}_n ]，每步系数矩阵 [ \mathbf{A}_r ] 相同（同一等效刚度矩阵），避免重复分解。

---

### 公式 (86)：单子步的隐式方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L70-L72)

```matlab
% src/TimeSolverRExpn.m，第 70–72 行（公式 (86) 的实现）
z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
% 这两行共同实现 (rI-A)*z^(ip) = z^(ip+1)（注意代码中 ip+1 对应 k-1，ip 对应 k）
% 等价于求解：A_r * z^(k) = z^(k-1)（公式 (86)）
```

**代码解读**：公式 (86) 给出单子步隐式方程 [ (r\mathbf{I}-\mathbf{A})\mathbf{z}^{(k)} = \mathbf{z}^{(k-1)} ]，其中 [ \mathbf{z}^{(k)} ] 为待求、[ \mathbf{z}^{(k-1)} ] 为已知。代码通过将 [ \mathbf{z} ] 分块（上半为 [ \dot{\mathbf{u}}^\circ ]，下半为 [ \mathbf{u} ]），利用矩阵 [ \mathbf{A} ]（公式 (8)）的分块结构，将 [ 2n\times 2n ] 的问题化简为 [ n\times n ] 的求解（公式 (88)）加一个显式更新（公式 (89)），大幅降低计算量。

**对应公式 (86)**：

\[
(r\mathbf{I} - \mathbf{A}) \mathbf{z}^{(k)} = \mathbf{z}^{(k-1)} \tag{86}
\]

单子步的隐式方程，[ \mathbf{A}_r = r\mathbf{I}-\mathbf{A} ] 为所有子步共用的系数矩阵，通过分块可化简为 [ n ] 维线性方程组（公式 (88)）。

---

### 公式 (87)：状态向量 z^(k) 的分块

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L70-L72)

```matlab
% src/TimeSolverRExpn.m，第 70–72 行（分块结构直接体现）
% 上半（1:ndof）：z_1^(k) = dot{u}^circle（无量纲速度）
% 下半（ndof+1:2*ndof）：z_2^(k) = u（位移）
z(1:ndof,ip)          % = z_1^(k)
z(ndof+1:2*ndof,ip+1) % = z_2^(k-1)
z(ndof+1:2*ndof,ip)   % = z_2^(k)
```

**代码解读**：公式 (87) 将 [ \mathbf{z}^{(k)} ] 分为上半 [ \mathbf{z}_1^{(k)} ]（无量纲速度 [ \dot{\mathbf{u}}^\circ ]）和下半 [ \mathbf{z}_2^{(k)} ]（位移 [ \mathbf{u} ]），与公式 (6) 的状态向量定义一致。代码中通过行索引 `1:ndof` 和 `ndof+1:2*ndof` 实现分块访问，直接对应公式 (87) 的分块符号表示，是将 [ 2n\times 2n ] 方程化简为 [ n\times n ] 方程（公式 (88)）的基础。

**对应公式 (87)**：

\[
\mathbf{z}^{(k)} = \begin{Bmatrix} \mathbf{z}_1^{(k)} \\ \mathbf{z}_2^{(k)} \end{Bmatrix}, \quad \mathbf{z}^{(k-1)} = \begin{Bmatrix} \mathbf{z}_1^{(k-1)} \\ \mathbf{z}_2^{(k-1)} \end{Bmatrix} \tag{87}
\]

状态向量的分块表示，上半为无量纲速度，下半为位移，直接对应代码中的行索引分割。

---

### 公式 (88)：子步上半分量 z_1^(k) 的方程

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L40-L42) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L70-L71)

```matlab
% src/TimeSolverRExpn.m，第 40–42 行（预分解等效刚度矩阵，公式 (88) 的左端）
Kd = sparse((r*r)*M + r*C + K);   % r²M + r(ΔtC) + Δt²K（对应 r²M+rΔtC+Δt²K）
dKd = decomposition(Kd);          % LU 分解，一次分解用于所有子步

% src/TimeSolverRExpn.m，第 70–71 行（公式 (88) 的实际求解）
z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
% 右端：r*M*z1^{k-1} - K*z2^{k-1}（与公式 (88) 右端 rM*z1^{k-1} - Δt²K*z2^{k-1} 对应）
% 注意代码中 M,C,K 已预乘 Δt 幂次（见公式 (5) 的缩放），故 r²M+rC+K = Δt²(r²/Δt²·M+r/Δt·C+K)
```

**代码解读**：公式 (88) 是子步隐式方程的关键：[ (r^2\mathbf{M}+r\Delta t\mathbf{C}+\Delta t^2\mathbf{K})\mathbf{z}_1^{(k)} = r\mathbf{M}\mathbf{z}_1^{(k-1)} - \Delta t^2\mathbf{K}\mathbf{z}_2^{(k-1)} ]。代码中 `M,C,K` 已在函数入口预乘 [ \Delta t ] 幂次（`K = dt²K`，`C = dt*C`，`M = M`），故 `Kd = r²*M + r*C + K` 对应 [ r^2\mathbf{M}+r\Delta t\mathbf{C}+\Delta t^2\mathbf{K} ]，右端 `r*(M*z1^{k-1}) - K*z2^{k-1}` 对应 [ r\mathbf{M}\mathbf{z}_1^{(k-1)} - \Delta t^2\mathbf{K}\mathbf{z}_2^{(k-1)} ]，与公式 (88) 完全一致。`dKd = decomposition(Kd)` 只分解一次（对所有子步有效），是算法效率的关键。

**对应公式 (88)**：

\[
(r^2\mathbf{M} + r\Delta t\mathbf{C} + \Delta t^2\mathbf{K}) \mathbf{z}_1^{(k)} = r\mathbf{M}\mathbf{z}_1^{(k-1)} - \Delta t^2\mathbf{K}\mathbf{z}_2^{(k-1)} \tag{88}
\]

子步 [ k ] 上半分量的 [ n\times n ] 线性方程，等效刚度矩阵 [ r^2\mathbf{M}+r\Delta t\mathbf{C}+\Delta t^2\mathbf{K} ] 对所有子步相同（只需一次矩阵分解）。

---

### 公式 (89)：子步下半分量 z_2^(k) 的显式更新

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L72)

```matlab
% src/TimeSolverRExpn.m，第 72 行（公式 (89) 的直接实现）
z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
% = (z1^(k) + z2^(k-1)) / r（公式 (89)）
```

**代码解读**：公式 (89) 给出下半分量 [ \mathbf{z}_2^{(k)} = (\mathbf{z}_1^{(k)} + \mathbf{z}_2^{(k-1)})/r ]，是一个纯显式的向量运算（无需求解线性方程组）。代码第 72 行直接实现此公式：`(z(1:ndof,ip) + z(ndof+1:2*ndof,ip+1))/r`，其中 `z(1:ndof,ip)` 是刚求得的 [ \mathbf{z}_1^{(k)} ]，`z(ndof+1:2*ndof,ip+1)` 是已知的 [ \mathbf{z}_2^{(k-1)} ]。这一显式更新与公式 (88) 的隐式方程结合，将 [ 2n\times 2n ] 问题化简为 [ n\times n ] 问题，计算量减半。

**对应公式 (89)**：

\[
\mathbf{z}_2^{(k)} = \frac{1}{r} \left(\mathbf{z}_1^{(k)} + \mathbf{z}_2^{(k-1)}\right) \tag{89}
\]

下半分量的显式更新公式，由上半分量 [ \mathbf{z}_1^{(k)} ] 和前一子步的下半分量 [ \mathbf{z}_2^{(k-1)} ] 直接计算，无需求解线性方程组。

---

### 公式 (90)：辅助变量与 A_r 幂次的关系

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L53-L56) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L68-L69)

```matlab
% 子步求解完成后：
% z(:,1) = z_n（新状态，即 A_r^0*z_n）
% z(:,2) = A_r^{M-1}*z_n（从子步求解过程中得到）
% ...
% z(:,p+1) = A_r^M 作用前的右端项

% 在初始化时（公式 (90) 的正向计算）：
for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
    % z(:,ii+1) = A_r * z(:,ii) = A_r^ii * z_{n-1}
end
% 因此 z(:,k+1) = A_r^k * z_{n-1}（对应公式 (90) 的正向版本）
```

**代码解读**：公式 (90) [ \mathbf{A}_r^k\mathbf{z}_n = \mathbf{z}^{(M-k)} ] 给出子步求解后辅助变量的含义：[ \mathbf{z}^{(M-k)} ] 等于 [ \mathbf{A}_r^k ] 作用于新状态 [ \mathbf{z}_n ] 的结果。代码中子步循环结束后，`z(:,1) = z_n`，`z(:,2) = A_r*z_n`（因为最后一步 `ip=1` 求得的是 [ \mathbf{z}^{(1)} = \mathbf{A}_r^{-1}\mathbf{z}^{(0)} ] 的反向，而存储的是 [ \mathbf{z}^{(M-1)} ]）。这些辅助变量在下一时步被 `prcoe` 加权求和，实现公式 (91) 的无矩阵乘法计算。

**对应公式 (90)**：

\[
\mathbf{A}_r^k \mathbf{z}_n = \mathbf{z}^{(M-k)} \tag{90}
\]

子步求解的中间结果 [ \mathbf{z}^{(M-k)} ] 等于 [ \mathbf{A}_r^k\mathbf{z}_n ]，可用于下一时步的 [ \mathbf{P}_r\mathbf{z}_n ] 计算（公式 (91)）。

---

### 公式 (91)：P_r z_n 的无矩阵-向量乘法计算

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L68-L69)

```matlab
% src/TimeSolverRExpn.m，第 68–69 行（公式 (91) 的直接实现）
z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
% (z(:,1:p+1))*prcoe' = p_{r0}*z(:,1) + p_{r1}*z(:,2) + ... + p_{rM}*z(:,p+1)
%                     = p_{r0}*z_n + p_{r1}*(A_r*z_n) + ... + p_{rM}*(A_r^M*z_n)
%                     = P_r * z_n（公式 (91)）
% z(:,k) 存储 A_r^{k-1} * z_n（来自上一时步子步求解后的存储）
```

**代码解读**：公式 (91) [ \mathbf{P}_r\mathbf{z}_n = p_{r0}\mathbf{z}_n + p_{r1}\mathbf{z}^{(M-1)} + \cdots + p_{rM}\mathbf{z}^{(0)} ] 利用公式 (90) 存储的辅助变量，无需任何矩阵-向量乘法（仅需向量的线性组合），大幅降低每时步计算量。代码中 `(z(:,1:p+1))*prcoe'` 是矩阵（其列为各 [ \mathbf{A}_r^k\mathbf{z}_{n-1} ]）与系数向量 `prcoe` 的点积，实现公式 (91) 的 [ \mathbf{P}_r\mathbf{z}_{n-1} ] 计算。这正是论文中所说「不需要矩阵-向量乘积（No matrix-vector product is needed）」的代码体现。

**对应公式 (91)**：

\[
\mathbf{P}_r \mathbf{z}_n = p_{r0}\mathbf{z}_n + p_{r1}\mathbf{z}^{(M-1)} + \cdots + p_{rM}\mathbf{z}^{(0)} \tag{91}
\]

利用存储的辅助变量 [ \mathbf{z}^{(M-k)} = \mathbf{A}_r^k\mathbf{z}_n ] 计算 [ \mathbf{P}_r\mathbf{z}_n ]，无需显式矩阵-向量乘法，每时步仅需 [ M+1 ] 次向量线性组合。

---

### 公式 (92)：时步末端加速度的计算式

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L58-L63) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L77-L78)

```matlab
% src/TimeSolverRExpn.m，第 58–63 行（初始时刻加速度）
fn = fHist(0); % 外力 f(t=0)
acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);
% = A*z_n 的上半（即 dot{z}_n 的上半部分）+ F_n(s=1)（见公式 (93)）

% src/TimeSolverRExpn.m，第 77–78 行（一般时步加速度）
fn = fHist(is*dt);  % f(t_{n-1})
acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);
% = [r*z_n - z^{(M-1)}] + F_n(s=1)（公式 (93)）
```

**代码解读**：公式 (92) [ \dot{\mathbf{z}}_n = \mathbf{A}\mathbf{z}_n + \mathbf{F}_n(s=1) ] 给出时步末端 [ t=t_n ] 处状态向量的时间导数，其上半部分为 [ \ddot{\mathbf{u}}^{\circ\circ}_n = \Delta t^2\ddot{\mathbf{u}}_n ]。代码中通过公式 (93) 的变换避免了直接计算 [ \mathbf{A}\mathbf{z}_n ]（该乘法代价较高），而是用已有辅助变量表达。`fn*Fb(pDOF)` 对应 [ \mathbf{F}_n(s=1) ] 的上半分量（已经过 [ \mathbf{M}^{-1} ] 作用）。

**对应公式 (92)**：

\[
\dot{\mathbf{z}}_n = \mathbf{A}\mathbf{z}_n + \mathbf{F}_n(s=1) \tag{92}
\]

时步末端一阶 ODE（公式 (7)）的直接应用，[ \dot{\mathbf{z}}_n ] 的上半即 [ \ddot{\mathbf{u}}^{\circ\circ}_n ]（经 [ \Delta t^2 ] 缩放后为物理加速度）。

---

### 公式 (93)：加速度的向量运算实现

**对应代码文件**：[`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L63) 与 [`src/TimeSolverRExpn.m`](https://github.com/xingzhiyuan1229/CompsiteTimeIntegration/blob/main/src/TimeSolverRExpn.m#L78)

```matlab
% src/TimeSolverRExpn.m，第 63 行和第 78 行（公式 (93) 的直接实现）
acc(it,:) = r*z(pDOF,1) - z(pDOF,2) + fn*Fb(pDOF);
% r*z_n - z^{(M-1)} = A*z_n（利用公式 (85) 最后一行：A_r*z_n = z^{(M-1)}）
% 即：A*z_n = r*z_n - A_r*z_n = r*z_n - z^{(M-1)}
% fn*Fb(pDOF) = F_n(s=1) 的上半分量（Δt²*M^{-1}*f(t_n)）
% 上半 = ddot{u}°°_n = Δt²*ddot{u}_n → acc 后面除以 dt^2 得物理加速度
acc = acc/(dt*dt);  % 恢复物理加速度（第 83 行）
```

**代码解读**：公式 (93) [ \dot{\mathbf{z}}_n = r\mathbf{z}_n - \mathbf{z}^{(M-1)} + \mathbf{F}_n(s=1) ] 利用公式 (80)/(85) 的关系（[ \mathbf{A}_r\mathbf{z}_n = \mathbf{z}^{(M-1)} ]）将矩阵-向量乘法 [ \mathbf{A}\mathbf{z}_n = (r\mathbf{I}-\mathbf{A}_r)\mathbf{z}_n = r\mathbf{z}_n - \mathbf{z}^{(M-1)} ] 转化为纯向量运算。代码中 `r*z(pDOF,1) - z(pDOF,2)` 正是 [ r\mathbf{z}_n - \mathbf{z}^{(M-1)} ] 的实现（`z(:,1) = z_n`，`z(:,2) = \mathbf{z}^{(M-1)}`），与公式 (93) 的向量运算形式完全一致，彻底避免了大型矩阵 [ \mathbf{A} ] 与 [ \mathbf{z}_n ] 的显式乘法。

**对应公式 (93)**：

\[
\dot{\mathbf{z}}_n = r\mathbf{z}_n - \mathbf{z}^{(M-1)} + \mathbf{F}_n(s=1) \tag{93}
\]

通过 [ \mathbf{A}\mathbf{z}_n = r\mathbf{z}_n - \mathbf{A}_r\mathbf{z}_n = r\mathbf{z}_n - \mathbf{z}^{(M-1)} ] 将矩阵-向量乘法替换为纯向量运算，[ \dot{\mathbf{z}}_n ] 的上半分量（[ \ddot{\mathbf{u}}^{\circ\circ}_n ]）除以 [ \Delta t^2 ] 即得物理加速度。

---

*本文档（Part 4）覆盖公式 (66)–(93)，共 28 个公式，对应论文第 3.2 节 (M+1)-格式（公式 (66)–(79)）及第 4 节算法实现（公式 (80)–(93)）。续篇见 [公式解读_Part5.md](./公式解读_Part5.md)，涵盖公式 (94)–(99)（第 4 节续及第 5 节算例）。*
