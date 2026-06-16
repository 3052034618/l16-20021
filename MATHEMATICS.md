# 曲线与插值引擎 - 数学原理说明

## 目录
1. [线性插值 (Linear Interpolation)
2. [三次样条插值 (Cubic Spline)
3. [三对角线性方程组求解]
4. [贝塞尔曲线与德卡斯特里奥算法]
5. [Catmull-Rom 样条]
6. [弧长均匀采样]

---

## 1. 线性插值 (Linear Interpolation)

### 基本原理
给定两个点 $P_0$ 和 $P_1$，线性插值在两点之间构造一条直线段：

$$
L(t) = (1-t)P_0 + tP_1, \quad t \in [0,1]
$$

也可写成：
$$
L(t) = P_0 + t(P_1 - P_0)
$$

### 多个控制点序列的分段线性插值：
对于 $n$ 个控制点 $P_0, P_1, \ldots, P_{n-1}$，在每个相邻点对之间使用线性插值，每段参数归一化到 $t \in [0, 1]$。

实现代码: [linear.ts](file:///d:/trae-bz/TraeProjects/20021/src/linear.ts)

---

## 2. 三次样条插值 (Cubic Spline)

### 2.1 每段三次多项式

给定 $n+1$ 个控制点 $(x_0,y_0), \ldots, (x_n,y_n)$，三次样条在每两个相邻点 $[x_i, x_{i+1}]$ 之间定义一段三次多项式 $S_i(x)$：

$$
S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
$$

或使用更一般的参数形式（对任意曲线）：

$$
S_i(t) = a_i + b_i t + c_i t^2 + d_i t^3, \quad t \in [0,1]
$$

其中 $t = \frac{x - x_i}{h_i}$，$h_i = x_{i+1} - x_i$。

### 2.2 连续性约束

为了使整条曲线整体光滑，需要相邻两段在连接点处满足：

1. **位置连续（C⁰）：$S_i(1) = S_{i+1}(0) = y_{i+1}$

2. **一阶导连续（C¹）：$S_i'(1) = S_{i+1}'(0)$

3. **二阶导连续（C²）：$S_i''(1) = S_{i+1}''(0)$

### 2.3 使用弯矩（M-值）公式

设 $M_i = S_i''(0)$，即每段起点的二阶导数。对第 $i$ 段的二阶导数线性变化：

$$
S_i''(t) = M_i + (M_{i+1} - M_i)t
$$

积分两次得到：

$$
S_i(t) = \frac{M_i (1-t)^3 + M_{i+1} t^3}{6 h_i} + \cdots
$$

### 2.4 边界条件

需要两个额外条件确定 $M_0$ 和 $M_n$：

- **自然样条 (Natural Spline）：$M_0 = 0, \quad M_n = 0

- **固定斜率 (Clamped Spline）：指定端点一阶导数 $S_0'(0) = y_0', \quad S_{n-1}'(1) = y_{n-1}'$

- **Not-a-Knot**：在 $x_1$ 和 $x_{n-1}$ 处三阶导连续

### 2.5 为什么归结为三对角方程组

对每个内部点 $i=1,\ldots,n-1$，由 C¹ 连续条件可推导得：

$$
h_{i-1} M_{i-1} + 2(h_{i-1}+h_i) M_i + h_i M_{i+1} = 6\left(\frac{y_{i+1}-y_i}{h_i} - \frac{y_i-y_{i-1}}{h_{i-1}}\right)
$$

这形成了系数矩阵为三对角带状的线性方程组。

实现代码: [cubic_spline.ts](file:///d:/trae-bz/TraeProjects/20021/src/cubic_spline.ts)

---

## 3. 三对角线性方程组求解

### 3.1 问题形式

求解形如：

$$
\begin{bmatrix}
b_0 & c_0 & & & \\
a_1 & b_1 & c_1 & & \\
& a_2 & b_2 & c_2 & \\
& & \ddots & \ddots & \ddots & \\
& & & a_{n-1} & b_{n-1} & c_{n-1} \\
& & & & a_n & b_n
\end{bmatrix}
\begin{bmatrix}
x_0 \\ x_1 \\ x_2 \\ \vdots \\ x_{n-1} \\ x_n
\end{bmatrix}
=
\begin{bmatrix}
d_0 \\ d_1 \\ d_2 \\ \vdots \\ d_{n-1} \\ d_n
\end{bmatrix}
$$

### 3.2 Thomas 算法（追赶法）

**第 1 步 - 前向消去（追）：
$$
c'_0 = c_0 / b_0 \\
d'_0 = d_0 / b_0 \\
c'_i = c_i / (b_i - a_i c'_{i-1}) \\
d'_i = (d_i - a_i d'_{i-1}) / (b_i - a_i c'_{i-1})
$$

**第 2 步 - 回代（赶）：
$$
x_n = d'_n \\
x_i = d'_i - c'_i x_{i+1}
$$

算法复杂度仅为 $O(n)$，远优于一般的高斯消元的 $O(n^3)$。

实现代码: [tridiagonal.ts](file:///d:/trae-bz/TraeProjects/20021/src/tridiagonal.ts)

---

## 4. 贝塞尔曲线 (Bezier Curve)

### 4.1 定义

给定 $n+1$ 个控制点 $P_0, \ldots, P_n$，$n$ 次贝塞尔曲线为：

$$
\mathbf{B}(t) = \sum_{i=0}^{n} B_{n,i}(t) P_i
$$

其中 $B_{n,i}(t)$ 为 Bernstein 基函数：

$$
B_{n,i}(t) = \binom{n}{i} t^i (1-t)^{n-i}
$$

Bernstein 多项式的性质：
- 非负性：$B_{n,i}(t) \ge 0$ 对 $t \in [0,1]$
- 单位分解：$\sum_{i=0}^{n} B_{n,i}(t) = 1$
- 对称性：$B_{n,i}(t) = B_{n,n-i}(1-t)$

### 4.2 德卡斯特里奥 (de Casteljau) 算法

通过反复线性插值求出曲线上的点：

**递归公式：**
$$
P_i^{(0)}(t) = P_i \\
P_i^{(k)}(t) = (1-t) P_i^{(k-1)}(t) + t P_{i+1}^{(k-1)}(t)
$$

最终结果：
$$
\mathbf{B}(t) = P_0^{(n)}(t)
$$

几何意义：从控制点开始，逐层做线性插值，每层减少一个点，最终剩下的一个点就是曲线上的点。

### 4.3 曲线细分（一分为二）

使用德卡斯特里奥算法在 $t = 0.5$ 时产生的整个三角形表格：

- **左半曲线控制点**：三角形的左边缘 $P_0^{(0)}, P_0^{(1)}, \ldots, P_0^{(n)}$
- **右半曲线控制点**：三角形的右边缘 $P_0^{(n)}, P_1^{(n-1)}, \ldots, P_n^{(0)}$

这样细分后的两段贝塞尔曲线拼接处保持 C⁰, C¹, C² 连续。

### 4.4 导数公式

贝塞尔曲线的导数仍是贝塞尔曲线（次数降 1）：

$$
\mathbf{B}'(t) = n \sum_{i=0}^{n-1} B_{n-1,i}(t) (P_{i+1} - P_i)
$$

实现代码: [bezier.ts](file:///d:/trae-bz/TraeProjects/20021/src/bezier.ts)

---

## 5. Catmull-Rom 样条

### 5.1 定义

给定控制点序列 $P_0, P_1, \ldots, P_n$，Catmull-Rom 样条在每段 $[P_i, P_{i+1}]$ 上使用相邻的四个点 $P_{i-1}, P_i, P_{i+1}, P_{i+2}$ 构造一段三次曲线。

**标准公式：

$$
\mathbf{CR}_i(t) = \tfrac{1}{2} \begin{bmatrix} 1 & t & t^2 & t^3 \end{bmatrix}
\begin{bmatrix}
0 & 2 & 0 & 0 \\
-1 & 0 & 1 & 0 \\
2 & -5 & 4 & -1 \\
-1 & 3 & -3 & 1
\end{bmatrix}
\begin{bmatrix}
P_{i-1} \\
P_i \\
P_{i+1} \\
P_{i+2}
\end{bmatrix}
$$

展开后：

$$
\mathbf{CR}_i(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3
$$

其中：
$$
a_0 = P_i \\
a_1 = \alpha (P_{i+1} - P_{i-1}) \\
a_2 = 2\alpha P_{i-1} + (\alpha-3) P_i + (3-2\alpha) P_{i+1} - \alpha P_{i+2} \\
a_3 = -\alpha P_{i-1} + (2-\alpha) P_i + (\alpha-2) P_{i+1} + \alpha P_{i+2}
$$

$\alpha = 0.5$ 为标准张力参数。

### 5.2 Catmull-Rom vs 贝塞尔的区别

| 特性 | Catmull-Rom | 贝塞尔 (Bézier) |
|------|-------------|-----------------|
| **是否过控制点 | ✅ 经过所有控制点 | ❌ 通常不过内部控制点（只过端点） |
| 连续阶 | C¹ 连续 | C∞（无限可导） |
| 控制点影响范围 | 局部（相邻4点） | 全局（所有控制点） |
| 控制点数量 | 任意多（分段） | 单段 n+1 个控制点决定 n 次曲线 |
| 凸包性 | 不保证 | ✅ 始终在凸包内 |
| 用途 | 相机路径、动画轨迹 | 字体设计、矢量图形 |

### 5.3 转换为贝塞尔形式

Catmull-Rom 每段可等价转换为三次贝塞尔：

$$
\mathbf{B}_0 = P_i \\
\mathbf{B}_1 = P_i + \frac{\alpha}{3}(P_{i+1} - P_{i-1}) \\
\mathbf{B}_2 = P_{i+1} - \frac{\alpha}{3}(P_{i+2} - P_i) \\
\mathbf{B}_3 = P_{i+1}
$$

实现代码: [catmull_rom.ts](file:///d:/trae-bz/TraeProjects/20021/src/catmull_rom.ts)

---

## 6. 弧长均匀采样

### 6.1 为什么需要按弧长采样？

参数均匀采样（$t$ 等间隔）得到的点，在曲线上的**空间距离往往不均匀。曲线弯曲部分点密集，伸直部分点稀疏。

动画中需要物体沿曲线**匀速运动**，必须按弧长等间隔采样。

### 6.2 弧长计算

参数曲线 $\mathbf{C}(t)$ 从 $t=a$ 到 $t=b$ 的弧长为：

$$
L(a,b) = \int_a^b |\mathbf{C}'(t)| \, dt
$$

总弧长：
$$
L_{\text{total}} = L(0,1) = \int_0^1 \sqrt{(x'(t))^2 + (y'(t))^2} \, dt
$$

### 6.3 数值方法

由于弧长函数 $s(t) = L(0,t)$ 通常没有解析反函数，使用数值方法：

**步骤 1 — 建立查找表**
离散采样 $N$ 个 $t$ 值，记录累积弧长 $(t_i, s_i)$：

$$
s_i = \sum_{k=1}^{i} |\mathbf{C}(t_k) - \mathbf{C}(t_{k-1})|
$$

**步骤 2 — 二分查找**
给定目标弧长 $s$，在表中用二分查找定位区间，线性插值求对应 $t$。

**步骤 3 — Newton-Raphson 精化（可选）
$$
t_{n+1} = t_n + \frac{s - s(t_n)}{|\mathbf{C}'(t_n)|}
$$

### 6.4 算法流程

```
输入: 目标弧长 s
输出: 对应参数 t 和点 C(t)

1. 二分查找弧长表，找到 s_i ≤ s ≤ s_{i+1}
2. 线性插值得初始 t ≈ t_i + (s-s_i)/(s_{i+1}-s_i) * (t_{i+1}-t_i)
3. 可选: Newton 迭代精化 t
4. 返回 C(t)
```

实现代码: [arc_length.ts](file:///d:/trae-bz/TraeProjects/20021/src/arc_length.ts)

---

## 项目文件结构

```
src/
├── types.ts              # 基础类型、工具函数（Point、lerp、距离等）
├── index.ts              # 统一导出入口
├── linear.ts             # 分段线性插值
├── tridiagonal.ts        # 三对角方程组 Thomas 算法
├── cubic_spline.ts        # 三次样条插值
├── bezier.ts             # 贝塞尔曲线（de Casteljau、细分、升阶、自适应细分
├── catmull_rom.ts        # Catmull-Rom 样条（过控制点）
└── arc_length.ts         # 弧长计算与均匀采样重参数化

demo/demo.ts              # 各类曲线功能演示
tests/test_curves.ts        # 单元测试（70项全部通过）
```

## 运行方式

```bash
# 安装依赖
npm install

# 运行单元测试
npm test

# 运行完整演示
npm run demo

# 编译 TypeScript
npm run build
```

## API 速查

```typescript
// 线性插值
const lin = new LinearInterpolation(points);
lin.evaluate(0.5);          // 参数求值
lin.sampleByArcLength(20);     // 弧长均匀采样

// 三次样条
const spline = new CubicSpline(points, { endCondition: 'natural' });
spline.checkContinuity();      // 检查 C0/C1/C2 连续性
spline.getSegments();          // 获取各段三次多项式系数

// 贝塞尔曲线
const bezier = new BezierCurve(controlPoints);
bezier.evaluate(0.3);        // de Casteljau 递归求值
bezier.evaluateBernstein(0.3);   // Bernstein 基函数求值
const { left, right } = bezier.subdivide(0.5);  // 曲线一分为二
bezier.adaptiveSubdivision(0.1); // 自适应细分

// Catmull-Rom 样条（过所有点！）
const cr = new CatmullRomSpline(points, { tension: 0.5, closed: false });
cr.verifyInterpolation();      // 验证确实经过控制点
cr.toBezierSegments();        // 转换为贝塞尔段列表

// 弧长采样（任意曲线都支持）
curve.arcLength(0, 0.5);     // 部分弧长
curve.sampleByArcLength(50);  // 按弧长采 50 个等距点
```
