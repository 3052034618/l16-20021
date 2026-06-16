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

## 7. 非均匀 Knots 与时间戳参数化

### 7.1 为什么需要非均匀参数化？

默认均匀分段假设每段控制点间隔时间相同。但在动画中，关键帧往往间隔不均匀（如 0s, 1s, 3s, 7s），如果仍按均匀分段插值，动画的速度会与真实时间脱节。

### 7.2 Knot 类型

- **Uniform (均匀)**：$t_i = i/n$，适用于控制点等时间间隔
- **Chordal (弦长)**：$t_i \propto \sum_{k=1}^{i} |P_k-P_{k-1}|$，按相邻点距离分配参数
- **Centripetal (向心)**：$t_i \propto \sum \sqrt{|P_k-P_{k-1}|}$，缓解曲线自交
- **自定义时间戳**：直接用时间戳数组 `[0, 1, 3, 7]`，自动归一化到 [0,1]

### 7.3 非均匀分段的样条推导

设全局参数 $t \in [0,1]$，段 $i$ 对应 $t \in [k_i, k_{i+1}]$，段宽 $\Delta_i = k_{i+1}-k_i$。定义局部参数 $s = (t - k_i) / \Delta_i \in [0,1]$。

每段多项式 $S_i(s) = a_i + b_i s + c_i s^2 + d_i s^3$，对全局 $t$ 的导数：

$$
\frac{dS}{dt} = \frac{dS}{ds} \cdot \frac{ds}{dt} = \frac{b_i + 2 c_i s + 3 d_i s^2}{\Delta_i}
$$

$$
\frac{d^2 S}{dt^2} = \frac{2 c_i + 6 d_i s}{\Delta_i^2}
$$

在段边界处 C¹ 和 C² 连续的约束方程推导后，仍归结为三对角方程组：

$$
\Delta_{i-1} M_{i-1} + 2(\Delta_{i-1} + \Delta_i) M_i + \Delta_i M_{i+1} = 6\left(\frac{y_{i+1}-y_i}{\Delta_i} - \frac{y_i-y_{i-1}}{\Delta_{i-1}}\right)
$$

系数为：
$$
a_i = y_i, \quad b_i = (y_{i+1}-y_i) - \Delta_i^2 (2 M_i + M_{i+1})/6
$$
$$
c_i = M_i \Delta_i^2 / 2, \quad d_i = (M_{i+1} - M_i) \Delta_i^2 / 6
$$

实现代码: [cubic_spline.ts](file:///d:/trae-bz/TraeProjects/20021/src/cubic_spline.ts) · [catmull_rom.ts](file:///d:/trae-bz/TraeProjects/20021/src/catmull_rom.ts)

---

## 8. 曲线分析

### 8.1 曲率

参数曲线 $\mathbf{r}(t) = (x(t), y(t))$ 的曲率为：

$$
\kappa(t) = \frac{x'(t) y''(t) - y'(t) x''(t)}{\left(x'(t)^2 + y'(t)^2\right)^{3/2}}
$$

曲率半径 $R(t) = 1/|\kappa(t)|$，表示曲线在该点的密切圆半径。

### 8.2 最近点 / 投影参数

给定点 $Q$，寻找 $t^*$ 使 $|\mathbf{r}(t^*) - Q|$ 最小。

**算法：牛顿迭代法**
目标：最小化 $f(t) = |\mathbf{r}(t) - Q|^2$。
$f'(t) = 2(\mathbf{r}(t) - Q) \cdot \mathbf{r}'(t)$，$f''(t) = 2|\mathbf{r}'(t)|^2 + 2(\mathbf{r}(t)-Q) \cdot \mathbf{r}''(t)$

迭代：
$$
t_{n+1} = t_n - \frac{(\mathbf{r}(t_n)-Q) \cdot \mathbf{r}'(t_n)}{|\mathbf{r}'(t_n)|^2 + (\mathbf{r}(t_n)-Q) \cdot \mathbf{r}''(t_n)}
$$

当分母为负（凹点）时退化为梯度下降步长。

实现代码: [curve_analyzer.ts](file:///d:/trae-bz/TraeProjects/20021/src/curve_analyzer.ts)

### 8.3 包围盒

采样曲线 N 个点取 min/max，得到：
$$
\text{BBox} = [\min x, \max x] \times [\min y, \max y]
$$

### 8.4 速度分析

$\mathbf{v}(t) = \mathbf{r}'(t)$ 为速度向量，$|\mathbf{v}(t)|$ 为速率。
加速度 $\mathbf{a}(t) = \mathbf{r}''(t)$。

可用于路径跟随的速度规划、碰撞预估等。

---

## 9. 曲线互导

### 9.1 SVG Path 导出

将曲线转换为 SVG `d` 属性字符串：
- `M x y` — 移动到起点
- `L x y` — 直线段（折线）
- `C x1 y1, x2 y2, x y` — 三次贝塞尔曲线段

### 9.2 Catmull-Rom → 贝塞尔段

每段 CR 的 4 控制点 $P_{i-1}, P_i, P_{i+1}, P_{i+2}$ 转成贝塞尔：
$$
B_0 = P_i, \quad B_1 = P_i + \frac{\alpha}{3}(P_{i+1}-P_{i-1})
$$
$$
B_2 = P_{i+1} - \frac{\alpha}{3}(P_{i+2}-P_i), \quad B_3 = P_{i+1}
$$

### 9.3 采样点 → 贝塞尔拟合（最小二乘）

给定采样点集 $\{Q_k\}$，拟合多段 CubicBezier。

**步骤：**
1. 用弦长参数化估计参数 $u_k$
2. 每段用端点条件和最小二乘求控制点 $B_1, B_2$
3. 若误差超限则在最大误差处递归细分

实现代码: [curve_io.ts](file:///d:/trae-bz/TraeProjects/20021/src/curve_io.ts)

---

## 10. 数值稳定性

### 10.1 常见退化场景

| 场景 | 处理方式 |
|------|---------|
| NaN / Infinity 点 | 自动过滤，保留有效点 |
| 连续重复点 | 去重（距离 < 1e-10） |
| 所有点重合 | 标记 `allCollapsed`，构造器抛出友好错误 |
| 极短段 $|P_{i+1}-P_i| < \epsilon$ | 归一化时用极小值替代，防止除零 |
| 大量控制点 (>10000) | 打印性能警告 |
| 段查找 | 二分查找 $O(\log n)$ 替代线性遍历 |

### 10.2 安全工具函数

- `safeNumber(v, fallback)` — 非有限数字替换
- `safePoint(p, fallback)` — 非有限坐标点替换
- `safeNormalize(v, fallback)` — 零向量替代
- `safeDivide(a, b, fallback)` — 防除零

实现代码: [stability.ts](file:///d:/trae-bz/TraeProjects/20021/src/stability.ts)

---

## 项目文件结构

```
src/
├── types.ts              # 基础类型、工具函数（Point、lerp、距离等）
├── index.ts              # 统一导出入口
├── linear.ts             # 分段线性插值
├── tridiagonal.ts        # 三对角方程组 Thomas 算法
├── cubic_spline.ts       # 三次样条插值（支持非均匀knots/时间戳）
├── bezier.ts             # 贝塞尔曲线（de Casteljau、细分、升阶、自适应细分）
├── catmull_rom.ts        # Catmull-Rom 样条（过控制点，支持时间戳）
├── arc_length.ts         # 弧长计算与均匀采样重参数化
├── stability.ts          # 数值稳定性、输入清洗、包围盒
├── curve_analyzer.ts     # 曲线分析（最近点、曲率、包围盒、速度）
└── curve_io.ts           # 互导（SVG path、贝塞尔拟合、Catmull-Rom转贝塞尔）

demo/demo.ts              # 各类曲线功能演示（11项演示）
tests/test_curves.ts      # 单元测试（119项全部通过）
```

## 运行方式

```bash
# 安装依赖
npm install

# 运行单元测试（119项）
npm test

# 运行完整演示
npm run demo

# 编译 TypeScript
npm run build
```

## API 速查

```typescript
// ============ 1. 各类曲线 ============
// 线性插值
const lin = new LinearInterpolation(points);

// 三次样条（支持时间戳/knots）
const spline = new CubicSpline(points, {
  endCondition: 'natural', // 'natural' | 'clamped' | 'not-a-knot'
  knotType: 'chordal',     // 'uniform' | 'chordal' | 'centripetal'
  timeStamps: [0, 1, 3, 7] // 直接指定关键帧时间
});
spline.sampleByTimestamps([0, 0.5, 1, 1.5, 2]); // 按真实时间采样

// 贝塞尔曲线
const bezier = new BezierCurve(controlPoints);
bezier.evaluate(0.3);
bezier.subdivide(0.5);   // 一分为二
bezier.adaptiveSubdivision(0.1);
const cubic = new CubicBezier(p0, c1, c2, p3);

// Catmull-Rom（经过所有控制点）
const cr = new CatmullRomSpline(points, {
  tension: 0.5, closed: false,
  timeStamps: [0, 0.5, 2, 5] // 支持非均匀时间
});

// ============ 2. 曲线分析 ============
const analyzer = spline.getAnalyzer();
analyzer.getBoundingBox();        // 包围盒
analyzer.curvature(0.5);          // t=0.5 处曲率
analyzer.speed(0.3);              // t=0.3 处速度大小
analyzer.velocity(0.3);           // 速度向量
analyzer.acceleration(0.3);       // 加速度向量
analyzer.tangent(0.3);            // 单位切向量
analyzer.normal(0.3);             // 单位法向量
analyzer.findNearestPoint(target);// 最近点 + 参数 + 距离
analyzer.isPointNearCurve(p, 5);  // 距曲线<5？
analyzer.speedProfile(100);       // 速度、曲率、加速度剖面
analyzer.analyze();               // 完整分析报告

// ============ 3. 曲线互导 ============
CurveIO.toSVGPath(curve);           // 导出 SVG path 字符串
CurveIO.toSVGPath(bezierSegments);  // CubicBezier[] → SVG
CurveIO.curveToFullSVG(curve, { width: 800, height: 600 }); // 完整 SVG 文件
CurveIO.fitCubicBezier(sampledPoints, { maxError: 0.5 });  // 采样点反拟合贝塞尔

// ============ 4. 弧长采样（所有曲线都支持） ============
curve.arcLength(0, 0.5);         // 部分弧长
curve.getTotalLength();          // 总弧长
curve.sampleByArcLength(50);     // 按弧长均匀采 50 点

// ============ 5. 数值稳定性工具 ============
cleanInputPoints(points);        // 去重 + 过滤NaN + 标记退化
autoKnots(n, 'chordal', pts);   // 自动生成 knots
validateKnots(knots, n);         // 验证 knots 合法性
computeBoundingBox(points);      // 包围盒
safePoint(p), safeNumber(v), safeDivide(a, b)
```

