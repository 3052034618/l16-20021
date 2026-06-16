import {
  Point, createPoint, LinearInterpolation, CubicSpline,
  BezierCurve, CubicBezier, CatmullRomSpline, solveTridiagonal,
  bernsteinPolynomial, ArcLengthSampler
} from '../src/index';

function logTitle(title: string): void {
  const separator = '='.repeat(60);
  console.log(`\n${separator}`);
  console.log(`  ${title}`);
  console.log(separator);
}

function logPoint(label: string, p: Point): void {
  console.log(`${label}: (${p.x.toFixed(4)}, ${p.y.toFixed(4)})`);
}

function testLinearInterpolation(): void {
  logTitle('1. 线性插值 (Linear Interpolation)');

  const points: Point[] = [
    createPoint(0, 0),
    createPoint(1, 1),
    createPoint(2, 0.5),
    createPoint(3, 2),
    createPoint(4, 1)
  ];

  console.log('控制点:');
  points.forEach((p, i) => logPoint(`  P${i}`, p));

  const curve = new LinearInterpolation(points);

  console.log('\n在 t=0, 0.25, 0.5, 0.75, 1.0 处求值:');
  [0, 0.25, 0.5, 0.75, 1].forEach(t => {
    logPoint(`  t=${t}`, curve.evaluate(t));
  });

  console.log('\n总弧长:', curve.getTotalLength().toFixed(4));
  console.log('从 t=0.25 到 t=0.75 的弧长:', curve.arcLength(0.25, 0.75).toFixed(4));

  const samples = curve.sampleByArcLength(5);
  console.log('\n按弧长均匀采样 (5个点):');
  samples.forEach((p, i) => logPoint(`  点${i}`, p));
}

function testCubicSpline(): void {
  logTitle('2. 三次样条插值 (Cubic Spline)');

  const points: Point[] = [
    createPoint(0, 0),
    createPoint(1, 2),
    createPoint(2, 1),
    createPoint(3, 3),
    createPoint(4, 0)
  ];

  console.log('控制点:');
  points.forEach((p, i) => logPoint(`  P${i}`, p));

  const natural = new CubicSpline(points, { endCondition: 'natural' });

  console.log('\n自然样条 - 在控制点处求值:');
  for (let i = 0; i < points.length; i++) {
    const t = i / (points.length - 1);
    const result = natural.evaluate(t);
    const expected = points[i];
    const err = Math.max(
      Math.abs(result.x - expected.x),
      Math.abs(result.y - expected.y)
    );
    console.log(`  t=${t.toFixed(2)}: 计算值(${result.x.toFixed(4)}, ${result.y.toFixed(4)}) 误差: ${err.toExponential(2)}`);
  }

  const continuity = natural.checkContinuity();
  console.log('\n连续性检查:');
  console.log(`  位置连续: ${continuity.position}`);
  console.log(`  一阶导连续: ${continuity.firstDerivative}`);
  console.log(`  二阶导连续: ${continuity.secondDerivative}`);
  if (continuity.errors.length > 0) {
    console.log('  误差详情:', continuity.errors);
  }

  console.log('\n样条段系数 (x方向三次多项式):');
  natural.getSegments().forEach((seg, i) => {
    const c = seg.xCoeff;
    console.log(`  段${i}: S(t) = ${c.a.toFixed(4)} + ${c.b.toFixed(4)}t + ${c.c.toFixed(4)}t² + ${c.d.toFixed(4)}t³`);
  });

  console.log('\n总弧长:', natural.getTotalLength().toFixed(4));

  const clamped = new CubicSpline(points, {
    endCondition: 'clamped',
    startTangent: createPoint(1, 1),
    endTangent: createPoint(1, -1)
  });
  console.log('\n固定端点斜率样条 - 端点导数值:');
  logPoint('  t=0 处导数', clamped.derivative(0));
  logPoint('  t=1 处导数', clamped.derivative(1));
}

function testTridiagonal(): void {
  logTitle('3. 三对角方程组求解 (Tridiagonal Solver)');

  const a = [0, 1, 1, 1];
  const b = [4, 4, 4, 4];
  const c = [1, 1, 1, 0];
  const d = [1, 2, 2, 1];

  console.log('三对角矩阵:');
  console.log('  对角线 a (下): [0, 1, 1, 1]');
  console.log('  对角线 b (主): [4, 4, 4, 4]');
  console.log('  对角线 c (上): [1, 1, 1, 0]');
  console.log('  右端项 d:     [1, 2, 2, 1]');

  const x = solveTridiagonal({ a, b, c, d });
  console.log('\n解 x:', x.map(v => v.toFixed(6)));

  console.log('\n验证 (计算 b[i]*x[i] + a[i]*x[i-1] + c[i]*x[i+1]):');
  for (let i = 0; i < 4; i++) {
    let sum = b[i] * x[i];
    if (i > 0) sum += a[i] * x[i - 1];
    if (i < 3) sum += c[i] * x[i + 1];
    console.log(`  i=${i}: ${sum.toFixed(6)} (预期: ${d[i]})`);
  }
}

function testBezier(): void {
  logTitle('4. 贝塞尔曲线 (Bezier Curve - de Casteljau算法)');

  const p0 = createPoint(0, 0);
  const p1 = createPoint(1, 3);
  const p2 = createPoint(3, 3);
  const p3 = createPoint(4, 0);

  console.log('三次贝塞尔控制点:');
  [p0, p1, p2, p3].forEach((p, i) => logPoint(`  P${i}`, p));

  const bezier = new CubicBezier(p0, p1, p2, p3);

  console.log('\n使用三种方法在 t=0.5 处求值:');
  logPoint('  de Casteljau (递归)', bezier.evaluate(0.5));
  logPoint('  de Casteljau (迭代)', bezier.evaluateIterative(0.5));
  logPoint('  Bernstein 基函数', bezier.evaluateBernstein(0.5));

  console.log('\nde Casteljau 三角形 (t=0.5):');
  const table = bezier.getDeCasteljauTable(0.5);
  table.forEach((row, level) => {
    const indent = '  '.repeat(level);
    const pts = row.map(p => `(${p.x.toFixed(2)},${p.y.toFixed(2)})`).join(', ');
    console.log(`  L${level}:${indent}${pts}`);
  });

  console.log('\n细分曲线 (t=0.5):');
  const { left, right } = bezier.subdivide(0.5);
  console.log('  左半部分控制点:');
  left.forEach((p, i) => logPoint(`    L${i}`, p));
  console.log('  右半部分控制点:');
  right.forEach((p, i) => logPoint(`    R${i}`, p));

  console.log('\n验证细分后曲线在 t=0.25:');
  logPoint('  原曲线 t=0.25', bezier.evaluate(0.25));
  logPoint('  左曲线 t=0.5', new BezierCurve(left).evaluate(0.5));

  console.log('\n导数 (速度向量) 在 t=0, 0.5, 1:');
  [0, 0.5, 1].forEach(t => {
    const d = bezier.derivative(t);
    console.log(`  t=${t}: (${d.x.toFixed(4)}, ${d.y.toFixed(4)}) |v|=${Math.sqrt(d.x * d.x + d.y * d.y).toFixed(4)}`);
  });

  console.log('\n切线 (单位向量) 在 t=0.5:');
  logPoint('  切线', bezier.tangent(0.5));
  logPoint('  法线', bezier.normal(0.5));

  console.log('\n自适应细分 (flatness=0.1):');
  const adaptive = bezier.adaptiveSubdivision(0.1);
  console.log(`  段数: ${adaptive.length}`);

  console.log('\nBernstein 基函数 (n=3, t=0.3):');
  for (let i = 0; i <= 3; i++) {
    console.log(`  B_3,${i}(0.3) = ${bernsteinPolynomial(3, i, 0.3).toFixed(6)}`);
  }
  console.log('  验证总和 =', [0, 1, 2, 3].map(i => bernsteinPolynomial(3, i, 0.3)).reduce((a, b) => a + b).toFixed(6));

  console.log('\n总弧长:', bezier.getTotalLength().toFixed(4));
}

function testCatmullRom(): void {
  logTitle('5. Catmull-Rom 样条 (经过所有控制点)');

  const points: Point[] = [
    createPoint(0, 0),
    createPoint(1, 2),
    createPoint(2, 0),
    createPoint(3, 3),
    createPoint(4, 1),
    createPoint(5, 2)
  ];

  console.log('控制点:');
  points.forEach((p, i) => logPoint(`  P${i}`, p));

  const spline = new CatmullRomSpline(points, { tension: 0.5 });

  console.log('\n验证经过所有控制点:');
  const n = points.length;
  for (let i = 0; i < n; i++) {
    const t = i / (n - 1);
    const result = spline.evaluate(t);
    const expected = points[i];
    const err = Math.max(
      Math.abs(result.x - expected.x),
      Math.abs(result.y - expected.y)
    );
    console.log(`  P${i} (t=${t.toFixed(2)}): 误差=${err.toExponential(2)}`);
  }
  console.log(`  全部通过? ${spline.verifyInterpolation()}`);

  console.log(`\nC¹ 连续? ${spline.checkC1Continuity()}`);

  console.log('\n张力 (tension) 对比:');
  const tensions = [0, 0.5, 1.0];
  tensions.forEach(tension => {
    const s = new CatmullRomSpline(points, { tension });
    console.log(`  tension=${tension}: t=0.5 处值 = (${s.evaluate(0.5).x.toFixed(4)}, ${s.evaluate(0.5).y.toFixed(4)})`);
  });

  console.log('\n转换为贝塞尔段:');
  const beziers = spline.toBezierSegments();
  beziers.forEach((b, i) => {
    const cps = b.getControlPoints();
    console.log(`  段${i}: P0(${cps[0].x.toFixed(2)},${cps[0].y.toFixed(2)}) C1(${cps[1].x.toFixed(2)},${cps[1].y.toFixed(2)}) C2(${cps[2].x.toFixed(2)},${cps[2].y.toFixed(2)}) P1(${cps[3].x.toFixed(2)},${cps[3].y.toFixed(2)})`);
  });

  console.log('\n总弧长:', spline.getTotalLength().toFixed(4));

  console.log('\n闭合 Catmull-Rom 样条:');
  const closedPoints = [
    createPoint(0, 0),
    createPoint(2, 1),
    createPoint(1, 3),
    createPoint(-1, 2)
  ];
  const closed = new CatmullRomSpline(closedPoints, { closed: true });
  console.log(`  闭合? ${closed.isClosed()}`);
  console.log(`  t=0 处 = (${closed.evaluate(0).x.toFixed(4)}, ${closed.evaluate(0).y.toFixed(4)})`);
  console.log(`  t=1 处 = (${closed.evaluate(1).x.toFixed(4)}, ${closed.evaluate(1).y.toFixed(4)})`);
  console.log(`  经过所有点? ${closed.verifyInterpolation()}`);
}

function testArcLengthSampling(): void {
  logTitle('6. 弧长均匀采样 vs 参数均匀采样');

  const p0 = createPoint(0, 0);
  const p1 = createPoint(0, 5);
  const p2 = createPoint(10, 5);
  const p3 = createPoint(10, 0);
  const bezier = new CubicBezier(p0, p1, p2, p3);

  console.log('贝塞尔曲线 (形成 S 形):');
  [p0, p1, p2, p3].forEach((p, i) => logPoint(`  P${i}`, p));

  console.log('\n总弧长:', bezier.getTotalLength().toFixed(4));

  const N = 10;

  console.log(`\n按参数 t 均匀采样 ${N} 个点 (相邻弧长):`);
  const byT = bezier.sample(N);
  for (let i = 1; i < N; i++) {
    const dx = byT[i].x - byT[i - 1].x;
    const dy = byT[i].y - byT[i - 1].y;
    const dist = Math.sqrt(dx * dx + dy * dy);
    console.log(`  段${i - 1}-${i}: ${dist.toFixed(4)}`);
  }

  console.log(`\n按弧长均匀采样 ${N} 个点 (相邻弧长):`);
  const byArc = bezier.sampleByArcLength(N);
  const lengths: number[] = [];
  for (let i = 1; i < N; i++) {
    const dx = byArc[i].x - byArc[i - 1].x;
    const dy = byArc[i].y - byArc[i - 1].y;
    const dist = Math.sqrt(dx * dx + dy * dy);
    lengths.push(dist);
    console.log(`  段${i - 1}-${i}: ${dist.toFixed(4)}`);
  }

  const avg = lengths.reduce((a, b) => a + b, 0) / lengths.length;
  const variance = lengths.reduce((s, l) => s + (l - avg) ** 2, 0) / lengths.length;
  console.log(`\n弧长均匀采样统计: 平均=${avg.toFixed(4)}, 标准差=${Math.sqrt(variance).toFixed(6)}`);

  console.log(`参数均匀采样相邻段长标准差: (应远大于弧长均匀采样)`);

  console.log('\n自定义弧长参数化示例 (走 1/3 弧长处的点): 找到对应 t 值:');
  const sampler = new ArcLengthSampler((t) => bezier.evaluate(t), 2000);
  const targetS = bezier.getTotalLength() / 3;
  const tAtS = sampler.arcLengthToT(targetS);
  console.log(`  s = ${targetS.toFixed(4)} → t = ${tAtS.toFixed(6)}`);
  console.log(`  验证: s(t) = ${sampler.tToArcLength(tAtS).toFixed(4)}`);
}

function testAnimationPath(): void {
  logTitle('7. 动画轨迹示例');

  const keyframes: Point[] = [
    createPoint(0, 0),
    createPoint(50, 100),
    createPoint(150, 50),
    createPoint(200, 200),
    createPoint(300, 100),
    createPoint(400, 0)
  ];

  console.log('关键帧 (动画路径):');
  keyframes.forEach((p, i) => logPoint(`  帧${i}`, p));

  const crSpline = new CatmullRomSpline(keyframes, { tension: 0.5 });

  console.log('\n10帧动画 (按弧长匀速):');
  const frames = crSpline.sampleByArcLength(11);
  frames.forEach((p, i) => {
    console.log(`  帧${i}: (${p.x.toFixed(1)}, ${p.y.toFixed(1)})`);
  });

  console.log('\n对比参数均匀 (前几段更快):');
  const framesT = crSpline.sample(11);
  framesT.forEach((p, i) => {
    console.log(`  帧${i}: (${p.x.toFixed(1)}, ${p.y.toFixed(1)})`);
  });
}

function main(): void {
  console.log('\n╔══════════════════════════════════════════════════════╗');
  console.log('║       曲线与插值引擎 - 完整演示                    ║');
  console.log('╚══════════════════════════════════════════════════════╝');

  testLinearInterpolation();
  testTridiagonal();
  testCubicSpline();
  testBezier();
  testCatmullRom();
  testArcLengthSampling();
  testAnimationPath();

  console.log('\n╔══════════════════════════════════════════════════════╗');
  console.log('║              演示完成！所有功能正常工作                 ║');
  console.log('╚══════════════════════════════════════════════════════╝\n');
}

main();
