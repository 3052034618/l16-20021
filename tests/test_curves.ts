import {
  Point, createPoint, LinearInterpolation, CubicSpline,
  BezierCurve, CubicBezier, CatmullRomSpline, CompositeBezierCurve, solveTridiagonal,
  verifyTridiagonalSolution, bernsteinPolynomial, binomialCoefficient,
  ArcLengthSampler, distance, CurveAnalyzer, CurveIO, cleanInputPoints,
  computeBoundingBox, isFinitePoint, autoKnots, bboxContainsPoint
} from '../src/index';

let passed = 0;
let failed = 0;

function assert(condition: boolean, message: string): void {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.log(`  ✗ ${message}`);
  }
}

function assertApprox(a: number, b: number, tolerance: number, message: string): void {
  const err = Math.abs(a - b);
  assert(err <= tolerance, `${message} (|${a} - ${b}| = ${err.toExponential(2)} <= ${tolerance})`);
}

function assertPointApprox(a: Point, b: Point, tolerance: number, message: string): void {
  const err = Math.max(Math.abs(a.x - b.x), Math.abs(a.y - b.y));
  assert(err <= tolerance, `${message} (err=${err.toExponential(2)})`);
}

function testSection(title: string, fn: () => void): void {
  console.log(`\n${title}`);
  console.log('-'.repeat(50));
  fn();
}

testSection('测试 1: 线性插值', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(1, 1),
    createPoint(2, 0),
    createPoint(3, 1)
  ];
  const lin = new LinearInterpolation(pts);

  assertPointApprox(lin.evaluate(0), pts[0], 1e-10, 't=0 等于第一个点');
  assertPointApprox(lin.evaluate(1), pts[3], 1e-10, 't=1 等于最后一个点');
  assertPointApprox(lin.evaluate(1 / 3), createPoint(1, 1), 1e-10, 't=1/3 等于 P1');
  assertPointApprox(lin.evaluate(0.5), createPoint(1.5, 0.5), 1e-10, 't=0.5 在中间');
  assert(lin.sample(5).length === 5, 'sample 返回正确数量的点');
  const manualLen = 3 * Math.sqrt(2);
  assertApprox(lin.getTotalLength(), manualLen, 1e-6, '总弧长正确');
});

testSection('测试 2: 三对角求解器', () => {
  const a = [0, 1, 1, 1];
  const b = [4, 4, 4, 4];
  const c = [1, 1, 1, 0];
  const d = [5, 10, 10, 5];

  const x = solveTridiagonal({ a, b, c, d });
  const residual = verifyTridiagonalSolution({ a, b, c, d }, x);
  const maxResidual = residual.reduce((m: number, r: number) => Math.max(m, Math.abs(r)), 0);

  assert(maxResidual < 1e-12, `残差足够小 (${maxResidual.toExponential(2)})`);

  const a2 = [0, 1, 1, 1];
  const b2 = [2, 2, 2, 2];
  const c2 = [1, 1, 1, 0];
  const d2 = [4, 7, 7, 4];
  const x2 = solveTridiagonal({ a: a2, b: b2, c: c2, d: d2 });
  assertApprox(x2[0], 1, 1e-10, '解 x0=1');
  assertApprox(x2[1], 2, 1e-10, '解 x1=2');
  assertApprox(x2[2], 2, 1e-10, '解 x2=2');
  assertApprox(x2[3], 1, 1e-10, '解 x3=1');
});

testSection('测试 3: 三次样条插值', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(1, 2),
    createPoint(2, 1),
    createPoint(3, 3),
    createPoint(4, 0)
  ];

  const natural = new CubicSpline(pts, { endCondition: 'natural' });

  for (let i = 0; i < pts.length; i++) {
    const t = i / (pts.length - 1);
    assertPointApprox(natural.evaluate(t), pts[i], 1e-8, `经过控制点 P${i}`);
  }

  const cont = natural.checkContinuity(1e-6);
  assert(cont.position, '位置连续');
  assert(cont.firstDerivative, '一阶导连续');
  assert(cont.secondDerivative, '二阶导连续');

  const samples = natural.sample(100);
  assert(samples.length === 100, '采样100个点');

  const clamped = new CubicSpline(pts, {
    endCondition: 'clamped',
    startTangent: createPoint(1, 0),
    endTangent: createPoint(1, 0)
  });
  const d0 = clamped.derivative(0);
  const d1 = clamped.derivative(1);
  assert(d0.x !== 0 || d0.y !== 0, '起点导数非零');
  assert(d1.x !== 0 || d1.y !== 0, '终点导数非零');
});

testSection('测试 4: 贝塞尔曲线与 de Casteljau', () => {
  const p0 = createPoint(0, 0);
  const p1 = createPoint(0, 2);
  const p2 = createPoint(2, 2);
  const p3 = createPoint(2, 0);
  const bz = new CubicBezier(p0, p1, p2, p3);

  assertPointApprox(bz.evaluate(0), p0, 1e-10, 't=0 等于 P0');
  assertPointApprox(bz.evaluate(1), p3, 1e-10, 't=1 等于 P3');

  const mid = createPoint(1, 1.5);
  assertPointApprox(bz.evaluate(0.5), mid, 1e-10, 't=0.5 正确');

  const pRecursive = bz.evaluate(0.3);
  const pIterative = bz.evaluateIterative(0.3);
  const pBernstein = bz.evaluateBernstein(0.3);
  assertPointApprox(pRecursive, pIterative, 1e-10, '递归=迭代');
  assertPointApprox(pRecursive, pBernstein, 1e-10, '递归=Bernstein');

  const { left, right } = bz.subdivide(0.5);
  const leftCurve = new BezierCurve(left);
  const rightCurve = new BezierCurve(right);
  assertPointApprox(leftCurve.evaluate(1), rightCurve.evaluate(0), 1e-10, '左右相接');
  assertPointApprox(leftCurve.evaluate(1), bz.evaluate(0.5), 1e-10, '中点正确');
  assertPointApprox(leftCurve.evaluate(0.5), bz.evaluate(0.25), 1e-8, '左曲线 t=0.5 = 原曲线 t=0.25');
  assertPointApprox(rightCurve.evaluate(0.5), bz.evaluate(0.75), 1e-8, '右曲线 t=0.5 = 原曲线 t=0.75');

  let total = 0;
  for (let i = 0; i <= 3; i++) {
    total += bernsteinPolynomial(3, i, 0.4);
  }
  assertApprox(total, 1, 1e-10, 'Bernstein 多项式和为1');

  const elevated = bz.elevateDegree();
  assert(elevated.getDegree() === 4, '升阶后 degree=4');
  assertPointApprox(elevated.evaluate(0.3), bz.evaluate(0.3), 1e-8, '升阶不改变曲线');

  const d = bz.derivative(0.5);
  assert(d.x !== 0 || d.y !== 0, '导数非零');

  const adaptive = bz.adaptiveSubdivision(0.01);
  assert(adaptive.length >= 2, '自适应细分至少2段');
});

testSection('测试 5: Catmull-Rom 样条', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(1, 3),
    createPoint(2, 1),
    createPoint(3, 4),
    createPoint(4, 0)
  ];

  const cr = new CatmullRomSpline(pts, { tension: 0.5 });

  for (let i = 0; i < pts.length; i++) {
    const t = i / (pts.length - 1);
    assertPointApprox(cr.evaluate(t), pts[i], 1e-8, `经过控制点 P${i}`);
  }

  assert(cr.verifyInterpolation(1e-6), '验证插值通过');
  assert(cr.checkC1Continuity(1e-5), 'C¹ 连续');

  const closed = new CatmullRomSpline(pts, { closed: true });
  assert(closed.isClosed(), '闭合标记正确');
  const closedPts = closed.getControlPoints();
  const nClosed = closedPts.length;
  let closedInterpOk = true;
  for (let i = 0; i < nClosed; i++) {
    const t = i / nClosed;
    const result = closed.evaluate(t);
    const expected = closedPts[i];
    const err = Math.max(Math.abs(result.x - expected.x), Math.abs(result.y - expected.y));
    if (err > 1e-6) {
      closedInterpOk = false;
      break;
    }
  }
  assert(closedInterpOk, '闭合曲线经过所有点');
  assertPointApprox(closed.evaluate(0), closed.evaluate(1), 1e-6, '首尾相接');

  const beziers = cr.toBezierSegments();
  assert(beziers.length === pts.length - 1, '贝塞尔段数量正确');
  for (let i = 0; i < beziers.length; i++) {
    const b = beziers[i];
    const t = (i + 0.5) / beziers.length;
    const fromBz = b.evaluate(0.5);
    const fromCr = cr.evaluate(t);
    assertPointApprox(fromBz, fromCr, 1e-5, `贝塞尔转换段${i}中点匹配`);
  }
});

testSection('测试 6: 弧长均匀采样', () => {
  const p0 = createPoint(0, 0);
  const p1 = createPoint(0, 10);
  const p2 = createPoint(10, 10);
  const p3 = createPoint(10, 0);
  const bz = new CubicBezier(p0, p1, p2, p3);

  const sampler = new ArcLengthSampler((t) => bz.evaluate(t), 2000);

  assertApprox(sampler.tToArcLength(0), 0, 1e-6, 't=0 对应 s=0');
  assertApprox(sampler.tToArcLength(1), sampler.getTotalLength(), 1e-6, 't=1 对应总弧长');

  const s = sampler.getTotalLength() / 2;
  const t = sampler.arcLengthToT(s);
  assertApprox(sampler.tToArcLength(t), s, 1e-4, '弧长映射互逆');

  const N = 20;
  const uniform = bz.sampleByArcLength(N);
  const lengths: number[] = [];
  for (let i = 1; i < N; i++) {
    lengths.push(distance(uniform[i - 1], uniform[i]));
  }
  const avg = lengths.reduce((a, b) => a + b, 0) / lengths.length;
  const maxDev = lengths.reduce((m, l) => Math.max(m, Math.abs(l - avg)), 0);
  assert(maxDev / avg < 0.02, `弧长均匀度偏差小 (${(maxDev / avg * 100).toFixed(2)}%)`);

  const tSamples = bz.sample(N);
  let maxDevT = 0;
  let avgT = 0;
  const tLengths: number[] = [];
  for (let i = 1; i < N; i++) {
    tLengths.push(distance(tSamples[i - 1], tSamples[i]));
  }
  avgT = tLengths.reduce((a, b) => a + b, 0) / tLengths.length;
  maxDevT = tLengths.reduce((m, l) => Math.max(m, Math.abs(l - avgT)), 0);
  assert(maxDevT / avgT > maxDev / avg, '参数均匀采样比弧长均匀采样不均匀度大');
});

testSection('测试 7: 数学函数', () => {
  assert(binomialCoefficient(0, 0) === 1, 'C(0,0)=1');
  assert(binomialCoefficient(5, 0) === 1, 'C(5,0)=1');
  assert(binomialCoefficient(5, 5) === 1, 'C(5,5)=1');
  assert(binomialCoefficient(5, 2) === 10, 'C(5,2)=10');
  assert(binomialCoefficient(6, 3) === 20, 'C(6,3)=20');

  assertApprox(bernsteinPolynomial(1, 0, 0.5), 0.5, 1e-10, 'B_1,0(0.5)=0.5');
  assertApprox(bernsteinPolynomial(1, 1, 0.5), 0.5, 1e-10, 'B_1,1(0.5)=0.5');
  assertApprox(bernsteinPolynomial(2, 1, 0.5), 0.5, 1e-10, 'B_2,1(0.5)=0.5');
});

testSection('测试 8: 边界与异常情况', () => {
  const twoPts = [createPoint(0, 0), createPoint(1, 1)];
  const lin2 = new LinearInterpolation(twoPts);
  assertPointApprox(lin2.evaluate(0.5), createPoint(0.5, 0.5), 1e-10, '两点线性插值');

  const cr2 = new CatmullRomSpline(twoPts, { tension: 0.5 });
  assertPointApprox(cr2.evaluate(0), twoPts[0], 1e-1, '两点CR起点');
  assertPointApprox(cr2.evaluate(1), twoPts[1], 1e-1, '两点CR终点');

  const quadBz = new BezierCurve([createPoint(0, 0), createPoint(1, 2), createPoint(2, 0)]);
  assert(quadBz.getDegree() === 2, '二次贝塞尔 degree=2');
  assertPointApprox(quadBz.evaluate(0.5), createPoint(1, 1), 1e-10, '二次贝塞尔中点');

  let errorThrown = false;
  try {
    new LinearInterpolation([createPoint(0, 0)]);
  } catch {
    errorThrown = true;
  }
  assert(errorThrown, '少于2点抛出异常');
});

testSection('测试 9: 数值稳定性与输入清洗', () => {
  const ptsWithDup: Point[] = [
    createPoint(0, 0),
    createPoint(0, 0),
    createPoint(1, 1),
    createPoint(1, 1),
    createPoint(2, 0)
  ];
  const cleaned = cleanInputPoints(ptsWithDup);
  assert(cleaned.removedDuplicates === true, '重复点被检测');
  assert(cleaned.points.length === 3, '去重后剩3点');
  assert(cleaned.originalIndices.length === 3, '原始索引对应');

  const allSame: Point[] = [
    createPoint(5, 5),
    createPoint(5, 5),
    createPoint(5, 5)
  ];
  const allSameCleaned = cleanInputPoints(allSame);
  assert(allSameCleaned.allCollapsed === true, '全重合点标记 collapsed');

  const nanPts: Point[] = [
    createPoint(0, 0),
    createPoint(NaN, NaN),
    createPoint(1, 1)
  ];
  const nanCleaned = cleanInputPoints(nanPts);
  assert(nanCleaned.points.length === 2, 'NaN 点被过滤');
  assert(isFinitePoint(nanCleaned.points[0]) === true, '剩余点为有限值');

  const bboxPts = [createPoint(0, 0), createPoint(3, 4)];
  const bbox = computeBoundingBox([createPoint(0, 0), createPoint(3, 4)]);
  assert(bbox.minX === 0 && bbox.maxX === 3, '包围盒 X 范围正确');
  assert(bbox.minY === 0 && bbox.maxY === 4, '包围盒 Y 范围正确');
  assert(bboxContainsPoint(bbox, createPoint(1.5, 2)), '包围盒包含内部点');
  assert(bboxContainsPoint(bbox, createPoint(10, 10)) === false, '包围盒不包含外部点');
  assert(bbox.width === 3 && bbox.height === 4, '包围盒宽高正确');

  const knots = autoKnots(5, 'uniform');
  assert(knots.length === 5, '均匀knots数量正确');
  assert(knots[0] === 0, '第一个knot=0');
  assert(Math.abs(knots[4] - 1) < 1e-15, '最后一个knot=1');
});

testSection('测试 10: 自定义 knots / 时间戳', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(1, 2),
    createPoint(3, 1),
    createPoint(6, 3)
  ];

  const uniform = new CubicSpline(pts, { knotType: 'uniform' });
  const uk = uniform.getKnots();
  assertApprox(uk[0], 0, 1e-12, 'uniform knots[0]=0');
  assertApprox(uk[uk.length - 1], 1, 1e-12, 'uniform knots[-1]=1');

  for (let i = 0; i < pts.length; i++) {
    const t = uk[i];
    assertPointApprox(uniform.evaluate(t), pts[i], 1e-6, 'custom-knot 经过 P' + i);
  }

  const chordal = new CubicSpline(pts, { knotType: 'chordal' });
  const ck = chordal.getKnots();
  assert(ck[0] === 0 && ck[ck.length - 1] === 1, 'chordal knots 归一化');

  const timestamps = [0, 10, 30, 60];
  const timeBased = new CubicSpline(pts, { timeStamps: timestamps });
  const tk = timeBased.getKnots();
  assertApprox(tk[0], 0, 1e-12, '时间戳归一化 t0=0');
  assertApprox(tk[1], 10 / 60, 1e-6, '时间戳归一化 t1=1/6');
  assertApprox(tk[tk.length - 1], 1, 1e-12, '时间戳归一化 tend=1');

  const crTimed = new CatmullRomSpline(pts, { timeStamps: [0, 1, 2, 5] });
  for (let i = 0; i < pts.length; i++) {
    const t = crTimed.getKnots()[i];
    assertPointApprox(crTimed.evaluate(t), pts[i], 1e-6, 'CR 时间戳 P' + i);
  }
});

testSection('测试 11: 曲线分析器', () => {
  const p0 = createPoint(0, 0);
  const p1 = createPoint(0, 10);
  const p2 = createPoint(10, 10);
  const p3 = createPoint(10, 0);
  const cr = new CatmullRomSpline([p0, p1, p2, p3], { tension: 0.5 });
  const analyzer = cr.getAnalyzer() as CurveAnalyzer;

  const analysis = analyzer.analyze();
  assert(analysis.boundingBox.width > 0, '包围盒有宽度');
  assert(analysis.boundingBox.height > 0, '包围盒有高度');
  assert(analysis.arcLength > 0, '弧长大于0');
  assert(analysis.estimatedMaxSpeed >= 0, '最大速度非负');

  const bbox = analyzer.getBoundingBox();
  assert(bbox.minX <= 0 && bbox.maxX >= 10, 'X范围覆盖');
  assert(bbox.minY <= 0 && bbox.maxY >= 10, 'Y范围覆盖');

  const speed = analyzer.speed(0.5);
  assert(speed > 0, 't=0.5 速度>0');

  const target = createPoint(5, 5);
  const proj = analyzer.findNearestPoint(target);
  assert(proj.distance < 10, '最近点距离有限');
  assert(proj.parameter >= 0 && proj.parameter <= 1, '最近点参数在[0,1]');

  const k = analyzer.curvature(0.5);
  assert(Number.isFinite(k), '曲率是有限值');

  const nearCheck = analyzer.isPointNearCurve(createPoint(5, 5), 20);
  assert(nearCheck.near === true, '附近点检测');
  const farCheck = analyzer.isPointNearCurve(createPoint(100, 100), 1);
  assert(farCheck.near === false, '远处点检测');
});

testSection('测试 12: SVG 导出与互转', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(1, 3),
    createPoint(2, 1),
    createPoint(3, 4),
    createPoint(4, 0)
  ];

  const cr = new CatmullRomSpline(pts);
  const svgPath = CurveIO.toSVGPath(cr);
  assert(svgPath.startsWith('M '), 'SVG path 以 M 开头');
  assert(svgPath.includes('C '), 'SVG path 包含 C 指令');

  const lin = new LinearInterpolation(pts);
  const polySvg = CurveIO.toSVGPath(lin.sample(20));
  assert(polySvg.startsWith('M '), '折线 SVG 以 M 开头');

  const bzSegments = cr.toBezierSegments();
  const compositeSvg = CurveIO.toSVGPath(bzSegments);
  assert(compositeSvg.includes('C '), '贝塞尔段 SVG 正确');

  const noisyPts = cr.sample(30);
  const fit = CurveIO.fitCubicBezier(noisyPts, { maxError: 0.5 });
  assert(fit.curves.length >= 1, '贝塞尔拟合至少1段');
  assert(fit.curves[0] instanceof CubicBezier, '拟合结果为 CubicBezier');

  const firstCP = fit.curves[0].getControlPoints();
  assertPointApprox(firstCP[0], noisyPts[0], 0.5, '拟合起点接近采样起点');

  const fullSvg = CurveIO.curveToFullSVG(cr, { width: 200, height: 150 });
  assert(fullSvg.includes('<svg'), '完整 SVG 含 svg 标签');
  assert(fullSvg.includes('<path'), '完整 SVG 含 path 标签');
});

testSection('测试 13: Catmull-Rom 非均匀时间戳速度 C¹ 连续', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(50, 100),
    createPoint(150, 80),
    createPoint(250, 200),
    createPoint(350, 0)
  ];
  const cr = new CatmullRomSpline(pts, {
    timeStamps: [0, 1, 5, 6, 10],
    tension: 0.5
  });

  assert(cr.verifyInterpolation(1e-6) === true, 'CR 经过所有控制点');
  assert(cr.checkC1Continuity(1e-5) === true, 'CR C¹ 连续通过');

  const analyzer = cr.getAnalyzer();
  const profile = analyzer.speedProfile(200);
  for (const p of profile) {
    assert(isFinite(p.speed) && p.speed >= 0, '速度值有效');
    assert(p.curvature === undefined || isFinite(p.curvature), '曲率值有效');
  }

  const knots = cr.getKnots();
  for (let i = 1; i < knots.length - 1; i++) {
    const k = knots[i];
    const dPrev = cr.derivative(Math.max(0, k - 1e-6));
    const dNext = cr.derivative(Math.min(1, k + 1e-6));
    const lenPrev = Math.sqrt(dPrev.x * dPrev.x + dPrev.y * dPrev.y);
    const lenNext = Math.sqrt(dNext.x * dNext.x + dNext.y * dNext.y);
    if (lenPrev > 1e-6 && lenNext > 1e-6) {
      const dot = (dPrev.x * dNext.x + dPrev.y * dNext.y) / (lenPrev * lenNext);
      assert(Math.abs(dot - 1) < 1e-4, `关键帧${i}处方向连续`);
      const speedRatio = Math.abs(lenPrev - lenNext) / Math.max(lenPrev, lenNext);
      assert(speedRatio < 0.05, `关键帧${i}处速度跳变<5%`);
    }
  }
});

testSection('测试 14: SVG 导出形状与 evaluate 严格重合', () => {
  const pts: Point[] = [
    createPoint(0, 0),
    createPoint(50, 100),
    createPoint(150, 80),
    createPoint(250, 200)
  ];

  const cr = new CatmullRomSpline(pts, { timeStamps: [0, 1, 5, 10] });
  const bzSegs = cr.toBezierSegments();
  const knots = cr.getKnots();

  assert(bzSegs.length === pts.length - 1, 'CR 段数与控制点数量一致');

  for (let i = 0; i < bzSegs.length; i++) {
    const cp = bzSegs[i].getControlPoints();
    assertPointApprox(cp[0], pts[i], 1e-6, `贝塞尔段${i}起点与CR控制点重合`);
    assertPointApprox(cp[3], pts[i + 1], 1e-6, `贝塞尔段${i}终点与CR控制点重合`);
  }

  const perSegChecks = 20;
  let maxErr = 0;
  for (let i = 0; i < bzSegs.length; i++) {
    const kStart = knots[i];
    const kEnd = knots[i + 1];
    for (let j = 0; j <= perSegChecks; j++) {
      const s = j / perSegChecks;
      const globalT = kStart + s * (kEnd - kStart);
      const crPt = cr.evaluate(globalT);
      const bzPt = bzSegs[i].evaluate(s);
      const err = Math.max(
        Math.abs(crPt.x - bzPt.x),
        Math.abs(crPt.y - bzPt.y)
      );
      if (err > maxErr) maxErr = err;
      assert(err < 1e-5, `段${i} s=${s.toFixed(3)} CR与贝塞尔形状一致, err=${err}`);
    }
  }

  const crSvg = CurveIO.toSVGPath(cr);
  const bzSvg = CurveIO.toSVGPath(bzSegs);
  assert(crSvg === bzSvg, 'CR SVG 与贝塞尔段 SVG 字符串完全相同');
});

testSection('测试 15: 大量控制点性能 / 稳定性', () => {
  const N_SMALL = 1000;
  const smallPts: Point[] = [];
  for (let i = 0; i < N_SMALL; i++) {
    const t = i / (N_SMALL - 1);
    smallPts.push(createPoint(t * 1000, Math.sin(t * Math.PI * 4) * 100 + t * 200));
  }

  const t1 = Date.now();
  const linBig = new LinearInterpolation(smallPts);
  const t2 = Date.now();
  const buildLin = t2 - t1;
  assert(buildLin < 2000, `LinearInterpolation(${N_SMALL}点) 构建<2s, 实际=${buildLin}ms`);

  const sample100 = linBig.sample(100);
  assert(sample100.length === 100, 'Linear 采样100个点正确');
  for (const p of sample100) {
    assert(Number.isFinite(p.x) && Number.isFinite(p.y), 'Linear 采样有限值');
  }

  const t3 = Date.now();
  const splineBig = new CubicSpline(smallPts, { endCondition: 'natural' });
  const t4 = Date.now();
  const buildSpline = t4 - t3;
  assert(buildSpline < 5000, `CubicSpline(${N_SMALL}点) 构建<5s, 实际=${buildSpline}ms`);

  const splineSample = splineBig.sample(50);
  assert(splineSample.length === 50, 'Spline 采样50点正确');
  for (const p of splineSample) {
    assert(Number.isFinite(p.x) && Number.isFinite(p.y), 'Spline 采样有限值');
  }

  const N_BEZIER = 200;
  const bzPts: Point[] = [];
  for (let i = 0; i < N_BEZIER; i++) {
    const t = i / (N_BEZIER - 1);
    bzPts.push(createPoint(t * 500, Math.sin(t * 8) * 50 + 100));
  }

  const t5 = Date.now();
  const bzBig = new BezierCurve(bzPts);
  const t6 = Date.now();
  const buildBz = t6 - t5;
  assert(buildBz < 5000, `Bezier(degree=${N_BEZIER - 1}) 构建<5s, 实际=${buildBz}ms`);

  for (let i = 0; i <= 20; i++) {
    const t = i / 20;
    const p = bzBig.evaluate(t);
    assert(Number.isFinite(p.x) && Number.isFinite(p.y), `Bezier evaluate(${t}) 有限值`);
  }

  const bzSample = bzBig.sample(20);
  assert(bzSample.length === 20, 'Bezier sample(20) 正确');

  const CR_N = 1000;
  const crPts: Point[] = [];
  for (let i = 0; i < CR_N; i++) {
    const tt = i / (CR_N - 1);
    crPts.push(createPoint(tt * 2000, Math.cos(tt * Math.PI * 6) * 80 + tt * 500));
  }
  const t7 = Date.now();
  const crBig = new CatmullRomSpline(crPts);
  const t8 = Date.now();
  const buildCR = t8 - t7;
  assert(buildCR < 5000, `CatmullRom(${CR_N}点) 构建<5s, 实际=${buildCR}ms`);
  const crMid = crBig.evaluate(0.5);
  assert(Number.isFinite(crMid.x) && Number.isFinite(crMid.y), 'CR 中点求值有限');
});

testSection('测试 16: SVG 导出坐标变换一致性', () => {
  const p0 = createPoint(10, 20);
  const p1 = createPoint(50, 100);
  const p2 = createPoint(150, 80);
  const p3 = createPoint(200, 200);
  const cubicBz = new CubicBezier(p0, p1, p2, p3);

  const pathStr = CurveIO.toSVGPath(cubicBz, { digits: 2 });
  assert(pathStr.startsWith('M 10 20'), 'SVG path 起点正确');
  assert(pathStr.includes('C 50 100, 150 80, 200 200'), '三次贝塞尔 C 指令控制点正确');
  console.log(`  → C 指令: ${pathStr}`);

  const fullSVG = CurveIO.curveToFullSVG(cubicBz, {
    width: 500, height: 300, padding: 20,
    stroke: '#3366ff', strokeWidth: 2
  });

  const fullSVG2 = CurveIO.curveToFullSVG(new CompositeBezierCurve([cubicBz]), {
    width: 500, height: 300, padding: 20,
    stroke: '#3366ff', strokeWidth: 2
  });

  const path1 = fullSVG.match(/d="([^"]+)"/)?.[1] || '';
  const path2 = fullSVG2.match(/d="([^"]+)"/)?.[1] || '';
  assert(path1 === path2, '单个CubicBezier与CompositeBezier导出的变换后path一致');

  const ptsMatch = fullSVG.match(/<circle cx="([\d.]+)" cy="([\d.]+)"/g);
  const controlPointCount = ptsMatch ? ptsMatch.length : 0;
  assert(controlPointCount === 4, `SVG 中显示了全部 4 个控制点 (实际: ${controlPointCount})`);

  const scaleX = 460 / 190;
  const scaleY = 260 / 180;
  const scale = Math.min(scaleX, scaleY);
  const offsetX = 20 - 10 * scale;
  const offsetY = 20 + 200 * scale;
  const tp1x = (50 * scale + offsetX).toFixed(2);
  const tp1y = (-100 * scale + offsetY).toFixed(2);
  assert(path1.includes(`${tp1x} ${tp1y}`), `变换后控制点 P1 位置正确 (${tp1x},${tp1y})`);
  assert(fullSVG.includes(`cx="${tp1x}" cy="${tp1y}"`), '控制点 marker 与 path 使用相同变换');

  const svgPts = (ptsMatch || []).map(m => {
    const mm = m.match(/cx="([\d.]+)" cy="([\d.]+)"/);
    return { x: parseFloat(mm![1]), y: parseFloat(mm![2]) };
  });
  const tP0x = (10 * scale + offsetX).toFixed(2);
  const tP3x = (200 * scale + offsetX).toFixed(2);
  const tP0y = (-20 * scale + offsetY).toFixed(2);
  const tP3y = (-200 * scale + offsetY).toFixed(2);
  assert(Math.abs(svgPts[0].x - parseFloat(tP0x)) < 0.01, '起点 marker 与计算一致');
  assert(Math.abs(svgPts[3].x - parseFloat(tP3x)) < 0.01, '终点 marker 与计算一致');
  assert(Math.abs(svgPts[0].y - parseFloat(tP0y)) < 0.01, '起点 marker Y 与计算一致');
  assert(Math.abs(svgPts[3].y - parseFloat(tP3y)) < 0.01, '终点 marker Y 与计算一致');

  const qbz = new BezierCurve([createPoint(0, 0), createPoint(50, 100), createPoint(100, 0)]);
  const qPath = CurveIO.toSVGPath(qbz, { digits: 2 });
  assert(qPath.startsWith('M 0 0'), '二次贝塞尔起点正确');
  assert(qPath.includes('Q 50 100, 100 0'), '二次贝塞尔 Q 指令正确');

  const lbz = new BezierCurve([createPoint(10, 20), createPoint(50, 60)]);
  const lPath = CurveIO.toSVGPath(lbz, { digits: 2 });
  assert(lPath.startsWith('M 10 20'), '一次贝塞尔起点正确');
  assert(lPath.includes('L 50 60'), '一次贝塞尔 L 指令正确');
});

testSection('测试 17: 高阶贝塞尔端点稳定性', () => {
  const Ns = [10, 50, 100, 200, 500, 1000];
  for (const N of Ns) {
    const pts: Point[] = [];
    for (let i = 0; i <= N; i++) {
      const angle = (i / N) * Math.PI * 2;
      pts.push(createPoint(Math.cos(angle) * 100, Math.sin(angle) * 100));
    }
    const bz = new BezierCurve(pts);

    const p0 = bz.evaluate(0);
    assertPointApprox(p0, pts[0], 1e-9, `degree=${N} evaluate(0) == 起点`);
    assert(Math.abs(p0.x - pts[0].x) < 1e-10, `evaluate(0).x 不跳变`);
    assert(Math.abs(p0.y - pts[0].y) < 1e-10, `evaluate(0).y 不跳变`);

    const p1 = bz.evaluate(1);
    assertPointApprox(p1, pts[N], 1e-9, `degree=${N} evaluate(1) == 终点`);
    assert(Math.abs(p1.x - pts[N].x) < 1e-10, `evaluate(1).x 不跳变`);
    assert(Math.abs(p1.y - pts[N].y) < 1e-10, `evaluate(1).y 不跳变`);

    let hasNaN = false, hasZeroJump = false;
    for (let i = 0; i <= 50; i++) {
      const t = i / 50;
      const p = bz.evaluate(t);
      if (!Number.isFinite(p.x) || !Number.isFinite(p.y)) {
        hasNaN = true;
      }
      if (t > 0 && t < 1 && Math.abs(p.x) < 1e-10 && Math.abs(p.y) < 1e-10 && pts[0].x !== 0) {
        hasZeroJump = true;
      }
    }
    assert(!hasNaN, `degree=${N} 无 NaN 采样点`);
    assert(!hasZeroJump, `degree=${N} 无跳变到 (0,0)`);

    const t05 = bz.evaluate(0.5);
    assert(Number.isFinite(t05.x) && Number.isFinite(t05.y), `degree=${N} t=0.5 有限`);

    const d0 = bz.derivative(0);
    const d1 = bz.derivative(1);
    assert(Number.isFinite(d0.x) && Number.isFinite(d0.y), `degree=${N} derivative(0) 有限`);
    assert(Number.isFinite(d1.x) && Number.isFinite(d1.y), `degree=${N} derivative(1) 有限`);

    const an = bz.getAnalyzer();
    const sp = an.speedProfile(50);
    const allFinite = sp.every(p => Number.isFinite(p.speed));
    assert(allFinite, `degree=${N} speedProfile 全部有限`);
    const speedAtEnd = sp[sp.length - 1].speed;
    const derivMag = Math.sqrt(d1.x ** 2 + d1.y ** 2);
    assertApprox(speedAtEnd, derivMag, 1e-6, `degree=${N} 末端 speedProfile 与 derivative(1) 一致`);
  }
});

testSection('测试 18: sample 末端与 evaluate(1) 一致性', () => {
  const pts: Point[] = [];
  for (let i = 0; i <= 500; i++) {
    pts.push(createPoint(i, Math.sin(i * 0.1) * 50 + i * 0.5));
  }

  const cr = new CatmullRomSpline(pts, {
    timeStamps: pts.map((_, i) => i * 0.1),
    tension: 0.5,
    closed: false
  });

  const bzHigh = new BezierCurve(pts.slice(0, 200));

  const testCurves = [
    { name: 'Linear(1000)', c: new LinearInterpolation(pts) },
    { name: 'CubicSpline(1000)', c: new CubicSpline(pts) },
    { name: 'CatmullRom(1000)', c: cr },
    { name: 'Bezier(degree=199)', c: bzHigh }
  ];

  for (const { name, c } of testCurves) {
    const samples = c.sample(200);
    const lastSample = samples[samples.length - 1];
    const eval1 = c.evaluate(1);
    assertPointApprox(lastSample, eval1, 1e-9, `${name} sample() 最后一点 == evaluate(1)`);

    const arcSamples = c.sampleByArcLength(200);
    const lastArc = arcSamples[arcSamples.length - 1];
    assertPointApprox(lastArc, eval1, 1e-6, `${name} sampleByArcLength() 最后一点 == evaluate(1)`);

    const totalLen = c.getTotalLength();
    const cAny = c as any;
    const evalByArc = cAny.getArcLengthSampler
      ? cAny.getArcLengthSampler().evaluateByArcLength(totalLen)
      : cAny.evaluateByArcLength
        ? cAny.evaluateByArcLength(totalLen)
        : eval1;
    assertPointApprox(evalByArc, eval1, 1e-6, `${name} evaluateByArcLength(total) == evaluate(1)`);

    const derivAt1 = c.derivative(1);
    const analyzer = c.getAnalyzer();
    const sp = analyzer.speedProfile(200);
    const endSpeed = sp[sp.length - 1].speed;
    const derivMag = Math.sqrt(derivAt1.x ** 2 + derivAt1.y ** 2);
    assertApprox(endSpeed, derivMag, 1e-6, `${name} 末端 speed 与 derivative(1) 一致`);

    const curvatures = analyzer.speedProfile(50);
    const lastCurv = curvatures[curvatures.length - 1];
    assert(
      lastCurv.curvature === undefined || Number.isFinite(lastCurv.curvature),
      `${name} 末端曲率有限`
    );

    const firstSample = samples[0];
    const eval0 = c.evaluate(0);
    assertPointApprox(firstSample, eval0, 1e-9, `${name} sample() 第一点 == evaluate(0)`);

    const arcFirst = arcSamples[0];
    assertPointApprox(arcFirst, eval0, 1e-6, `${name} sampleByArcLength() 第一点 == evaluate(0)`);
  }

  const pts2: Point[] = [];
  for (let i = 0; i <= 5000; i++) {
    pts2.push(createPoint(i * 0.1, Math.cos(i * 0.05) * 30 + Math.sin(i * 0.02) * 20));
  }
  const cr5k = new CatmullRomSpline(pts2, { closed: false, tension: 0.5 });
  const s5k = cr5k.sample(500);
  const endPt = pts2[pts2.length - 1];
  assertPointApprox(s5k[s5k.length - 1], cr5k.evaluate(1), 1e-9, 'CatmullRom(5000) sample 末端 == evaluate(1)');
  assertPointApprox(cr5k.evaluate(0), pts2[0], 1e-9, 'CatmullRom(5000) evaluate(0) == 起点');
  assertPointApprox(cr5k.evaluate(1), endPt, 1e-9, 'CatmullRom(5000) evaluate(1) == 终点');
});

console.log('\n' + '='.repeat(50));
console.log(` 测试结果: ${passed} 通过, ${failed} 失败`);
console.log('='.repeat(50));

if (failed > 0) {
  process.exit(1);
}
