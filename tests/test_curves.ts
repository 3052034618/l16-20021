import {
  Point, createPoint, LinearInterpolation, CubicSpline,
  BezierCurve, CubicBezier, CatmullRomSpline, solveTridiagonal,
  verifyTridiagonalSolution, bernsteinPolynomial, binomialCoefficient,
  ArcLengthSampler, distance
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

console.log('\n' + '='.repeat(50));
console.log(` 测试结果: ${passed} 通过, ${failed} 失败`);
console.log('='.repeat(50));

if (failed > 0) {
  process.exit(1);
}
