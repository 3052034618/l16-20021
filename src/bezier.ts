import {
  Point, Curve, BezierSubdivision, lerpPoint, createPoint, clamp
} from './types';
import { ArcLengthSampler } from './arc_length';
import { CurveAnalyzer } from './curve_analyzer';
import { cleanInputPoints, safePoint, safeNumber, STABILITY_EPSILON } from './stability';

export function binomialCoefficient(n: number, k: number): number {
  if (k < 0 || k > n) return 0;
  if (k === 0 || k === n) return 1;
  k = Math.min(k, n - k);
  let result = 1;
  for (let i = 1; i <= k; i++) {
    result = result * (n - k + i) / i;
  }
  return Math.round(result);
}

export function bernsteinPolynomial(n: number, i: number, t: number): number {
  return binomialCoefficient(n, i) * Math.pow(t, i) * Math.pow(1 - t, n - i);
}

function evaluateBernsteinStable(points: Point[], t: number): Point {
  const n = points.length - 1;
  if (n < 0) return createPoint(0, 0);
  if (n === 0) return { ...points[0] };

  if (t <= 0) return { ...points[0] };
  if (t >= 1) return { ...points[n] };

  if (n <= 500) {
    let pts = points.slice();
    for (let r = 1; r <= n; r++) {
      for (let i = 0; i <= n - r; i++) {
        pts[i] = {
          x: (1 - t) * pts[i].x + t * pts[i + 1].x,
          y: (1 - t) * pts[i].y + t * pts[i + 1].y
        };
      }
    }
    return createPoint(
      Number.isFinite(pts[0].x) ? pts[0].x : points[n].x,
      Number.isFinite(pts[0].y) ? pts[0].y : points[n].y
    );
  }

  const logT = Math.log(t);
  const logOneMinusT = Math.log(1 - t);

  const logBasis: number[] = new Array(n + 1);
  let logBC = 0;
  for (let i = 1; i <= n; i++) {
    logBC += Math.log(n - i + 1) - Math.log(i);
  }
  let maxLogBasis = -Infinity;
  for (let i = 0; i <= n; i++) {
    let lb = 0;
    if (i === 0) {
      lb = n * logOneMinusT;
    } else if (i === n) {
      lb = n * logT;
    } else {
      let lbc = 0;
      for (let k = 1; k <= i; k++) {
        lbc += Math.log(n - k + 1) - Math.log(k);
      }
      lb = lbc + i * logT + (n - i) * logOneMinusT;
    }
    logBasis[i] = lb;
    if (lb > maxLogBasis) maxLogBasis = lb;
  }

  const weights: number[] = new Array(n + 1);
  let sumWeights = 0;
  for (let i = 0; i <= n; i++) {
    const w = Math.exp(logBasis[i] - maxLogBasis);
    weights[i] = Number.isFinite(w) ? w : 0;
    sumWeights += weights[i];
  }

  if (sumWeights < 1e-300) {
    return { ...points[n] };
  }

  let x = 0, y = 0;
  for (let i = 0; i <= n; i++) {
    x += points[i].x * weights[i];
    y += points[i].y * weights[i];
  }

  return createPoint(
    Number.isFinite(x) ? x / sumWeights : points[n].x,
    Number.isFinite(y) ? y / sumWeights : points[n].y
  );
}

export class BezierCurve implements Curve {
  private controlPoints: Point[];
  private degree: number;
  private arcLengthSampler: ArcLengthSampler;
  private analyzer: CurveAnalyzer;

  constructor(controlPoints: Point[]) {
    const cleaned = cleanInputPoints(controlPoints);
    if (cleaned.points.length < 2) {
      throw new Error('BezierCurve requires at least 2 valid distinct control points');
    }
    this.controlPoints = cleaned.points;
    this.degree = this.controlPoints.length - 1;

    if (this.degree >= 2000) {
      console.warn(
        `[BezierCurve] degree=${this.degree} is extremely high, performance may be limited. ` +
        `Consider using CompositeBezierCurve or CubicSpline for many data points.`
      );
    } else if (this.degree >= 500) {
      console.warn(
        `[BezierCurve] degree=${this.degree} is high, consider using CompositeBezierCurve instead.`
      );
    }

    let samplerRes = 2000;
    let analyzerRes = 200;
    if (this.degree > 2000) {
      samplerRes = 200;
      analyzerRes = 40;
    } else if (this.degree > 500) {
      samplerRes = 500;
      analyzerRes = 60;
    } else if (this.degree > 100) {
      samplerRes = 1000;
      analyzerRes = 100;
    }

    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), samplerRes);
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      analyzerRes
    );
  }

  getDegree(): number {
    return this.degree;
  }

  getControlPoints(): Point[] {
    return [...this.controlPoints];
  }

  evaluate(t: number): Point {
    t = clamp(t, 0, 1);
    if (t <= 0) return { ...this.controlPoints[0] };
    if (t >= 1) return { ...this.controlPoints[this.degree] };
    if (this.degree <= 200) {
      return this.deCasteljauIterative(this.controlPoints, t);
    }
    return safePoint(evaluateBernsteinStable(this.controlPoints, t));
  }

  private deCasteljau(points: Point[], t: number): Point {
    const n = points.length;
    if (n === 1) return { ...points[0] };

    const next: Point[] = [];
    for (let i = 0; i < n - 1; i++) {
      next.push(lerpPoint(points[i], points[i + 1], t));
    }
    return this.deCasteljau(next, t);
  }

  evaluateIterative(t: number): Point {
    t = clamp(t, 0, 1);
    return this.deCasteljauIterative(this.controlPoints, t);
  }

  private deCasteljauIterative(points: Point[], t: number): Point {
    let current = [...points];
    while (current.length > 1) {
      const next: Point[] = [];
      const nCur = current.length;
      for (let i = 0; i < nCur - 1; i++) {
        next.push(lerpPoint(current[i], current[i + 1], t));
      }
      current = next;
    }
    return safePoint({ ...current[0] });
  }

  evaluateBernstein(t: number): Point {
    t = clamp(t, 0, 1);
    return safePoint(evaluateBernsteinStable(this.controlPoints, t));
  }

  getDeCasteljauTable(t: number): Point[][] {
    t = clamp(t, 0, 1);
    if (this.degree > 100) {
      console.warn(`[BezierCurve] getDeCasteljauTable degree=${this.degree} is high, table may be very large.`);
    }
    const table: Point[][] = [];
    table.push([...this.controlPoints]);

    const maxLevels = Math.min(this.degree, 200);
    for (let level = 1; level <= maxLevels; level++) {
      const prev = table[level - 1];
      const current: Point[] = [];
      for (let i = 0; i < prev.length - 1; i++) {
        current.push(lerpPoint(prev[i], prev[i + 1], t));
      }
      table.push(current);
    }
    return table;
  }

  subdivide(t: number): BezierSubdivision {
    t = clamp(t, 0, 1);
    if (this.degree > 100) {
      console.warn(`[BezierCurve] subdivide with degree=${this.degree} may be slow, use adaptiveSubdivision instead.`);
    }
    const table = this.getDeCasteljauTable(t);
    const left: Point[] = [];
    const right: Point[] = [];

    const n = Math.min(this.degree + 1, table.length);
    for (let i = 0; i < n; i++) {
      left.push(table[i][0]);
    }
    for (let i = 0; i < n; i++) {
      right.push(table[n - 1 - i][i]);
    }
    return { left, right };
  }

  split(t: number): [BezierCurve, BezierCurve] {
    if (this.degree > 100) {
      throw new Error(`[BezierCurve] split not supported for degree=${this.degree}, use sampling instead.`);
    }
    const { left, right } = this.subdivide(t);
    return [new BezierCurve(left), new BezierCurve(right)];
  }

  derivative(t: number): Point {
    t = clamp(t, 0, 1);
    if (this.degree === 1) {
      const dx = this.controlPoints[1].x - this.controlPoints[0].x;
      const dy = this.controlPoints[1].y - this.controlPoints[0].y;
      return createPoint(dx, dy);
    }

    if (this.degree <= 50) {
      const derivativePoints: Point[] = [];
      const n = this.degree;
      for (let i = 0; i < n; i++) {
        derivativePoints.push({
          x: n * (this.controlPoints[i + 1].x - this.controlPoints[i].x),
          y: n * (this.controlPoints[i + 1].y - this.controlPoints[i].y)
        });
      }
      const dc = new BezierCurve(derivativePoints);
      return dc.evaluate(t);
    }

    const h = 1e-6;
    const t1 = Math.max(0, t - h);
    const t2 = Math.min(1, t + h);
    const p1 = this.evaluate(t1);
    const p2 = this.evaluate(t2);
    const denom = t2 - t1;
    if (denom < 1e-15) return createPoint(0, 0);
    return safePoint({
      x: (p2.x - p1.x) / denom,
      y: (p2.y - p1.y) / denom
    });
  }

  tangent(t: number): Point {
    const d = this.derivative(t);
    const len = Math.sqrt(d.x * d.x + d.y * d.y);
    if (len < 1e-15) return createPoint(0, 0);
    return createPoint(d.x / len, d.y / len);
  }

  normal(t: number): Point {
    const tan = this.tangent(t);
    return createPoint(-tan.y, tan.x);
  }

  elevateDegree(): BezierCurve {
    const n = this.degree;
    const newPoints: Point[] = [];
    newPoints.push({ ...this.controlPoints[0] });
    for (let i = 1; i <= n; i++) {
      const alpha = i / (n + 1);
      newPoints.push({
        x: alpha * this.controlPoints[i - 1].x + (1 - alpha) * this.controlPoints[i].x,
        y: alpha * this.controlPoints[i - 1].y + (1 - alpha) * this.controlPoints[i].y
      });
    }
    newPoints.push({ ...this.controlPoints[n] });
    return new BezierCurve(newPoints);
  }

  adaptiveSubdivision(flatnessTolerance: number = 0.5, maxDepth: number = 10): Point[][] {
    if (this.degree > 200) {
      console.warn(`[BezierCurve] adaptiveSubdivision: degree=${this.degree}, falling back to sample-based subdivision.`);
      const samples = this.sample(Math.max(10, Math.min(1000, this.degree / 2)));
      const result: Point[][] = [];
      for (let i = 0; i < samples.length - 1; i++) {
        result.push([samples[i], samples[i + 1]]);
      }
      return result;
    }
    const segments: Point[][] = [];
    this.adaptiveSubdivideRecursive(this.controlPoints, segments, flatnessTolerance, 0, Math.min(maxDepth, 15));
    return segments;
  }

  private adaptiveSubdivideRecursive(
    points: Point[],
    segments: Point[][],
    tolerance: number,
    depth: number,
    maxDepth: number
  ): void {
    if (this.isFlatEnough(points, tolerance) || depth >= maxDepth) {
      segments.push([...points]);
      return;
    }

    const curve = new BezierCurve(points);
    const { left, right } = curve.subdivide(0.5);
    this.adaptiveSubdivideRecursive(left, segments, tolerance, depth + 1, maxDepth);
    this.adaptiveSubdivideRecursive(right, segments, tolerance, depth + 1, maxDepth);
  }

  private isFlatEnough(points: Point[], tolerance: number): boolean {
    if (points.length <= 2) return true;

    const p0 = points[0];
    const pn = points[points.length - 1];
    let maxDist = 0;

    for (let i = 1; i < points.length - 1; i++) {
      const d = this.pointToLineDistance(points[i], p0, pn);
      if (d > maxDist) maxDist = d;
    }

    return maxDist <= tolerance;
  }

  private pointToLineDistance(p: Point, a: Point, b: Point): number {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const lenSq = dx * dx + dy * dy;
    if (lenSq < 1e-15) {
      return Math.sqrt((p.x - a.x) ** 2 + (p.y - a.y) ** 2);
    }
    let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / lenSq;
    t = clamp(t, 0, 1);
    const projX = a.x + t * dx;
    const projY = a.y + t * dy;
    return Math.sqrt((p.x - projX) ** 2 + (p.y - projY) ** 2);
  }

  sample(count: number): Point[] {
    if (count < 2) {
      throw new Error('Sample count must be at least 2');
    }
    const n = Math.min(count, 100000);
    const result: Point[] = [];
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      result.push(this.evaluate(t));
    }
    return result;
  }

  arcLength(t0: number = 0, t1: number = 1): number {
    return this.arcLengthSampler.arcLength(t0, t1);
  }

  getTotalLength(): number {
    return this.arcLengthSampler.getTotalLength();
  }

  sampleByArcLength(count: number): Point[] {
    return this.arcLengthSampler.sampleByArcLength(count);
  }

  getAnalyzer(): CurveAnalyzer {
    return this.analyzer;
  }
}

export class QuadraticBezier extends BezierCurve {
  constructor(p0: Point, p1: Point, p2: Point) {
    super([p0, p1, p2]);
  }
}

export class CubicBezier extends BezierCurve {
  constructor(p0: Point, p1: Point, p2: Point, p3: Point) {
    super([p0, p1, p2, p3]);
  }
}

export class CompositeBezierCurve implements Curve {
  private curves: BezierCurve[];
  private totalLength: number;
  private arcLengthSampler: ArcLengthSampler;
  private segmentStartParams: number[] = [];
  private closed: boolean;
  private analyzer: CurveAnalyzer;

  constructor(curves: BezierCurve[], closed: boolean = false) {
    if (curves.length === 0) {
      throw new Error('CompositeBezierCurve requires at least one Bezier curve');
    }
    this.curves = [...curves];
    this.closed = closed;
    this.segmentStartParams = [0];
    let total = 0;
    for (const curve of curves) {
      total += 1;
      this.segmentStartParams.push(total);
    }

    const nCurves = this.curves.length;
    if (nCurves >= 2000) {
      console.warn(`[CompositeBezierCurve] ${nCurves} segments, performance protection enabled.`);
    }

    let samplerRes = 3000;
    let analyzerRes = 200;
    if (nCurves > 10000) {
      samplerRes = 500;
      analyzerRes = 40;
    } else if (nCurves > 2000) {
      samplerRes = 1000;
      analyzerRes = 80;
    } else if (nCurves > 500) {
      samplerRes = 2000;
      analyzerRes = 120;
    }

    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), samplerRes);
    this.totalLength = this.arcLengthSampler.getTotalLength();
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      analyzerRes
    );
  }

  private findSegment(t: number): { curveIndex: number; localT: number } {
    t = clamp(t, 0, 1);
    const scaledT = t * this.curves.length;
    const index = Math.min(Math.floor(scaledT), this.curves.length - 1);
    const localT = scaledT - index;
    return { curveIndex: index, localT };
  }

  evaluate(t: number): Point {
    const { curveIndex, localT } = this.findSegment(t);
    return this.curves[curveIndex].evaluate(localT);
  }

  derivative(t: number): Point {
    const { curveIndex, localT } = this.findSegment(clamp(t, 0, 1));
    const d = this.curves[curveIndex].derivative(localT);
    const scale = this.curves.length;
    return createPoint(d.x * scale, d.y * scale);
  }

  sample(count: number): Point[] {
    if (count < 2) {
      throw new Error('Sample count must be at least 2');
    }
    const n = Math.min(count, 100000);
    const result: Point[] = [];
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      result.push(this.evaluate(t));
    }
    return result;
  }

  arcLength(t0: number = 0, t1: number = 1): number {
    return this.arcLengthSampler.arcLength(t0, t1);
  }

  getTotalLength(): number {
    return this.totalLength;
  }

  sampleByArcLength(count: number): Point[] {
    return this.arcLengthSampler.sampleByArcLength(count);
  }

  getCurveCount(): number {
    return this.curves.length;
  }

  getCurve(index: number): BezierCurve {
    return this.curves[index];
  }

  isClosed(): boolean {
    return this.closed;
  }

  getAnalyzer(): CurveAnalyzer {
    return this.analyzer;
  }
}
