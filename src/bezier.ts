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
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 2000);
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      200
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
    return this.deCasteljau(this.controlPoints, t);
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
      for (let i = 0; i < current.length - 1; i++) {
        next.push(lerpPoint(current[i], current[i + 1], t));
      }
      current = next;
    }
    return { ...current[0] };
  }

  evaluateBernstein(t: number): Point {
    t = clamp(t, 0, 1);
    const n = this.degree;
    let x = 0, y = 0;
    for (let i = 0; i <= n; i++) {
      const b = bernsteinPolynomial(n, i, t);
      x += this.controlPoints[i].x * b;
      y += this.controlPoints[i].y * b;
    }
    return createPoint(x, y);
  }

  getDeCasteljauTable(t: number): Point[][] {
    t = clamp(t, 0, 1);
    const table: Point[][] = [];
    table.push([...this.controlPoints]);

    for (let level = 1; level <= this.degree; level++) {
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
    const table = this.getDeCasteljauTable(t);
    const left: Point[] = [];
    const right: Point[] = [];

    for (let i = 0; i <= this.degree; i++) {
      left.push(table[i][0]);
    }
    for (let i = 0; i <= this.degree; i++) {
      right.push(table[this.degree - i][i]);
    }
    return { left, right };
  }

  split(t: number): [BezierCurve, BezierCurve] {
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

    const derivativePoints: Point[] = [];
    const n = this.degree;
    for (let i = 0; i < n; i++) {
      derivativePoints.push({
        x: n * (this.controlPoints[i + 1].x - this.controlPoints[i].x),
        y: n * (this.controlPoints[i + 1].y - this.controlPoints[i].y)
      });
    }

    const derivativeCurve = new BezierCurve(derivativePoints);
    return derivativeCurve.evaluate(t);
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
    const segments: Point[][] = [];
    this.adaptiveSubdivideRecursive(this.controlPoints, segments, flatnessTolerance, 0, maxDepth);
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
    const result: Point[] = [];
    for (let i = 0; i < count; i++) {
      const t = i / (count - 1);
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
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 3000);
    this.totalLength = this.arcLengthSampler.getTotalLength();
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      200
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
    const result: Point[] = [];
    for (let i = 0; i < count; i++) {
      const t = i / (count - 1);
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
