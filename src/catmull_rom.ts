import {
  Point, Curve, createPoint, clamp, addPoints, subPoints, scalePoint, distance
} from './types';
import { ArcLengthSampler } from './arc_length';
import { CubicBezier, BezierCurve } from './bezier';
import { CurveAnalyzer } from './curve_analyzer';
import {
  cleanInputPoints, validateKnots, autoKnots, safePoint, safeNumber, STABILITY_EPSILON
} from './stability';

export interface CatmullRomOptions {
  tension?: number;
  closed?: boolean;
  knots?: number[];
  knotType?: 'uniform' | 'chordal' | 'centripetal';
  timeStamps?: number[];
}

export class CatmullRomSpline implements Curve {
  private controlPoints: Point[];
  private tension: number;
  private arcLengthSampler: ArcLengthSampler;
  private analyzer: CurveAnalyzer;
  private totalLength: number;
  private closed: boolean;
  private knots: number[] = [];
  private options: CatmullRomOptions;

  constructor(points: Point[], options: CatmullRomOptions = {}) {
    this.options = { ...options };
    const cleaned = cleanInputPoints(points);

    if (cleaned.points.length < 2) {
      throw new Error('CatmullRomSpline requires at least 2 valid distinct control points');
    }

    this.controlPoints = cleaned.points;
    this.tension = options.tension ?? 0.5;
    this.closed = options.closed ?? false;

    if (this.closed && cleaned.points.length < 3) {
      throw new Error('Closed CatmullRomSpline requires at least 3 control points');
    }

    this.setupKnots();

    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 2000);
    this.totalLength = this.arcLengthSampler.getTotalLength();
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      200
    );
  }

  private setupKnots(): void {
    const n = this.getSegmentCount() + 1;

    if (this.options.timeStamps) {
      const ts = this.options.timeStamps;
      const nPts = this.controlPoints.length;
      const validation = validateKnots(ts, nPts);
      if (!validation.valid) {
        console.warn(`[CatmullRom] Invalid timeStamps: ${validation.reason}, falling back.`);
        this.knots = this.buildKnotsFromSegmentCount(this.getSegmentCount());
      } else {
        let minT = Infinity, maxT = -Infinity;
        for (const t of ts) {
          if (t < minT) minT = t;
          if (t > maxT) maxT = t;
        }
        const range = maxT - minT;
        if (range < STABILITY_EPSILON) {
          this.knots = this.buildKnotsFromSegmentCount(this.getSegmentCount());
        } else {
          this.knots = ts.map(t => (t - minT) / range);
        }
      }
      return;
    }

    if (this.options.knots) {
      const kts = this.options.knots;
      const nPts = this.controlPoints.length;
      const validation = validateKnots(kts, nPts);
      if (!validation.valid) {
        console.warn(`[CatmullRom] Invalid knots: ${validation.reason}, falling back.`);
        this.knots = this.buildKnotsFromSegmentCount(this.getSegmentCount());
      } else {
        this.knots = kts;
      }
      return;
    }

    this.knots = autoKnots(
      this.getSegmentCount() + 1,
      this.options.knotType ?? 'uniform',
      this.closed
        ? [...this.controlPoints, this.controlPoints[0]]
        : this.controlPoints
    );
  }

  private buildKnotsFromSegmentCount(nSegments: number): number[] {
    const result: number[] = [];
    for (let i = 0; i <= nSegments; i++) {
      result.push(i / nSegments);
    }
    return result;
  }

  getKnots(): number[] {
    return [...this.knots];
  }

  private getPoint(index: number): Point {
    const n = this.controlPoints.length;
    if (this.closed) {
      return this.controlPoints[((index % n) + n) % n];
    } else {
      if (index < 0) {
        const first = this.controlPoints[0];
        const second = this.controlPoints[1];
        return {
          x: 2 * first.x - second.x,
          y: 2 * first.y - second.y
        };
      }
      if (index >= n) {
        const last = this.controlPoints[n - 1];
        const secondLast = this.controlPoints[n - 2];
        return {
          x: 2 * last.x - secondLast.x,
          y: 2 * last.y - secondLast.y
        };
      }
      return this.controlPoints[index];
    }
  }

  private getSegmentCount(): number {
    return this.closed ? this.controlPoints.length : this.controlPoints.length - 1;
  }

  private findSegment(t: number): { segmentIndex: number; localT: number } {
    t = clamp(t, 0, 1);
    const n = this.getSegmentCount();
    if (n === 0) return { segmentIndex: 0, localT: 0 };

    if (t >= 1) {
      return { segmentIndex: n - 1, localT: 1 };
    }

    let lo = 0, hi = n - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.knots[mid + 1] < t) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }

    const idx = lo;
    const kStart = this.knots[idx] ?? 0;
    const kEnd = this.knots[idx + 1] ?? 1;
    const span = kEnd - kStart;
    const localT = span > STABILITY_EPSILON ? (t - kStart) / span : 0;

    return { segmentIndex: idx, localT: safeNumber(localT, 0) };
  }

  evaluate(t: number): Point {
    const { segmentIndex, localT } = this.findSegment(t);
    const n = this.controlPoints.length;

    let p0Idx: number, p1Idx: number, p2Idx: number, p3Idx: number;

    if (this.closed) {
      p0Idx = (segmentIndex - 1 + n) % n;
      p1Idx = segmentIndex;
      p2Idx = (segmentIndex + 1) % n;
      p3Idx = (segmentIndex + 2) % n;
    } else {
      p0Idx = segmentIndex - 1;
      p1Idx = segmentIndex;
      p2Idx = segmentIndex + 1;
      p3Idx = segmentIndex + 2;
    }

    const p0 = this.getPoint(p0Idx);
    const p1 = this.controlPoints[this.closed ? p1Idx : Math.min(Math.max(p1Idx, 0), n - 1)];
    const p2 = this.controlPoints[this.closed ? p2Idx : Math.min(Math.max(p2Idx, 0), n - 1)];
    const p3 = this.getPoint(p3Idx);

    return safePoint(this.catmullRom(p0, p1, p2, p3, localT, this.tension));
  }

  private catmullRom(
    p0: Point, p1: Point, p2: Point, p3: Point,
    t: number, alpha: number = 0.5
  ): Point {
    const t2 = t * t;
    const t3 = t2 * t;

    const a0 = p1;
    const a1 = {
      x: alpha * (p2.x - p0.x),
      y: alpha * (p2.y - p0.y)
    };
    const a2 = {
      x: 2 * alpha * p0.x + (alpha - 3) * p1.x + (3 - 2 * alpha) * p2.x - alpha * p3.x,
      y: 2 * alpha * p0.y + (alpha - 3) * p1.y + (3 - 2 * alpha) * p2.y - alpha * p3.y
    };
    const a3 = {
      x: -alpha * p0.x + (2 - alpha) * p1.x + (alpha - 2) * p2.x + alpha * p3.x,
      y: -alpha * p0.y + (2 - alpha) * p1.y + (alpha - 2) * p2.y + alpha * p3.y
    };

    return {
      x: a0.x + a1.x * t + a2.x * t2 + a3.x * t3,
      y: a0.y + a1.y * t + a2.y * t2 + a3.y * t3
    };
  }

  private catmullRomDerivative(
    p0: Point, p1: Point, p2: Point, p3: Point,
    t: number, alpha: number = 0.5
  ): Point {
    const t2 = t * t;

    const a1 = {
      x: alpha * (p2.x - p0.x),
      y: alpha * (p2.y - p0.y)
    };
    const a2 = {
      x: 2 * alpha * p0.x + (alpha - 3) * p1.x + (3 - 2 * alpha) * p2.x - alpha * p3.x,
      y: 2 * alpha * p0.y + (alpha - 3) * p1.y + (3 - 2 * alpha) * p2.y - alpha * p3.y
    };
    const a3 = {
      x: -alpha * p0.x + (2 - alpha) * p1.x + (alpha - 2) * p2.x + alpha * p3.x,
      y: -alpha * p0.y + (2 - alpha) * p1.y + (alpha - 2) * p2.y + alpha * p3.y
    };

    return {
      x: a1.x + 2 * a2.x * t + 3 * a3.x * t2,
      y: a1.y + 2 * a2.y * t + 3 * a3.y * t2
    };
  }

  derivative(t: number): Point {
    const { segmentIndex, localT } = this.findSegment(clamp(t, 0, 1));
    const n = this.controlPoints.length;

    let p0Idx: number, p1Idx: number, p2Idx: number, p3Idx: number;

    if (this.closed) {
      p0Idx = (segmentIndex - 1 + n) % n;
      p1Idx = segmentIndex;
      p2Idx = (segmentIndex + 1) % n;
      p3Idx = (segmentIndex + 2) % n;
    } else {
      p0Idx = segmentIndex - 1;
      p1Idx = segmentIndex;
      p2Idx = segmentIndex + 1;
      p3Idx = segmentIndex + 2;
    }

    const p0 = this.getPoint(p0Idx);
    const p1 = this.controlPoints[this.closed ? p1Idx : Math.min(Math.max(p1Idx, 0), n - 1)];
    const p2 = this.controlPoints[this.closed ? p2Idx : Math.min(Math.max(p2Idx, 0), n - 1)];
    const p3 = this.getPoint(p3Idx);

    const nSegments = this.getSegmentCount();
    const kStart = this.knots[segmentIndex] ?? 0;
    const kEnd = this.knots[segmentIndex + 1] ?? 1;
    const span = kEnd - kStart;
    const spanFactor = span > STABILITY_EPSILON ? 1 / span : nSegments;

    const d = this.catmullRomDerivative(p0, p1, p2, p3, localT, this.tension);
    return safePoint({
      x: d.x * spanFactor,
      y: d.y * spanFactor
    });
  }

  tangent(t: number): Point {
    const d = this.derivative(t);
    const len = Math.sqrt(d.x * d.x + d.y * d.y);
    if (len < STABILITY_EPSILON) return createPoint(0, 0);
    return createPoint(d.x / len, d.y / len);
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

  sampleByTimestamps(timeStamps: number[]): Point[] {
    if (!this.options.timeStamps) {
      return this.sample(timeStamps.length);
    }
    const ts = this.options.timeStamps;
    let minT = Infinity, maxT = -Infinity;
    for (const t of ts) {
      if (t < minT) minT = t;
      if (t > maxT) maxT = t;
    }
    const range = maxT - minT;
    if (range < STABILITY_EPSILON) {
      return timeStamps.map(() => this.evaluate(0));
    }
    return timeStamps.map(t => {
      const normT = (t - minT) / range;
      return this.evaluate(normT);
    });
  }

  arcLength(t0: number = 0, t1: number = 1): number {
    return safeNumber(this.arcLengthSampler.arcLength(t0, t1), 0);
  }

  getTotalLength(): number {
    return safeNumber(this.totalLength, 0);
  }

  sampleByArcLength(count: number): Point[] {
    return this.arcLengthSampler.sampleByArcLength(count);
  }

  toBezierSegments(): CubicBezier[] {
    const segments: CubicBezier[] = [];
    const n = this.getSegmentCount();
    const alpha = this.tension;
    const cp = this.controlPoints;
    const numPoints = cp.length;

    for (let i = 0; i < n; i++) {
      let p0Idx: number, p1Idx: number, p2Idx: number, p3Idx: number;

      if (this.closed) {
        p0Idx = (i - 1 + numPoints) % numPoints;
        p1Idx = i;
        p2Idx = (i + 1) % numPoints;
        p3Idx = (i + 2) % numPoints;
      } else {
        p0Idx = i - 1;
        p1Idx = i;
        p2Idx = i + 1;
        p3Idx = i + 2;
      }

      const p0 = this.getPoint(p0Idx);
      const p1 = cp[this.closed ? p1Idx : Math.min(Math.max(p1Idx, 0), numPoints - 1)];
      const p2 = cp[this.closed ? p2Idx : Math.min(Math.max(p2Idx, 0), numPoints - 1)];
      const p3 = this.getPoint(p3Idx);

      const tangent1 = scalePoint(subPoints(p2, p0), alpha);
      const tangent2 = scalePoint(subPoints(p3, p1), alpha);

      const c1 = addPoints(p1, scalePoint(tangent1, 1 / 3));
      const c2 = subPoints(p2, scalePoint(tangent2, 1 / 3));

      segments.push(new CubicBezier(p1, c1, c2, p2));
    }

    return segments;
  }

  getControlPoints(): Point[] {
    return [...this.controlPoints];
  }

  getTension(): number {
    return this.tension;
  }

  isClosed(): boolean {
    return this.closed;
  }

  getAnalyzer(): CurveAnalyzer {
    return this.analyzer;
  }

  verifyInterpolation(tolerance: number = 1e-9): boolean {
    const n = this.controlPoints.length;
    const nSegments = this.getSegmentCount();

    const checkPoints = this.closed ? n : n;
    for (let i = 0; i < checkPoints; i++) {
      const t = this.knots[i] ?? (i / nSegments);
      const expected = this.controlPoints[this.closed ? i % n : i];

      const result = this.evaluate(t);
      const err = Math.max(
        Math.abs(result.x - expected.x),
        Math.abs(result.y - expected.y)
      );
      if (err > tolerance) return false;
    }

    return true;
  }

  checkC1Continuity(tolerance: number = 1e-6): boolean {
    const nSegments = this.getSegmentCount();
    const limit = this.closed ? nSegments : nSegments - 1;

    for (let i = 1; i <= limit; i++) {
      const t = this.knots[i] ?? (i / nSegments);
      const tBefore = Math.max(0, t - 1e-8);
      const tAfter = Math.min(1, t + 1e-8);

      const dBefore = this.tangent(tBefore);
      const dAfter = this.tangent(tAfter);

      const dotP = dBefore.x * dAfter.x + dBefore.y * dAfter.y;
      if (Math.abs(dotP - 1) > tolerance) return false;
    }

    return true;
  }
}

export function catmullRomToBezier(
  points: Point[],
  tension: number = 0.5,
  closed: boolean = false
): CubicBezier[] {
  const spline = new CatmullRomSpline(points, { tension, closed });
  return spline.toBezierSegments();
}
