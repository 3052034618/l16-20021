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

interface CRSegment {
  a: Point;
  b: Point;
  c: Point;
  d: Point;
  knotStart: number;
  knotEnd: number;
  span: number;
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
  private segments: CRSegment[] = [];

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
    this.buildSegments();

    const nPts = this.controlPoints.length;
    const samplerRes = nPts > 2000 ? 400 : nPts > 500 ? 800 : 2000;
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), samplerRes);
    this.totalLength = this.arcLengthSampler.getTotalLength();
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      nPts > 2000 ? 50 : nPts > 500 ? 100 : 200
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

  private getKnot(index: number): number {
    const n = this.knots.length;
    if (this.closed) {
      if (index < 0) return this.knots[index + n] - 1;
      if (index >= n) return this.knots[index - n] + 1;
      return this.knots[index];
    } else {
      if (index < 0) {
        const h = this.knots[1] - this.knots[0];
        return this.knots[0] - h;
      }
      if (index >= n) {
        const h = this.knots[n - 1] - this.knots[n - 2];
        return this.knots[n - 1] + h;
      }
      return this.knots[index];
    }
  }

  private getSegmentCount(): number {
    return this.closed ? this.controlPoints.length : this.controlPoints.length - 1;
  }

  private computeTangent(
    pPrev: Point, pCurr: Point, pNext: Point,
    hPrev: number, hCurr: number
  ): Point {
    const hp = hPrev < STABILITY_EPSILON ? STABILITY_EPSILON : hPrev;
    const hc = hCurr < STABILITY_EPSILON ? STABILITY_EPSILON : hCurr;
    const denom = hp + hc;

    if (denom < STABILITY_EPSILON) {
      return createPoint(0, 0);
    }

    const term1 = scalePoint(subPoints(pCurr, pPrev), hc / (hp * denom));
    const term2 = scalePoint(subPoints(pNext, pCurr), hp / (hc * denom));

    return addPoints(term1, term2);
  }

  private buildSegments(): void {
    this.segments = [];
    const n = this.getSegmentCount();
    const nPts = this.controlPoints.length;
    const alpha = this.tension;

    for (let i = 0; i < n; i++) {
      let p0Idx: number, p1Idx: number, p2Idx: number, p3Idx: number;

      if (this.closed) {
        p0Idx = (i - 1 + nPts) % nPts;
        p1Idx = i;
        p2Idx = (i + 1) % nPts;
        p3Idx = (i + 2) % nPts;
      } else {
        p0Idx = i - 1;
        p1Idx = i;
        p2Idx = i + 1;
        p3Idx = i + 2;
      }

      const p0 = this.getPoint(p0Idx);
      const p1 = this.controlPoints[this.closed ? p1Idx : clamp(p1Idx, 0, nPts - 1)];
      const p2 = this.controlPoints[this.closed ? p2Idx : clamp(p2Idx, 0, nPts - 1)];
      const p3 = this.getPoint(p3Idx);

      const knotStart = this.knots[i] ?? 0;
      const knotEnd = this.knots[i + 1] ?? 1;
      let hi = knotEnd - knotStart;
      if (hi < STABILITY_EPSILON) hi = STABILITY_EPSILON;

      const hPrev = this.getKnot(i) - this.getKnot(i - 1);
      const hNext = this.getKnot(i + 2) - this.getKnot(i + 1);

      const mCurr = this.computeTangent(p0, p1, p2, hPrev, hi);
      const mNext = this.computeTangent(p1, p2, p3, hi, hNext);

      const deltaP = subPoints(p2, p1);

      const a = { ...p1 };
      const b = scalePoint(mCurr, alpha * hi);
      const c = subPoints(
        scalePoint(deltaP, 3),
        addPoints(
          scalePoint(mCurr, alpha * hi * 2),
          scalePoint(mNext, alpha * hi)
        )
      );
      const d = addPoints(
        scalePoint(deltaP, -2),
        addPoints(
          scalePoint(mCurr, alpha * hi),
          scalePoint(mNext, alpha * hi)
        )
      );

      this.segments.push({
        a, b, c, d,
        knotStart,
        knotEnd,
        span: hi
      });
    }
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
    const seg = this.segments[idx];
    const span = seg?.span ?? 1;
    const kStart = seg?.knotStart ?? 0;
    const localT = span > STABILITY_EPSILON ? (t - kStart) / span : 0;

    return { segmentIndex: idx, localT: safeNumber(localT, 0) };
  }

  evaluate(t: number): Point {
    const { segmentIndex, localT } = this.findSegment(t);
    const seg = this.segments[segmentIndex];
    if (!seg) return safePoint(this.controlPoints[0] ?? createPoint(0, 0));
    const s = localT;
    const s2 = s * s;
    const s3 = s2 * s;
    return safePoint({
      x: seg.a.x + seg.b.x * s + seg.c.x * s2 + seg.d.x * s3,
      y: seg.a.y + seg.b.y * s + seg.c.y * s2 + seg.d.y * s3
    });
  }

  derivative(t: number): Point {
    const { segmentIndex, localT } = this.findSegment(clamp(t, 0, 1));
    const seg = this.segments[segmentIndex];
    if (!seg) return createPoint(0, 0);
    const s = localT;
    const s2 = s * s;
    const hi = seg.span > STABILITY_EPSILON ? seg.span : STABILITY_EPSILON;
    return safePoint({
      x: (seg.b.x + 2 * seg.c.x * s + 3 * seg.d.x * s2) / hi,
      y: (seg.b.y + 2 * seg.c.y * s + 3 * seg.d.y * s2) / hi
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
    const n = Math.min(count, 100000);
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
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

    for (const seg of this.segments) {
      const a = seg.a, b = seg.b, c = seg.c, d = seg.d;

      const p0 = { ...a };
      const c1 = addPoints(a, scalePoint(b, 1 / 3));
      const c2 = addPoints(addPoints(a, scalePoint(b, 2 / 3)), scalePoint(c, 1 / 3));
      const p3 = addPoints(addPoints(addPoints(a, b), c), d);

      segments.push(new CubicBezier(p0, c1, c2, p3));
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

      const dBefore = this.derivative(tBefore);
      const dAfter = this.derivative(tAfter);
      const lenBefore = Math.sqrt(dBefore.x * dBefore.x + dBefore.y * dBefore.y);
      const lenAfter = Math.sqrt(dAfter.x * dAfter.x + dAfter.y * dAfter.y);

      if (lenBefore < STABILITY_EPSILON && lenAfter < STABILITY_EPSILON) {
        continue;
      }

      if (lenBefore > STABILITY_EPSILON && lenAfter > STABILITY_EPSILON) {
        const dotP = (dBefore.x * dAfter.x + dBefore.y * dAfter.y) / (lenBefore * lenAfter);
        if (Math.abs(dotP - 1) > tolerance) return false;
      }

      if (lenBefore > STABILITY_EPSILON && lenAfter > STABILITY_EPSILON) {
        const speedRatio = Math.abs(lenBefore - lenAfter) / Math.max(lenBefore, lenAfter);
        if (speedRatio > tolerance * 1000) {
          return false;
        }
      }
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
