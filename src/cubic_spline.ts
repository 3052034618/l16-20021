import {
  Point, Curve, CubicCoefficients, SplineSegment,
  SplineOptions, EndCondition, createPoint, clamp, EPSILON
} from './types';
import { solveTridiagonal, createTridiagonalSystem } from './tridiagonal';
import { ArcLengthSampler } from './arc_length';
import { CurveAnalyzer } from './curve_analyzer';
import {
  cleanInputPoints, validateKnots, autoKnots, safePoint, safeNumber, STABILITY_EPSILON
} from './stability';

export interface CubicSplineExtendedOptions extends SplineOptions {
  knots?: number[];
  knotType?: 'uniform' | 'chordal' | 'centripetal';
  timeStamps?: number[];
}

export class CubicSpline implements Curve {
  private controlPoints: Point[];
  private segments: SplineSegment[] = [];
  private arcLengthSampler: ArcLengthSampler;
  private analyzer: CurveAnalyzer;
  private options: Required<SplineOptions>;
  private extendedOptions: CubicSplineExtendedOptions;
  private totalLength: number = 0;
  private knots: number[] = [];

  constructor(points: Point[], options: CubicSplineExtendedOptions = {}) {
    this.extendedOptions = { ...options };
    const cleaned = cleanInputPoints(points);

    if (cleaned.points.length < 2) {
      throw new Error('CubicSpline requires at least 2 valid distinct control points');
    }

    this.options = {
      endCondition: options.endCondition || 'natural',
      startTangent: options.startTangent || { x: 0, y: 0 },
      endTangent: options.endTangent || { x: 0, y: 0 }
    };

    this.controlPoints = cleaned.points;

    this.setupKnots();
    this.buildSpline();
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 2000);
    this.totalLength = this.arcLengthSampler.getTotalLength();
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      200
    );
  }

  private setupKnots(): void {
    const n = this.controlPoints.length;

    if (this.extendedOptions.timeStamps) {
      const ts = this.extendedOptions.timeStamps;
      const validation = validateKnots(ts, n);
      if (!validation.valid) {
        console.warn(`[CubicSpline] Invalid timeStamps: ${validation.reason}, falling back to uniform.`);
        this.knots = autoKnots(n, 'uniform');
      } else {
        let minT = Infinity, maxT = -Infinity;
        for (const t of ts) {
          if (t < minT) minT = t;
          if (t > maxT) maxT = t;
        }
        const range = maxT - minT;
        if (range < STABILITY_EPSILON) {
          this.knots = autoKnots(n, 'uniform');
        } else {
          this.knots = ts.map(t => (t - minT) / range);
        }
      }
      return;
    }

    if (this.extendedOptions.knots) {
      const kts = this.extendedOptions.knots;
      const validation = validateKnots(kts, n);
      if (!validation.valid) {
        console.warn(`[CubicSpline] Invalid knots: ${validation.reason}, falling back to uniform.`);
        this.knots = autoKnots(n, 'uniform');
      } else {
        this.knots = kts;
      }
      return;
    }

    this.knots = autoKnots(n, this.extendedOptions.knotType ?? 'uniform', this.controlPoints);
  }

  getKnots(): number[] {
    return [...this.knots];
  }

  private buildSpline(): void {
    const n = this.controlPoints.length - 1;
    if (n === 0) return;

    const h: number[] = [];
    for (let i = 0; i < n; i++) {
      const diff = this.knots[i + 1] - this.knots[i];
      h.push(diff < STABILITY_EPSILON ? STABILITY_EPSILON : diff);
    }

    const xValues = this.controlPoints.map(p => p.x);
    const yValues = this.controlPoints.map(p => p.y);

    const endCond = this.options.endCondition;
    let xStartTan = 0, xEndTan = 0, yStartTan = 0, yEndTan = 0;
    if (endCond === 'clamped') {
      xStartTan = this.options.startTangent.x;
      xEndTan = this.options.endTangent.x;
      yStartTan = this.options.startTangent.y;
      yEndTan = this.options.endTangent.y;
    }

    const xCoeffs = this.computeSplineCoefficients(xValues, h, xStartTan, xEndTan);
    const yCoeffs = this.computeSplineCoefficients(yValues, h, yStartTan, yEndTan);

    this.segments = [];
    for (let i = 0; i < n; i++) {
      this.segments.push({
        xCoeff: xCoeffs[i],
        yCoeff: yCoeffs[i],
        tStart: this.knots[i],
        tEnd: this.knots[i + 1]
      });
    }
  }

  private computeSplineCoefficients(
    y: number[], h: number[],
    startTangent: number = 0, endTangent: number = 0
  ): CubicCoefficients[] {
    const n = y.length - 1;
    const coeffs: CubicCoefficients[] = [];

    if (n === 0) {
      return [{ a: y[0], b: 0, c: 0, d: 0 }];
    }

    if (n === 1) {
      return [{
        a: y[0],
        b: y[1] - y[0],
        c: 0,
        d: 0
      }];
    }

    const system = createTridiagonalSystem(n + 1);
    const { a, b, c, d } = system;

    for (let i = 1; i < n; i++) {
      a[i] = h[i - 1];
      b[i] = 2 * (h[i - 1] + h[i]);
      c[i] = h[i];
      const dy1 = h[i] > STABILITY_EPSILON ? (y[i + 1] - y[i]) / h[i] : 0;
      const dy2 = h[i - 1] > STABILITY_EPSILON ? (y[i] - y[i - 1]) / h[i - 1] : 0;
      d[i] = 6 * (dy1 - dy2);
    }

    this.applyBoundaryConditions(a, b, c, d, h, y, startTangent, endTangent);

    const M = solveTridiagonal(system);

    for (let i = 0; i < n; i++) {
      const hi = h[i] > STABILITY_EPSILON ? h[i] : STABILITY_EPSILON;
      const hi2 = hi * hi;
      const a_coeff = y[i];
      const b_coeff = (y[i + 1] - y[i]) - hi2 * (2 * M[i] + M[i + 1]) / 6;
      const c_coeff = M[i] * hi2 / 2;
      const d_coeff = (M[i + 1] - M[i]) * hi2 / 6;
      coeffs.push({ a: a_coeff, b: b_coeff, c: c_coeff, d: d_coeff });
    }

    return coeffs;
  }

  private applyBoundaryConditions(
    a: number[], b: number[], c: number[], d: number[],
    h: number[], y: number[],
    startTangent: number, endTangent: number
  ): void {
    const n = y.length - 1;
    const endCond = this.options.endCondition;

    switch (endCond) {
      case 'natural':
        b[0] = 1;
        c[0] = 0;
        d[0] = 0;
        a[n] = 0;
        b[n] = 1;
        d[n] = 0;
        break;

      case 'clamped':
        const h0 = h[0] > STABILITY_EPSILON ? h[0] : STABILITY_EPSILON;
        const hn = h[n - 1] > STABILITY_EPSILON ? h[n - 1] : STABILITY_EPSILON;
        b[0] = 2 * h0;
        c[0] = h0;
        d[0] = 6 * ((y[1] - y[0]) / h0 - startTangent);
        a[n] = hn;
        b[n] = 2 * hn;
        d[n] = 6 * (endTangent - (y[n] - y[n - 1]) / hn);
        break;

      case 'not-a-knot':
        if (n >= 2) {
          b[0] = h[1];
          c[0] = -(h[0] + h[1]);
          const tmp = h[0];
          d[0] = tmp;
          a[n] = h[n - 2];
          b[n] = -(h[n - 2] + h[n - 1]);
          d[n] = h[n - 1];
          Object.assign(d, { [0]: 0, [n]: 0 });
        } else {
          b[0] = 1;
          c[0] = 0;
          d[0] = 0;
          a[n] = 0;
          b[n] = 1;
          d[n] = 0;
        }
        break;

      default:
        b[0] = 1;
        c[0] = 0;
        d[0] = 0;
        a[n] = 0;
        b[n] = 1;
        d[n] = 0;
    }
  }

  private findSegment(t: number): { segment: SplineSegment; localT: number; index: number } {
    t = clamp(t, 0, 1);
    if (this.segments.length === 0) {
      throw new Error('No segments built');
    }

    let lo = 0, hi = this.segments.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.segments[mid].tEnd < t) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }

    const idx = lo;
    const seg = this.segments[idx];
    const segSpan = seg.tEnd - seg.tStart;
    const localT = segSpan > STABILITY_EPSILON ? (t - seg.tStart) / segSpan : 0;
    return { segment: seg, localT: safeNumber(localT, 0), index: idx };
  }

  private evaluatePolynomial(coeff: CubicCoefficients, t: number): number {
    return coeff.a + coeff.b * t + coeff.c * t * t + coeff.d * t * t * t;
  }

  private evaluateDerivative(coeff: CubicCoefficients, t: number): number {
    return coeff.b + 2 * coeff.c * t + 3 * coeff.d * t * t;
  }

  private evaluateSecondDerivative(coeff: CubicCoefficients, t: number): number {
    return 2 * coeff.c + 6 * coeff.d * t;
  }

  evaluate(t: number): Point {
    const { segment, localT } = this.findSegment(t);
    return safePoint(createPoint(
      this.evaluatePolynomial(segment.xCoeff, localT),
      this.evaluatePolynomial(segment.yCoeff, localT)
    ));
  }

  derivative(t: number): Point {
    const { segment, localT } = this.findSegment(clamp(t, 0, 1));
    const segSpan = segment.tEnd - segment.tStart;
    const dt = segSpan > STABILITY_EPSILON ? segSpan : 1;
    return safePoint(createPoint(
      this.evaluateDerivative(segment.xCoeff, localT) / dt,
      this.evaluateDerivative(segment.yCoeff, localT) / dt
    ));
  }

  secondDerivative(t: number): Point {
    const { segment, localT } = this.findSegment(clamp(t, 0, 1));
    const segSpan = segment.tEnd - segment.tStart;
    const dt = segSpan > STABILITY_EPSILON ? segSpan : 1;
    return safePoint(createPoint(
      this.evaluateSecondDerivative(segment.xCoeff, localT) / (dt * dt),
      this.evaluateSecondDerivative(segment.yCoeff, localT) / (dt * dt)
    ));
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
    if (!this.extendedOptions.timeStamps) {
      return this.sample(timeStamps.length);
    }
    const ts = this.extendedOptions.timeStamps;
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

  getSegments(): SplineSegment[] {
    return this.segments.map(s => ({
      xCoeff: { ...s.xCoeff },
      yCoeff: { ...s.yCoeff },
      tStart: s.tStart,
      tEnd: s.tEnd
    }));
  }

  getControlPoints(): Point[] {
    return [...this.controlPoints];
  }

  getSegmentCount(): number {
    return this.segments.length;
  }

  getAnalyzer(): CurveAnalyzer {
    return this.analyzer;
  }

  checkContinuity(tolerance: number = 1e-9): {
    position: boolean;
    firstDerivative: boolean;
    secondDerivative: boolean;
    errors: { t: number; type: string; error: number }[];
  } {
    const errors: { t: number; type: string; error: number }[] = [];
    let posOk = true;
    let firstOk = true;
    let secondOk = true;

    for (let i = 1; i < this.knots.length - 1; i++) {
      const t = this.knots[i];
      const tMinus = Math.max(0, t - EPSILON);
      const tPlus = Math.min(1, t + EPSILON);

      const pLeft = this.evaluate(tMinus);
      const pRight = this.evaluate(tPlus);
      const posErr = Math.max(Math.abs(pLeft.x - pRight.x), Math.abs(pLeft.y - pRight.y));
      if (posErr > tolerance) {
        posOk = false;
        errors.push({ t, type: 'position', error: posErr });
      }

      const dLeft = this.derivative(tMinus);
      const dRight = this.derivative(tPlus);
      const dErr = Math.max(Math.abs(dLeft.x - dRight.x), Math.abs(dLeft.y - dRight.y));
      if (dErr > tolerance) {
        firstOk = false;
        errors.push({ t, type: 'first derivative', error: dErr });
      }

      const sLeft = this.secondDerivative(tMinus);
      const sRight = this.secondDerivative(tPlus);
      const sErr = Math.max(Math.abs(sLeft.x - sRight.x), Math.abs(sLeft.y - sRight.y));
      if (sErr > tolerance) {
        secondOk = false;
        errors.push({ t, type: 'second derivative', error: sErr });
      }
    }

    return { position: posOk, firstDerivative: firstOk, secondDerivative: secondOk, errors };
  }
}
