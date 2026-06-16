import {
  Point, Curve, CubicCoefficients, SplineSegment,
  SplineOptions, EndCondition, createPoint, clamp, EPSILON
} from './types';
import { solveTridiagonal, createTridiagonalSystem } from './tridiagonal';
import { ArcLengthSampler } from './arc_length';

export class CubicSpline implements Curve {
  private controlPoints: Point[];
  private segments: SplineSegment[] = [];
  private arcLengthSampler: ArcLengthSampler;
  private options: Required<SplineOptions>;
  private totalLength: number = 0;

  constructor(points: Point[], options: SplineOptions = {}) {
    if (points.length < 2) {
      throw new Error('CubicSpline requires at least 2 control points');
    }

    this.options = {
      endCondition: options.endCondition || 'natural',
      startTangent: options.startTangent || { x: 0, y: 0 },
      endTangent: options.endTangent || { x: 0, y: 0 }
    };

    this.controlPoints = [...points];
    this.buildSpline();
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 2000);
    this.totalLength = this.arcLengthSampler.getTotalLength();
  }

  private buildSpline(): void {
    const n = this.controlPoints.length - 1;
    if (n === 0) return;

    const h: number[] = [];
    for (let i = 0; i < n; i++) {
      h.push(1);
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
        tStart: i / n,
        tEnd: (i + 1) / n
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
      const slope = (y[1] - y[0]) / h[0];
      return [{
        a: y[0],
        b: slope,
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
      d[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }

    this.applyBoundaryConditions(a, b, c, d, h, y, startTangent, endTangent);

    const M = solveTridiagonal(system);

    for (let i = 0; i < n; i++) {
      const a_coeff = y[i];
      const b_coeff = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * M[i] + M[i + 1]) / 6;
      const c_coeff = M[i] / 2;
      const d_coeff = (M[i + 1] - M[i]) / (6 * h[i]);
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
        b[0] = 2 * h[0];
        c[0] = h[0];
        d[0] = 6 * ((y[1] - y[0]) / h[0] - startTangent);
        a[n] = h[n - 1];
        b[n] = 2 * h[n - 1];
        d[n] = 6 * (endTangent - (y[n] - y[n - 1]) / h[n - 1]);
        break;

      case 'not-a-knot':
        if (n >= 2) {
          b[0] = h[1];
          c[0] = -(h[0] + h[1]);
          d[0] = h[0];
          const d0 = 0;

          a[n] = h[n - 2];
          b[n] = -(h[n - 2] + h[n - 1]);
          d[n] = h[n - 1];
          const dn = 0;
          Object.assign(d, { [0]: d0, [n]: dn });
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
    if (t >= 1) {
      const idx = this.segments.length - 1;
      const seg = this.segments[idx];
      return { segment: seg, localT: (t - seg.tStart) / (seg.tEnd - seg.tStart), index: idx };
    }

    const n = this.segments.length;
    const index = Math.min(Math.floor(t * n), n - 1);
    const seg = this.segments[index];
    const localT = (t - seg.tStart) / (seg.tEnd - seg.tStart);
    return { segment: seg, localT, index };
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
    const n = this.segments.length;
    const scaledLocalT = localT;
    return createPoint(
      this.evaluatePolynomial(segment.xCoeff, scaledLocalT),
      this.evaluatePolynomial(segment.yCoeff, scaledLocalT)
    );
  }

  derivative(t: number): Point {
    const { segment, localT } = this.findSegment(clamp(t, 0, 1));
    const n = this.segments.length;
    const scaledLocalT = localT;
    const dt = 1.0 / n;
    return createPoint(
      this.evaluateDerivative(segment.xCoeff, scaledLocalT) / dt,
      this.evaluateDerivative(segment.yCoeff, scaledLocalT) / dt
    );
  }

  secondDerivative(t: number): Point {
    const { segment, localT } = this.findSegment(clamp(t, 0, 1));
    const n = this.segments.length;
    const scaledLocalT = localT;
    const dt = 1.0 / n;
    return createPoint(
      this.evaluateSecondDerivative(segment.xCoeff, scaledLocalT) / (dt * dt),
      this.evaluateSecondDerivative(segment.yCoeff, scaledLocalT) / (dt * dt)
    );
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

    for (let i = 1; i < this.controlPoints.length - 1; i++) {
      const t = i / this.segments.length;
      const tMinus = t - EPSILON;
      const tPlus = t + EPSILON;

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
