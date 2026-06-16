import { Point, clamp, distance } from './types';

export class ArcLengthSampler {
  private evaluateFn: (t: number) => Point;
  private cumulativeLengths: number[] = [];
  private parameterValues: number[] = [];
  private totalArcLength: number = 0;
  private resolution: number;

  constructor(evaluateFn: (t: number) => Point, resolution: number = 1000) {
    this.evaluateFn = evaluateFn;
    this.resolution = Math.max(10, resolution);
    this.buildLookupTable();
  }

  private buildLookupTable(): void {
    this.cumulativeLengths = [];
    this.parameterValues = [];
    this.cumulativeLengths.push(0);
    this.parameterValues.push(0);

    let prevPoint = this.evaluateFn(0);
    for (let i = 1; i <= this.resolution; i++) {
      const t = i / this.resolution;
      const point = this.evaluateFn(t);
      const segLen = distance(prevPoint, point);
      const cumLen = this.cumulativeLengths[i - 1] + segLen;
      this.cumulativeLengths.push(cumLen);
      this.parameterValues.push(t);
      prevPoint = point;
    }
    this.totalArcLength = this.cumulativeLengths[this.resolution];
  }

  refresh(): void {
    this.buildLookupTable();
  }

  getTotalLength(): number {
    return this.totalArcLength;
  }

  arcLength(t0: number, t1: number): number {
    t0 = clamp(t0, 0, 1);
    t1 = clamp(t1, 0, 1);
    if (t0 === t1) return 0;
    if (t0 > t1) [t0, t1] = [t1, t0];

    const s0 = this.tToArcLength(t0);
    const s1 = this.tToArcLength(t1);
    return s1 - s0;
  }

  tToArcLength(t: number): number {
    t = clamp(t, 0, 1);
    if (t <= 0) return 0;
    if (t >= 1) return this.totalArcLength;

    const idx = t * this.resolution;
    const i = Math.floor(idx);
    const frac = idx - i;

    if (i >= this.resolution) return this.totalArcLength;

    return this.cumulativeLengths[i] +
      frac * (this.cumulativeLengths[i + 1] - this.cumulativeLengths[i]);
  }

  arcLengthToT(s: number): number {
    s = clamp(s, 0, this.totalArcLength);
    if (s <= 0) return 0;
    if (s >= this.totalArcLength) return 1;

    let lo = 0;
    let hi = this.resolution;
    while (lo < hi) {
      const mid = Math.floor((lo + hi) / 2);
      if (this.cumulativeLengths[mid] < s) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }

    const i = Math.max(0, lo - 1);
    const s_i = this.cumulativeLengths[i];
    const s_next = this.cumulativeLengths[Math.min(i + 1, this.resolution)];
    const t_i = this.parameterValues[i];
    const t_next = this.parameterValues[Math.min(i + 1, this.resolution)];

    if (Math.abs(s_next - s_i) < 1e-15) {
      return t_i;
    }

    const frac = (s - s_i) / (s_next - s_i);
    return t_i + frac * (t_next - t_i);
  }

  private arcLengthToTRefined(s: number): number {
    const tApprox = this.arcLengthToT(s);
    const epsilon = 1e-10;
    const maxIter = 10;

    let t = tApprox;
    for (let i = 0; i < maxIter; i++) {
      const currentS = this.tToArcLength(t);
      const error = s - currentS;
      if (Math.abs(error) < epsilon) break;

      const dt = 1e-7;
      const ds = this.tToArcLength(Math.min(1, t + dt)) - currentS;
      if (Math.abs(ds) < 1e-15) break;

      t = clamp(t + (error / ds) * dt, 0, 1);
    }
    return t;
  }

  sampleByArcLength(count: number): Point[] {
    if (count < 2) {
      throw new Error('Sample count must be at least 2');
    }

    const result: Point[] = [];
    for (let i = 0; i < count; i++) {
      const s = (i / (count - 1)) * this.totalArcLength;
      const t = this.arcLengthToTRefined(s);
      result.push(this.evaluateFn(t));
    }
    return result;
  }

  evaluateByArcLength(s: number): Point {
    const t = this.arcLengthToTRefined(s);
    return this.evaluateFn(t);
  }

  numericalDerivative(t: number, h: number = 1e-6): Point {
    t = clamp(t, 0, 1);
    const t1 = clamp(t - h, 0, 1);
    const t2 = clamp(t + h, 0, 1);
    const p1 = this.evaluateFn(t1);
    const p2 = this.evaluateFn(t2);
    const denom = t2 - t1;
    if (denom < 1e-15) {
      return { x: 0, y: 0 };
    }
    return {
      x: (p2.x - p1.x) / denom,
      y: (p2.y - p1.y) / denom
    };
  }
}

export function numericalIntegration(
  f: (x: number) => number,
  a: number,
  b: number,
  n: number = 1000
): number {
  const h = (b - a) / n;
  let sum = 0.5 * (f(a) + f(b));
  for (let i = 1; i < n; i++) {
    sum += f(a + i * h);
  }
  return sum * h;
}
