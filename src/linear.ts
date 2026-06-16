import { Point, Curve, lerpPoint, distance, createPoint, clamp, EPSILON } from './types';
import { ArcLengthSampler } from './arc_length';
import { CurveAnalyzer } from './curve_analyzer';
import { cleanInputPoints, safePoint, safeNumber, STABILITY_EPSILON } from './stability';

export class LinearInterpolation implements Curve {
  private points: Point[];
  private segmentLengths: number[];
  private totalLength: number;
  private arcLengthSampler: ArcLengthSampler;
  private analyzer: CurveAnalyzer;

  constructor(points: Point[]) {
    const cleaned = cleanInputPoints(points);
    if (cleaned.points.length < 2) {
      throw new Error('LinearInterpolation requires at least 2 valid distinct points');
    }
    this.points = cleaned.points;
    this.segmentLengths = this.computeSegmentLengths();
    this.totalLength = this.segmentLengths.reduce((a, b) => a + b, 0);
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 1000);
    this.analyzer = new CurveAnalyzer(
      (t) => this.evaluate(t),
      (t) => this.derivative(t),
      200
    );
  }

  private computeSegmentLengths(): number[] {
    const lengths: number[] = [];
    for (let i = 0; i < this.points.length - 1; i++) {
      lengths.push(distance(this.points[i], this.points[i + 1]));
    }
    return lengths;
  }

  private findSegment(t: number): { index: number; localT: number } {
    t = clamp(t, 0, 1);
    const n = this.points.length - 1;
    if (t >= 1) {
      return { index: n - 1, localT: 1 };
    }
    const scaledT = t * n;
    const index = Math.min(Math.floor(scaledT), n - 1);
    const localT = scaledT - index;
    return { index, localT };
  }

  evaluate(t: number): Point {
    const { index, localT } = this.findSegment(t);
    return safePoint(lerpPoint(this.points[index], this.points[index + 1], localT));
  }

  derivative(t: number): Point {
    const { index } = this.findSegment(Math.max(0, Math.min(t, 1 - EPSILON)));
    const p0 = this.points[index];
    const p1 = this.points[index + 1];
    const n = this.points.length - 1;
    return safePoint(createPoint((p1.x - p0.x) * n, (p1.y - p0.y) * n));
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
    t0 = clamp(t0, 0, 1);
    t1 = clamp(t1, 0, 1);
    if (t0 > t1) [t0, t1] = [t1, t0];
    return safeNumber(this.arcLengthSampler.arcLength(t0, t1), 0);
  }

  getTotalLength(): number {
    return safeNumber(this.totalLength, 0);
  }

  sampleByArcLength(count: number): Point[] {
    return this.arcLengthSampler.sampleByArcLength(count);
  }

  getPoints(): Point[] {
    return [...this.points];
  }

  getSegmentCount(): number {
    return this.points.length - 1;
  }

  getAnalyzer(): CurveAnalyzer {
    return this.analyzer;
  }
}
