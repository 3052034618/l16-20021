import { Point, Curve, lerpPoint, distance, createPoint, clamp, EPSILON } from './types';
import { ArcLengthSampler } from './arc_length';

export class LinearInterpolation implements Curve {
  private points: Point[];
  private segmentLengths: number[];
  private totalLength: number;
  private arcLengthSampler: ArcLengthSampler;

  constructor(points: Point[]) {
    if (points.length < 2) {
      throw new Error('LinearInterpolation requires at least 2 points');
    }
    this.points = [...points];
    this.segmentLengths = this.computeSegmentLengths();
    this.totalLength = this.segmentLengths.reduce((a, b) => a + b, 0);
    this.arcLengthSampler = new ArcLengthSampler((t) => this.evaluate(t), 1000);
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
    const scaledT = t * n;
    const index = Math.min(Math.floor(scaledT), n - 1);
    const localT = scaledT - index;
    return { index, localT };
  }

  evaluate(t: number): Point {
    const { index, localT } = this.findSegment(t);
    return lerpPoint(this.points[index], this.points[index + 1], localT);
  }

  derivative(t: number): Point {
    const { index } = this.findSegment(Math.max(0, Math.min(t, 1 - EPSILON)));
    const p0 = this.points[index];
    const p1 = this.points[index + 1];
    const dx = (p1.x - p0.x) * (this.points.length - 1);
    const dy = (p1.y - p0.y) * (this.points.length - 1);
    return createPoint(dx, dy);
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
    return this.arcLengthSampler.arcLength(t0, t1);
  }

  getTotalLength(): number {
    return this.totalLength;
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
}
