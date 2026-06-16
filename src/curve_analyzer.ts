import {
  Point, Curve, createPoint, distance, subPoints, dot, scalePoint, addPoints
} from './types';
import {
  BoundingBox, computeBoundingBox, bboxContainsPoint,
  clampT, safeNumber, safePoint, STABILITY_EPSILON, safeDivide, safeNormalize
} from './stability';

export interface CurveAnalysis {
  boundingBox: BoundingBox;
  arcLength: number;
  startPoint: Point;
  endPoint: Point;
  estimatedMinSpeed: number;
  estimatedMaxSpeed: number;
  estimatedAvgSpeed: number;
  estimatedMinCurvature: number;
  estimatedMaxCurvature: number;
}

export interface ProjectionResult {
  point: Point;
  parameter: number;
  distance: number;
  error: number;
}

export interface SpeedAnalysis {
  t: number;
  position: Point;
  speed: number;
  velocity: Point;
  acceleration?: Point;
  curvature?: number;
  radiusOfCurvature?: number;
}

export class CurveAnalyzer {
  private evaluateFn: (t: number) => Point;
  private derivativeFn?: (t: number) => Point;
  private sampleResolution: number;

  constructor(
    evaluateFn: (t: number) => Point,
    derivativeFn?: (t: number) => Point,
    sampleResolution: number = 200
  ) {
    this.evaluateFn = evaluateFn;
    this.derivativeFn = derivativeFn;
    this.sampleResolution = sampleResolution;
  }

  numericalDerivative(t: number, h: number = 1e-6): Point {
    if (this.derivativeFn) {
      return safePoint(this.derivativeFn(t));
    }
    const t1 = Math.max(0, t - h);
    const t2 = Math.min(1, t + h);
    const p1 = safePoint(this.evaluateFn(t1));
    const p2 = safePoint(this.evaluateFn(t2));
    const dt = t2 - t1;
    if (dt < STABILITY_EPSILON) {
      return createPoint(0, 0);
    }
    return createPoint(
      (p2.x - p1.x) / dt,
      (p2.y - p1.y) / dt
    );
  }

  numericalSecondDerivative(t: number, h: number = 1e-5): Point {
    const t1 = Math.max(0, t - h);
    const t2 = Math.min(1, t + h);
    const tm = t;
    const p1 = safePoint(this.evaluateFn(t1));
    const pm = safePoint(this.evaluateFn(tm));
    const p2 = safePoint(this.evaluateFn(t2));
    const dt2 = h * h;
    if (dt2 < STABILITY_EPSILON * STABILITY_EPSILON) {
      return createPoint(0, 0);
    }
    return createPoint(
      (p2.x - 2 * pm.x + p1.x) / dt2,
      (p2.y - 2 * pm.y + p1.y) / dt2
    );
  }

  speed(t: number): number {
    const v = this.numericalDerivative(t);
    return Math.sqrt(v.x * v.x + v.y * v.y);
  }

  velocity(t: number): Point {
    return this.numericalDerivative(t);
  }

  acceleration(t: number): Point {
    return this.numericalSecondDerivative(t);
  }

  curvature(t: number): number {
    const v = this.numericalDerivative(t);
    const a = this.numericalSecondDerivative(t);
    const cross = v.x * a.y - v.y * a.x;
    const speedCubed = Math.pow(v.x * v.x + v.y * v.y, 1.5);
    return safeDivide(cross, speedCubed, 0);
  }

  radiusOfCurvature(t: number): number {
    const k = this.curvature(t);
    if (Math.abs(k) < STABILITY_EPSILON) {
      return Infinity;
    }
    return 1 / Math.abs(k);
  }

  tangent(t: number): Point {
    return safeNormalize(this.numericalDerivative(t));
  }

  normal(t: number): Point {
    const tan = this.tangent(t);
    return createPoint(-tan.y, tan.x);
  }

  getBoundingBox(): BoundingBox {
    const samples: Point[] = [];
    for (let i = 0; i <= this.sampleResolution; i++) {
      const t = i / this.sampleResolution;
      samples.push(safePoint(this.evaluateFn(t)));
    }
    return computeBoundingBox(samples);
  }

  findNearestPoint(
    target: Point,
    options: { tolerance?: number; maxIterations?: number; initialSamples?: number } = {}
  ): ProjectionResult {
    const tolerance = options.tolerance ?? 1e-10;
    const maxIterations = options.maxIterations ?? 50;
    const initialSamples = options.initialSamples ?? 100;

    let bestT = 0;
    let bestDist = Infinity;
    let bestPoint: Point = safePoint(this.evaluateFn(0));

    for (let i = 0; i <= initialSamples; i++) {
      const t = i / initialSamples;
      const p = safePoint(this.evaluateFn(t));
      const d = distance(p, target);
      if (d < bestDist) {
        bestDist = d;
        bestT = t;
        bestPoint = p;
      }
    }

    let t = bestT;
    for (let iter = 0; iter < maxIterations; iter++) {
      const dt = 1e-6;
      const tm = clampT(t - dt);
      const tp = clampT(t + dt);
      const pm = safePoint(this.evaluateFn(tm));
      const p0 = safePoint(this.evaluateFn(t));
      const pp = safePoint(this.evaluateFn(tp));

      const dm = distance(pm, target);
      const d0 = distance(p0, target);
      const dp = distance(pp, target);

      if (d0 <= dm && d0 <= dp && Math.min(t - 0, 1 - t) < STABILITY_EPSILON) {
        break;
      }

      const dt_actual = tp - tm;
      if (dt_actual < STABILITY_EPSILON) break;

      const dDeriv = (dp - dm) / (2 * dt_actual);
      const dSecond = (dp - 2 * d0 + dm) / (dt_actual * dt_actual);

      let delta: number;
      if (Math.abs(dSecond) > STABILITY_EPSILON && dSecond > 0) {
        delta = -dDeriv / dSecond;
      } else {
        delta = -dDeriv * 0.01;
      }

      const newT = clampT(t + delta);
      if (Math.abs(newT - t) < tolerance) {
        t = newT;
        break;
      }
      t = newT;
    }

    const finalPoint = safePoint(this.evaluateFn(t));
    const finalDist = distance(finalPoint, target);

    return {
      point: finalPoint,
      parameter: t,
      distance: finalDist,
      error: Math.abs(finalDist - Math.min(finalDist, bestDist))
    };
  }

  projectPoint(target: Point, tolerance?: number): ProjectionResult {
    return this.findNearestPoint(target, { tolerance });
  }

  speedProfile(sampleCount: number = 100): SpeedAnalysis[] {
    const result: SpeedAnalysis[] = [];
    for (let i = 0; i < sampleCount; i++) {
      const t = i / (sampleCount - 1);
      const pos = safePoint(this.evaluateFn(t));
      const vel = this.numericalDerivative(t);
      const spd = Math.sqrt(vel.x * vel.x + vel.y * vel.y);
      const acc = this.numericalSecondDerivative(t);
      const k = this.curvature(t);
      result.push({
        t,
        position: pos,
        speed: spd,
        velocity: vel,
        acceleration: acc,
        curvature: k,
        radiusOfCurvature: Math.abs(k) < STABILITY_EPSILON ? Infinity : 1 / Math.abs(k)
      });
    }
    return result;
  }

  analyze(): CurveAnalysis {
    const bbox = this.getBoundingBox();
    const profile = this.speedProfile(Math.max(50, this.sampleResolution));

    let minSpeed = Infinity, maxSpeed = -Infinity, totalSpeed = 0;
    let minCurv = Infinity, maxCurv = -Infinity;

    for (const p of profile) {
      if (p.speed < minSpeed) minSpeed = p.speed;
      if (p.speed > maxSpeed) maxSpeed = p.speed;
      totalSpeed += p.speed;
      if (p.curvature !== undefined) {
        if (p.curvature < minCurv) minCurv = p.curvature;
        if (p.curvature > maxCurv) maxCurv = p.curvature;
      }
    }

    const startPoint = safePoint(this.evaluateFn(0));
    const endPoint = safePoint(this.evaluateFn(1));

    let arcLength = 0;
    for (let i = 1; i < profile.length; i++) {
      arcLength += distance(profile[i - 1].position, profile[i].position);
    }

    return {
      boundingBox: bbox,
      arcLength: safeNumber(arcLength),
      startPoint,
      endPoint,
      estimatedMinSpeed: safeNumber(minSpeed === Infinity ? 0 : minSpeed),
      estimatedMaxSpeed: safeNumber(maxSpeed === -Infinity ? 0 : maxSpeed),
      estimatedAvgSpeed: safeNumber(totalSpeed / profile.length),
      estimatedMinCurvature: safeNumber(minCurv === Infinity ? 0 : minCurv),
      estimatedMaxCurvature: safeNumber(maxCurv === -Infinity ? 0 : maxCurv)
    };
  }

  detectCorners(angleThresholdDegrees: number = 30, sampleCount: number = 200): number[] {
    const threshold = Math.cos(angleThresholdDegrees * Math.PI / 180);
    const corners: number[] = [];
    const tangents: Point[] = [];

    for (let i = 0; i < sampleCount; i++) {
      const t = i / (sampleCount - 1);
      tangents.push(this.tangent(t));
    }

    for (let i = 1; i < sampleCount - 1; i++) {
      const t1 = tangents[i - 1];
      const t2 = tangents[i + 1];
      const d = dot(t1, t2);
      if (d < threshold) {
        corners.push(i / (sampleCount - 1));
      }
    }

    if (corners.length > 1) {
      const merged: number[] = [corners[0]];
      for (let i = 1; i < corners.length; i++) {
        if (corners[i] - merged[merged.length - 1] > 5 / sampleCount) {
          merged.push(corners[i]);
        }
      }
      return merged;
    }

    return corners;
  }

  static pointToSegmentDistance(p: Point, a: Point, b: Point): { distance: number; parameter: number; point: Point } {
    const ab = subPoints(b, a);
    const abLenSq = dot(ab, ab);

    if (abLenSq < STABILITY_EPSILON) {
      return { distance: distance(p, a), parameter: 0, point: a };
    }

    const ap = subPoints(p, a);
    let t = dot(ap, ab) / abLenSq;
    t = clampT(t);

    const proj = addPoints(a, scalePoint(ab, t));
    return {
      distance: distance(p, proj),
      parameter: t,
      point: proj
    };
  }

  isPointNearCurve(
    point: Point,
    maxDistance: number
  ): { near: boolean; projection?: ProjectionResult } {
    const bbox = this.getBoundingBox();
    const margin = maxDistance;
    if (!bboxContainsPoint(bbox, point, margin)) {
      return { near: false };
    }
    const proj = this.findNearestPoint(point);
    return { near: proj.distance <= maxDistance, projection: proj };
  }
}
