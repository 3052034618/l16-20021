import { Point, createPoint, distance, EPSILON } from './types';

export const STABILITY_EPSILON = 1e-10;
export const VERY_SMALL = 1e-15;
export const MAX_POINTS_WARNING = 10000;

export interface CleanedPoints {
  points: Point[];
  originalIndices: number[];
  removedDuplicates: boolean;
  allCollapsed: boolean;
}

export function isFinitePoint(p: Point): boolean {
  return Number.isFinite(p.x) && Number.isFinite(p.y);
}

export function safePoint(p: Point, fallback?: Point): Point {
  if (isFinitePoint(p)) return p;
  return fallback || createPoint(0, 0);
}

export function safeNumber(v: number, fallback: number = 0): number {
  if (Number.isFinite(v)) return v;
  return fallback;
}

export function pointsEqual(a: Point, b: Point, tolerance: number = STABILITY_EPSILON): boolean {
  return Math.abs(a.x - b.x) <= tolerance && Math.abs(a.y - b.y) <= tolerance;
}

export function cleanInputPoints(
  points: Point[],
  options: { minDistance?: number; warnOnLargeCount?: boolean } = {}
): CleanedPoints {
  const minDist = options.minDistance ?? STABILITY_EPSILON;

  if (!Array.isArray(points) || points.length === 0) {
    return {
      points: [],
      originalIndices: [],
      removedDuplicates: false,
      allCollapsed: true
    };
  }

  const validPoints: Point[] = [];
  const originalIndices: number[] = [];

  points.forEach((p, idx) => {
    if (isFinitePoint(p)) {
      validPoints.push({ x: p.x, y: p.y });
      originalIndices.push(idx);
    }
  });

  if (validPoints.length === 0) {
    return {
      points: [],
      originalIndices: [],
      removedDuplicates: false,
      allCollapsed: true
    };
  }

  if (points.length > MAX_POINTS_WARNING && options.warnOnLargeCount !== false) {
    console.warn(`[CurveEngine] Warning: received ${points.length} points, performance may be affected.`);
  }

  if (validPoints.length === 1) {
    return {
      points: validPoints,
      originalIndices: [originalIndices[0]],
      removedDuplicates: false,
      allCollapsed: false
    };
  }

  const deduped: Point[] = [validPoints[0]];
  const dedupedIndices: number[] = [originalIndices[0]];
  let removed = false;

  for (let i = 1; i < validPoints.length; i++) {
    const last = deduped[deduped.length - 1];
    const curr = validPoints[i];
    if (distance(last, curr) >= minDist) {
      deduped.push(curr);
      dedupedIndices.push(originalIndices[i]);
    } else {
      removed = true;
    }
  }

  if (deduped.length === 1) {
    return {
      points: deduped,
      originalIndices: dedupedIndices,
      removedDuplicates: removed,
      allCollapsed: true
    };
  }

  return {
    points: deduped,
    originalIndices: dedupedIndices,
    removedDuplicates: removed,
    allCollapsed: false
  };
}

export function safeNormalize(v: Point, fallback: Point = createPoint(1, 0)): Point {
  const lenSq = v.x * v.x + v.y * v.y;
  if (lenSq < VERY_SMALL) {
    return fallback;
  }
  const len = Math.sqrt(lenSq);
  return createPoint(v.x / len, v.y / len);
}

export function safeDivide(a: number, b: number, fallback: number = 0): number {
  if (Math.abs(b) < VERY_SMALL) {
    return fallback;
  }
  return a / b;
}

export function validateKnots(knots: number[], nPoints: number): { valid: boolean; knots: number[]; reason?: string } {
  if (knots.length !== nPoints) {
    return { valid: false, knots: [], reason: `knots length ${knots.length} != points length ${nPoints}` };
  }

  const cleanKnots = knots.map(k => safeNumber(k, NaN));

  for (const k of cleanKnots) {
    if (!Number.isFinite(k)) {
      return { valid: false, knots: [], reason: 'knot contains non-finite value' };
    }
  }

  for (let i = 1; i < cleanKnots.length; i++) {
    if (cleanKnots[i] < cleanKnots[i - 1]) {
      return { valid: false, knots: [], reason: 'knots not monotonically increasing' };
    }
  }

  let minDiff = Infinity;
  for (let i = 1; i < cleanKnots.length; i++) {
    const diff = cleanKnots[i] - cleanKnots[i - 1];
    if (diff < minDiff) minDiff = diff;
  }

  if (minDiff < VERY_SMALL) {
    return { valid: false, knots: [], reason: 'consecutive knots are too close (distance < 1e-15)' };
  }

  return { valid: true, knots: cleanKnots };
}

export function autoKnots(n: number, type: 'uniform' | 'chordal' | 'centripetal' = 'uniform', points?: Point[]): number[] {
  const knots: number[] = new Array(n);

  if (type === 'uniform' || !points) {
    for (let i = 0; i < n; i++) {
      knots[i] = i / (n - 1);
    }
    return knots;
  }

  if (n < 2) {
    return n === 1 ? [0] : [];
  }

  const dists: number[] = [];
  let totalDist = 0;

  for (let i = 1; i < n; i++) {
    let d: number;
    if (type === 'chordal') {
      d = Math.max(distance(points[i - 1], points[i]), VERY_SMALL);
    } else {
      d = Math.max(Math.sqrt(distance(points[i - 1], points[i])), VERY_SMALL);
    }
    dists.push(d);
    totalDist += d;
  }

  if (totalDist < VERY_SMALL) {
    for (let i = 0; i < n; i++) {
      knots[i] = i / (n - 1);
    }
    return knots;
  }

  knots[0] = 0;
  let accum = 0;
  for (let i = 1; i < n; i++) {
    accum += dists[i - 1];
    knots[i] = accum / totalDist;
  }

  return knots;
}

export interface BoundingBox {
  minX: number;
  maxX: number;
  minY: number;
  maxY: number;
  width: number;
  height: number;
  center: Point;
}

export function computeBoundingBox(points: Point[]): BoundingBox {
  if (points.length === 0) {
    return { minX: 0, maxX: 0, minY: 0, maxY: 0, width: 0, height: 0, center: createPoint(0, 0) };
  }

  let minX = Infinity, minY = Infinity;
  let maxX = -Infinity, maxY = -Infinity;

  for (const p of points) {
    if (!isFinitePoint(p)) continue;
    if (p.x < minX) minX = p.x;
    if (p.x > maxX) maxX = p.x;
    if (p.y < minY) minY = p.y;
    if (p.y > maxY) maxY = p.y;
  }

  if (!Number.isFinite(minX)) {
    minX = 0; maxX = 0; minY = 0; maxY = 0;
  }

  const width = maxX - minX;
  const height = maxY - minY;

  return {
    minX, maxX, minY, maxY,
    width, height,
    center: createPoint((minX + maxX) / 2, (minY + maxY) / 2)
  };
}

export function bboxContainsPoint(bbox: BoundingBox, p: Point, margin: number = 0): boolean {
  return p.x >= bbox.minX - margin && p.x <= bbox.maxX + margin &&
         p.y >= bbox.minY - margin && p.y <= bbox.maxY + margin;
}

export function bboxesIntersect(a: BoundingBox, b: BoundingBox, margin: number = 0): boolean {
  return a.minX - margin <= b.maxX + margin &&
         a.maxX + margin >= b.minX - margin &&
         a.minY - margin <= b.maxY + margin &&
         a.maxY + margin >= b.minY - margin;
}

export function clampT(t: number): number {
  if (t < 0) return 0;
  if (t > 1) return 1;
  return t;
}

export function lerpUnclamped(a: number, b: number, t: number): number {
  return a + (b - a) * t;
}
