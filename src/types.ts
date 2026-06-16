export interface Point {
  x: number;
  y: number;
}

export interface Point3D extends Point {
  z: number;
}

export type Vec2 = [number, number];
export type Vec3 = [number, number, number];
export type Vector = Vec2 | Vec3;

export interface Curve {
  evaluate(t: number): Point;
  sample(count: number): Point[];
  derivative(t: number): Point;
  arcLength(t0?: number, t1?: number): number;
  getTotalLength(): number;
  sampleByArcLength(count: number): Point[];
  getAnalyzer?(): any;
}

export interface CubicCoefficients {
  a: number;
  b: number;
  c: number;
  d: number;
}

export interface SplineSegment {
  xCoeff: CubicCoefficients;
  yCoeff: CubicCoefficients;
  tStart: number;
  tEnd: number;
}

export type EndCondition = 'natural' | 'clamped' | 'not-a-knot';

export interface SplineOptions {
  endCondition?: EndCondition;
  startTangent?: Point;
  endTangent?: Point;
}

export interface BezierSubdivision {
  left: Point[];
  right: Point[];
}

export const EPSILON = 1e-12;

export function createPoint(x: number, y: number): Point {
  return { x, y };
}

export function addPoints(a: Point, b: Point): Point {
  return { x: a.x + b.x, y: a.y + b.y };
}

export function subPoints(a: Point, b: Point): Point {
  return { x: a.x - b.x, y: a.y - b.y };
}

export function scalePoint(p: Point, s: number): Point {
  return { x: p.x * s, y: p.y * s };
}

export function dot(a: Point, b: Point): number {
  return a.x * b.x + a.y * b.y;
}

export function length(p: Point): number {
  return Math.sqrt(dot(p, p));
}

export function distance(a: Point, b: Point): number {
  return length(subPoints(a, b));
}

export function normalize(p: Point): Point {
  const len = length(p);
  if (len < EPSILON) return { x: 0, y: 0 };
  return scalePoint(p, 1 / len);
}

export function lerpNumber(a: number, b: number, t: number): number {
  return a + (b - a) * t;
}

export function lerpPoint(a: Point, b: Point, t: number): Point {
  return {
    x: lerpNumber(a.x, b.x, t),
    y: lerpNumber(a.y, b.y, t)
  };
}

export function clamp(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}
