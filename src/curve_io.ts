import { Point, Curve, createPoint, distance } from './types';
import { BezierCurve, CubicBezier, CompositeBezierCurve } from './bezier';
import { CatmullRomSpline } from './catmull_rom';
import { CubicSpline } from './cubic_spline';
import { LinearInterpolation } from './linear';
import { safePoint, cleanInputPoints, STABILITY_EPSILON } from './stability';

export interface SVGOptions {
  digits?: number;
  relative?: boolean;
  closed?: boolean;
}

export interface BezierFitOptions {
  maxError?: number;
  maxSegments?: number;
  reparameterize?: boolean;
  maxIterations?: number;
}

export interface FitResult {
  curves: CubicBezier[];
  error: number;
  iterations: number;
}

export class CurveIO {
  static toSVGPath(
    curve: Curve | CubicBezier[] | Point[],
    options: SVGOptions = {}
  ): string {
    const digits = options.digits ?? 3;

    if (Array.isArray(curve) && curve.length > 0) {
      if ('getControlPoints' in (curve[0] as any)) {
        return this.bezierSegmentsToSVG(curve as CubicBezier[], options);
      }
      return this.pointsToSVGPath(curve as Point[], options);
    }

    const c = curve as Curve;
    if (c instanceof CompositeBezierCurve) {
      const segments: CubicBezier[] = [];
      for (let i = 0; i < c.getCurveCount(); i++) {
        segments.push(c.getCurve(i) as CubicBezier);
      }
      return this.bezierSegmentsToSVG(segments, options);
    }

    if (c instanceof CatmullRomSpline) {
      const segments = c.toBezierSegments();
      return this.bezierSegmentsToSVG(segments, options);
    }

    const samples = c.sample(100);
    return this.pointsToSVGPath(samples, options);
  }

  static pointsToSVGPath(points: Point[], options: SVGOptions = {}): string {
    const digits = options.digits ?? 3;
    const closed = options.closed ?? false;

    if (points.length === 0) return '';

    const fmt = (n: number) => Number(n.toFixed(digits)).toString();

    let path = `M ${fmt(points[0].x)} ${fmt(points[0].y)}`;

    for (let i = 1; i < points.length; i++) {
      path += ` L ${fmt(points[i].x)} ${fmt(points[i].y)}`;
    }

    if (closed) path += ' Z';

    return path;
  }

  static bezierSegmentsToSVG(segments: CubicBezier[], options: SVGOptions = {}): string {
    const digits = options.digits ?? 3;
    const closed = options.closed ?? false;

    if (segments.length === 0) return '';

    const fmt = (n: number) => Number(n.toFixed(digits)).toString();
    const firstCP = segments[0].getControlPoints();
    let path = `M ${fmt(firstCP[0].x)} ${fmt(firstCP[0].y)}`;

    for (let i = 0; i < segments.length; i++) {
      const cps = segments[i].getControlPoints();
      path += ` C ${fmt(cps[1].x)} ${fmt(cps[1].y)}, ${fmt(cps[2].x)} ${fmt(cps[2].y)}, ${fmt(cps[3].x)} ${fmt(cps[3].y)}`;
    }

    if (closed) path += ' Z';

    return path;
  }

  static pointsToPolylineSVG(points: Point[], strokeColor: string = '#000', strokeWidth: number = 2): string {
    const coords = points.map(p => `${p.x.toFixed(2)},${p.y.toFixed(2)}`).join(' ');
    return `<polyline points="${coords}" fill="none" stroke="${strokeColor}" stroke-width="${strokeWidth}" stroke-linecap="round" stroke-linejoin="round"/>`;
  }

  static controlPointsToSvgMarkers(
    points: Point[],
    fill: string = '#f00',
    radius: number = 3
  ): string {
    return points.map(p =>
      `<circle cx="${p.x.toFixed(2)}" cy="${p.y.toFixed(2)}" r="${radius}" fill="${fill}"/>`
    ).join('\n');
  }

  static curveToFullSVG(
    curve: Curve,
    options: { width?: number; height?: number; padding?: number; stroke?: string; strokeWidth?: number; showControlPoints?: boolean } = {}
  ): string {
    const analyzer = (curve as any).getAnalyzer
      ? (curve as any).getAnalyzer()
      : null;

    let width = options.width ?? 500;
    let height = options.height ?? 300;
    const padding = options.padding ?? 20;
    const stroke = options.stroke ?? '#3366ff';
    const strokeWidth = options.strokeWidth ?? 2;

    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    const samples = curve.sample(200);
    for (const p of samples) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.y > maxY) maxY = p.y;
    }

    if (!Number.isFinite(minX)) { minX = 0; maxX = 100; minY = 0; maxY = 100; }

    const dataW = maxX - minX || 1;
    const dataH = maxY - minY || 1;
    const scaleX = (width - 2 * padding) / dataW;
    const scaleY = (height - 2 * padding) / dataH;
    const scale = Math.min(scaleX, scaleY);
    const offsetX = padding - minX * scale;
    const offsetY = padding + maxY * scale;

    const transform = (p: Point): Point => createPoint(p.x * scale + offsetX, -p.y * scale + offsetY);

    const path = this.toSVGPath(curve, { digits: 2 });

    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
    svg += `<g stroke="${stroke}" stroke-width="${strokeWidth}" fill="none">`;

    const transformed = path.replace(/([MLQC])\s+(-?[\d.]+)\s+(-?[\d.]+)/g,
      (_m, cmd, x, y) => {
        const tx = (parseFloat(x) * scale + offsetX).toFixed(2);
        const ty = (-parseFloat(y) * scale + offsetY).toFixed(2);
        return `${cmd} ${tx} ${ty}`;
      }).replace(/C\s+(-?[\d.]+)\s+(-?[\d.]+),\s*(-?[\d.]+)\s+(-?[\d.]+),\s*(-?[\d.]+)\s+(-?[\d.]+)/g,
        (_m, x1, y1, x2, y2, x, y) => {
          const px1 = (parseFloat(x1) * scale + offsetX).toFixed(2);
          const py1 = (-parseFloat(y1) * scale + offsetY).toFixed(2);
          const px2 = (parseFloat(x2) * scale + offsetX).toFixed(2);
          const py2 = (-parseFloat(y2) * scale + offsetY).toFixed(2);
          const px = (parseFloat(x) * scale + offsetX).toFixed(2);
          const py = (-parseFloat(y) * scale + offsetY).toFixed(2);
          return `C ${px1} ${py1}, ${px2} ${py2}, ${px} ${py}`;
        });

    svg += `<path d="${transformed}"/>`;
    svg += `</g>`;

    if (options.showControlPoints !== false) {
      if ('getControlPoints' in (curve as any)) {
        const pts = (curve as any).getControlPoints() as Point[];
        svg += `<g fill="#ff3333">`;
        for (const p of pts) {
          const tp = transform(p);
          svg += `<circle cx="${tp.x.toFixed(2)}" cy="${tp.y.toFixed(2)}" r="4"/>`;
        }
        svg += `</g>`;
      }
    }

    svg += `</svg>`;
    return svg;
  }

  static catmullRomToCompositeBezier(spline: CatmullRomSpline): CompositeBezierCurve {
    return new CompositeBezierCurve(spline.toBezierSegments(), spline.isClosed());
  }

  static fitCubicBezier(
    points: Point[],
    options: BezierFitOptions = {}
  ): FitResult {
    const cleaned = cleanInputPoints(points);
    if (cleaned.points.length < 2) {
      return { curves: [], error: 0, iterations: 0 };
    }
    if (cleaned.points.length === 2) {
      const [p0, p1] = cleaned.points;
      const mid1 = createPoint(p0.x + (p1.x - p0.x) / 3, p0.y + (p1.y - p0.y) / 3);
      const mid2 = createPoint(p0.x + 2 * (p1.x - p0.x) / 3, p0.y + 2 * (p1.y - p0.y) / 3);
      return { curves: [new CubicBezier(p0, mid1, mid2, p1)], error: 0, iterations: 1 };
    }

    const maxError = options.maxError ?? 1.0;
    const maxSegments = options.maxSegments ?? 32;
    const maxIter = options.maxIterations ?? 4;
    const doReparam = options.reparameterize ?? true;

    const n = cleaned.points.length;
    const params = this.computeChordLengthParams(cleaned.points);

    const result: CubicBezier[] = [];
    let totalError = 0;
    let totalIter = 0;

    this.fitRecursive(
      cleaned.points, 0, n - 1, params, 0, 1,
      result, maxError, maxIter, doReparam, maxSegments
    );

    for (const bezier of result) {
      let maxDist = 0;
      const samples = bezier.sample(50);
      for (let i = 0; i < cleaned.points.length; i++) {
        let best = Infinity;
        for (const s of samples) {
          const d = distance(s, cleaned.points[i]);
          if (d < best) best = d;
        }
        if (best > maxDist) maxDist = best;
      }
      totalError = Math.max(totalError, maxDist);
    }

    return { curves: result, error: totalError, iterations: totalIter };
  }

  private static computeChordLengthParams(points: Point[]): number[] {
    const n = points.length;
    const params: number[] = new Array(n);
    params[0] = 0;

    let total = 0;
    const dists: number[] = [0];

    for (let i = 1; i < n; i++) {
      const d = distance(points[i - 1], points[i]);
      total += d;
      dists.push(total);
    }

    if (total < STABILITY_EPSILON) {
      for (let i = 0; i < n; i++) params[i] = i / (n - 1);
      return params;
    }

    for (let i = 0; i < n; i++) {
      params[i] = dists[i] / total;
    }

    return params;
  }

  private static fitRecursive(
    points: Point[],
    first: number, last: number,
    params: number[], tStart: number, tEnd: number,
    result: CubicBezier[],
    maxError: number, maxIter: number, reparam: boolean,
    maxSegments: number
  ): void {
    const nPts = last - first + 1;

    if (nPts === 2) {
      const p0 = points[first];
      const p1 = points[last];
      const d = distance(p0, p1);
      const tDir = createPoint((p1.x - p0.x) / 3, (p1.y - p0.y) / 3);
      const c1 = createPoint(p0.x + tDir.x, p0.y + tDir.y);
      const c2 = createPoint(p1.x - tDir.x, p1.y - tDir.y);
      result.push(new CubicBezier(p0, c1, c2, p1));
      return;
    }

    if (result.length >= maxSegments) {
      const p0 = points[first];
      const p1 = points[last];
      const d = distance(p0, p1);
      const tDir = createPoint((p1.x - p0.x) / 3, (p1.y - p0.y) / 3);
      const c1 = createPoint(p0.x + tDir.x, p0.y + tDir.y);
      const c2 = createPoint(p1.x - tDir.x, p1.y - tDir.y);
      result.push(new CubicBezier(p0, c1, c2, p1));
      return;
    }

    const localParams: number[] = [];
    for (let i = first; i <= last; i++) {
      const t = (params[i] - tStart) / (tEnd - tStart || 1);
      localParams.push(Math.max(0, Math.min(1, t)));
    }

    const bezier = this.fitSingleSegment(points.slice(first, last + 1), localParams, maxIter, reparam);

    let maxDist = 0;
    let splitIdx = Math.floor((last - first) / 2);
    const sampleCount = Math.max(20, (last - first) * 2);

    for (let i = 0; i <= sampleCount; i++) {
      const tt = i / sampleCount;
      const bp = bezier.evaluate(tt);
      let best = Infinity;
      for (let j = first; j <= last; j++) {
        const d = distance(bp, points[j]);
        if (d < best) best = d;
      }
      if (best > maxDist) {
        maxDist = best;
      }
    }

    if (maxDist <= maxError) {
      result.push(bezier);
      return;
    }

    let maxErr = 0;
    for (let i = first + 1; i < last; i++) {
      let bestErr = Infinity;
      for (let s = 0; s <= 20; s++) {
        const tt = s / 20;
        const bp = bezier.evaluate(tt);
        const err = distance(bp, points[i]);
        if (err < bestErr) bestErr = err;
      }
      if (bestErr > maxErr) {
        maxErr = bestErr;
        splitIdx = i;
      }
    }

    if (splitIdx <= first || splitIdx >= last) {
      splitIdx = first + Math.floor((last - first) / 2);
    }

    const tMid = params[splitIdx];

    this.fitRecursive(
      points, first, splitIdx, params, tStart, tMid,
      result, maxError, maxIter, reparam, maxSegments
    );
    this.fitRecursive(
      points, splitIdx, last, params, tMid, tEnd,
      result, maxError, maxIter, reparam, maxSegments
    );
  }

  private static fitSingleSegment(
    points: Point[],
    params: number[],
    maxIter: number,
    reparam: boolean
  ): CubicBezier {
    const n = points.length;
    const P0 = points[0];
    const Pn = points[n - 1];

    let currentParams = [...params];

    for (let iter = 0; iter < Math.max(1, maxIter); iter++) {
      let A: [number, number][] = [];
      let cC1 = createPoint(0, 0);
      let cC2 = createPoint(0, 0);
      let X1 = createPoint(0, 0);
      let X2 = createPoint(0, 0);

      for (let i = 0; i < n; i++) {
        const t = currentParams[i];
        const t2 = t * t;
        const t3 = t2 * t;
        const b0 = (1 - t) * (1 - t) * (1 - t);
        const b1 = 3 * t * (1 - t) * (1 - t);
        const b2 = 3 * t2 * (1 - t);
        const b3 = t3;

        A.push([b1, b2]);

        const tmp = createPoint(
          points[i].x - b0 * P0.x - b3 * Pn.x,
          points[i].y - b0 * P0.y - b3 * Pn.y
        );

        X1 = createPoint(X1.x + b1 * tmp.x, X1.y + b1 * tmp.y);
        X2 = createPoint(X2.x + b2 * tmp.x, X2.y + b2 * tmp.y);
      }

      let C00 = 0, C01 = 0, C11 = 0;
      for (const [a, b] of A) {
        C00 += a * a;
        C01 += a * b;
        C11 += b * b;
      }

      const det = C00 * C11 - C01 * C01;
      let alphaL: number, alphaR: number;

      if (Math.abs(det) < STABILITY_EPSILON) {
        const dist = distance(P0, Pn);
        alphaL = dist / 3;
        alphaR = dist / 3;
      } else {
        alphaL = (C11 * (X1.x * 1) - C01 * (X2.x * 1)) / det;
        alphaR = (-C01 * (X1.x * 1) + C00 * (X2.x * 1)) / det;

        const denom = det;
        alphaL = Math.max(0, (C11 * X1.x - C01 * X2.x) / denom);
        alphaR = Math.max(0, (-C01 * X1.x + C00 * X2.x) / denom);

        const alphaL2 = (C11 * X1.y - C01 * X2.y) / denom;
        const alphaR2 = (-C01 * X1.y + C00 * X2.y) / denom;
        alphaL = (alphaL + alphaL2) / 2;
        alphaR = (alphaR + alphaR2) / 2;
      }

      if (!Number.isFinite(alphaL) || alphaL < STABILITY_EPSILON) alphaL = distance(P0, Pn) / 6;
      if (!Number.isFinite(alphaR) || alphaR < STABILITY_EPSILON) alphaR = distance(P0, Pn) / 6;

      let tanL: Point, tanR: Point;
      if (n >= 3) {
        tanL = createPoint(points[1].x - P0.x, points[1].y - P0.y);
        const lenL = Math.sqrt(tanL.x * tanL.x + tanL.y * tanL.y) || 1;
        tanL = createPoint(tanL.x / lenL, tanL.y / lenL);

        tanR = createPoint(Pn.x - points[n - 2].x, Pn.y - points[n - 2].y);
        const lenR = Math.sqrt(tanR.x * tanR.x + tanR.y * tanR.y) || 1;
        tanR = createPoint(tanR.x / lenR, tanR.y / lenR);
      } else {
        tanL = createPoint(1, 0);
        tanR = createPoint(1, 0);
      }

      const ctrl1 = createPoint(P0.x + alphaL * tanL.x, P0.y + alphaL * tanL.y);
      const ctrl2 = createPoint(Pn.x - alphaR * tanR.x, Pn.y - alphaR * tanR.y);

      if (!reparam || iter === maxIter - 1) {
        return new CubicBezier(P0, ctrl1, ctrl2, Pn);
      }

      const bezier = new CubicBezier(P0, ctrl1, ctrl2, Pn);

      for (let i = 1; i < n - 1; i++) {
        let bestT = currentParams[i];
        let bestD = Infinity;
        for (let s = 0; s <= 10; s++) {
          const tt = bestT + (s - 5) * 0.01;
          const t = Math.max(0, Math.min(1, tt));
          const bp = bezier.evaluate(t);
          const d = distance(bp, points[i]);
          if (d < bestD) {
            bestD = d;
            bestT = t;
          }
        }
        currentParams[i] = bestT;
      }
    }

    const bezier = new CubicBezier(P0, P0, Pn, Pn);
    return bezier;
  }
}
