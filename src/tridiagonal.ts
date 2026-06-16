export interface TridiagonalSystem {
  a: number[];
  b: number[];
  c: number[];
  d: number[];
}

export function solveTridiagonal(system: TridiagonalSystem): number[] {
  const { a, b, c, d } = system;
  const n = d.length;

  if (a.length !== n || b.length !== n || c.length !== n) {
    throw new Error('Tridiagonal matrix arrays must have the same length');
  }

  const cPrime: number[] = new Array(n).fill(0);
  const dPrime: number[] = new Array(n).fill(0);
  const x: number[] = new Array(n).fill(0);

  if (Math.abs(b[0]) < 1e-15) {
    throw new Error('First diagonal element is zero, cannot solve');
  }

  cPrime[0] = c[0] / b[0];
  dPrime[0] = d[0] / b[0];

  for (let i = 1; i < n; i++) {
    const denom = b[i] - a[i] * cPrime[i - 1];
    if (Math.abs(denom) < 1e-15) {
      throw new Error(`Zero pivot at row ${i}, system may be singular`);
    }
    cPrime[i] = c[i] / denom;
    dPrime[i] = (d[i] - a[i] * dPrime[i - 1]) / denom;
  }

  x[n - 1] = dPrime[n - 1];
  for (let i = n - 2; i >= 0; i--) {
    x[i] = dPrime[i] - cPrime[i] * x[i + 1];
  }

  return x;
}

export function verifyTridiagonalSolution(system: TridiagonalSystem, x: number[]): number[] {
  const { a, b, c, d } = system;
  const n = d.length;
  const residual: number[] = new Array(n).fill(0);

  for (let i = 0; i < n; i++) {
    let sum = b[i] * x[i];
    if (i > 0) sum += a[i] * x[i - 1];
    if (i < n - 1) sum += c[i] * x[i + 1];
    residual[i] = sum - d[i];
  }

  return residual;
}

export function createTridiagonalSystem(n: number): TridiagonalSystem {
  return {
    a: new Array(n).fill(0),
    b: new Array(n).fill(0),
    c: new Array(n).fill(0),
    d: new Array(n).fill(0)
  };
}
