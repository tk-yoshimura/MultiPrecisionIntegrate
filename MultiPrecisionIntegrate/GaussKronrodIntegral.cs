﻿using MultiPrecision;
using System.Collections.ObjectModel;

namespace MultiPrecisionIntegrate {
    public enum GaussKronrodOrder : int {
        G3K7 = 3,
        G4K9 = 4,
        G7K15 = 7,
        G8K17 = 8,
        G15K31 = 15,
        G16K33 = 16,
        G31K63 = 31,
        G32K65 = 32
    }

    public static class GaussKronrodIntegral<N> where N : struct, IConstant {

        public static (MultiPrecision<N> value, MultiPrecision<N> error) Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, GaussKronrodOrder order = GaussKronrodOrder.G7K15) {
            ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)> ps = GaussKronrodPoints<N>.Table[order];

            MultiPrecision<N> sg = MultiPrecision<N>.Zero, sk = MultiPrecision<N>.Zero;
            MultiPrecision<N> r = b - a;

            if (!r.IsFinite) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            for (int i = 0; i < ps.Count; i++) {
                MultiPrecision<N> x = ps[i].x;
                MultiPrecision<N> x_shifted = x * r + a;

                MultiPrecision<N> y = f(x_shifted);

                sk += ps[i].wk * y;

                if ((i & 1) == 1) {
                    sg += ps[i].wg * y;
                }
            }

            sk *= r;
            sg *= r;

            MultiPrecision<N> error = MultiPrecision<N>.Abs(sk - sg);

            return (sk, error);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error) AdaptiveIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

            if (!(error > eps) || depth <= 0) {
                return (value, error);
            }

            MultiPrecision<N> c = (a + b) / 2, eps_half = eps / 2;

            (MultiPrecision<N> value1, MultiPrecision<N> error1) = AdaptiveIntegrateFiniteInterval(f, a, c, eps_half, order, depth - 1);
            (MultiPrecision<N> value2, MultiPrecision<N> error2) = AdaptiveIntegrateFiniteInterval(f, c, b, eps_half, order, depth - 1);

            return (value1 + value2, error1 + error2);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error) AdaptiveIntegrateInfiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            if (a.IsNaN || b.IsNaN) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (MultiPrecision<N> value, MultiPrecision<N> error) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, order, depth);

                return (-value, error);
            }

            if (a.IsInfinity && b.IsInfinity) {
                if (a.Sign == b.Sign) {
                    throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
                }

                MultiPrecision<N> g(MultiPrecision<N> t) {
                    if (t.IsZero) {
                        return MultiPrecision<N>.Zero;
                    }

                    MultiPrecision<N> u = (1 - t) / t;

                    return (f(u) + f(-u)) / (t * t);
                }

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
            }

            if (a.IsFinite && b.IsInfinity) {
                if (a.IsZero) {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (t.IsZero) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (1 - t) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (t.IsZero) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (1 - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
            }

            if (a.IsInfinity && b.IsFinite) {
                if (b.IsZero) {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (t.IsZero) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (t - 1) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (t.IsZero) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (t - 1) / t;

                        return f(b + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
            }

            throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
        }

        public static (MultiPrecision<N> value, MultiPrecision<N> error) AdaptiveIntegrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            if (a.IsFinite && b.IsFinite) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, order, depth);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, order, depth);
            }
        }
    }
}
