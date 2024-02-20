using MultiPrecision;
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
        G32K65 = 32,
        G63K127 = 63,
        G64K129 = 64,
        G127K255 = 127,
        G128K257 = 128,
    }

    public static class GaussKronrodIntegral<N> where N : struct, IConstant {

        public static (MultiPrecision<N> value, MultiPrecision<N> error) Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, GaussKronrodOrder order = GaussKronrodOrder.G31K63) {
            ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)> ps = GaussKronrodPoints<N>.Table[order];

            MultiPrecision<N> sg = MultiPrecision<N>.Zero, sk = MultiPrecision<N>.Zero;
            MultiPrecision<N> r = b - a;

            if (!MultiPrecision<N>.IsFinite(r)) {
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

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int depth) {
            (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

            long eval_points = 1 + 2 * (int)order;
            if (!(error > eps) || depth == 0) {
                return (value, error, eval_points);
            }

            MultiPrecision<N> c = (a + b) / 2, eps_half = eps / 2;

            (MultiPrecision<N> value1, MultiPrecision<N> error1, long eval_points1) = AdaptiveIntegrateFiniteInterval(f, a, c, eps_half, order, depth > 0 ? depth - 1 : -1);
            (MultiPrecision<N> value2, MultiPrecision<N> error2, long eval_points2) = AdaptiveIntegrateFiniteInterval(f, c, b, eps_half, order, depth > 0 ? depth - 1 : -1);

            return (value1 + value2, error1 + error2, eval_points + eval_points1 + eval_points2);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrateInfiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int depth) {
            if (MultiPrecision<N>.IsNaN(a) || MultiPrecision<N>.IsNaN(b)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, order, depth);

                return (-value, error, eval_points);
            }

            if (MultiPrecision<N>.IsInfinity(a) && MultiPrecision<N>.IsInfinity(b)) {
                if (a.Sign == b.Sign) {
                    throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
                }

                MultiPrecision<N> g(MultiPrecision<N> t) {
                    if (MultiPrecision<N>.IsZero(t)) {
                        return MultiPrecision<N>.Zero;
                    }

                    MultiPrecision<N> u = (1 - t) / t;

                    return (f(u) + f(-u)) / (t * t);
                }

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
            }

            if (MultiPrecision<N>.IsFinite(a) && MultiPrecision<N>.IsInfinity(b)) {
                if (MultiPrecision<N>.IsZero(a)) {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (1 - t) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (1 - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
            }

            if (MultiPrecision<N>.IsInfinity(a) && MultiPrecision<N>.IsFinite(b)) {
                if (MultiPrecision<N>.IsZero(b)) {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (t - 1) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
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

        public static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order = GaussKronrodOrder.G31K63, int maxdepth = -1) {
            if (maxdepth < -1) {
                throw new ArgumentOutOfRangeException(nameof(maxdepth), "Invalid param. maxdepth=-1: infinite, maxdepth>=0: finite");
            }

            if (MultiPrecision<N>.IsFinite(a) && MultiPrecision<N>.IsFinite(b)) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, order, maxdepth);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, order, maxdepth);
            }
        }
    }
}
