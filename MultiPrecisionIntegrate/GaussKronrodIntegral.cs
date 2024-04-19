using MultiPrecision;
using System.Collections.ObjectModel;
using System.Diagnostics;

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
        G255K511 = 255,
        G256K513 = 256,
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

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) UnlimitedIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order) {
            Stack<(MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps)> stack = new();
            stack.Push((a, b, eps));

            long eval_points_sum = 0;
            MultiPrecision<N> value_sum = 0d, error_sum = 0d;

            while (stack.Count > 0) {
                (a, b, eps) = stack.Pop();

                (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps)) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                MultiPrecision<N> c = MultiPrecision<N>.Ldexp(a + b, -1), eps_half = MultiPrecision<N>.Ldexp(eps, -1);
                stack.Push((a, c, eps_half));
                stack.Push((c, b, eps_half));
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) LimitedDepthIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int maxdepth) {
            Debug.Assert(maxdepth >= 0);

            Stack<(MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, int depth)> stack = new();
            stack.Push((a, b, eps, maxdepth));

            long eval_points_sum = 0;
            MultiPrecision<N> value_sum = 0d, error_sum = 0d;

            while (stack.Count > 0) {
                (a, b, eps, int depth) = stack.Pop();

                (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || depth <= 0) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                MultiPrecision<N> c = MultiPrecision<N>.Ldexp(a + b, -1), eps_half = MultiPrecision<N>.Ldexp(eps, -1);
                depth -= 1;
                stack.Push((a, c, eps_half, depth));
                stack.Push((c, b, eps_half, depth));
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) LimitedEvalIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, long discontinue_eval_points) {
            Debug.Assert(discontinue_eval_points >= 0);

            PriorityQueue<(MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps), long> queue = new();
            queue.Enqueue((a, b, eps), 0);

            long eval_points_sum = 0;
            MultiPrecision<N> value_sum = 0d, error_sum = 0d;

            while (queue.Count > 0) {
                (a, b, eps) = queue.Dequeue();

                (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || eval_points_sum > discontinue_eval_points) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                MultiPrecision<N> c = MultiPrecision<N>.Ldexp(a + b, -1), eps_half = MultiPrecision<N>.Ldexp(eps, -1);
                long priority = double.ILogB((double)error);
                queue.Enqueue((a, c, eps_half), -priority);
                queue.Enqueue((c, b, eps_half), -priority);
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) LimitedDepthAndEvalIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int maxdepth, long discontinue_eval_points) {
            Debug.Assert(maxdepth >= 0);
            Debug.Assert(discontinue_eval_points >= 0);

            PriorityQueue<(MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, int depth), long> queue = new();
            queue.Enqueue((a, b, eps, maxdepth), 0);

            long eval_points_sum = 0;
            MultiPrecision<N> value_sum = 0d, error_sum = 0d;

            while (queue.Count > 0) {
                (a, b, eps, int depth) = queue.Dequeue();

                (MultiPrecision<N> value, MultiPrecision<N> error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || depth <= 0 || eval_points_sum > discontinue_eval_points) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                MultiPrecision<N> c = MultiPrecision<N>.Ldexp(a + b, -1), eps_half = MultiPrecision<N>.Ldexp(eps, -1);
                long priority = double.ILogB((double)error);
                depth -= 1;
                queue.Enqueue((a, c, eps_half, depth), -priority);
                queue.Enqueue((c, b, eps_half, depth), -priority);
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrateFiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int maxdepth, long discontinue_eval_points) {
            if (maxdepth >= 0 && discontinue_eval_points >= 0) {
                return LimitedDepthAndEvalIntegrateFiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
            if (maxdepth >= 0) {
                return LimitedDepthIntegrateFiniteInterval(f, a, b, eps, order, maxdepth);
            }
            if (discontinue_eval_points >= 0) {
                return LimitedEvalIntegrateFiniteInterval(f, a, b, eps, order, discontinue_eval_points);
            }

            return UnlimitedIntegrateFiniteInterval(f, a, b, eps, order);
        }

        private static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrateInfiniteInterval(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order, int depth, long discontinue_eval_points) {
            if (MultiPrecision<N>.IsNaN(a) || MultiPrecision<N>.IsNaN(b)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, order, depth, discontinue_eval_points);

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

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
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

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (1 - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
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

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
                else {
                    MultiPrecision<N> g(MultiPrecision<N> t) {
                        if (MultiPrecision<N>.IsZero(t)) {
                            return MultiPrecision<N>.Zero;
                        }

                        MultiPrecision<N> u = (t - 1) / t;

                        return f(b + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
            }

            throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
        }

        public static (MultiPrecision<N> value, MultiPrecision<N> error, long eval_points) AdaptiveIntegrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, MultiPrecision<N> eps, GaussKronrodOrder order = GaussKronrodOrder.G31K63, int maxdepth = -1, long discontinue_eval_points = -1) {
            if (maxdepth < -1) {
                throw new ArgumentOutOfRangeException(nameof(maxdepth), "Invalid param. maxdepth=-1: infinite, maxdepth>=0: finite");
            }
            if (discontinue_eval_points < -1) {
                throw new ArgumentOutOfRangeException(nameof(discontinue_eval_points), "Invalid param. discontinue_eval_points=-1: infinite, discontinue_eval_points>=0: finite");
            }
            if (!(MultiPrecision<N>.IsFinite(eps) && MultiPrecision<N>.IsPositive(eps))) {
                throw new ArgumentOutOfRangeException(nameof(eps), "Invalid param. eps must be nonnegative value");
            }

            if (MultiPrecision<N>.IsFinite(a) && MultiPrecision<N>.IsFinite(b)) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
        }
    }
}
