using MultiPrecision;
using System.Collections.ObjectModel;

namespace MultiPrecisionIntegrate {
    public static class GaussHermiteIntegral<N> where N : struct, IConstant {
        public static MultiPrecision<N> Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, int n, bool f_expscaled = false) {
            if (n < GaussHermitePoints<N>.MinPoints || n > GaussHermitePoints<N>.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)> ps = GaussHermitePoints<N>.Table[n];

            MultiPrecision<N> s = MultiPrecision<N>.Zero;

            if ((n & 1) == 0) {
                if (f_expscaled) {
                    foreach ((MultiPrecision<N> x, MultiPrecision<N> w, _) in ps) {
                        s += w * (f(x) + f(-x));
                    }
                }
                else {
                    foreach ((MultiPrecision<N> x, _, MultiPrecision<N> w) in ps) {
                        s += w * (f(x) + f(-x));
                    }
                }
            }
            else {
                if (f_expscaled) {
                    s += ps[0].w * f(0);

                    foreach ((MultiPrecision<N> x, MultiPrecision<N> w, _) in ps.Skip(1)) {
                        s += w * (f(x) + f(-x));
                    }
                }
                else {
                    s += ps[0].wexp * f(0);

                    foreach ((MultiPrecision<N> x, _, MultiPrecision<N> w) in ps.Skip(1)) {
                        s += w * (f(x) + f(-x));
                    }
                }
            }

            return s;
        }
    }
}
