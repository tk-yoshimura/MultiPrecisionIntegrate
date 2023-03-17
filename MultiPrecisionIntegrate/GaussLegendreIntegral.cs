using MultiPrecision;
using System.Collections.ObjectModel;

namespace MultiPrecisionIntegrate {
    public static class GaussLegendreIntegral<N> where N : struct, IConstant {
        public static MultiPrecision<N> Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, int n) {
            MultiPrecision<N> r = b - a;

            if (!r.IsFinite) {
                throw new ArgumentOutOfRangeException($"{nameof(a)},{nameof(b)}");
            }

            if (n < GaussLegendrePoints<N>.MinPoints || n > GaussLegendrePoints<N>.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w)> ps = GaussLegendrePoints<N>.Table[n];

            MultiPrecision<N> s = MultiPrecision<N>.Zero;

            foreach ((MultiPrecision<N> x, MultiPrecision<N> w) in ps) {
                MultiPrecision<N> x_shifted = x * r + a;

                s += w * f(x_shifted);
            }

            s *= r;

            return s;
        }
    }
}
