using MultiPrecision;
using System.Collections.ObjectModel;

namespace MultiPrecisionIntegrate {
    public static class GaussLaguerreIntegral<N> where N : struct, IConstant {
        public static MultiPrecision<N> Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, int n, bool f_expscaled = false) {
            if (n < GaussLegendrePoints<N>.MinPoints || n > GaussLegendrePoints<N>.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)> ps = GaussLaguerrePoints<N>.Table[n];
            MultiPrecision<N> s = MultiPrecision<N>.Zero;

            if (f_expscaled) {
                foreach ((MultiPrecision<N> x, MultiPrecision<N> w, _) in ps) {
                    s += w * f(x);
                }
            }
            else {
                foreach ((MultiPrecision<N> x, _, MultiPrecision<N> w) in ps) {
                    s += w * f(x);
                }
            }

            return s;
        }
    }
}
