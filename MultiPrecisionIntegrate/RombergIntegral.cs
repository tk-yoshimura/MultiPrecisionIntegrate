using MultiPrecision;

namespace MultiPrecisionIntegrate {

    public static class RombergIntegral<N> where N : struct, IConstant {

        public static MultiPrecision<N> Integrate(Func<MultiPrecision<N>, MultiPrecision<N>> f, MultiPrecision<N> a, MultiPrecision<N> b, int precision_level = 16) {
            MultiPrecision<N> h = b - a;

            if (!h.IsFinite) {
                throw new ArgumentOutOfRangeException($"{nameof(a)},{nameof(b)}");
            }

            if (precision_level < 1 || precision_level > 24) {
                throw new ArgumentOutOfRangeException(nameof(precision_level));
            }

            int max_div = 1 << precision_level;
            MultiPrecision<N> min_h = h / max_div;
            MultiPrecision<N>[] v = new MultiPrecision<N>[max_div + 1];
            RichardsonExtrapolation<N> conv = new();

            for (int i = 0; i <= max_div; i++) {
                v[i] = f(a + i * min_h);
            }

            MultiPrecision<N> t = h * (v[0] + v[max_div]) / 2, new_t;
            conv.Inject(t);

            for (int s = max_div; s > 1; s /= 2) {
                new_t = 0;
                for (int i = s / 2; i < max_div; i += s) {
                    new_t += v[i];
                }

                h /= 2;
                t = t / 2 + h * new_t;

                conv.Inject(t);
            }

            return conv.ConvergenceValue;
        }
    }
}
