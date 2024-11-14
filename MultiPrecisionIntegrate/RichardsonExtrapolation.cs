using MultiPrecision;

namespace MultiPrecisionIntegrate {

    internal class RichardsonExtrapolation<N> where N : struct, IConstant {
        static readonly List<MultiPrecision<N>> rs = new() { MultiPrecision<N>.NaN };
        readonly List<MultiPrecision<N>[]> values = new();

        public void Inject(MultiPrecision<N> new_value) {
            lock (values) {
                if (SeriesCount <= 0) {
                    values.Add([new_value]);
                    return;
                }

                MultiPrecision<N>[] t = values[SeriesCount - 1], t_next = new MultiPrecision<N>[SeriesCount + 1];

                t_next[0] = new_value;

                for (int i = 1; i <= SeriesCount; i++) {
                    t_next[i] = t_next[i - 1] + (t_next[i - 1] - t[i - 1]) * R(i);
                }

                values.Add(t_next);
            }
        }

        private static MultiPrecision<N> R(int i) {
            if (i < rs.Count) { 
                return rs[i];
            }

            lock (rs) {
                for (int k = rs.Count; k <= i; k++) {
                    MultiPrecision<N> r = 1d / (MultiPrecision<N>.Ldexp(1d, k * 2) - 1);

                    rs.Add(r);
                }

                return rs[i];
            }
        }

        public IEnumerable<MultiPrecision<N>> Series {
            get {
                for (int i = 0; i < values.Count; i++) {
                    yield return values[i][i];
                }
            }
        }

        public MultiPrecision<N> ConvergenceValue {
            get {
                if (SeriesCount <= 0) {
                    throw new InvalidOperationException();
                }

                return values[SeriesCount - 1][SeriesCount - 1];
            }
        }

        public int SeriesCount => values.Count;

        public MultiPrecision<N> Epsilon {
            get {
                if (SeriesCount <= 1) {
                    throw new InvalidOperationException();
                }

                return MultiPrecision<N>.Abs(values[SeriesCount - 1][SeriesCount - 1] - values[SeriesCount - 2][SeriesCount - 2]);
            }
        }
    }
}
