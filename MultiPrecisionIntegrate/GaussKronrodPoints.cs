using MultiPrecision;
using System.Collections.ObjectModel;
using System.IO.Compression;

namespace MultiPrecisionIntegrate {

    public static class GaussKronrodPoints {
        public const int MaxLength = 64;
    }

    public static class GaussKronrodPoints<N> where N : struct, IConstant {
        public static ReadOnlyDictionary<GaussKronrodOrder, ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)>> Table { private set; get; }

        static GaussKronrodPoints() {
            if (MultiPrecision<N>.Length > GaussKronrodPoints.MaxLength) {
                throw new NotSupportedException($"Too large multi precision. (N <= {GaussKronrodPoints.MaxLength})");
            }

            Table = new(new Dictionary<GaussKronrodOrder, ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)>>(){
                { GaussKronrodOrder.G3K7, LoadTable(Resource.G3K7_n64, n: 3) },
                { GaussKronrodOrder.G4K9, LoadTable(Resource.G4K9_n64, n: 4) },
                { GaussKronrodOrder.G7K15, LoadTable(Resource.G7K15_n64, n: 7) },
                { GaussKronrodOrder.G8K17, LoadTable(Resource.G8K17_n64, n: 8) },
                { GaussKronrodOrder.G15K31, LoadTable(Resource.G15K31_n64, n: 15) },
                { GaussKronrodOrder.G16K33, LoadTable(Resource.G16K33_n64, n: 16) },
                { GaussKronrodOrder.G31K63, LoadTable(Resource.G31K63_n64, n: 31) },
                { GaussKronrodOrder.G32K65, LoadTable(Resource.G32K65_n64, n: 32) },
            });
        }

        static ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)> LoadTable(byte[] resource_bytes, int n) {
            using MemoryStream stream = new();
            using (MemoryStream resorce = new(resource_bytes)) {
                using GZipStream decompressor = new(resorce, CompressionMode.Decompress);

                decompressor.CopyTo(stream);
            }

            stream.Position = 0;

            (MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)[] ps = new (MultiPrecision<N> x, MultiPrecision<N> wg, MultiPrecision<N> wk)[checked(2 * n + 1)];

            for (int i = 0; i <= 2 * n; i += 2) {
                ps[i].wg = MultiPrecision<N>.Zero;
            }

            using (BinaryReader sr = new(stream)) {
                for (int i = 0; i <= n; i++) {
                    MultiPrecision<Pow2.N64> v = sr.ReadMultiPrecision<Pow2.N64>();

                    ps[i].x = ((1 - v) / 2).Convert<N>();
                    ps[^(i + 1)].x = ((1 + v) / 2).Convert<N>();
                }
                for (int i = 0; i < (n + 1) / 2; i++) {
                    MultiPrecision<Pow2.N64> v = sr.ReadMultiPrecision<Pow2.N64>();
                    MultiPrecision<N> w = (v / 2).Convert<N>();

                    ps[i * 2 + 1].wg = w;
                    ps[^(i * 2 + 2)].wg = w;
                }
                for (int i = 0; i < n + 1; i++) {
                    MultiPrecision<Pow2.N64> v = sr.ReadMultiPrecision<Pow2.N64>();
                    MultiPrecision<N> w = (v / 2).Convert<N>();

                    ps[i].wk = w;
                    ps[^(i + 1)].wk = w;
                }
            }

            return Array.AsReadOnly(ps);
        }
    }
}
