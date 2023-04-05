using MultiPrecision;
using System.Collections.ObjectModel;
using System.IO.Compression;

namespace MultiPrecisionIntegrate {

    public static class GaussHermitePoints {
        public const int MinPoints = 4, MaxPoints = 64, MaxLength = 64;
    }

    public static class GaussHermitePoints<N> where N : struct, IConstant {
        public static ReadOnlyDictionary<int, ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)>> Table { private set; get; }

        static GaussHermitePoints() {
            if (MultiPrecision<N>.Length > GaussHermitePoints.MaxLength) {
                throw new NotSupportedException($"Too large multi precision. (N <= {GaussHermitePoints.MaxLength})");
            }

            using MemoryStream stream = new();
            using (MemoryStream resorce = new(Resource.hermite_n64)) {
                using GZipStream decompressor = new(resorce, CompressionMode.Decompress);

                decompressor.CopyTo(stream);
            }

            stream.Position = 0;

            Dictionary<int, ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)>> table = new();

            using (BinaryReader sr = new(stream)) {
                for (int n = GaussHermitePoints.MinPoints; n <= GaussHermitePoints.MaxPoints; n++) {
                    int ns = sr.ReadInt32();
                    if (n != ns) {
                        throw new IOException("The format of resource file is invalid.");
                    }

                    List<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)> vals = new();

                    for (int i = 0; i < (n + 1) / 2; i++) {
                        MultiPrecision<Pow2.N64> x = sr.ReadMultiPrecision<Pow2.N64>();
                        MultiPrecision<Pow2.N64> w = sr.ReadMultiPrecision<Pow2.N64>();
                        MultiPrecision<Pow2.N64> wexp = sr.ReadMultiPrecision<Pow2.N64>();

                        vals.Add((x.Convert<N>(), w.Convert<N>(), wexp.Convert<N>()));
                    }

                    table.Add(n, Array.AsReadOnly(vals.ToArray()));
                }
            }

            Table = new ReadOnlyDictionary<int, ReadOnlyCollection<(MultiPrecision<N> x, MultiPrecision<N> w, MultiPrecision<N> wexp)>>(table);
        }
    }
}
