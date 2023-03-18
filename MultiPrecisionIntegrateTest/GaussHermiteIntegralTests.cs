using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass]
    public class GaussHermiteIntegralTests {
        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = 4; n <= 128; n++) {
                MultiPrecision<Pow2.N8> y = GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x), n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x), 32, f_expscaled: false)), 1e-29);
            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x), 33, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 128; n++) {
                MultiPrecision<Pow2.N8> y = GaussHermiteIntegral<Pow2.N8>.Integrate((x) => 1, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => 1, 32, f_expscaled: true)), 1e-29);
            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => 1, 33, f_expscaled: true)), 1e-29);
        }

        [TestMethod]
        public void IntegrateXExpTest() {
            for (int n = 4; n <= 128; n++) {
                MultiPrecision<Pow2.N8> y = GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) / 2 - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x) * x * x, 32, f_expscaled: false)), 1e-29);
            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) / 2 - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x) * x * x, 33, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 128; n++) {
                MultiPrecision<Pow2.N8> y = GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x / 2) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI * 2) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x / 2) * x * x, 63, f_expscaled: false)), 1e-25);
            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI * 2) - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x * x / 2) * x * x, 64, f_expscaled: false)), 1e-25);

            for (int n = 4; n <= 128; n++) {
                MultiPrecision<Pow2.N8> y = GaussHermiteIntegral<Pow2.N8>.Integrate((x) => x * x, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) / 2 - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => x * x, 32, f_expscaled: true)), 1e-29);
            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.PI) / 2 - GaussHermiteIntegral<Pow2.N8>.Integrate((x) => x * x, 33, f_expscaled: true)), 1e-29);
        }
    }
}
