using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass]
    public class GaussLaguerreIntegralTests {
        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = GaussLaguerrePoints.MinPoints; n <= GaussLaguerrePoints.MaxPoints; n++) {
                MultiPrecision<Pow2.N8> y = GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x), n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x), 32, f_expscaled: false)), 1e-29);

            for (int n = GaussLaguerrePoints.MinPoints; n <= GaussLaguerrePoints.MaxPoints; n++) {
                MultiPrecision<Pow2.N8> y = GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => 1, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => 1, 32, f_expscaled: true)), 1e-29);
        }

        [TestMethod]
        public void IntegrateXExpTest() {
            for (int n = GaussLaguerrePoints.MinPoints; n <= GaussLaguerrePoints.MaxPoints; n++) {
                MultiPrecision<Pow2.N8> y = GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x) * x * x, 32, f_expscaled: false)), 1e-29);

            for (int n = GaussLaguerrePoints.MinPoints; n <= GaussLaguerrePoints.MaxPoints; n++) {
                MultiPrecision<Pow2.N8> y = GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x / 2) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(16 - GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => MultiPrecision<Pow2.N8>.Exp(-x / 2) * x * x, 40, f_expscaled: false)), 1e-29);

            for (int n = GaussLaguerrePoints.MinPoints; n <= GaussLaguerrePoints.MaxPoints; n++) {
                MultiPrecision<Pow2.N8> y = GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => x * x, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLaguerreIntegral<Pow2.N8>.Integrate((x) => x * x, 32, f_expscaled: true)), 1e-29);
        }
    }
}
