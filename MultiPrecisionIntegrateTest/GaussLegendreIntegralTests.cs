using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass]
    public class GaussLegendreIntegralTests {
        [TestMethod]
        public void IntegrateSinTest() {
            for (int n = 4; n <= 256; n++) {
                MultiPrecision<Pow2.N8> y = GaussLegendreIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, n);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLegendreIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, 32)), 1e-30);
        }

        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = 4; n <= 256; n++) {
                MultiPrecision<Pow2.N8> y = GaussLegendreIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, n);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.E * (MultiPrecision<Pow2.N8>.Cube(MultiPrecision<Pow2.N8>.E) - 1) - GaussLegendreIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, 32)), 1e-29);
        }
    }
}
