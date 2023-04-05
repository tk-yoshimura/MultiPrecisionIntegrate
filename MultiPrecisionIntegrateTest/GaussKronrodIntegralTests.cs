using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass]
    public class GaussKronrodIntegralTests {
        [TestMethod]
        public void IntegrateSinTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err) = GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, order);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, GaussKronrodOrder.G32K65).value), 1e-30);
        }

        [TestMethod]
        public void IntegrateExpTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err) = GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, order);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.E * (MultiPrecision<Pow2.N8>.Cube(MultiPrecision<Pow2.N8>.E) - 1) - GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, GaussKronrodOrder.G32K65).value), 1e-29);
        }

        [TestMethod]
        public void AdaptiveIntegrateExpTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, 1e-40, order, depth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }
        }
    }
}
