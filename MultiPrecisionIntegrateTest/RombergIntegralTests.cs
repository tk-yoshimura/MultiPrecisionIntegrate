using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass()]
    public class RombergIntegralTests {
        [TestMethod()]
        public void IntegrateTest() {
            static MultiPrecision<Pow2.N8> f(MultiPrecision<Pow2.N8> x) => MultiPrecision<Pow2.N8>.Sqrt(1 - x * x);

            for (int level = 1; level <= 16; level++) {
                MultiPrecision<Pow2.N8> v = RombergIntegral<Pow2.N8>.Integrate(f, 0, MultiPrecision<Pow2.N8>.Sqrt(2) / 2, level);

                Console.WriteLine($"{level}\t{v}");
            }

            {
                MultiPrecision<Pow2.N8> v = RombergIntegral<Pow2.N8>.Integrate(f, 0, MultiPrecision<Pow2.N8>.Sqrt(2) / 2, 20);
                Assert.AreEqual(0, (double)((MultiPrecision<Pow2.N8>.PI + 2) / 8 - v), 1e-20);
            }
        }
    }
}