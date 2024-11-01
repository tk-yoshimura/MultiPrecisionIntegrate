using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using MultiPrecisionIntegrate;

namespace MultiPrecisionIntegrateTest {
    [TestClass]
    public class GaussKronrodIntegralTests {
        [TestMethod]
        public void IntegrateSinTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err) = GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.Pi, order);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussKronrodIntegral<Pow2.N8>.Integrate(MultiPrecision<Pow2.N8>.Sin, MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.Pi, GaussKronrodOrder.G32K65).value), 1e-30);
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
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, 1e-40, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            Assert.AreEqual(0d, (double)(MultiPrecision<Pow2.N8>.E * (MultiPrecision<Pow2.N8>.Cube(MultiPrecision<Pow2.N8>.E) - 1) - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, 1, 4, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest1() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), 0, MultiPrecision<Pow2.N8>.PositiveInfinity, 1e-40, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), 0, MultiPrecision<Pow2.N8>.PositiveInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest2() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), 1, MultiPrecision<Pow2.N8>.PositiveInfinity, 1e-40, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            MultiPrecision<Pow2.N8> expected = "1.353352832366126918939994949724844034076315459095758814681588724e-1";

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), 2, MultiPrecision<Pow2.N8>.PositiveInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest3() {
            Assert.AreEqual(0d, (double)(-1 - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), MultiPrecision<Pow2.N8>.PositiveInfinity, 0, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest4() {
            MultiPrecision<Pow2.N8> expected = "-1.353352832366126918939994949724844034076315459095758814681588724e-1";

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x), MultiPrecision<Pow2.N8>.PositiveInfinity, 2, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest5() {
            Assert.AreEqual(0d, (double)(1 - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, MultiPrecision<Pow2.N8>.NegativeInfinity, 0, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest6() {
            MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.E * MultiPrecision<Pow2.N8>.E;

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, MultiPrecision<Pow2.N8>.NegativeInfinity, 2, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest7() {
            Assert.AreEqual(0d, (double)(-1 - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, 0, MultiPrecision<Pow2.N8>.NegativeInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest8() {
            MultiPrecision<Pow2.N8> expected = -MultiPrecision<Pow2.N8>.E * MultiPrecision<Pow2.N8>.E;

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(MultiPrecision<Pow2.N8>.Exp, 2, MultiPrecision<Pow2.N8>.NegativeInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest9() {
            MultiPrecision<Pow2.N8> expected = MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x * x), MultiPrecision<Pow2.N8>.NegativeInfinity, MultiPrecision<Pow2.N8>.PositiveInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest10() {
            MultiPrecision<Pow2.N8> expected = -MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x * x), MultiPrecision<Pow2.N8>.PositiveInfinity, MultiPrecision<Pow2.N8>.NegativeInfinity, 1e-40, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-40);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest11() {
            MultiPrecision<Pow2.N8> expected = -MultiPrecision<Pow2.N8>.Sqrt(MultiPrecision<Pow2.N8>.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(-x * x), MultiPrecision<Pow2.N8>.PositiveInfinity, MultiPrecision<Pow2.N8>.NegativeInfinity, eps: 0, order: GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-70);
        }

        [TestMethod]
        public void LimitedEvalIntegrateExpTest() {
            (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(MultiPrecision<Pow2.N8>.Sin(16 * x)), 0, 4, 1e-25, discontinue_eval_points: 1024);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedDepthIntegrateExpTest() {
            (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(MultiPrecision<Pow2.N8>.Sin(16 * x)), 0, 4, 1e-25, maxdepth: 10);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedEvalAndDepthIntegrateExpTest() {
            (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(MultiPrecision<Pow2.N8>.Sin(16 * x)), 0, 4, 1e-25, maxdepth: 10, discontinue_eval_points: 1024);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void UnlimitedIntegrateExpTest() {
            (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(MultiPrecision<Pow2.N8>.Sin(16 * x)), 0, 4, 1e-25);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(MultiPrecision<Pow2.N8>.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedEvalIntegrateExpIncrementTest() {
            for (int discontinue_eval_points = 16; discontinue_eval_points <= 8192; discontinue_eval_points *= 2) {
                (MultiPrecision<Pow2.N8> y, MultiPrecision<Pow2.N8> err, long eval_points) = GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(x => MultiPrecision<Pow2.N8>.Exp(MultiPrecision<Pow2.N8>.Sin(16 * x)), 0, 4, eps: 0, discontinue_eval_points: discontinue_eval_points);

                Console.WriteLine(y);
                Console.WriteLine(err);
                Console.WriteLine(eval_points);
            }
        }
    }
}
