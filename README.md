# MultiPrecisionIntegrate
 MultiPrecision Numerical Integration Implements

## Requirement
.NET 8.0  
AVX2 suppoted CPU. (Intel:Haswell(2013)-, AMD:Excavator(2015)-)  
[MultiPrecision](https://github.com/tk-yoshimura/MultiPrecision)

## Install

[Download DLL](https://github.com/tk-yoshimura/MultiPrecisionIntegrate/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.multiprecision.integrate/)  

## Usage
```csharp
// Gauss-Legendre Integrate 32 Points: sin(t) t=0 to pi
GaussLegendreIntegral<Pow2.N8>.Integrate(
    MultiPrecision<Pow2.N8>.Sin, 
    MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, 
    n: 32
);

// Gauss-Kronrod Adaptive Integrate 7-15: exp(t) t=1 to 4
GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(
    MultiPrecision<Pow2.N8>.Exp, 
    1, 4, 
    eps: 1e-40, 
    order: GaussKronrodOrder.G7K15, 
    maxdepth: 10
);

// Gauss-Kronrod Adaptive Integrate 32-65: exp(-t^2) t=-inf to +inf
GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(
    x => MultiPrecision<Pow2.N8>.Exp(-x * x), 
    MultiPrecision<Pow2.N8>.NegativeInfinity, MultiPrecision<Pow2.N8>.PositiveInfinity, 
    eps: 1e-40, 
    order: GaussKronrodOrder.G32K65, 
    maxdepth: 10
);

// Gauss-Kronrod Adaptive Integrate 32-65: exp(-t^2) t=-inf to +inf
GaussKronrodIntegral<Pow2.N8>.AdaptiveIntegrate(
    x => MultiPrecision<Pow2.N8>.Exp(-x * x), 
    MultiPrecision<Pow2.N8>.NegativeInfinity, MultiPrecision<Pow2.N8>.PositiveInfinity, 
    eps: 0, // Auto epsilon
    order: GaussKronrodOrder.G32K65, 
    maxdepth: 10
);
```

## Licence
[MIT](https://github.com/tk-yoshimura/MultiPrecisionIntegrate/blob/main/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
