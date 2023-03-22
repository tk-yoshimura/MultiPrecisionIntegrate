# MultiPrecisionIntegrate
 Float MultiPrecision Numerical Integration Implements

## Requirement
.NET 6.0

## Install

[Download DLL](https://github.com/tk-yoshimura/MultiPrecisionIntegrate/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.multiprecision.integrate/)  

- Import MultiPrecision(https://github.com/tk-yoshimura/MultiPrecision)

## Usage
```csharp
// Gauss-Legendre Integrate 32 Points: sin(t) t=0 to pi

GaussLegendreIntegral<Pow2.N8>.Integrate(
    MultiPrecision<Pow2.N8>.Sin, 
    MultiPrecision<Pow2.N8>.Zero, MultiPrecision<Pow2.N8>.PI, 
    n: 32
);
```

## Licence
[MIT](https://github.com/tk-yoshimura/MultiPrecisionIntegrate/blob/main/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
