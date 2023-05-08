# Paralelismo, modelos, límites y problemas

## Resumen

* Límites.
* Roofline Model.
* Moore's Law
	* Dennard scaling.
* The free lunch is over.
* Bill Dally quote, Pollack's Rule.
* Dark silicon.
* Taxonomía Flynn: SISD, SIMD, MIMD.
* Amdahl's Law vs. Gustafson-Barsis' Law
	* Strong vs. weak scaling.


## Problemas y Limites

Obtener lo máximo de un núcleo.
* Operaciones de punto flotante (GFLOPS).
* Ancho de banda de memoria (GBps).

Dependiendo de la aplicación los límites pueden ser:

* CPU-bound.
* Memory-bound.

(usualmente memory-bound, excepto tests para lucirse: TOP500)

### Performance pico

Usualmente referido a GFLOPS.

`#FPU x Clock x #Core`

`Clock x FMA (vale 2) x #FPU x vectores x #Cores`


> Supuestos:
> Ejecución OoO es capaz de tirar 1op por clock.
> Hay #FPU unidades de punto flotante independientes. ¿Punto flotante simple fp32 o doble fp64?
>	TOP500 usa fp64.


* CPU: BLAS3
* Memoria: STREAM
* Nomenclatura
	* Max: real, medida con un programa.
	* Peak: teórica.
* Eficiencia = Max/Peak
* Eficiencias por debajo del 80% marcan una arquitectura desbalanceada.
