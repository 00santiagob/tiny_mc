# Laboratorio 1 - Optimización secuencial (ILP)

## Cosas para hacer

1. Encontrar una métrica de performance del problema (Que sea comparable para cualquier tamaño del problema).
    - K photons per second: 0.001 * #photons / seconds
    - FLOPS/IPS.

2. (Hint) Tomar decisiones sobre dónde mirar primero en el código haciendo profiling. (perf, VTune)

3. Mejorar la performance cambiando cosas, por ejemplo:
    - Compiladores (dpkg -l | grep -w "compiler" | grep -w "C"):
        CC = gcc clang icc nvcc pgi
    - Opciones de compilación.
        
    - Mejoras algorítmicas y/o numéricas (si las hubiera, e.g. RNG).
    - Optimizaciones de cálculos (que no haga el compilador).
    - Unrolling de loops y otras fuentes de ILP (que no haga el compilador).
    - Sistema de memoria: Hugepages y estrategias cache-aware (altamente probable que no rindan hasta agregar paralelismo, ni para sistemas pequeños) --> No nos hace falta.

## Hints

- Automatizar TODO, es una inversión para todo el cuatrimestre:
    - Compilación.
    - Tests para detectar rápidamente problemas en el código.
    - Ejecución y medición de performance.
    - Procesamiento de la salida del programa (salida en CSV es fácil de ingerir).
    - Generación de gráficas.

## Entrega

- Características del hardware y del software:
    - CPU: modelo y velocidad.
        - Poder de cómputo de un core medido con Empirical Roofline Toolkit o LINPACK.
    - Memoria: capacidad, velocidad, cantidad de canales ocupados.
        - Ancho de banda para un core medido con Empirical Roofline Toolkit o STREAM.
    - Compiladores: nombres y versiones.
    - Sistema Operativo: nombre, versión, arquitectura.

- Gráficas de scaling para la versión más rápida obtenida.
    - Performance vs. tamaño del problema (usualmente lin-log).
    - No va a dar scaling lineal, hay que explorar tamaños encontrando relaciones con la jerarquía de memoria.
    - Considerar la calidad estadística de los resultados.

- Optimizaciones probadas y sus resultados.
    - Explicación y mediciones que validen la explicación.
    - Intentar medir las causas además de la performance.