# Laboratorio 2 - Vectorización

## Cosas para hacer

- Explorar (poco) que se puede hacer respecto a autovectorización:
    - Opciones de compilación.
    - Cambios de código en líneas que previenen la vectorización para que el autovectorizador funcione.
    - Intentar en gcc, icc, clang y otros compiladores.
- Vectorizar lo que tenga sentido a mano con ISPC y/o intrinsics. ([ISPC](https://ispc.github.io/))
    - Leer de a vectores desde memoria.
    - Procesar múltiples elementos simultáneamente.
    - Escribir de a vectores hacia memoria.

## Entrega

- Presentación de los resultados en clase (10 minutos) e informe breve.
    - Explicación de lo que se intentó a nivel de autovectorización por parte del compilador.
    - Detalle de la vectorización realizada.
    - Comparativa contra la mejor versión CPU obtenida en el Laboratorio 1.
    - Gráficas de scaling para distintos tamaños de problema.
    - Potenciales mejoras en la vectorización.