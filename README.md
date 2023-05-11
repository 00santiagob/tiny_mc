# Tiny Monte Carlo

### Alumnos:
  * __Balog Santiago Alberto__ (@00santiagob) ~ _santiagobalog1998@gmail.com_
  * __Maddalozzo Nicola__

## Contexto

- [Página en Wikipedia sobre el problema](https://en.wikipedia.org/wiki/Monte_Carlo_method_for_photon_transport)
- [Código original](https://omlc.org/software/mc/) de [Scott Prahl](https://omlc.org/~prahl/)

### Datos sobre el código

* No es memory bank: No esta limitado por el ancho de banda de memoria.
* Toda su limitación esta dada por la potencia de calculo de los procesadores.
* IMPORTA mucho el generador de numeros aleatorios, se pueden conseguir buenos resultados. Este proyecto requiere MUCHOS números al azar.
* No es requerido, pero se pueden hacer visualizaciones.
* Tiene desbalance de carga: los fotones que caen, algunos llegan hasta el fondo y otros quedan ahi nomas.
* Problematica en comun: computar las diferentes partes de algo, no todas tardan lo mismo.

## Instrucciónes

Compilar el programa haciendo:

```bash
❯ make clean
❯ make
```

Se generara un archivo __a.out__ el cual se debe ejecutar en la terminal: `make run`.

## Estructura del Proyecto

      tiny_mc/
        |
        |__ computer_for_test/
        |    |__ asus_zenbook.md
        |    |__ atom.md
		    |    |__ zx81.md
        |
        |__ img/
        |
        |__ src/
        |    |__ params.h
        |    |__ tiny_mc.c
        |    |__ wtime.c
        |    |__ wtime.h
        |
        |__ .clang-format
        |__ .gitignore
        |__ Makefile
        |__ meson.build
        |__ README.md

## Laboratorios

Al finalizar los Labs, escribir sobre cada uno de ellos.

  * Lab 1: Optimización de código, Paralelismo de Instrucciones

  Se realizaron analisis respecto de la ejecución del programa con diferentes compiladores y flags de compilación, además se probó con cantidades mayores de fotones (32768, 131072, 262144, 524288).

```bash
make perf PFLAG="stat -r 16" CFLAG="-O1 -march=native -std=c17 -ffast-math -finline-functions"

make perf PFLAG="stat -r 16" CFLAG="-O1 -march=native -std=c17 -fguess-branch-probability -funroll-loops"

make perf PFLAG="stat -r 16" CFLAG="-O1 -march=native -std=c17 -ffast-math -fguess-branch-probability -funroll-loops -finline-functions"

make perf PFLAG="stat -r 16" CFLAG="-O3 -march=native -std=c17 -ffast-math -finline-functions"

make perf PFLAG="stat -r 16" CFLAG="-O3 -march=native -std=c17 -fguess-branch-probability -funroll-loops"

make perf PFLAG="stat -r 16" CFLAG="-O3 -march=native -std=c17 -ffast-math -fguess-branch-probability -funroll-loops -finline-functions"

make perf PFLAG="stat -r 16" CFLAG="-Ofast -march=native -std=c17 -ffast-math -finline-functions"

make perf PFLAG="stat -r 16" CFLAG="-Ofast -march=native -std=c17 -fguess-branch-probability -funroll-loops"

make perf PFLAG="stat -r 16" CFLAG="-Ofast -march=native -std=c17 -ffast-math -fguess-branch-probability -funroll-loops -finline-functions"
```

  * Lab 2: Limites (CPU & Memory Bound), Paralelismo Vectorial
  * Lab 3:
  * Lab 4: