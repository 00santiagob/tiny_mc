# Tiny Monte Carlo

### Alumnos:
  * __Balog Santiago Alberto__ (@00santiagob) ~ _santiagobalog1998@gmail.com_
  * __Maddalozzo Nicola__

## Resumen

Al finalizar los Labs, escribir sobre cada uno de ellos

  * Lab 1:
  * Lab 2:
  * Lab 3:
  * Lab 4:

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

Se generara un archivo __tiny_mc__ el cual se debe ejecutar en la terminal: `./tiny_mc`.

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

-std=gnu11 -Wall -Wextra -DVERBOSE=0 -O3 -march=native --fast-math