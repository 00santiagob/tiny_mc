#---------------#
# INSTRUCCIONES #
#---------------#

# "make"							--> Compila todos los archivos necesarios para generar el binario.
# "make SRC=X"						--> Compila todos los archivos relacionado al archivo principal.
# "make COMPILER=compilador"		--> Compila todos los archivos con el compilador especificado.
# "make CFLAG='cflags'"				--> Compila todos los archivos con los flags de compilación proporcionados, puede ser mas de uno siempre que se lo encapsule con " ó '.
# "make run"						--> Ejecuta el programa.
# "make perf"
# "make valgrind"					--> Ejecuta los test con valgrind.
# "make valgrind VFLAG='vflags'"	--> Ejecuta los test con valgrind con los flags de valgrind proporcionados, puede ser mas de uno siempre que se lo encapsule con " ó '.
# "make debug"						--> Ejecuta el programa con gdb.
# "make clean"						--> Elimina los archivos y las carpetas generadas.

#-----------#
# VARIABLES #
#-----------#

# Variables de directorios
# Directorio de código fuente
SRC_DIR = src
# Directorio de objetos
OBJ_DIR = obj
# Directorio de ejecutables
BIN_DIR = bin
# Directorio de headers
INCLUDE_DIR = include

# CC es la lista de compiladores que vamos a usar. Por ejemplo gcc clang icc
CC = gcc
# Sobrescribir CC si se especifica como argumento
ifdef COMPILER
	CC = $(COMPILER)
endif

# CCFLAGS son los flags/opciones que usaremos para compilar el programa.
CFLAGS = -I$(INCLUDE_DIR) -Wall -Wextra -Werror -Wundef -Wshadow -g
# Agregar flags si se especifican como argumento
ifdef CFLAG
	CFLAGS += $(CFLAG)
endif

########################### FLAGS PARA CONTROL DE WARNINGs & ERRORs

# -Wall
# -Wextra
# -Werror
# -Wundef
# -Wshadow

########################### FLAGS PARA INCLUIR LIBRERIAS

# -I$(INCLUDE_DIR)			Ayuda a encontrar los HEADERS en el directorio INCLUDE_DIR.

########################### FLAGS PARA GENERAR LOS ARCHIVOS INTERMEDIOS DEL COMPILADOR

# -save-temps
# -fverbose-asm				Is useful if you're compiling with -S to examine the assembly output - it adds some informative comments.
# -fdump-tree-all
# -fdump-rtl-all
# -fdump-ipa-all

########################### FLAGS PARA DEBUG

# -g
# -DNDEBUG					Para eliminar los mensajes de debug.

########################### FLAGS PARA OPTIMIZAR

# -O1, -O2, -O3, -Os			Estas flags activan distintos niveles de optimización. -O1 realiza optimizaciones simples, mientras que -O3 realiza optimizaciones más agresivas y puede aumentar significativamente el tiempo de compilación. -Os optimiza para el tamaño del archivo ejecutable.
# -Ofast						Este flag habilita una optimización agresiva del código, lo que puede generar código más rápido pero también menos preciso. Es útil para maximizar el rendimiento en entornos en los que la precisión no es crítica.
# -ffast-math					Esta flag habilita las optimizaciones de matemáticas rápidas, como la reasociación de operaciones y la eliminación de comprobaciones de error. Puede ser útil en aplicaciones que realizan cálculos intensivos.
# -funsafe-math-optimizations

########################### FLAGS PARA MEMORIA

# -falign-functions			Esta flag indica al compilador que alinee las funciones en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -falign-jumps				Esta flag indica al compilador que alinee las etiquetas de salto en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -falign-loops				Esta flag indica al compilador que alinee los bucles en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -finline-functions		Esta flag indica al compilador que intente alinear las funciones pequeñas en lugar de llamarlas. Esto puede mejorar el rendimiento en funciones que se llaman con frecuencia.

# -fstrict-aliasing			Esta flag activa un estricto control de aliasing, lo que puede permitir al compilador realizar optimizaciones más agresivas. Puede ser útil en aplicaciones que acceden a datos mediante punteros.
# -fno-strict-aliasing		Esta flag desactiva las optimizaciones que dependen de un estricto control de aliasing. Puede ser útil en aplicaciones que acceden a datos mediante punteros.

# -freorder-blocks			Esta flag indica al compilador que reordene los bloques básicos del programa para mejorar la localidad de referencia de la caché.

# -mcmodel=medium 			Este flag especifica que el modelo de memoria utilizado durante la compilación debe ser mediano. Esto significa que el código generado se puede utilizar para programas grandes, pero no tan grandes como los que se podrían generar con el modelo de memoria "large".
# -mcmodel=large 			Este flag especifica que el modelo de memoria utilizado durante la compilación debe ser grande. Esto significa que el código generado se puede utilizar para programas muy grandes.

########################### FLAGS PARA PUNTEROS

# -fipa-pta					Esta flag habilita la propagación de análisis de punteros interprocedural, que puede mejorar la precisión de las optimizaciones de punteros.
# -fipa-profile				Esta flag habilita la optimización interprocedural basada en perfiles, que utiliza información de ejecución para optimizar el código.
# -fomit-frame-pointer		Esta flag indica al compilador que omita el puntero de marco en la pila, lo que puede aumentar la velocidad de las llamadas a función, pero también hace que el depurador sea menos útil.

########################### FLAGS PARA VECTORIZACIÓN

# -fvectorize
# -ftree-vectorize	 		Este flag habilita la optimización de vectorización de bucles en el árbol de sintaxis abstracta del compilador.
# -ftree-loop-vectorize
# -ftree-slp-vectorize

# -fopt-info-vec			Este flag habilita la impresión de información de optimización para vectores durante la compilación. Es útil para diagnosticar problemas con la vectorización del código.
# -fopt-info-vec-missed 	Este flag habilita la impresión de información detallada sobre los bucles que no se pudieron vectorizar durante la compilación.

########################### FLAGS PARA PARALELISMO

# -msse4 					Este flag habilita la generación de código que utiliza instrucciones SSE4, que son extensiones del conjunto de instrucciones de procesadores x86 para realizar operaciones de punto flotante y enteros más rápido.
# -mavx2 					Este flag habilita la generación de código que utiliza instrucciones AVX2, que son extensiones del conjunto de instrucciones de procesadores x86 para realizar cálculos en paralelo.
# -mavx2-512 				Este flag habilita la generación de código que utiliza instrucciones AVX-512, que son extensiones más nuevas del conjunto de instrucciones de procesadores x86 para realizar cálculos en paralelo.
# -mavx2-512f

# -fopenmp					Este flag habilita la compilación de código que utiliza OpenMP, una API de programación en paralelo para sistemas multiprocesador.

########################### FLAGS EXTRAS

# -funroll-loops			Esta flag activa la optimización que desenrolla los bucles for. Esto puede mejorar el rendimiento en bucles que se ejecutan muchas veces.

# -fprofile-use				Esta flag habilita la optimización basada en perfiles, que utiliza información de ejecución para optimizar el código.
# -fprofile-generate		Esta flag habilita la generación de perfiles de ejecución, que puede ser útil para identificar cuellos de botella en el código.

########################### FLAGS PARA DIFERENTES ARQUITECTURAS
# Esto permite al compilador generar código específico para cada arquitectura.

# -march=native				Esta flag indica al compilador que genere código para la arquitectura del procesador en el que se está compilando. Esto puede mejorar el rendimiento al aprovechar las características específicas del procesador.
# -march=knl: 				Este flag especifica la arquitectura de destino para la compilación como Intel Knights Landing.
# -march=znver3 			Este flag especifica la arquitectura de destino para la compilación como AMD Zen 3.
# -march=znver4 			Este flag especifica la arquitectura de destino para la compilación como AMD Zen 4.

########################### FLAGS PARA DIFERENTES ESTANDARES

# -std=c89					Para especificar la versión de C de 1989. Se debe compilar según el estándar C89.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c99					Para especificar la versión de C de 1999. Se debe compilar según el estándar C99.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c11					Para especificar la versión de C de 2011. Se debe compilar según el estándar C11.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c17					Para especificar la versión de C de 2017. Se debe compilar según el estándar C17.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c23					Para especificar la versión de C de 2023. Se debe compilar según el estándar C23.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.



# LDFLAGS son las opciones del enlazador.
LDFLAGS = -lm
# Agregar flags si se especifican como argumento
ifdef LDFLAG
	LDFLAGS += $(LDFLAG)
endif

# Flags especificas de Perf
PFLAGS =
# Agregar flags si se especifican como argumento
ifdef PFLAG
	PFLAGS += $(PFLAG)
endif
# stat
# top
# -e <event>		Este flag permite especificar el evento de hardware a monitorear. Los eventos incluyen cosas como ciclos de CPU, instrucciones de ejecución, acceso a la caché, y mucho más.
# -p <pid>		Este flag indica que Perf debe monitorear un proceso específico. El PID se puede especificar como un número o como el nombre del proceso.
# -t <tid>		Este flag indica que Perf debe monitorear un hilo específico. El TID se puede especificar como un número.
# -c <freq>		Este flag indica la frecuencia de muestreo que debe utilizar Perf. Por defecto, Perf muestrea una vez por segundo.
# -a			Este flag indica que Perf debe monitorear todos los procesos del sistema.
# -o <file>		Este flag indica que Perf debe guardar los resultados en un archivo específico.
# -F <format>		Este flag permite especificar el formato de salida de los resultados.
# --call-graph <type>	Este flag indica que Perf debe generar un gráfico de llamadas para el perfilado. El tipo de gráfico de llamadas se puede especificar como "dwarf", "fp", "lbr", entre otros.

# Flags especificas de Valgrind
VFLAGS = --show-reachable=yes --leak-check=full --tool=memcheck --track-origins=yes
# Agregar flags si se especifican como argumento
ifdef VFLAG
	VFLAGS += $(VFLAG)
endif
# --leak-check=full		Este flag indica que Valgrind debe buscar y reportar cualquier memoria que se haya asignado a una aplicación y que no se haya liberado antes de que la aplicación finalice.
# --show-reachable=yes		Este flag indica que Valgrind debe buscar y reportar cualquier memoria que aún sea accesible por la aplicación, aunque se haya asignado pero no liberado.
# --suppressions=filename	Este flag indica que Valgrind debe utilizar un archivo de supresión de errores específico que puede ayudar a reducir la cantidad de mensajes de error que se generan.
# --tool=memcheck		Este flag indica que Valgrind debe utilizar la herramienta de memoria predeterminada "memcheck" para detectar errores de memoria en la aplicación.
# --track-origins=yes		Este flag indica que Valgrind debe rastrear la fuente de cualquier valor no inicializado que se utilice en la aplicación, para ayudar a identificar la causa de los errores.


# Nombre del ejecutable
TARGET = a.out

BINARY = $(BIN_DIR)/$(TARGET)

FILES = $(BINARY)
ifdef FILE
	FILES = $(FILE)
endif

# Estos HEADERS son los que se buscaran con la flag -I$(SRC_DIR).
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)

# Lista de fuentes
#SOURCES = $(wildcard $(SRC_DIR)/*.c)
SRC_ORIGINAL = $(filter-out $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*twister.c) $(wildcard $(SRC_DIR)/*rng*.c) $(wildcard $(SRC_DIR)/*lab*.c), $(wildcard $(SRC_DIR)/*.c))
RNG_XOR = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*twister.c) $(wildcard $(SRC_DIR)/*lab*.c), $(wildcard $(SRC_DIR)/*.c))
RNG_MT = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*lab*.c), $(wildcard $(SRC_DIR)/*.c))
LAB1 = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*twister.c), $(wildcard $(SRC_DIR)/*.c))
LAB1_MT = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*rng*.c) $(SRC_DIR)/tiny_mc_lab1.c, $(wildcard $(SRC_DIR)/*.c))
LAB2_MT = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*rng*.c) $(wildcard $(SRC_DIR)/*lab1*.c), $(wildcard $(SRC_DIR)/*.c))
LAB2_MT_ALL_IN_ONE = $(SRC_DIR)/tiny_mc_lab2_mtwister_all_in_one.c $(SRC_DIR)/wtime.c
SOURCES = $(SRC_ORIGINAL)
ifdef SRC
	SOURCES = $($(SRC))
endif

# Lista de objetos
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

#-----------#
# OBJETIVOS #
#-----------#

# Regla principal
# Si uno solo escribe "make", se ejecuta all, ya que es el primer objetivo que aparece en el Makefile.
all: $(BINARY)

# Regla de compilación
# Se va a generar la carpeta BIN_DIR si no existe
# $(CC) $(CFLAGS) -MD -o $@ $? $(LDFLAGS)
$(BINARY): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# Regla de generación de objetos
# Todos los archivos *.o se crearan en la carpeta OBJ_DIR.
# Los archivos *.c los buscara en el directorio SRC_DIR.
# % es la forma que tiene el Makefile de hacer *, p.e. %.c es igual a *.c
# Se va a generar la carpeta OBJ_DIR si no existe
# -MD
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -MD -c $< -o $@

# Incluimos los archivos de dependencia generados automáticamente
-include $(OBJ_DIR)/*.d

run: $(BINARY)
	./$(BINARY)

perf: $(BINARY)
	perf $(PFLAGS) ./$(BINARY)

valgrind: $(BINARY)
	valgrind $(VFLAGS) ./$(BINARY)

hex: $(BINARY)
	hexdump -C $(FILES)

file: $(BINARY)
	file $(FILES)

obj: $(BINARY)
	objdump -d $(FILES)

# Regla para ejecutar el debuger
debug: $(BINARY)
	gdb ./$(BINARY)

# Regla para limpiar los archivos generados
clean:
	rm -rf $(OBJECTS)
	rm -rf $(OBJ_DIR)
	rm -rf $(TARGET)
	rm -rf $(BIN_DIR)
	rm -rf $(SRC_DIR)/*.gch
	rm -rf *.gdb_history
	rm -rf *.clang-format
	rm -rf *.bin
	rm -rf *.bc
	rm -rf *.i
	rm -rf *.s


# Esta regla indica que si ocurre un error en la compilación de algún archivo, se borren los archivos generados y se detenga el proceso.
.DELETE_ON_ERROR:

# Esta regla se utiliza para indicar que estos objetivos no son archivos.
.PHONY: all clean run perf hex file obj valgrind debug
