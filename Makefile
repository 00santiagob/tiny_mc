#---------------#
# INSTRUCCIONES #
#---------------#

# "make"			--> Compila todos los archivos necesarios para generar el binario.
# "make COMPILER=compilador"	--> Compila todos los archivos con el compilador especificado.
# "make CFLAG='cflags'"		--> Compila todos los archivos con los flags de compilación proporcionados, puede ser mas de uno siempre que se lo encapsule con " ó '.
# "make run"			--> Ejecuta el programa.
# "make perf"
# "make valgrind"		--> Ejecuta los test con valgrind.
# "make valgrind VFLAG='vflags'"--> Ejecuta los test con valgrind con los flags de valgrind proporcionados, puede ser mas de uno siempre que se lo encapsule con " ó '.
# "make debug"			--> Ejecuta el programa con gdb.
# "make clean"			--> Elimina los archivos y las carpetas generadas.

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
CFLAGS = -I$(INCLUDE_DIR) -Wall -Wextra -Werror -g
# Agregar flags si se especifican como argumento
ifdef CFLAG
	CFLAGS += $(CFLAG)
endif
# -I$(INCLUDE_DIR)	Ayuda a encontrar los HEADERS en el directorio INCLUDE_DIR.
# -DNDEBUG		Para eliminar los mensajes de debug.
# -O1, -O2, -O3, -Os	Estas flags activan distintos niveles de optimización. -O1 realiza optimizaciones simples, mientras que -O3 realiza optimizaciones más agresivas y puede aumentar significativamente el tiempo de compilación. -Os optimiza para el tamaño del archivo ejecutable.
# -Ofast
# -funroll-loops	Esta flag activa la optimización que desenrolla los bucles for. Esto puede mejorar el rendimiento en bucles que se ejecutan muchas veces.
# -finline-functions	Esta flag indica al compilador que intente alinear las funciones pequeñas en lugar de llamarlas. Esto puede mejorar el rendimiento en funciones que se llaman con frecuencia.
# -fomit-frame-pointer	Esta flag indica al compilador que omita el puntero de marco en la pila, lo que puede aumentar la velocidad de las llamadas a función, pero también hace que el depurador sea menos útil.
# -march=native		Esta flag indica al compilador que genere código para la arquitectura del procesador en el que se está compilando. Esto puede mejorar el rendimiento al aprovechar las características específicas del procesador.
# -ffast-math		Esta flag habilita las optimizaciones de matemáticas rápidas, como la reasociación de operaciones y la eliminación de comprobaciones de error. Puede ser útil en aplicaciones que realizan cálculos intensivos.
# -fno-strict-aliasing	Esta flag desactiva las optimizaciones que dependen de un estricto control de aliasing. Puede ser útil en aplicaciones que acceden a datos mediante punteros.
# -fstrict-aliasing	Esta flag activa un estricto control de aliasing, lo que puede permitir al compilador realizar optimizaciones más agresivas. Puede ser útil en aplicaciones que acceden a datos mediante punteros.
# -falign-functions	Esta flag indica al compilador que alinee las funciones en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -falign-loops		Esta flag indica al compilador que alinee los bucles en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -falign-jumps		Esta flag indica al compilador que alinee las etiquetas de salto en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
# -fprofile-use		Esta flag habilita la optimización basada en perfiles, que utiliza información de ejecución para optimizar el código.
# -fprofile-generate	Esta flag habilita la generación de perfiles de ejecución, que puede ser útil para identificar cuellos de botella en el código.
# -freorder-blocks	Esta flag indica al compilador que reordene los bloques básicos del programa para mejorar la localidad de referencia de la caché.
# -fipa-pta		Esta flag habilita la propagación de análisis de punteros interprocedural, que puede mejorar la precisión de las optimizaciones de punteros.
# -fipa-profile		Esta flag habilita la optimización interprocedural basada en perfiles, que utiliza información de ejecución para optimizar el código.
# -std=c89		Para especificar la versión de C de 1989.
# -std=c99		Para especificar la versión de C de 1999.
# -std=c11		Para especificar la versión de C de 2011.
# -std=c17		Para especificar la versión de C de 2017.


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
SRC_ORIGINAL = $(filter-out $(wildcard $(SRC_DIR)/*xorshift.c) $(wildcard $(SRC_DIR)/*rng*.c), $(wildcard $(SRC_DIR)/*.c))
RNG_XOR = $(filter-out $(SRC_DIR)/tiny_mc.c $(SRC_DIR)/tiny_mc_rng_mersenne_twister.c, $(wildcard $(SRC_DIR)/*.c))
RNG_MT = $(filter-out $(SRC_DIR)/tiny_mc.c $(wildcard $(SRC_DIR)/*xorshift.c), $(wildcard $(SRC_DIR)/*.c))
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

# Esta regla indica que si ocurre un error en la compilación de algún archivo, se borren los archivos generados y se detenga el proceso.
.DELETE_ON_ERROR:

# Esta regla se utiliza para indicar que estos objetivos no son archivos.
.PHONY: all clean run perf hex file obj valgrind debug
