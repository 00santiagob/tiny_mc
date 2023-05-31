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
CC_LIST = gcc clang
# Sobrescribir CC si se especifica como argumento
ifdef COMPILER
	CC = $(COMPILER)
endif



# CCFLAGS son los flags/opciones que usaremos para compilar el programa.
CFLAGS = -I$(INCLUDE_DIR) -Wall -Wextra -Werror -Wundef -g
# Agregar flags si se especifican como argumento
ifdef CFLAG
	CFLAGS += $(CFLAG)
endif

########################### FLAGS PARA CONTROL DE WARNINGs & ERRORs

# -Wall			This enables all the warnings about constructions that some users consider questionable, and that are easy to avoid (or modify to prevent the warning), even in conjunction with macros.
# -Wextra		This enables some extra warning flags that are not enabled by -Wall.
# -Werror		Make all warnings into errors.
# -Wundef		Warn if an undefined identifier is evaluated in an "#if" directive.  Such identifiers are replaced with zero.

########################### FLAGS PARA INCLUIR LIBRERIAS

# -I$(INCLUDE_DIR)			Ayuda a encontrar los HEADERS en el directorio INCLUDE_DIR.

########################### FLAGS PARA DIFERENTES ESTANDARES

# -std=c89					Para especificar la versión de C de 1989. Se debe compilar según el estándar C89.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c99					Para especificar la versión de C de 1999. Se debe compilar según el estándar C99.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c11					Para especificar la versión de C de 2011. Se debe compilar según el estándar C11.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c17					Para especificar la versión de C de 2017. Se debe compilar según el estándar C17.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.
# -std=c23					Para especificar la versión de C de 2023. Se debe compilar según el estándar C23.  Esto puede cambiar la forma en que el compilador interpreta ciertas construcciones del lenguaje, lo que puede afectar el comportamiento del programa generado.

########################### FLAGS PARA GENERAR LOS ARCHIVOS INTERMEDIOS DEL COMPILADOR

# -save-temps				Store the usual "temporary" intermediate files permanently; place them in the current directory and name them based on the source file.
# -fverbose-asm				Put extra commentary information in the generated assembly code to make it more readable.  This option is generally only of use to those who actually need to read the generated assembly code (perhaps while debugging the compiler itself).

# -fdump-tree-all			Control the dumping at various stages of processing the intermediate language tree to a file.
# -fdump-rtl-all			Produce all the dumps listed above.
# -fdump-ipa-all

########################### FLAGS PARA DEBUG

# -g

# -fopt-info | -fopt-info-options	Controls optimization dumps from various optimization passes. If the -options form is used, options is a list of - separated option keywords to select the dump details and optimizations.
	# -fopt-info-optall-all

	# The following options control which kinds of messages should be emitted:
		# optimized			Print information when an optimization is successfully applied. It is up to a pass to decide which information is relevant.
		# missed			Print information about missed optimizations. Individual passes control which information to include in the output.
		# note 				Print verbose information about optimizations, such as certain transformations, more detailed messages about decisions etc.

		# all 				Print detailed optimization information. This includes optimized, missed, and note.

	# The following option controls the dump verbosity:
		# internals			This option enables additional, more detailed, messages, which are likely to only be of interest to GCC developers.
		
		# ipa 				Enable dumps from all interprocedural optimizations.
		# loop 				Enable dumps from all loop optimizations.
		# inline			Enable dumps from all inlining optimizations.
		# omp 				Enable dumps from all OMP (Offloading and Multi Processing) optimizations.
		# vec 				Enable dumps from all vectorization optimizations.

		# optall			Enable dumps from all optimizations. This is a superset of the optimization groups listed above (menos insternals).

########################### FLAGS PARA OPTIMIZAR

# -Ofast								Disregard strict standards compliance.  -Ofast enables all -O3 optimizations.  It also enables optimizations that are not valid for all standard-compliant programs.  It turns on -ffast-math and the Fortran-specific -fstack-arrays
	# -ffast-math						This option causes the preprocessor macro "__FAST_MATH__" to be defined.  This option is not turned on by any -O option besides -Ofast since it can result in incorrect output for programs that depend on an exact implementation of IEEE or ISO rules/specifications for math functions.
	# -O3								Estas flags activan distintos niveles de optimización. -O1 realiza optimizaciones simples, mientras que -O3 realiza optimizaciones más agresivas y puede aumentar significativamente el tiempo de compilación.
		# -O2
			# -Os						Optimize for size.  -Os enables all -O2 optimizations except those that often increase code size.
				# -O1
					# -ftree-loop-optimize	Perform loop optimizations on trees.  This flag is enabled by default at -O and higher.
					# -fipa-profile		Perform interprocedural profile propagation.  The functions called only from cold functions are marked as cold. Also functions executed once (such as "cold", "noreturn", static constructors or destructors) are identified. Cold functions and loop less parts of functions executed once are then optimized for size.  Enabled by default at -O and higher.
					# -fomit-frame-pointer	Omit the frame pointer in functions that don't need one. This avoids the instructions to save, set up and restore the frame pointer; on many targets it also makes an extra register available.  On some targets this flag has no effect because the standard calling sequence always uses a frame pointer, so it cannot be omitted.
					# -O0
				# -freorder-blocks		Reorder basic blocks in the compiled function in order to reduce number of taken branches and improve code locality.
				# -freorder-functions	Reorder functions in the object file in order to improve code locality.
				# -finline-functions	Consider all functions for inlining, even if they are not declared inline.  The compiler heuristically decides which functions are worth integrating in this way.  If all calls to a given function are integrated, and the function is declared "static", then the function is normally not output as assembler code in its own right.  Enabled at levels -O3, -Os.  Also enabled by -fprofile-use and -fauto-profile.
				# -fstrict-aliasing		Esta flag activa un estricto control de aliasing, lo que puede permitir al compilador realizar optimizaciones más agresivas. Puede ser útil en aplicaciones que acceden a datos mediante punteros.
			# -falign-functions			Esta flag indica al compilador que alinee las funciones en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
			# -falign-jumps				Esta flag indica al compilador que alinee las etiquetas de salto en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
			# -falign-loops				Esta flag indica al compilador que alinee los bucles en la memoria, lo que puede mejorar el rendimiento al aprovechar la caché.
		# -floop-unroll-and-jam			Apply unroll and jam transformations on feasible loops.  In a loop nest this unrolls the outer loop by some factor and fuses the resulting multiple inner loops.  This flag is enabled by default at -O3.
# -funsafe-math-optimizations			Allow optimizations for floating-point arithmetic that assume that arguments and results are valid and may violate IEEE or ANSI standards.  When used at link time, it may include libraries or startup files that change the default FPU control word or other similar optimizations.  This option is not turned on by any -O option
	# -fassociative-math				Allow re-association of operands in series of floating-point operations. NOTE: re-ordering may change the sign of zero as well as ignore NaNs and inhibit or create underflow or overflow
# -fbranch-probabilities				After running a program compiled with -fprofile-arcs, you can compile it a second time using -fbranch-probabilities, to improve optimizations based on the number of times each branch was taken.
# -floop-block							Perform loop nest optimizations.  Same as -floop-nest-optimize.  To use this code transformation, GCC has to be configured with --with-isl to enable the Graphite loop transformation infrastructure.
# -floop-parallelize-all				Use the Graphite data dependence analysis to identify loops that can be parallelized.  Parallelize all the loops that can be analyzed to not contain loop carried dependences without checking that it is profitable to parallelize the loops.
# -faggressive-loop-optimizations		This option tells the loop optimizer to use language constraints to derive bounds for the number of iterations of a loop.  This assumes that loop code does not invoke undefined behavior by for example causing signed integer overflows or out-of-bound array accesses.  The bounds for the number of iterations of a loop are used to guide loop unrolling and peeling and loop exit test optimizations.  This option is enabled by default.
# -fipa-pta								Perform interprocedural pointer analysis and interprocedural modification and reference analysis.  This option can cause excessive memory and compile-time usage on large compilation units.  It is not enabled by default at any optimization level.
# -fvectorize							NO ES de GCC
# -ftree-vectorize	 					Este flag habilita la optimización de vectorización en el árbol de sintaxis abstracta del compilador.
	# -ftree-loop-vectorize
	# -ftree-slp-vectorize
# -fopenmp								Este flag habilita la compilación de código que utiliza OpenMP, una API de programación en paralelo para sistemas multiprocesador.

########################### FLAGS PARA MEMORIA

# -mcmodel=medium 			Este flag especifica que el modelo de memoria utilizado durante la compilación debe ser mediano. Esto significa que el código generado se puede utilizar para programas grandes, pero no tan grandes como los que se podrían generar con el modelo de memoria "large".
# -mcmodel=large 			Este flag especifica que el modelo de memoria utilizado durante la compilación debe ser grande. Esto significa que el código generado se puede utilizar para programas muy grandes.

########################### FLAGS QUE DEPENDEN DE LA MAQUINA

# -mbig-endian			Generate big-endian code.
# -mlittle-endian		Generate little-endian code.
# -mgeneral-regs-only	Generate code which uses only the general-purpose registers. This will prevent the compiler from using floating-point and Advanced SIMD registers but will not impose any restrictions on the assembler.
# -mcmodel=tiny			Generate code for the tiny code model.  The program and its statically defined symbols must be within 1MB of each other. Programs can be statically or dynamically linked.
# -march=cpu-type 		Generate instructions for the machine type.
	# native 			This selects the CPU to generate code for at compilation time by determining the processor type of the compiling machine.
	# x86-64			A generic CPU with 64-bit extensions.
	# i386				Original Intel i386 CPU.
	# etc...

# -mfpu=fpu				Enables support for specific floating-point hardware extensions.
	# fpus_all			Enables support for all single-precision floating-point hardware extensions.
	# fpud_all			Enables support for all single- and double-precision floating-point hardware extensions.
# -mfpu=name 			This specifies what floating-point hardware (or hardware emulation) is available on the target.
	# auto 				The setting auto is the default and is special.  It causes the compiler to select the floating-point and Advanced SIMD instructions based on the settings of -mcpu and -march.
# -mandroid				Compile code compatible with Android platform.

# -m32 | -m64			Generate code for 32-bit or 64-bit.

# -msimd

# -mfma4
# -mprefer-avx128
# -mprefer-vector-width=opt
# -msse
# -msse2
# -msse3
# -mssse3
# -msse4 					Este flag habilita la generación de código que utiliza instrucciones SSE4, que son extensiones del conjunto de instrucciones de procesadores x86 para realizar operaciones de punto flotante y enteros más rápido.
# -msse4a
# -msse4.1
# -msse4.2
# -msse2avx
# -mavx
# -mavx2 					Este flag habilita la generación de código que utiliza instrucciones AVX2, que son extensiones del conjunto de instrucciones de procesadores x86 para realizar cálculos en paralelo.
# -mavx2-512 				Este flag habilita la generación de código que utiliza instrucciones AVX-512, que son extensiones más nuevas del conjunto de instrucciones de procesadores x86 para realizar cálculos en paralelo.
# -mavx2-512f
# -mavx512f
# -mavx512pf
# -mavx512er
# -mavx512cd
# -mavx512vl
# -mavx512bw
# -mavx512dq
# -mavx512ifma
# -mavx512vbmi
# -mavx512vbmi2
# -mavx512bitalg
# -mavx512vpopcntdq
# -mavx5124fmaps
# -mavx512vnni
# -mavx5124vnniw

########################### FLAGS EXTRAS

# -funroll-loops			Unroll loops whose number of iterations can be determined at compile time or upon entry to the loop.

# -fauto-profile			Enable sampling-based feedback-directed optimizations
# -fprofile-use				Esta flag habilita la optimización basada en perfiles, que utiliza información de ejecución para optimizar el código.
# -fprofile-generate		Esta flag habilita la generación de perfiles de ejecución, que puede ser útil para identificar cuellos de botella en el código.


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
SRC_ORIGINAL = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc.c

RNG_MT = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_mtwister.c $(SRC_DIR)/mtwister.c

LAB1 = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab1.c
LAB1_MT = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab1_mtwister.c $(SRC_DIR)/mtwister.c
LAB1_XS = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab1_xoshiro.c

LAB2_MT = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_mtwister.c $(SRC_DIR)/mtwister.c
LAB2_XS = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_xoshiro.c
LAB2_MT_ARR = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_mtwister_arr.c $(SRC_DIR)/mtwister.c
LAB2_XS_ARR = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_xoshiro_arr.c
LAB2_MT_INTR = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_mtwister_intr.c $(SRC_DIR)/mtwister.c
LAB2_XS_INTR = $(SRC_DIR)/wtime.c $(SRC_DIR)/tiny_mc_lab2_xoshiro_intr.c
LAB2_XOR128 = $(SRC_DIR)/wtime.c $(SRC_DIR)/xor128.c $(SRC_DIR)/simdxorshift128plus.c


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
	rm -rf perf.*
	rm -rf *.bin
	rm -rf *.bc
	rm -rf *.i
	rm -rf *.s


# Esta regla indica que si ocurre un error en la compilación de algún archivo, se borren los archivos generados y se detenga el proceso.
.DELETE_ON_ERROR:

# Esta regla se utiliza para indicar que estos objetivos no son archivos.
.PHONY: all clean run perf hex file obj valgrind debug
