#---------------#
# INSTRUCCIONES #
#---------------#

# "make run" ejecuta el programa.
# "make valgrind" ejecuta los test con valgrind.
# "make debug" ejecuta el programa con gdb.
# "make clean" elimina los archivos generados.

#-----------#
# VARIABLES #
#-----------#

# Variables de directorios
SRC_DIR = src               # Directorio de código fuente
OBJ_DIR = obj               # Directorio de objetos
BIN_DIR = bin               # Directorio de ejecutables

# Nombre del ejecutable
TARGET = a.out			# Nombre del ejecutable a generar

# CC es el compilador que vamos a usar.
# Tambien se podria usar el compilador "cclang" o algun otro.
CC = gcc

# CCFLAGS son los flags/opciones que usaremos para compilar el programa.
# -I$(SRC_DIR) ayuda a encontrar el HEADERS en el directorio SRC_DIR.
# Para eliminar los mensajes de debug usar "-DNDEBUG".
CFLAGS = -I$(SRC_DIR) -Wall -Wextra -Werror -g
# Otros flags:
# -std=c11
# -std=c99
# -O0
# -O1
# -O2
# -O3

# LDFLAGS son las opciones del enlazador.
LDFLAGS = -lm

# Estos HEADERS son los que se buscaran con la flag -I$(SRC_DIR).
HEADERS = $(wildcard $(SRC_DIR)/*.h)

# Lista de fuentes
SOURCES = $(wildcard $(SRC_DIR)/*.c)

# Lista de objetos
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

#-----------#
# OBJETIVOS #
#-----------#

# Regla principal
# Si uno solo escribe "make", se ejecuta all, ya que es el primer objetivo que aparece en el Makefile.
all: clean $(BIN_DIR)/$(TARGET)

# Regla de compilación
$(BIN_DIR)/$(TARGET): $(OBJECTS)
		$(CC) $(LDFLAGS) $^ -o $@

# Regla de generación de objetos
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
		$(CC) $(CFLAGS) -c $< -o $@

# Regla para crear los directorios necesarios
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

# Dependencias adicionales para las reglas principales
all: | $(BIN_DIR)

$(OBJECTS): | $(OBJ_DIR)

main: tiny_mc.c $(OBJECTS)
        $(CC) $(CFLAGS) -o $(BIN_DIR)/$(TARGET) $? $(LDFLAGS)

run: $(BIN_DIR)/$(TARGET)
        ./$(BIN_DIR)/$(TARGET)

valgrind: $(BIN_DIR)/$(TARGET)
        valgrind --show-reachable=yes --leak-check=full --tool=memcheck --track-origins=yes ./$(BIN_DIR)/$(TARGET)

# Regla para ejecutar el debuger
debug: $(BIN_DIR)/$(TARGET)
        gdb ./$(BIN_DIR)/$(TARGET)

# Regla para limpiar los archivos generados
clean:
		rm -rf $(OBJ_DIR) $(BIN_DIR) *.gdb_history%

# Esta regla indica que si ocurre un error en la compilación de algún archivo, se borren los archivos generados y se detenga el proceso.
.DELETE_ON_ERROR:

# Esta regla se utiliza para indicar que estos objetivos no son archivos.
.PHONY: all clean run valgrind debug