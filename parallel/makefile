
OBJ_DIR = debug

SOURCE_DIR = source

SOURCE_FILES = $(SOURCE_DIR)/bcond.c \
	$(SOURCE_DIR)/geocalc.c \
	$(SOURCE_DIR)/ioutils.c \
	$(SOURCE_DIR)/main.c \
	$(SOURCE_DIR)/material.c \
	$(SOURCE_DIR)/mesh.c \
	$(SOURCE_DIR)/octree.c \
	$(SOURCE_DIR)/param.c \
	$(SOURCE_DIR)/post.c \
	$(SOURCE_DIR)/decomp.c \
	$(SOURCE_DIR)/reorder.c \
	$(SOURCE_DIR)/rcm.c \
	$(SOURCE_DIR)/parser.c \
	$(SOURCE_DIR)/gradient.c \
	$(SOURCE_DIR)/restart.c \
	$(SOURCE_DIR)/vector.c \
	$(SOURCE_DIR)/matrix.c \
	$(SOURCE_DIR)/itersolv.c \
	$(SOURCE_DIR)/setup.c \
	$(SOURCE_DIR)/solve.c \
	$(SOURCE_DIR)/gamma.c \
	$(SOURCE_DIR)/velocity.c \
	$(SOURCE_DIR)/pressure.c \
	$(SOURCE_DIR)/temperature.c \
	$(SOURCE_DIR)/msolver.c \
	$(SOURCE_DIR)/fill.c 

SOURCE_OBJ = $(OBJ_DIR)/bcond.o \
	$(OBJ_DIR)/geocalc.o \
	$(OBJ_DIR)/ioutils.o \
	$(OBJ_DIR)/main.o \
	$(OBJ_DIR)/material.o \
	$(OBJ_DIR)/mesh.o \
	$(OBJ_DIR)/octree.o \
	$(OBJ_DIR)/param.o \
	$(OBJ_DIR)/post.o \
	$(OBJ_DIR)/decomp.o \
	$(OBJ_DIR)/reorder.o \
	$(OBJ_DIR)/rcm.o \
	$(OBJ_DIR)/parser.o \
	$(OBJ_DIR)/gradient.o \
	$(OBJ_DIR)/restart.o \
	$(OBJ_DIR)/vector.o \
	$(OBJ_DIR)/matrix.o \
	$(OBJ_DIR)/itersolv.o \
	$(OBJ_DIR)/setup.o \
	$(OBJ_DIR)/solve.o \
	$(OBJ_DIR)/gamma.o \
	$(OBJ_DIR)/velocity.o \
	$(OBJ_DIR)/pressure.o \
	$(OBJ_DIR)/temperature.o \
	$(OBJ_DIR)/msolver.o \
	$(OBJ_DIR)/fill.o

# METIS
METIS_LIB = $(METIS_DIR)/$(METIS_ARCH)/libmetis.a

CLEANFILES = $(SOURCE_OBJ) $(LASPACK_OBJ) $(EXECUTABLE) 
include ${PETSC_DIR}/bmake/common/base

# Debug: -g Performance: -O2
CFLAGS = $(PETSC_INCLUDE) -O2

CC = $(CLINKER)

EXECUTABLE = ../examples/OpenFVM

all: $(SOURCE_OBJ) chkopts

	@echo ""
	@echo "***** Building OpenFVM-Flow parallel version *****"
	@echo ""

	$(CC) $(SOURCE_OBJ) $(LASPACK_OBJ) $(METIS_LIB) ${PETSC_LIB} -o $(EXECUTABLE)

$(OBJ_DIR)/bcond.o: $(SOURCE_DIR)/bcond.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/bcond.c -o $(OBJ_DIR)/bcond.o

$(OBJ_DIR)/geocalc.o: $(SOURCE_DIR)/geocalc.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/geocalc.c -o $(OBJ_DIR)/geocalc.o

$(OBJ_DIR)/ioutils.o: $(SOURCE_DIR)/ioutils.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/ioutils.c -o $(OBJ_DIR)/ioutils.o

$(OBJ_DIR)/main.o: $(SOURCE_DIR)/main.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/main.c -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/material.o: $(SOURCE_DIR)/material.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/material.c -o $(OBJ_DIR)/material.o

$(OBJ_DIR)/mesh.o: $(SOURCE_DIR)/mesh.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/mesh.c -o $(OBJ_DIR)/mesh.o

$(OBJ_DIR)/octree.o: $(SOURCE_DIR)/octree.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/octree.c -o $(OBJ_DIR)/octree.o

$(OBJ_DIR)/param.o: $(SOURCE_DIR)/param.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/param.c -o $(OBJ_DIR)/param.o

$(OBJ_DIR)/post.o: $(SOURCE_DIR)/post.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/post.c -o $(OBJ_DIR)/post.o

$(OBJ_DIR)/decomp.o: $(SOURCE_DIR)/decomp.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/decomp.c -o $(OBJ_DIR)/decomp.o

$(OBJ_DIR)/parallel.o: $(SOURCE_DIR)/parallel.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/parallel.c -o $(OBJ_DIR)/parallel.o

$(OBJ_DIR)/reorder.o: $(SOURCE_DIR)/reorder.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/reorder.c -o $(OBJ_DIR)/reorder.o

$(OBJ_DIR)/rcm.o: $(SOURCE_DIR)/rcm.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/rcm.c -o $(OBJ_DIR)/rcm.o

$(OBJ_DIR)/parser.o: $(SOURCE_DIR)/parser.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/parser.c -o $(OBJ_DIR)/parser.o

$(OBJ_DIR)/gradient.o: $(SOURCE_DIR)/gradient.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/gradient.c -o $(OBJ_DIR)/gradient.o

$(OBJ_DIR)/restart.o: $(SOURCE_DIR)/restart.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/restart.c -o $(OBJ_DIR)/restart.o

$(OBJ_DIR)/vector.o: $(SOURCE_DIR)/vector.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/vector.c -o $(OBJ_DIR)/vector.o

$(OBJ_DIR)/matrix.o: $(SOURCE_DIR)/matrix.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/matrix.c -o $(OBJ_DIR)/matrix.o

$(OBJ_DIR)/itersolv.o: $(SOURCE_DIR)/itersolv.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/itersolv.c -o $(OBJ_DIR)/itersolv.o

$(OBJ_DIR)/setup.o: $(SOURCE_DIR)/setup.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/setup.c -o $(OBJ_DIR)/setup.o

$(OBJ_DIR)/solve.o: $(SOURCE_DIR)/solve.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/solve.c -o $(OBJ_DIR)/solve.o

$(OBJ_DIR)/gamma.o: $(SOURCE_DIR)/gamma.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/gamma.c -o $(OBJ_DIR)/gamma.o

$(OBJ_DIR)/velocity.o: $(SOURCE_DIR)/velocity.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/velocity.c -o $(OBJ_DIR)/velocity.o

$(OBJ_DIR)/pressure.o: $(SOURCE_DIR)/pressure.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/pressure.c -o $(OBJ_DIR)/pressure.o

$(OBJ_DIR)/temperature.o: $(SOURCE_DIR)/temperature.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/temperature.c -o $(OBJ_DIR)/temperature.o

$(OBJ_DIR)/msolver.o: $(SOURCE_DIR)/msolver.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/msolver.c -o $(OBJ_DIR)/msolver.o

$(OBJ_DIR)/fill.o: $(SOURCE_DIR)/fill.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/fill.c -o $(OBJ_DIR)/fill.o

rm: 
	rm $(SOURCE_OBJ) $(LASPACK_OBJ)

