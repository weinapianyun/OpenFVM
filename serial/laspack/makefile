# Modify the "-O" optimization option for best performance (-O3 on Linux):
CC = gcc

# Debug: -g Performance: O2
CFLAGS = -O2

AR = ar cr
RANLIB = ranlib

all: liblaspack.a

LASPACK = eigenval.o errhandl.o factor.o itersolv.o matrix.o \
	  operats.o precond.o qmatrix.o rtc.o vector.o

$(LASPACK): laspack.h makefile

liblaspack.a: $(LASPACK)
	@echo ""
	@echo "***** Building LASPACK *****"
	@echo ""
	$(AR) liblaspack.a $(LASPACK)
	$(RANLIB) liblaspack.a

rm:
	rm -f *.o

rmlib:
	rm -f *.a

