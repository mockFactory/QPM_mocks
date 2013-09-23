# this is the directory where Numerical Recipes code exists, which is not 
# provided in this codebase.
NRLIB= nrlib

CC=gcc
CFLAGS=-O2 -I $(NRLIB)
LINK= -lm 
EXEC=QPM.mock

# the list of NR routines that are required by this code
OBJ_NR = nrutil.o qromo.o midpnt.o  \
	midinf.o polint.o splint.o spline.o \
	zbrent.o qtrap.o trapzd.o ran1.o gasdev.o 
OBJ_NR := $(patsubst %,$(NRLIB)/%,$(OBJ_NR))

# just compile all source files into objects
OBJ := $(patsubst %.c,%.o,$(wildcard src/*.c))

$(EXEC): $(OBJ) $(OBJ_NR)
	$(CC) -o $@ $^ $(LINK)

clean:
	rm -f *.o
	rm -f src/*.o
	rm -f $(NRLIB)/*.o

real-clean: clean 
	rm $(EXEC)

