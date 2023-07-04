# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# List of executables...
EXC = t8_compress

# Include directories...
INCDIR += -I /home/lars/wrk/t8code/build/local/include

# Library directories...
LIBDIR += -L /home/lars/wrk/t8code/build/local/lib

# -----------------------------------------------------------------------------
# Set compiler flags...
# -----------------------------------------------------------------------------

# Set CC and CFLAGS...
CC = g++
CFLAGS = $(INCDIR) -W -Wall -pedantic -O3 -g -static

# Set LDFLAGS...
LDFLAGS = $(LIBDIR) -lt8 -lp4est -lsc -lm -lz

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

.PHONY : all

all: $(EXC)
	rm -f *~

$(EXC): %: %.cxx
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

check: all
	./t8_compress data/temp.xyz 5

clean:
	rm -rf $(EXC) *.o *~ forest*

cppcheck:
	cppcheck --enable=all ./

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.cxx *.h
