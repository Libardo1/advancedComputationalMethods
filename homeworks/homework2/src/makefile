CC = gcc
CC_OPTIONS =   -O2 -Wall
EXEC = shock.x
LIBS = -lm

OBJS = struct.o func.o io.o calibration.o shock.o 
INCL = Makefile struct.h io.h func.h calibration.h


all: $(EXEC)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
clean:
	rm -f $(OBJS) *~ core* ${EXEC} *.dat
