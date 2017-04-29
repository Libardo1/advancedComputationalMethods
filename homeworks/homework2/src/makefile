CC = gcc
CC_OPTIONS =   -O2 -Wall
LIBS = -lm

OBJS_SHOCK = structShock.o funcShock.o ioShock.o calibrationShock.o shock.o
INCL_SHOCK = Makefile structShock.h ioShock.h funcShock.h calibrationShock.h

OBJS_SEDOV = structSedov.o Zedov.o
INCL_SEDOV = Makefile structSedov.h

plotshock: shock
	python plotshock.py

plotsedov: sedov
		python plotsedov.py

shock: shock.x
	./shock.x

sedov: sedov.x
	./sedov.x

exec: shock.x sedov.x

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

shock.x:$(OBJS_SHOCK)
	$(CC) $(OBJS_SHOCK) $(LIBS) -o shock.x

sedov.x:$(OBJS_SEDOV)
	$(CC) $(OBJS_SEDOV) $(LIBS) -o sedov.x

clean:
	rm -f $(OBJS_SHOCK) $(OBJS_SEDOV) *~ core* *.x *.dat *.pdf *.txt
