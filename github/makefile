all: energy.png

energy.png: energy.dat
	python energyPlot.py
evolution.png: 2Dplot.py position.dat
	python 2Dplot.py
energy.dat: solid.out
	time ./solid.out 1 > time1.dat
	time ./solid.out 2 > time2.dat
	time ./solid.out 4 > time4.dat
solid.out: solid.c
	gfortran solid.c -fopenmp -o solid.out
clean:
	rm *.png
	rm *.out
	rm *.dat
