all: energy.png

energy.png: energy.dat
	python energyPlot.py
evolution.png: 2Dplot.py position.dat
	python 2Dplot.py
energy.dat: serial.out
	./serial.out
serial.out: solid_serial.c
	cc solid_serial.c -o serial.out
clean:
	rm *.png
	rm *.out
	rm *.dat
