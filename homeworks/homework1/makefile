serial.png: plot.py serial.txt
	python plot.py
serial.txt: serial
	./serial >serial.txt

serial: solid_serial.c
	cc solid_serial.c -o serial 
clean:
	rm *.png
	rm serial
	rm serial.txt