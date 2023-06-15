all:plot

naim.x: main.cpp
	g++ $< -o naim.x -O3
	./naim.x


plot: plotsFiles/a.gp naim.x
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	make contour

replot:
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	make contour

contour:
	gnuplot plotsFiles/BL.gp

clean:
	rm *.txt *.x Datos/*.txt Results/*.pdf
