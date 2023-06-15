all: plot

%.x: main.cpp
	g++ $< -o $@ -O1 -O2 -O3
	./$@

plot: main.x GnuplotFiles/plot.gp main.cpp
	gnuplot GnuplotFiles/plot.gp
	xdg-open Graphics/data.jpg

clean:
	rm -rf data.dat *.jpg *.x a.out Graphics Results
