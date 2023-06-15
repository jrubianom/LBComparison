all: plot

%.x: main.cpp
	g++ $< -o $@ -O3 -O2 -O1
	./$@

plot: main.x GnuplotFiles/plot.gp GnuplotFiles/animation.gp main.cpp AuxiliarLibraries/Parameters.h
	gnuplot GnuplotFiles/plot.gp
	xdg-open Graphics/datas.jpg
	gnuplot GnuplotFiles/animation.gp
	xdg-open Graphics/Animation.gif

clean:
	rm -rf *.dat *.jpg *.x a.out *.pdf *.log *.gif animation.gp Animation Graphics Results
