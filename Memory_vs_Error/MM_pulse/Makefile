all: plot

%.x: main.cpp
	g++ $< -o $@
	./$@

plot: main.x GnuplotFiles/plot.gp GnuplotFiles/animation.gp main.cpp AuxiliarLibraries/Parameters.h
	gnuplot GnuplotFiles/plot.gp
	gnuplot GnuplotFiles/errors.gp
	xdg-open Graphics/datas.jpg
	xdg-open Graphics/data5.jpg
	xdg-open Graphics/Errors.jpg
	gnuplot GnuplotFiles/animation.gp
	xdg-open Graphics/Animation.gif

clean:
	rm -rf *.dat *.jpg *.x a.out *.pdf *.log *.gif animation.gp Animation/ Graphics/ Results/
