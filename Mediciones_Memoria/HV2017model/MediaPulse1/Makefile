all: plot

%.x: Dielecttric_EMLBM.cpp
	g++ $< -o $@
	./$@

plot: Dielecttric_EMLBM.x plot.gp animation.gp Dielecttric_EMLBM.cpp Parameters.h
	gnuplot plot.gp
	gnuplot errors.gp
	xdg-open datas.jpg
	xdg-open data5.jpg
	xdg-open Errors.jpg
	gnuplot animation.gp
	xdg-open Animation.gif

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif animation.gp Animation/*
