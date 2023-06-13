all: plot

%.x: Dielecttric_EMLBM.cpp
	g++ $< -o $@ -O3 -O2 -O1
	./$@

plot: Dielecttric_EMLBM.x plot.gp animation.gp Dielecttric_EMLBM.cpp Parameters.h
	gnuplot plot.gp
	xdg-open datas.jpg
	gnuplot animation.gp
	xdg-open Animation.gif

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif animation.gp Animation/*
