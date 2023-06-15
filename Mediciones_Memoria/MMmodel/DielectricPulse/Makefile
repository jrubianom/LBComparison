all: plot

%.x: Dielecttric_EMLBM.cpp
	g++ $< -o $@
	./$@

plot: Dielecttric_EMLBM.x plot.gp Dielecttric_EMLBM.cpp Parameters.h
	gnuplot plot.gp
	gnuplot errors.gp
	xdg-open datas.jpg
	xdg-open data5.jpg
	xdg-open Errors.jpg

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log
