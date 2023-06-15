all: plot

%.x: Dielecttric_EMLBM.cpp
	g++ $< -o $@
	./$@

plot: Dielecttric_EMLBM.x plot.gp Dielecttric_EMLBM.cpp Parameters.h
	gnuplot plot.gp
	xdg-open datas.jpg

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log
