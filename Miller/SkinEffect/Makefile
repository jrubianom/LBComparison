all: plot

%.x: Skin_EMLBM.cpp
	g++ $< -o $@
	./$@

plot: Skin_EMLBM.x plot.gp Skin_EMLBM.cpp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.jpg *.x a.out
