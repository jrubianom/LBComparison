all: plot

%.x: LB_skin.cpp
	g++ $< -o $@ -O1 -O2 -O3
	./$@

plot: LB_skin.x plot.gp LB_skin.cpp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.jpg *.x a.out
