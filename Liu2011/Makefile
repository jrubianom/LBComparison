all: plot

%.x: LiuPulse.cpp
	g++ $< -o $@
	./$@ > data.dat

plot: LiuPulse.x plot.gp LiuPulse.cpp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.jpg *.x 
