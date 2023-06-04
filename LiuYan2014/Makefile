all: plot

%.x: LiuYanInterface.cpp
	g++ $< -o $@
	./$@ > data.dat

plot: LiuYanInterface.x plot.gp LiuYanInterface.cpp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.jpg *.x 
