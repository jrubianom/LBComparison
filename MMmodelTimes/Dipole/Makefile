all: plot

%.x: DipoleMM.cpp
	g++ $< -o $@ -O3
	./$@

plot: DipoleMM.x plot_DipoleMM.py DipoleMM.cpp Parameters.h Theory.h
	python plot_DipoleMM.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
