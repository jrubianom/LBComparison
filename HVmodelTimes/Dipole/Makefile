all: plot

%.x: DipoleHV.cpp
	g++ $< -o $@ -O3
	./$@

plot: DipoleHV.x plot_DipoleHV.py DipoleHV.cpp Parameters.h Theory.h
	python plot_DipoleHV.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
