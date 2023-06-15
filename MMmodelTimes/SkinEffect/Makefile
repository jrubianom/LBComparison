all: plot

%.x: SkinMM.cpp
	g++ $< -o $@ -O3
	./$@

plot: SkinMM.x plot_SkinMM.py SkinMM.cpp Parameters.h
	python plot_SkinMM.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
