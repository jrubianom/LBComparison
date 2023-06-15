all: plot

%.x: SkinEffectHV.cpp
	g++ $< -o $@ -O3
	./$@

plot: SkinEffectHV.x plot_SkinEffectHV.py SkinEffectHV.cpp Parameters.h
	python plot_SkinEffectHV.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
