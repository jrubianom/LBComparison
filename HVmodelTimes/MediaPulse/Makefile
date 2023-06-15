all: plot

%.x: MediaPulseHV.cpp
	g++ $< -o $@ -O3
	./$@

plot: MediaPulseHV.x plot_MediaPulseHV.py MediaPulseHV.cpp Parameters.h
	python plot_MediaPulseHV.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
