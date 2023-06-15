all: plot

%.x: MediaPulseMM.cpp
	g++ $< -o $@ -O3
	./$@

plot: MediaPulseMM.x plot_MediaPulse.py MediaPulseMM.cpp Parameters.h
	python plot_MediaPulse.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
