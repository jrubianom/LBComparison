all: plot

%.x: time_error.cpp
	g++ $< -o $@
	./$@

plot: time_error.x plots.py time_error.cpp Parameters.h
	python plots.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
