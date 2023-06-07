all: plot

%.x: time_cost_SE.cpp
	g++ $< -o $@
	./$@

plot: time_cost_SE.x plot_cost.py time_cost_SE.cpp Parameters.h
	python plot_cost.py

clean:
	rm *.dat *.jpg *.x a.out *.pdf *.log *.gif
