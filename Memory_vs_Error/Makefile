all: plot

mm_pulse.x:MM_pulse/main.cpp
	g++ $< -o mm_pulse.x -O3
	./mm_pulse.x


hv_pulse.x:HV_pulse/main.cpp
	g++ $< -o hv_pulse.x -O3
	./hv_pulse.x

plot: hv_pulse.x mm_pulse.x
	gnuplot mem_refs.gp
	gnuplot mem_trans.gp
	xdg-open Mems_Pulse_Ref.jpg
	xdg-open Mems_Pulse_Trans.jpg

clean:
	rm -rf *.txt *.x Datos *.txt Results Animation Graphics *.pdf *.dat  *.jpg



