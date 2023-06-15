all: plot

hv_pulse.x: HV2017model/MediaPulse1/Dielecttric_EMLBM.cpp
	g++ $< -o hv_pulse.x -O3
	./hv_pulse.x

mm_pulse.x: MMmodel/DielectricPulse/Dielecttric_EMLBM.cpp
	g++ $< -o mm_pulse.x -O3
	./mm_pulse.x

hv_dipole.x: HV2017model/HertzDipoleRad/main.cpp
	g++ $< -o hv_dipole.x -O3
	./hv_dipole.x

mm_dipole.x: MMmodel/HertzDipoleRad/main.cpp
	g++ $< -o mm_dipole.x -O3
	./mm_dipole.x

plot: hv_pulse.x mm_pulse.x hv_dipole.x mm_dipole.x
	gnuplot mems_pulse.gp
	gnuplot mems_dipole.gp
	xdg-open Mems_Dipole.jpg
	xdg-open Mems_Pulse.jpg

clean:
	rm -rf *.txt *.x Datos *.txt Results *.pdf *.dat  *.jpg

