set term jpeg enhanced
set o "Mems_Pulse.jpg"
set g
set xlabel sprintf("Refinement level")
set ylabel sprintf("Memory Used by Distribution Functions (kB)")
set xrange [0:6]
set size ratio 1
set key top left
set key box width 0.5 height 1

plot "Memory_Pulse_Miller.dat" using 1:2 pt 7 ps 2 title "MM Model", "Memory_Pulse_HV.dat" using 1:2 title "HV Model" pt 7 ps 2
