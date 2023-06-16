set term jpeg enhanced
set o "Mems_Pulse_Ref.jpg"
set g
set xlabel sprintf("Relative Error")
set ylabel sprintf("Memory Used by Distribution Functions (kB)")
set size ratio 1
set key top right
set key box width 0.5 height 1
set title sprintf("Reflected Pulse")

plot "MM_memory.dat" every ::2 using 1:3 pt 7 ps 2 title "MM Model", "HV_memory.dat" every ::2 using 1:3  title "HV Model" pt 7 ps 2 