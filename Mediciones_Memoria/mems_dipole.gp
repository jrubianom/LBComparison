set term jpeg enhanced
set o "Mems_Dipole.jpg"
set g
set xlabel sprintf("Refinement level")
set ylabel sprintf("Memory Used by Distribution Functions (MB)")
set xrange [0:3]
set size ratio 1
set key top left
set key box width 0.5 height 1

plot "Datos/Memory_Dipole_Miller.dat" using 1:2 pt 7 ps 2 title "MM Model", "Datos/Memory_Dipole_HV.dat" using 1:2 title "HV Model" pt 7 ps 2
