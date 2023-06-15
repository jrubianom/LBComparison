set term jpeg enhanced
set o "Graphics/data.jpg"

set g
set size ratio 1
set xrange [0.25:1]
set xlabel "z/L_z"
set ylabel "E/E_0"
plot "Results/data.dat" u 1:2 w l title "E_{max}",'' u 1:(-1*$2) w l title "E_{max}",'' u 1:3 w l title "Theoretical amplitud",'' u 1:(-1*$3) w l title "Theoretical amplitud",'' u 1:4 w l title "E/E_0 HV"
