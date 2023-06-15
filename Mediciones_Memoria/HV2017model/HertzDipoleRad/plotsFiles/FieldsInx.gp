set term pdf
set out "Results/Comparison.pdf"
set g

set xlabel "x"
set ylabel "B_y / J_0"

plot "Datos/dataBy.txt" u 1:3 w l title "teorico","" u 1:2 title "LB"
set yrange[-0.1:0.1]

set ylabel "-E_z / (J_0 Z_0)"
plot "Datos/dataEz.txt" u 1:3 w l title "teorico","" u 1:2 title "LB"
