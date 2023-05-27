set term jpeg enhanced
set o "data.jpg"
set g

set size ratio 1
set xlabel "z"
set ylabel "E/E_0"
plot "data.dat" u 1:2 w l title "amplitud"
