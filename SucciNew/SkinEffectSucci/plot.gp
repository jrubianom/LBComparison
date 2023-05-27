set term jpeg enhanced
set o "data.jpg"
set g

set size ratio 1
plot "data.dat" u 1:2 w l title "E_max",'' u 1:(-1*$2) w l title "E_max",'' u 1:3 w l title "Theoretical amplitud",'' u 1:(-1*$3) w l title "Theoretical amplitud",'' u 1:4 w l title "E/E_0"
