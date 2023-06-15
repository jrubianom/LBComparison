set term pdf
set out "Results/Power.pdf"

dataTeo = "Datos/TeoEPlane.txt"
dataLB = "Datos/EPlane.txt"

stats dataTeo u 2 nooutput
a = STATS_max
stats dataLB u 2 nooutput
b = STATS_max

set polar

set g

f(x) = x/a
g(x) = x/b

plot "Datos/EPlane.txt" u 1:(g($2)) title "S_E LB","Datos/TeoEPlane.txt" u 1:(f($2)) w l title "S_E Teorico","Datos/BPlane.txt" u 1:(g($2)) title "S_B LB","Datos/TeoBPlane.txt" u 1:(f($2)) w l title "S_B Teorico"
plot "Datos/EPlane.txt" u 1:(g($2)) title "S_E LB","Datos/BPlane.txt" u 1:(g($2)) title "S_B LB"
plot "Datos/TeoEPlane.txt" u 1:(f($2)) w l title "S_E Teorico","Datos/TeoBPlane.txt" u 1:(f($2)) w l title "S_B Teorico"
