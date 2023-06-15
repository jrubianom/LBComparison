set term jpeg enhanced
set o "Graphics/datas.jpg"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("E/E_0 Electric field")

set key top left
set size ratio 1
plot "Results/data4.dat" u 1:2 w l title "Electric Field"
