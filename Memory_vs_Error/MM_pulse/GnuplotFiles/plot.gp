set term jpeg enhanced
set o "Graphics/datas.jpg"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("A/A_0 Electric field")

set key top left
set size ratio 1
plot "Results/data1.dat" u 1:2 title "Nz = 100" w l,"Results/data2.dat" u 1:2 title "Nz = 200" w l,"Results/data3.dat" u 1:2 title "Nz = 300" w l,"Results/data4.dat" u 1:2 title "Nz = 400" w l,"Results/data5.dat" u 1:2 title "Nz = 500" w l lc "black"
set term jpeg enhanced
set o "Graphics/data5.jpg"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("A/A_0 Electric field")

set size ratio 1
plot "Results/data5.dat" u 1:2 title "Nz = 500" w l
