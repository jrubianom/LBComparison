set term jpeg enhanced
set o "datas.jpg"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("E/E_0 Electric field")

set key top left
set size ratio 1
plot "data1.dat" u 1:2 title "Nz = 100" w l,"data2.dat" u 1:2 title "Nz = 200" w l,"data3.dat" u 1:2 title "Nz = 300" w l,"data4.dat" u 1:2 title "Nz = 400" w l,"data5.dat" u 1:2 title "Nz = 500" w l lc "black"
set term jpeg enhanced
set o "data5.jpg"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("E/E_0 Electric field")

set size ratio 1
plot "data5.dat" u 1:2 title "Nz = 500" w l


set term pdf enhanced
set o "data5.pdf"
set g
set xlabel sprintf("z/L_z position")
set ylabel sprintf("E/E_0 Electric field")

set size ratio 1
plot "data5.dat" u 1:2 title "Nz = 500" w l
