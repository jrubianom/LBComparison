set term jpeg enhanced
set o "Errors.jpg"
set g
set xlabel sprintf("R Refinement level")
set ylabel sprintf("Relative error")

set logscale
f(x) = m*x+b
g(x) = m2*x+b2
fit f(x) "Errors.dat" every ::1 using (log($1)):(log($4)) via m,b
fit g(x) "Errors.dat" every ::1 using (log($1)):(log($7)) via m2,b2

ff(x) = exp(b)*x**m
gg(x) = exp(b2)*x**m2

Result = sprintf("{/Symbol D} E_r =%.2f*R^{%.2f} \n\n{/Symbol D} E_t=%.2f*R^{%.2f}",exp(b),m,exp(b2),m2)
set obj 1 rect from graph 0, 1 to graph 0.33, 0.85 fc rgb "white" front
set lab 1 Result at graph 0.01, 0.96 front

set size ratio 1

plot "Errors.dat" using 1:4 every ::1 title "E_r",ff(x) title "Fit E_r","Errors.dat" using 1:7 every ::1 title "E_t",gg(x) title "Fit E_t"
