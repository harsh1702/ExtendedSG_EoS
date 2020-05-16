# Gnuplot plot file to plot ek(T) with k as Liq, Vap

reset
# set style data lp
set style line 1 lw 3 lc rgb "blue"

set key bottom right
set xlabel "Temperature T (K)"
set ylabel "Internal Energy (kJ/kg)"

set title "e_g"
plot "eG_th.txt" u 1:($2*0.001) title "Theoric" w l ls 1, "../input/expData.txt" u 1:($8*0.001) title "Experimental" w points pt 1 ps 1.2 lc rgb "red"
pause(-1)

set title "e_l"
plot "eL_th.txt" u 1:($2*0.001) title "Theoric" w l ls 1, "../input/expData.txt" u 1:($9*0.001) title "Experimental" w points pt 1 ps 1.2 lc rgb "red"
pause(-1)
