unset label
unset terminal
set pointsize 0.75
set xrange [0:360]
set yrange [0:5]
set lmargin 15
set bmargin 7
set tics font ", 15"
set xlabel '$/phi$' font ",20" offset 0,0,0
set ylabel 'h^U' font ",20" offset -2.0,0,0


plot "Fig4_Testing.txt" with lines title ""
set out