set terminal wxt font 'times new roman,14'
file='p.dat'
set xlabel 't/s'
set ylabel 'pressure/pa'
plot file w l title "pressure"
