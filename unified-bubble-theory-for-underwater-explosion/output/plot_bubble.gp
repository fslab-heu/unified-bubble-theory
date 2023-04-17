set terminal wxt font 'times new roman,14'
file='bubble.dat'
set xlabel 't/s'
set ylabel 'radius/m'
set y2label 'migration/m'
plot file w l axis x1y1 title "radius", \
    file using 1:3 w l axis x1y2 title "migration"
