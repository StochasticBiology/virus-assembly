reset
set multiplot
set size 1,0.33
set yrange [-0.1:1.1]
set xrange [0.135: 0.235]
xoff = 0.002
sxoff = 0.0005


# plot capsid yield with temperature and crowding agents
set origin 0,0.66
plot "sim-summary.txt" u ($1 == 120 && $5 == 2.5 && $4 == 65 ? $6+xoff*-1+$10*sxoff : 1/0):($13/12.) title "0 crowders", "" u ($1 == 140 && $5 == 2.5 && $4 == 65 ? $6+xoff*0+$10*sxoff : 1/0):($13/12.) title "20 crowders", "" u ($1 == 160 && $5 == 2.5 && $4 == 65 ? $6+xoff*1+$10*sxoff : 1/0):($13/12.) title "40 crowders"

# plot capsid yield with temperature and capsomer height
set origin 0,0.33
plot "sim-summary.txt" u ($1 == 120 && $5 == 1.5 && $4 == 65 ? $6+xoff*-1+$10*sxoff : 1/0):($13/12.) title "h = 1.5", "" u ($1 == 120 && $5 == 2.5 && $4 == 65 ? $6+xoff*0+$10*sxoff : 1/0):($13/12.) title "h = 2.5", "" u ($1 == 120 && $5 == 3.5 && $4 == 65 ? $6+xoff*1+$10*sxoff : 1/0):($13/12.) title "h = 3.5"

# plot capsid yield with temperature and box dimension (density)
set origin 0,0.
plot "sim-summary.txt" u ($1 == 120 && $5 == 2.5 && $4 == 55 ? $6+xoff*-1+$10*sxoff : 1/0):($13/12.) title "L = 55", "" u ($1 == 120 && $5 == 2.5 && $4 == 65 ? $6+xoff*0+$10*sxoff : 1/0):($13/12.) title "L = 65", "" u ($1 == 120 && $5 == 2.5 && $4 == 75 ? $6+xoff*1+$10*sxoff : 1/0):($13/12.) title "L = 75"
