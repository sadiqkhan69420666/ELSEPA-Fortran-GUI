# Gnuplot version 5.0 patchlevel 3.
# Plots results from 'elscatm.f'

reset

# Define the appropriate terminal
 set terminal wxt size 900,600 enhanced font 'Arial,12' dashlength 2.0
#set terminal windows size 900,600 enhanced font 'Arial,12'
#set terminal qt size 900,600 font 'Arial,12' enhanced

set encoding utf8  # Character encoding
set mouse
set zero 1.0e-99

# Soft colors
#set linetype  1 lw 2 lc rgb 29390     # UB blue
#set linetype  2 lw 2 lc rgb 13382400  # soft red
#set linetype  3 lw 2 lc rgb 39168     # dark green
#set linetype  4 lw 2 lc rgb 10053171  # brown
#set linetype  5 lw 2 lc rgb 16724940  # pink

set linetype  1 lw 2 lc rgb 255       # blue
set linetype  2 lw 2 lc rgb 16711680  # red
set linetype  3 lw 2 lc rgb 65280     # green
set linetype  4 lw 2 lc rgb 10053171  # brown
set linetype  5 lw 2 lc rgb 16711935  # magenta

set linetype 10 lw 2 lc rgb 0         # black
set linetype 11 lw 2 lc rgb 8421504   # 50 percent black
set xzeroaxis linestyle 11
set yzeroaxis linestyle 11

#set palette rgbformulae 7,5,15
set palette defined ( 1 "blue", 2 "green", 3 "yellow" , 4 "orange", 5 "red")
#set palette defined ( 1 "blue", 2 "green", 3 "yellow" , 4 "white")
#set palette defined ( 1 "red", 2 "orange", 3 "yellow" , 4 "white")
set linetype 100 lw 1 lc rgb 0         # black

set dashtype 1 (8,3)       # dashed
set dashtype 2 (12,3,1,3)  # dot-dashed
set dashtype 3 (4,3)       # dotted

set tmargin at screen 0.90
set bmargin at screen 0.13
set rmargin at screen 0.94
set lmargin at screen 0.15

set size 1,1  # The plot fills the available canvas
set border lw 2.0
set pointsize 0.8
set bars 0.5

#set format xy "10^{%T}"
#set format xy "%.1te%+-3T"
#set format cb "%.1te%+-3T"

#-----------------------------------------------------------------------

set mxtics
set mytics

set grid xtics ytics mxtics mytics ls 11 lw 1, ls 11 lw 1

set format xy "%.1te%+-3T"

set xrange [1e-4:*]
set nologscale x
set logscale y; set mytics 10
set xlabel "{/Symbol q}  (deg)"
set ylabel "DCS  (cm^2)"
plot 'dcsm.dat' u 1:3 t 'coherent' w l ls 1, \
     'dcsm.dat' u 1:4 t 'incoherent' w l ls 2
pause -1  "Type enter to continue..."

set xrange [1e-6:180]
set logscale x; set mxtics 10
set logscale y
plot 'dcsm.dat' u 1:3 t 'coherent' w l ls 1, \
     'dcsm.dat' u 1:4 t 'incoherent' w l ls 2
pause -1  "Type enter to continue..."

set xrange [1e-6:1]
set logscale x
set logscale y
set xlabel "{/Symbol m} = sin^2({/Symbol q}/2)"
plot 'dcsm.dat' u 2:3 t 'coherent' w l ls 1, \
     'dcsm.dat' u 3:4 t 'incoherent' w l ls 2
pause -1  "Type enter to continue..."

set xrange [0:180]
set yrange [*:*]
set nologscale x; set mxtics
set nologscale y; set mytics
set xlabel "{/Symbol q}  (deg)"
set ylabel "spin polarization"
plot 'dcsm.dat' u 1:5 not w l ls 1
pause -1  "Type enter to continue..."

