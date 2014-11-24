set terminal gif size 800,600 animate delay 10
set output 'output.gif'
stats 'datafile' nooutput
set xrange [-2000:2001]
set yrange [-2000:2001]
do for [i=1:int(STATS_blocks)] {
	plot 'datafile' index (i-1) using 2:3 with circle lc rgb 'red'
}