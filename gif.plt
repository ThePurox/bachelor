reset
set terminal gif animate delay 5 size 1080,720 
file='dcurrent.dat'
set size square
set output 'output.gif'
stats file nooutput
if(STATS_columns==3){set pm3d map
	}
set xrange [] writeback
set yrange [] writeback
set zrange [] writeback

#unset autoscale 

#set zrange [0:2e-1]
#set xrange [0:10]
#set yrange [0:0.6]
#if(STATS_columns==3){set yrange [0:10]}


set nokey
do for [i=1:int(STATS_blocks)] {
if(STATS_columns==3){splot file index (i-1) w pm3d}
if(STATS_columns==2){plot file index (i-1) lc 2}
if(STATS_columns==4){plot file index (i-1) u 1:($2+$3+$4)}
set xrange restore
set yrange restore
set zrange restore
}
