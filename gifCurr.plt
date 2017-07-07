reset
set terminal gif animate delay 5 size 800,400 
file='./data/superCurr408.dat'
file2='./data/density408.dat'
set style data line
set output 'superCurr.gif'
stats file nooutput
if(STATS_columns==3){
#	set pm3d # map
#	set hidden3d 
	}
set xrange [] writeback
set yrange [] writeback
set zrange [] writeback
#set yrange [-0.0255:0.0255]
#unset autoscale 

#set zrange [-2e-4:2e-4]
#set xrange [0:10]
#set yrange [0:0.6]


#set nokey
do for [i=1:int(STATS_blocks)] { 
	plot file index (i) u 1:3 title 'current', file2  index i u 1:(0.01*$2) title '0.01*density'
	#splot file index (i) u 1:2:4, file index (i) u 1:2:5
	set xrange restore
set yrange restore
#set zrange restore
 }
