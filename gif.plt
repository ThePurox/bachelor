reset
set terminal gif animate delay 5 size 1920,1080 
file='./data/adiabPotential421.dat'
file2='./data/density184.dat'#'konst.dat'
set size square
set output 'output.gif'
stats file nooutput
if(STATS_columns==3){
	set pm3d  map
	set hidden3d 
	}
set xrange [] writeback
set yrange [-1:1]# writeback
set zrange [] writeback
#set yrange [-0.0255:0.0255]
#unset autoscale 

#set zrange [-2e-4:2e-4]
#set xrange [0:10]
#set yrange [0:0.6]
#if(STATS_columns==3){set yrange [0:10]}


set nokey
do for [i=1:int(STATS_blocks)] {
if(STATS_columns==3){splot file index (i-1)#, file2 index (i-1)#,file2 index (20*(i-1)) u ($1+4.95):($2+4.95):3 # w pm3d
	}
if(STATS_columns==2){plot file index (i-1) lc 2}
if(STATS_columns==4){
#	plot file index (i-1) u 1:2, file index (i-1) u 1:4
	plot file index (i-1) u 1:(1*$2+1*$3+1*$4),file2 index (i-1) u 1:(1*$2+1*$3+1*$4)
	}
if(STATS_columns==6){splot file index (i-1) u 1:2:($3-$6)}
if(STATS_columns==8){splot file index (i-1) u 1:2:7}
set xrange restore
#set yrange restore
set zrange restore
print i*100/STATS_blocks
}
reset
set output
set terminal x11
