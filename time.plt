
reset
version = 419	
data = './data/'
plots = '~/thesis/plots/'
suf = ''
set terminal cairolatex pdf 
set style data line
set xlabel '$t/\tau$'
set key spacing 1.25
set linetype 1 lc 'black' lw 3 dt 1
set linetype 2 lc 'royalblue' lw 3 dt 2
set linetype 4 lc '#006400' lw 3 dt 5
set linetype 3 lc 'coral' lw 3 dt 4
set linetype 5 lc '#993333' lw 3 dt 3


f = 'freeEnergy'
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'.tex'

set key bottom
set y2tics 
set y2tics nomirror
set ytics nomirror
set x2zeroaxis
set ylabel '$\left[\epsilon\right]$'
set y2label '$\left[\epsilon/\tau\right]$'
plot file u 1:2 title '$F$',file u 1:5 title '$\dot F_{num}$' axes x1y2, file u 1:4 title '$F_\mathrm{tot}$', file u 1:3 title '$\dot F$' axes x1y2
unset y2tics
unset x2zeroaxis
unset y2label
set key top
###############
f = 'power'
set xzeroaxis 
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'.tex'
f = 'dani'
file2 = sprintf('%s%s%03d.dat',data,f,version)
set ylabel '$\left[\epsilon/\tau\right]$'
set key bottom
plot file u 1:($6) title '$-\frac{\gamma}{2}\langle \vec v^2 \rangle$', file u 1:2 title '$\frac{\gamma}{2}\langle\vec v\cdot\nabla\Phi\rangle$', file u 1:3 title '$\frac{\gamma}{2}\langle\vec v\cdot\nabla V_\mathrm{ext}\rangle$', file u 1:4 title '$\frac{\gamma}{2}k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$', file u 1:5 title '$\langle\dot{V}_\mathrm{ext}\rangle$', file2 u 1:(-0.5*$2) title '-1/2*dani'
#plot  file u 1:4 title '$\frac{\gamma}{2}k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$',  file2 u 1:(-0.5*$2) title '-1/2*dani'
set key top
###############
f = 'power'
set xzeroaxis 
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'log.tex'
f = 'dani'
file2 = sprintf('%s%s%03d.dat',data,f,version)
set logscale y
set format y '$10^%T$'
plot file u 1:(abs($6)) title '$-\frac{\gamma}{2}\langle \vec v^2 \rangle$', file u 1:(abs($2)) title '$\frac{\gamma}{2}\langle\vec v\cdot\nabla\Phi\rangle$', file u 1:(abs($3)) title '$\frac{\gamma}{2}\langle\vec v\cdot\nabla V_\mathrm{ext}\rangle$', file u 1:(abs($4)) title '$\frac{\gamma}{2}k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$', file u 1:(abs($5)) title '$\langle\dot{V}_\mathrm{ext}\rangle$', file2 u 1:(abs(-0.5*$2)) title '-1/2*dani'

#plot file u 1:(abs($6)) title '$\langle \vec v^2 \rangle$', file u 1:(abs($2)) title '$-\langle\vec v\cdot\nabla\Phi\rangle$', file u 1:(abs($3)) title '$-\langle\vec v\cdot\nabla V_\mathrm{ext}\rangle$', file u 1:(abs($4)) title '$-k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$'#, file2 u 1:(abs(-0.5*$2)) title '$-k_B T\langle \vec{v} \cdot \nabla \ln \rho \rangle$
unset logscale
unset format
#############
f = 'dissipation'
set key bottom
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'.tex'
f = 'extPower'
file2 = sprintf('%s%s%03d.dat',data,f,version)
f = 'freeEnergy'
file3 = sprintf('%s%s%03d.dat',data,f,version)
f = 'superPower'
file4 = sprintf('%s%s%03d.dat',data,f,version)
plot file u 1:2 title '$P_t$', file u 1:3 title '$P_t^\mathrm{id}$', file u 1:($2-$3) title '$P_t^\mathrm{exc}$', file2 u 1:($2-$3) title '$\chi_t$', file3 u 1:3 title '$\dot F$'	, file4 title '$\int \vec{J}_{\mathrm{sup}}\cdot \vec{J}/\rho d^n \vec{r}$' 
#############
f = 'dissipation'
set key bottom
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'log.tex'
f = 'extPower'
file2 = sprintf('%s%s%03d.dat',data,f,version)
f = 'freeEnergy'
file3 = sprintf('%s%s%03d.dat',data,f,version)
f = 'superPower'
file4 = sprintf('%s%s%03d.dat',data,f,version)
set logscale y
set format y '$10^%T$'
plot file u 1:2 title '$P_t$', file u 1:3 title '$P_t^\mathrm{id}$', file u 1:(abs($2-$3)) title '$P_t^\mathrm{exc}$', file2 u 1:($2-$3) title '$\chi_t$', file3 u 1:(abs($3)) title '$\dot F$', file4 u 1:(abs($2)) title '$\int \vec{J}_{\mathrm{sup}}\cdot \vec{J} d^n \vec{r}$' 

unset logscale
unset format 
###############
f = 'powerAd'
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'.tex'
unset xzeroaxis
plot file u 1:5 title '$\langle \vec v^2 \rangle$', file u 1:2 title '$-\langle\vec v\cdot\nabla\Phi\rangle$', file u 1:3 title '$-\langle\vec v\cdot\nabla V_{ext}\rangle$', file u 1:4 title '$-k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$'
###############
f = 'squares'
set key top
file = sprintf('%s%s%03d.dat',data,f,version)
set output plots.f.suf.'.tex'
f = 'dani'
file2 = sprintf('%s%s%03d.dat',data,f,version)
unset xzeroaxis
#plot for[i=2:10]file u 1:i# title '$\langle \vec v^2 \rangle$', file u 1:2 title '$-\langle\vec v\cdot\nabla\Phi\rangle$', file u 1:3 title '$-\langle\vec v\cdot\nabla V_{ext}\rangle$', file u 1:4 title '$-k_BT \langle\vec v\cdot\nabla \ln \psi \rangle$'
col = "1 2 3 5 6 9"
names = "$(\\nabla\\phi)^2$ $\\nabla\\phi\\cdot\\nabla\\,V_e$ $k_BT\\nabla\\phi\\cdot\\nabla\\ln\\psi$ $(\\nabla\\,V_e)^2$ $k_BT\\nabla\\,V_e\\cdot\\nabla\\ln\\psi$ $(k_BT\\nabla\\ln\\psi)^2$" 

plot for[i=1:words(col)] file u 1:(column(word(col,i)+1)) title word(names,i),file2 u 1:(3*$3) title '3*dani real' , file2 u 1:(8*$4) title '8*fantasy'
###############
f = 'superCurr'
set ylabel '$1/\tau$'
set output plots.f.suf.'.tex'
file = sprintf('%s%s%03d.dat',data,f,version)
times = "1 3 5 10 30 50 99"
set xlabel '$x/\sigma$'
plot for [i=1:words(times)] file u 1:3 index (word(times,i)+0) title '$\vec{J}_\mathrm{super}(t='.word(times,i).'\Delta t)$' 
#plot  file u 1:3 index 1 , file u 1:3 index 30 , file u 1:3 index 99
###############


set output
set terminal x11
reset
