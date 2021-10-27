set terminal epslatex standalone size 4.8in,3.6in color header \
   "\\usepackage{bm}\n"

set output "plot.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set border 31 linewidth 2

#set logscale x 10

set xlabel "inner iteration"
set ylabel ''
set xrange [:]
#set xtics ( "256" 256, "1024" 1024, "2048" 2048, "3072" 3072)
set yrange [:]
set ytics nomirror
#set y2range [0.3:0.5]
#set y2label 'Fraction of peak performance'
#set y2tics ( "0.3" 0.3, "0.4" 0.4, "0.5" 0.5, "1.0" 1.0)

#set key width -10 vertical maxrows 3 tmargin center
set key right top

#set logscale y

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 1 pointsize default
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 2 pointsize default
set style line 3 linecolor rgb col3 linewidth 5 dashtype 1 pointtype 3 pointsize default
set style line 4 linecolor rgb col4 linewidth 5 dashtype 1 pointtype 4 pointsize default
set style line 5 linecolor rgb col2 linewidth 5 dashtype 2 pointtype 2 pointsize default
set style line 6 linecolor rgb col3 linewidth 5 dashtype 2 pointtype 3 pointsize default
set style line 7 linecolor rgb col4 linewidth 5 dashtype 2 pointtype 4 pointsize default



set arrow 1 from 77,0 to 77,0.5 nohead dt "-"
set arrow 2 from 151,0 to 151,0.5 nohead dt "-"
set arrow 3 from 233,0 to 233,0.5 nohead dt "-"
set arrow 4 from 316,0 to 316,0.5 nohead dt "-"

F=0.18436311216650209
plot [:403][0:0.5] "< grep rdmf_obj_eval sim_from_noiseless.out" u ($2):($6) with l ls 3 t "$\\langle \\hat W_1 \\rangle$/t", \
"< grep rdmf_obj_eval sim_from_noiseless.out" u ($2):($8) with l ls 4 t "$\\sum_i c_i^2$", \
F with l ls 6 t "$F_\\mathrm{exact}/t$"

set output "bla"
!latex plot.tex
!dvips plot.dvi
!ps2eps -f plot.ps
!ps2pdf plot.ps
!mv plot.pdf conv_noise_sim_from_noiseless.pdf
!rm plot.eps *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
