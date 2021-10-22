set terminal epslatex standalone size 4.5in,3.2in color

set output "plot.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set multiplot
set border 31 linewidth 2

#set logscale x 10

set xlabel "Interaction strength $U/t$"
set ylabel 'density-matrix functional/t'
set xrange [:]
#set xtics ( "256" 256, "1024" 1024, "2048" 2048, "3072" 3072)
set yrange [:]
set ytics nomirror
#set y2range [0.3:0.5]
#set y2label 'Fraction of peak performance'
#set y2tics ( "0.3" 0.3, "0.4" 0.4, "0.5" 0.5, "1.0" 1.0)

#set key width -7 vertical maxrows 2 tmargin right
set key right top

#set logscale y

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 1 pointsize default
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 2 pointsize default
set style line 3 linecolor rgb col3 linewidth 5 dashtype 1 pointtype 3 pointsize default
set style line 4 linecolor rgb col4 linewidth 5 dashtype 1 pointtype 4 pointsize default
set style line 5 linecolor rgb col2 linewidth 5 dashtype 2 pointtype 2 pointsize default
set style line 6 linecolor rgb col3 linewidth 5 dashtype 2 pointtype 3 pointsize default
set style line 7 linecolor rgb col4 linewidth 5 dashtype 2 pointtype 4 pointsize default

plot [0:20][:0.27] "U.dat"  u ($1):($3) with lp ls 1 t "$F_\\mathrm{qc, no\\ noise}/t$", \
"U.dat"  u ($1):($2) with lp ls 5 t "$F_\\mathrm{exact}/t$"

set size 0.5,0.4
set origin 0.2,0.15
set xlabel ""
set ylabel "$10^5 |F_\\mathrm{qc, no\\ noise}-F_\\mathrm{exact}|/t$" offset 1,1.2
set xtics ( "0" 0, "5" 5, "10" 10, "15" 15, "20" 20)
plot [0:20][0:] "U.dat"  u ($1):(100000.0*abs($2-$3)) with lp ls 3 t ""

unset multiplot


set output "bla"
!latex plot.tex
!dvips plot.dvi
!ps2eps -f plot.ps
!ps2pdf plot.ps
!mv plot.pdf U.pdf
!rm plot.eps *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
