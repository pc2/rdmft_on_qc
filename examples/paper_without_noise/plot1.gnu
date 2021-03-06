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

set xlabel "outer iteration of the augmented Lagrangian"
set ylabel ''
set xrange [:]
#set xtics ( "256" 256, "1024" 1024, "2048" 2048, "3072" 3072)
set yrange [:]
set ytics nomirror
#set y2range [0.3:0.5]
#set y2label 'Fraction of peak performance'
#set y2tics ( "0.3" 0.3, "0.4" 0.4, "0.5" 0.5, "1.0" 1.0)

set key width -10 vertical maxrows 3 tmargin center
#set key left top

set logscale y

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 1 pointsize default
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 2 pointsize default
set style line 3 linecolor rgb col3 linewidth 5 dashtype 1 pointtype 3 pointsize default
set style line 4 linecolor rgb col4 linewidth 5 dashtype 1 pointtype 4 pointsize default
set style line 5 linecolor rgb col2 linewidth 5 dashtype 2 pointtype 2 pointsize default
set style line 6 linecolor rgb col3 linewidth 5 dashtype 2 pointtype 3 pointsize default
set style line 7 linecolor rgb col4 linewidth 5 dashtype 2 pointtype 4 pointsize default

F=0.18436311216650209
plot [:][:0.1] "U_1/Uopt_True/Mueller_False/conv"  u ($0+1):(abs(F-$1)) with  lp ls 2 t "$|L-F_\\mathrm{exact}|$, $\\boldsymbol \\lambda_0=0$", \
"U_1/Uopt_True/Mueller_False/conv" u ($0+1):(abs(F-$2)) with lp ls 3 t "$|\\langle\\hat W_1\\rangle-F_\\mathrm{exact}|$, $\\boldsymbol \\lambda_0=0$", \
"U_1/Uopt_True/Mueller_False/conv" u ($0+1):(($3)) with  lp ls 4 t "$\\sum_i c_i^2$, $\\boldsymbol\\lambda_0=0$", \
"U_1/Uopt_True/Mueller_True/conv"  u ($0+1):(abs(F-$1)) with  lp ls 5 t "$|L-F_\\mathrm{exact}|$, $\\boldsymbol\\lambda_0=\\boldsymbol \\lambda_\\mathrm{M}$", \
"U_1/Uopt_True/Mueller_True/conv" u ($0+1):(abs(F-$2)) with  lp ls 6 t "$|\\langle \\hat W_1\\rangle-F_\\mathrm{exact}|$, $\\boldsymbol \\lambda_0=\\boldsymbol\\lambda_\\mathrm{M}$", \
"U_1/Uopt_True/Mueller_True/conv" u ($0+1):(($3)) with  lp ls 7 t "$\\sum_i c_i^2$, $\\boldsymbol \\lambda_0=\\boldsymbol\\lambda_\\mathrm{M}$"


set output "bla"
!latex plot.tex
!dvips plot.dvi
!ps2eps -f plot.ps
!ps2pdf plot.ps
!mv plot.pdf conv_U_1_Urot_True_Mueller_True_False.pdf
!rm plot.eps *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
