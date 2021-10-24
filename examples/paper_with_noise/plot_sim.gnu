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


#rdmf_obj_eval= 169 L= 0.6002007803755358 W= 0.20257568359375 sum(c^2)= 0.03463847431361082 t= 4.099917650222778 exact= L= 0.20734663972788656 W= 0.19720937418067067 sum(c^2)= 0.001102458554183269
#rdmf_obj_eval= 285 L= 0.9936461013631458 W= 0.22625732421875 sum(c^2)= 0.030680875861275173 t= 3.4839389324188232 exact= L= 0.20649711369849583 W= 0.22308050001331037 sum(c^2)= 0.0015340045045842752
#rdmf_obj_eval= 396 L= 1.5799549279716987 W= 0.23583984375 sum(c^2)= 0.03112278839461208 t= 3.8507864475250244 exact= L= 0.20155163374607357 W= 0.22422909622724635 sum(c^2)= 0.0015778840744654723
#rdmf_obj_eval= 495 L= 2.3138562402425915 W= 0.232421875 sum(c^2)= 0.028361071694584017 t= 3.6494455337524414 exact= L= 0.19430779425507377 W= 0.2249150805198196 sum(c^2)= 0.0015475815635530112
#rdmf_obj_eval= 588 L= 3.790901532901434 W= 0.2247314453125 sum(c^2)= 0.03218765637753422 t= 3.5248022079467773 exact= L= 0.19811726235834667 W= 0.22606542382157607 sum(c^2)= 0.00162522688290656
#rdmf_obj_eval= 680 L= 5.469314015727477 W= 0.2257080078125 sum(c^2)= 0.029369917019805958 t= 3.707080841064453 exact= L= 0.1765229709588972 W= 0.22656533915473653 sum(c^2)= 0.0015814471421286029
#rdmf_obj_eval= 771 L= 8.423981800655982 W= 0.23089599609375 sum(c^2)= 0.030697200575166658 t= 3.103827714920044 exact= L= 0.08222230194514524 W= 0.22585213505109947 sum(c^2)= 0.0015734803449400385
#rdmf_obj_eval= 896 L= 12.342349105789038 W= 0.231201171875 sum(c^2)= 0.030411253691196464 t= 3.582214593887329 exact= L= -0.761183625448103 W= 0.2224018693961612 sum(c^2)= 0.0021230607893963143
#rdmf_obj_eval= 993 L= 17.76972204393441 W= 0.23394775390625 sum(c^2)= 0.02761154036787614 t= 3.5253713130950928 exact= L= -0.9181643769957026 W= 0.22341248799262975 sum(c^2)= 0.0019440417485916091
#rdmf_obj_eval= 1083 L= 27.067197785112985 W= 0.2257080078125 sum(c^2)= 0.028882236498291507 t= 3.2499608993530273 exact= L= -1.0254623657904651 W= 0.22315579870157087 sum(c^2)= 0.002137001630305005


set arrow 1 from 169,0 to 169,0.5 nohead dt "-"
set arrow 2 from 285,0 to 285,0.5 nohead dt "-"
set arrow 3 from 396,0 to 396,0.5 nohead dt "-"
set arrow 4 from 495,0 to 495,0.5 nohead dt "-"

F=0.18436311216650209
plot [:495][:] "< grep rdmf_obj_eval sim_santiago.out" u ($2):($6) with l ls 3 t "$\\langle \\hat W_1 \\rangle$/t", \
"< grep rdmf_obj_eval sim_santiago.out" u ($2):($8) with l ls 4 t "$\\sum_i c_i^2$", \
F with l ls 6 t "$F_\\mathrm{exact}/t$"

set output "bla"
!latex plot.tex
!dvips plot.dvi
!ps2eps -f plot.ps
!ps2pdf plot.ps
!mv plot.pdf conv_noise_sim.pdf
!rm plot.eps *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
