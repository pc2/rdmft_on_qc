#F=0.26695476622333536
#F=0.18625440268126447
F=0.1087059437137474

#set terminal png
#set output "min_without_noise.png"
#plot [:5000][:1] "< grep \"rdmf_obj_eval=\" min_without_noise.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise.out" u 2:6 w l t "W",F t "F exact"
#
#set terminal png
#set output "min_without_noise2.png"
#plot [:][:1] "< grep \"rdmf_obj_eval=\" min_without_noise2.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise2.out" u 2:6 w l t "W",F t "F exact"
#
#set terminal png
#set output "min_without_noise3.png"
#plot [:][:1] "< grep \"rdmf_obj_eval=\" min_without_noise3.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise3.out" u 2:6 w l t "W",F t "F exact"
#
#set terminal png
#set output "min_with_noise.png"
#plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise.out" u 2:6 w l t "W",F t "F exact"
#
#set terminal png
#set output "min_with_noise2.png"
#plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise2.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise2.out" u 2:6 w l t "W",F t "F exact"
#
#set terminal png
#set output "min_with_noise3.png"
#plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise3.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise3.out" u 2:6 w l t "W",F t "F exact"
#
#F=0.18625440295141568
#set terminal png
#set output "min_with_noise4.png"
#plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise4.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise4.out" u 2:6 w l t "W",F t "F exact"

#F=0.2651807032752629
F=0.1862544028341545
#F=0.2669547662377094
set terminal png size 1920,1080
set output "qc.png"
ds=0.3
plot [:1000][:3] "< grep \"rdmf_obj_eval=\" qc.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" qc.out" u 2:6 w l t "W",F t "F exact","< grep \"rdmf_obj_eval=\" qc.out" u 2:(10*$8) w l t "10*sum(c^2)","< grep \"before f\" qc.out" u 0:3 w p ps ds pt 1 lc rgb "black" t "f","< grep \"before f\" qc.out" u 0:4 w p ps ds pt 1 lc rgb "black" t "f","< grep \"before f\" qc.out" u 0:5 w p ps ds pt 1 lc rgb "black" t "f","< grep \"before f\" qc.out" u 0:6 w p ps ds pt 1 lc rgb "black" t "f","< grep \"before f\" qc.out" u 0:($3+$4+$5+$6) w p ps ds pt 1 lc rgb "red" t "N",2 lc rgb "red"

