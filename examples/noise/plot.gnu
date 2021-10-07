#F=0.26695476622333536
#F=0.18625440268126447
F=0.1087059437137474

set terminal png
set output "min_without_noise.png"
plot [:5000][:1] "< grep \"rdmf_obj_eval=\" min_without_noise.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise.out" u 2:6 w l t "W",F t "F exact"

set terminal png
set output "min_without_noise2.png"
plot [:][:1] "< grep \"rdmf_obj_eval=\" min_without_noise2.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise2.out" u 2:6 w l t "W",F t "F exact"

set terminal png
set output "min_without_noise3.png"
plot [:][:1] "< grep \"rdmf_obj_eval=\" min_without_noise3.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_without_noise3.out" u 2:6 w l t "W",F t "F exact"

set terminal png
set output "min_with_noise.png"
plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise.out" u 2:6 w l t "W",F t "F exact"

set terminal png
set output "min_with_noise2.png"
plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise2.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise2.out" u 2:6 w l t "W",F t "F exact"

set terminal png
set output "min_with_noise3.png"
plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise3.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise3.out" u 2:6 w l t "W",F t "F exact"

F=0.18625440295141568
set terminal png
set output "min_with_noise4.png"
plot [:500][:1] "< grep \"rdmf_obj_eval=\" min_with_noise4.out" u 2:4 w l t "L","< grep \"rdmf_obj_eval=\" min_with_noise4.out" u 2:6 w l t "W",F t "F exact"
