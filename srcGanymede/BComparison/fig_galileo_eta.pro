;filename='../run_hall_conduct_13*_steady/GM/box*.out'
filename='../run_id*G2/GM/box*.outs'

!p.charsize=1.8
!x.thick=5
!y.thick=5
!p.thick=5

log_spacex=7
log_spacey=2

npict=51

;set_device,'fig_galileo_G2.eps',/land

;colors=[255,100,250]
colors=[255,100]
extract_galileo_times, firstpict=npict, lastpict=npict, stretch=1

;close_device,/pdf

;set_default_values
