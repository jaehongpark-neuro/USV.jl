using GLMakie

include("usv_temp.jl")

filename = "F:\\USV\\20211230\\PAG-RAm-ChR2(R)-headfix-sqeuak-ketamine-after_2021-12-30_0002.h5"


sp,laser,air_down,time_down = loading(filename);

#####
# within a function
# plot_usv(sp,laser,air_down,time_down)

# within a function
f = plot_usv(sp,laser,air_down,time_down,3)