using GLMakie


include("usv_temp.jl")


filename = "F:\\USV\\20211221\\NTS-RAm\\NTS-RAm-ChR2(L)-headfix-cont-2s_2021-12-21_0001.h5"


sp,laser,air_down,time_down = loading(filename);

#####
# within a function
# plot_usv(sp,laser,air_down,time_down)

# within a function
f = plot_usv(sp,laser,air_down,time_down,10)
save("F:\\labmeeting\\2022\\2022-1st\\20211227\\NTS-RAm-2s-cont.png",f)