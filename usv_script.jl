using GLMakie

include("usv_temp.jl")



filename = "F:\\USV\\20211209\\PAG-RAm-ChR2-1\\PAG-RAm-ChR2-1-100hz_2021-12-09_0001.h5";


@time d,t = loadNIDAQdataFromh5(filename);

original = d[:,1];
air_down,time_down = analog_down(original,250_000,1000);

figure = Figure(resolution = (1000,600))

spcolor = :black

axiss = Axis[];
push!(axiss,axis_object(figure,spcolor,i,1)) 
for i in 1:3
    push!(axiss,axis_object(figure,spcolor,i,1)) 
end

# lines!(ax,range(0,step=1/250_000,length=length(original)),original,color=(:blue,0.2))
lines!(axiss[1],time_down,air_down,color=:red,linewidth=2)
lines!(axiss[2],range(0,step=1/250_000,length=length(original)),d[:,2],color=:black)
heatmap!(ax3,sng[:,10000:30000]')


sng = spectrogram2(d[:,2],512);

# plot all spectrogram is not good, so I need to 

lsgrid = labelslidergrid!(figure,
    ["slope"],
    Ref(LinRange(0:1000:length(d[:,1])));
    labelkw = Dict([(:textsize, 30)]),
    sliderkw = Dict([(:linewidth, 24)]),
    valuekw = Dict([(:textsize, 30)])
)

set_close_to!(lsgrid.sliders[1], 1000)

sl_sublayout = GridLayout(height = 150)
figure[2, 1] = sl_sublayout
figure[2, 1] = lsgrid.layout


span = lsgrid.sliders[1].value



# initialize plot

fig = Figure(resolution = (840, 460))

# add axis

ax1 = fig[1, 1] = Axis(fig,
    # borders
    aspect = 1, targetlimits = BBox(-10, 10, -10, 10),
    # title
    title = "Sliders Tutorial",
    titlegap = 48, titlesize = 60,
    # x-axis
    xautolimitmargin = (0, 0), xgridwidth = 2, xticklabelsize = 36,
    xticks = LinearTicks(20), xticksize = 18,
    # y-axis
    yautolimitmargin = (0, 0), ygridwidth = 2, yticklabelpad = 14,
    yticklabelsize = 36, yticks = LinearTicks(20), yticksize = 18
)

# darken axes

vlines!(ax1, [0], linewidth = 2)
hlines!(ax1, [0], linewidth = 2)

# create sliders

lsgrid = labelslidergrid!(fig,
    ["slope", "y-intercept"],
    Ref(LinRange(-10:0.01:10));
    labelkw = Dict([(:textsize, 30)]),
    sliderkw = Dict([(:linewidth, 24)]),
    valuekw = Dict([(:textsize, 30)])
)
lsgrid = labelslidergrid!(fig,
    ["slope", "y-intercept"],
    Ref(LinRange(-10:0.01:10));
    formats = [x -> "$(round(x, digits = 2)) V"],
    labelkw = Dict([(:textsize, 30)]),
    sliderkw = Dict([(:linewidth, 24)]),
    valuekw = Dict([(:textsize, 30)])
)
# set starting position for slope

set_close_to!(lsgrid.sliders[1], 1.0)

# layout sliders

sl_sublayout = GridLayout(height = 150)
fig[2, 1] = sl_sublayout
fig[2, 1] = lsgrid.layout

# create listener

slope = lsgrid.sliders[1].value

intercept = lsgrid.sliders[2].value

x = -10:0.01:10

y = @lift($slope .* x .+ $intercept)

# add line plot

line1 = lines!(ax1, x, y, color = :blue, linewidth = 5)

# reset axes limits, if necessary

xlims!(ax1, -10, 10)
ylims!(ax1, -10, 10)

# add scatter plot

rx = -10:0.5:10
ry = rand(length(rx)) .+ -rx * 0.5 .+ 3
scatter1 = scatter!(ax1, rx, ry, color = :red, markersize = 15)

# reset axes limits, if necessary

xlims!(ax1, -10, 10)
ylims!(ax1, -10, 10)

using Meshes
fig = Figure()

ax = Axis(fig[1, 1])

lsgrid = labelslidergrid!(
    fig,
    ["Voltage", "Current", "Resistance"],
    [0:0.1:10, 0:0.1:20, 0:0.1:30];
    formats = "{:.1f}" .* ["V", "A", "Ω"],
    width = 350,
    tellheight = false)

fig[1, 2] = lsgrid.layout

sliderobservables = [s.value for s in lsgrid.sliders]
bars = lift(sliderobservables...) do slvalues...
    [slvalues...]
end

barplot!(ax, bars, color = [:yellow, :orange, :red])
ylims!(ax, 0, 30)

set_close_to!(lsgrid.sliders[1], 5.3)
set_close_to!(lsgrid.sliders[2], 10.2)
set_close_to!(lsgrid.sliders[3], 15.9)



fig = Figure()

ax = Axis(fig[1, 1])

sl_x = Slider(fig[2, 1], range = 0:0.01:10, startvalue = 3)
sl_y = Slider(fig[1, 2], range = 0:0.01:10, horizontal = false, startvalue = 6)

scatter!([sl_x.value],[sl_y.value], color = :red, markersize = 20)
scatter!([1],[2], color = :red, markersize = 20)
limits!(ax, 0, 10, 0, 10)