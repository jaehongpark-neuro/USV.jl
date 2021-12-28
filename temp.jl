using GLMakie, ModernGL
using SparseArrays
include("usv_temp.jl")

filename = "F:\\USV\\20211221\\PAG-RAm-ChR2(left)_2021-12-21_0002.h5";

@time d,t = loadNIDAQdataFromh5(filename);

## TTL downsampling

laser = ttl_down(t[:,1],250_000,1000);

original = d[:,1];
air_down,time_down = analog_down(original,250_000,1000);

audio_data = d[:,2];
@time sp = stft_params(audio=audio_data);
@time spectrogram!(sp)

figure = Figure(resolution = (1000,600))
spcolor = :black

sub_grid = GridLayout()
sub_grid[1:2,1] = [Axis(figure) for i in 1:2]
figure.layout[1,1] = sub_grid

lines!(sub_grid[1,1],time_down,air_down,color=:red,linewidth=2)
xlims!(sub_grid[1,1],0,20)
heatmap!(sub_grid[2,1],sp.t_sng[1:32000],sp.f,sp.fout[:,1:32000]')
xlims!(0,20)


#################
figure = Figure(resolution = (1000,600))
spcolor = :black


axes = [Axis(figure[i,1]) for i in 1:3]

sparse_sng = sparse(sp.fout)
lines!(axes[1],time_down,air_down,color=:red,linewidth=2)
lines!(axes[1],time_down,laser,color=:cyan,linewidth=2)
xlims!(axes[1],0,30)

lines!(axes[2],sp.t,sp.audio)

xlims!(axes[2],0,30)
heatmap!(axes[3],sp.t_sng,sp.f,sp.fout')
heatmap!(axes[3],rand(Float64,10000,10000))
xlims!(axes[3],0,30)


c =rand(5,10)
c[1] = 10.0
heatmap!(axes[2],1:5,1:100,c')

xlims!(axes[2],0,20)
sp.fout[1:10]

heatmap!(axes[2],sp.fout[1:257,1:32000]')

heatmap(sp.fout')

a = rand(Float64,2,40_000)
heatmap(a)


glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE)