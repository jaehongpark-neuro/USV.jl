# plot function

function axis_object(fig,backbone,row,column)
    ax = Axis(fig[row,column],backgroundcolor= :transparent,yticklabelcolor = backbone, 
    ytickcolor = backbone,xticklabelcolor = backbone, 
    xtickcolor = backbone, bottomspinecolor=backbone,topspinecolor=backbone,rightspinecolor=backbone,
    leftspinecolor=backbone)
end

## update Points in a vector array of Point2f0
function update_points(x::Vector{T},y::Vector{R}) where {T,R}
    point = @inbounds [Point2(i,j) for (i,j) in zip(x,y)]
end

function update_points(x::UnitRange{T},y::Vector{R}) where {T,R}
    point = @inbounds [Point2(i,j) for (i,j) in zip(x,y)]
end

function update_points(x::Observable{UnitRange{T}},y::Observable{Vector{R}}) where {T,R}
    point = update_points(x.val,y.val)
end

function update_points(x::Observable{StepRangeLen{T, Base.TwicePrecision{T}, Base.TwicePrecision{T}}},y::Observable{Vector{R}}) where {T,R}
    point = update_points(x.val,y.val)
end

function update_points(x::StepRangeLen{T, Base.TwicePrecision{T}, Base.TwicePrecision{T}},y::Vector{R}) where {T,R}
    point = @inbounds [Point2(i,j) for (i,j) in zip(x,y)]
end

#####
function plot_usv(sp,laser,air_down,time_down)
    figure = Figure(resolution = (1000,600))
    axes = [Axis(figure[i,1]) for i in 1:3]
    sl_x = labelslider!(figure,"current position",1:1:25000,telheight = false,startvalue=3)
    figure[4,1] = sl_x.layout

    # air flow
    x = 1 # as integer index
    air_points = Observable(update_points(time_down[x:x+2000],air_down[x:x+2000]))
    audio_points = Observable(update_points(sp.t[1:250*2000],sp.audio[1:250*2000]))
    sng = Observable(sp.fout[:,x:x+4000]')

    air_points = lift(sl_x.slider.value) do x
        x_range = x:x+5000
        
        ap = update_points(time_down[x_range],air_down[x_range])
        xlims!(axes[1],time_down[x_range[1]],time_down[x_range[end]])
    
        return ap
    end
    
    audio_points = lift(sl_x.slider.value) do x
        
        y = 250*(x-1) +1
        y_range = y:y+5000*250
        
        ap2 = update_points(sp.t[y_range],sp.audio[y_range])
        xlims!(axes[2],sp.t[y_range[1]],sp.t[y_range[end]])
    
        return ap2
    end
    
    ## Node for sng plots
    sng = lift(sl_x.slider.value) do x

        time_z = clamp(time_down[x],0,24.9)
        sind = findfirst(x->x>=time_z,sp.t_sng)

        sg = sp.fout[:,sind:sind+9765]'
        # xlims!(axes[3],time_z,time_z+5)
        return sg
    end
    
    
    lines!(axes[1],air_points,color=:red,linewidth=2)
    lines!(axes[2],audio_points,color=:blue,linewidth=2)
    heatmap!(axes[3],sng,colormap=:hot,colorrange=(0.0, 0.7))

    figure
end


#####
# plot airflow, raw audio, and the spectrogram; specifies the timespan for visualization
function plot_usv(sp,laser,air_down,time_down,timespan::Int64)
    figure = Figure(resolution = (1000,600))
    axes = [Axis(figure[i,1]) for i in 1:3]

    num = Int(sp.fs/1000)
    sl_length = Int((sp.npts - sp.fs*timespan)/(num)) # 1000, 1ms resolution, timespan as second scale; pull back the limit
    sl_x = labelslider!(figure,"current position",1:1:sl_length,telheight = false,startvalue=3)
    figure[4,1] = sl_x.layout

    # air flow
    x = 1 # as integer index
    # initialize
    air_points = Observable(update_points(time_down[x:x+2000],air_down[x:x+2000]))
    laser_points = Observable(update_points(time_down[x:x+2000],laser[x:x+2000]))
    audio_points = Observable(update_points(sp.t[1:250*2000],sp.audio[1:250*2000]))
    sng = Observable(sp.fout[:,x:x+4000]')

    air_points = lift(sl_x.slider.value) do x
        x_range = x:x+timespan*1000
        
        ap = update_points(time_down[x_range],air_down[x_range])
        xlims!(axes[1],time_down[x_range[1]],time_down[x_range[end]])
    
        return ap
    end
    
    laser_points = lift(sl_x.slider.value) do x
        x_range = x:x+timespan*1000
        
        lp = update_points(time_down[x_range],laser[x_range].*4)
        xlims!(axes[1],time_down[x_range[1]],time_down[x_range[end]])
    
        return lp
    end

    audio_points = lift(sl_x.slider.value) do x
        
        y = 250*(x-1) +1
        y_range = y:y+timespan*1000*num
        
        ap2 = update_points(sp.t[y_range],sp.audio[y_range])
        xlims!(axes[2],sp.t[y_range[1]],sp.t[y_range[end]])
    
        return ap2
    end
    
    time_limit = (sl_length/1000) - 0.004 # pull back 4ms to ensure the boundary
    num_col = Int(floor(timespan/sp.dt_sng))
    ## Node for sng plots
    sng = lift(sl_x.slider.value) do x

        time_z = clamp(time_down[x],0,time_limit)
        sind = findfirst(x->x>=time_z,sp.t_sng)

        sg = sp.fout[:,sind:sind+num_col]'
        # xlims!(axes[3],time_z,time_z+5)
        return sg
    end
    
    
    lines!(axes[1],air_points,color=:red,linewidth=2)
    lines!(axes[1],laser_points,color=:cyan,linewidth=1)
    lines!(axes[2],audio_points,color=:blue,linewidth=2)
    heatmap!(axes[3],sng,colormap=:hot,colorrange=(0.0, 0.7))

    figure
end
