
# apply hanning window filter to smoothe out the ends
function sig!(data::Array{Float64,1},f::Vector{Float64},ui::UnitRange{Int64},wa::Array{Float64,1})
    of =1
    for i = ui
        @inbounds f[of] = data[i]*wa[of]
        of+=1
    end
    nothing
end

# convert into real number
function fout!(fout::Array{Float64,2},tmp::Array{ComplexF64,1},offset::Int64)
    for i = 1:length(tmp)
        @inbounds fout[offset+i] = abs2(tmp[i])
    end
    nothing
end

@with_kw mutable struct stft_params

    audio::Vector{Float64} # raw audio data
    npts::Int64 = length(audio) # the number of data

    nfft::Int64 = 512 
    noverlap::Int64 = Int(ceil(0.75*nfft))
    skip::Int64 = nfft-noverlap
    nblocks::Int64 = Int(floor((npts-nfft)/skip) + 1) # repeat number of fft
    wa::Vector{Float64} = hanning(nfft) # hanning window for smoothing the ends
    nout::Int64 = (nfft >> 1)+1 # frequency resolution, 257

    fs::Int64 = 250_000 # sampling rate
    f = range(0,fs/2,length=nout)  # frequency range 0 to 125 khz

    # for 1d audio signal
    dt::Float64 = 1.0/fs
    t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = range(0,step=dt,length=npts)
    
    # for 2d spectrogram
    dt_sng::Float64 = skip/fs
    t_sng::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = range(0,step=dt_sng,length=nblocks)
    
    # fft params
    fin::Vector{Float64} = zeros(Float64,nfft)
    tmp = zeros(ComplexF64,nout)
    plan::FFTW.rFFTWPlan{Float64, -1, false, 1, UnitRange{Int64}} = plan_rfft(fin)

    # spectrogram
    fout::Matrix{Float64} = zeros(Float64,nout,nblocks)
    fout_denoised::Matrix{Float64} = zeros(Float64,nout,nblocks)
    
    offset::Int64 = 0

    # usv inds
    usvON::Vector{Int64} = Int64[]
    usvOFF::Vector{Int64} = Int64[]

end

function spectrogram!(param::stft_params)

    @inbounds begin for k = 1:param.nblocks
        period = (1:param.nfft) .+ (k-1) *param.skip
        sig!(param.audio,param.fin,period,param.wa)
        mul!(param.tmp,param.plan,param.fin)
        fout!(param.fout,param.tmp,param.offset)
        param.offset+=param.nout
        end
    end
    param.offset = 0
    nothing
end


function denoise_usv(sng::Matrix{Float64},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # container for denoised spectrogram
    usv2 = deepcopy(sng)

    tind1 = findlast(x->x<15000,f)
    tind2 = findfirst(x->x>40000,f)

    noise_vals = usv2[tind1:tind2,:]
    noise_vals_vec = reshape(noise_vals,(length(noise_vals),))
    noise_sorted = sort(noise_vals_vec)

    thresh_val = noise_sorted[Int(round(length(noise_sorted)*0.95))]

    a =findall(x->x<thresh_val,usv2)

    usv2[a] .= 0

    # % Get rid of energy outside of the USV frequency range

    find1 = findall(x->x<15_000,f)
    find2 = findall(x->x>110_000,f)
    usv2[find1,:] .=0
    usv2[find2,:] .=0
   return usv2
end

# normalize a 2d array in column; sum of each column is one
function normalize_sng!(sng::Matrix{Float64},pow::BitArray{1})
    for i in 1:size(sng,2)
        sum = 0
        for j in 1:size(sng,1)
            @inbounds sum += sng[j,i]
        end
        div,id = !iszero(sum) ? (sum,false) : (1,true)
        for j in 1:size(sng,1)
            @inbounds sng[j,i] /= div
        end
        @inbounds pow[i] = id # zero
    end
    nothing
end

# ith row, jth column sdsng
function dense_sng_diff(sng::Matrix{Float64},frange::UnitRange{Int64},offset::UnitRange{Int64},i::Int64,j::Int64)
    sdsng = 0
    for n = 1:length(frange) # cacluate all in a column
        sdsng  += abs(sng[frange[n],j+1] - sng[frange[n]+offset[i],j])
    end
    return sdsng
end

function specdiscont(sng::Matrix{Float64})
    zero_ind = falses(size(sng,2))
    normalize_sng!(sng,zero_ind)

    moff = 3
    offset = -moff:moff
    frange = 1+moff:size(sng,1)-moff

    sdsng = zeros(length(offset)) # not neccessary...
    dch = zeros(Float64,size(sng,2))
    dch[zero_ind] .= 2.0
    dch[end] = 2.0

    zero_ind .= .!zero_ind # change to non-zero inds
    zero_ind[end] = 0 # make sure no calculation happens beyond the array
    zz = sparse(zero_ind)

    @inbounds for j in zz.nzind # indices for non-zero values
        for i in 1:length(offset)
             sdsng[i] = dense_sng_diff(sng,frange,offset,i,j)
        end
        dch[j] = minimum(sdsng)
    end

    return dch
end

function spectralpurity(sng::Matrix{Float64})
    specpur = zeros(size(sng,2))
    for col = 1:size(sng,2)
        temp = @view sng[:,col]
        specpur[col] = maximum(temp)/sum(temp)
    end
    nanind = findall(x->x==1,isnan.(specpur))
    specpur[nanind] .= 0
    return specpur
end

function meanfrequency(sng::Matrix{Float64},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # sng: spectrogram
    # f: frequency range 0 to 125 khz
    mf = zeros(size(sng,2))
    @inbounds begin for col = 1:size(sng,2)
        up = 0.0
        down = 0.0
        for row = 1:size(sng,1)
             up += f[row]*sng[row,col]
             down += sng[row,col]
        end
        mf[col] = up/down
    end
    end
    nanind = findall(x->x==1,isnan.(mf))
    mf[nanind] .= 0

    return mf
end

function usv_preprocess(sng::Array{Float64,2},dt_sng::Float64,f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})

    mf = meanfrequency(sng,f)
    sp = spectralpurity(sng)
    sd = specdiscont(sng)

    mf[2:end] .= running_median(mf,Int(round(0.005/dt_sng)))
    sp[2:end] .= running_median(sp,Int(round(0.005/dt_sng)))
    sd[2:end] .= running_median(sd,Int(round(0.005/dt_sng)))

    return mf,sp,sd
end

function usv_ind(mf::Vector{Float64},sp::Vector{Float64},sd::Vector{Float64})
    mfthresh=20000;
    spthresh = 0.3;
    sdthresh = 0.85;

    spind = findall(x->x>spthresh,sp)
    mfind = findall(x->x>mfthresh,mf)
    sdind = findall(x->x<sdthresh,sd)

    usvind = intersect(mfind,spind,sdind)
    return usvind
end

function putative_usv(uind::Vector{Int64},t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    usv_putative = zeros(length(t))
    usv_putative[uind] .= 1
    usv_on_off = diff(usv_putative)
    usvON = findall(x->x==1,usv_on_off).-1

    if usv_putative[1] == 1
        prepend!(usvON,[1])
    end
    usvOFF = findall(x->x==-1,usv_on_off).-1
    if usv_putative[end] ==1
        append!(usvOFF,[length(t)])
    end

    return usvON, usvOFF
end

function usv_duration_check(usvON::Vector{Int64},usvOFF::Vector{Int64},dt_sng)
    keepind = Int64[];
    for z = 1:length(usvON)
        if usvOFF[z] > usvON[z] + Int(round(0.005/dt_sng))
            push!(keepind,z)
        end
    end
    usvON = usvON[keepind]
    usvOFF = usvOFF[keepind]
    return usvON, usvOFF
end

function usv_merge(usvON::Vector{Int64},usvOFF::Vector{Int64},dt_sng)
    mergeind = Int64[];
    for p = 1:length(usvOFF)-1
        if usvON[p+1] - usvOFF[p] < Int(round(0.03/dt_sng))
            push!(mergeind,p)
        end
    end
    deleteat!(usvON,mergeind.+1)
    deleteat!(usvOFF,mergeind)
    return usvON, usvOFF
end

function false_positive_usv(usvON::Vector{Int64},usvOFF::Vector{Int64},usv_putative,dt_sng)
    gapThreshold=Int(round(0.25/dt_sng))
    distanceBefore = 0
    distanceAfter = 0

    if length(usvON) > 1
        indicesToThrowOut = Int64[];
        for z = 1:length(usvON)
            if z >1
                distanceBefore=usvON[z]-usvOFF[z-1]
            else
                distanceBefore=usvON[z]
            end
            if z < length(usvON)
                distanceAfter=usvON[z+1]-usvOFF[z]
            else
                distanceAfter=length(usv_putative)-usvOFF[z]
            end
            if distanceAfter>gapThreshold;
                if distanceBefore>gapThreshold;
                    push!(indicesToThrowOut,z)
                end
            end
        end
        deleteat!(usvON,indicesToThrowOut)
        deleteat!(usvOFF,indicesToThrowOut)
    end
    return usvON, usvOFF
end

function usv_processing(uind::Vector{Int64},t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},dt_sng::Float64)
    on, off = putative_usv(uind,t)
    on,off = usv_duration_check(on,off,dt_sng)
    on,off = usv_merge(on,off,dt_sng)
    tt = Array(t)
    on,off = false_positive_usv(on,off,tt,dt_sng)

    return on,off
end

function usv_whole(sng::Matrix{Float64},dt_sng::Float64,t_sng::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},f::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    # sng
    # dt_sng
    # t_sng
    # f : frequency range
    mf2,sp2,sd2 = usv_preprocess(sng,dt_sng,f)
    ui =  usv_ind(mf2,sp2,sd2)
    ON,OFF = usv_processing(ui,t_sng,dt_sng)
    return ON,OFF
end



## loading 
function loading(filename)
    @time d,t = loadNIDAQdataFromh5(filename);
    ## TTL downsampling
    ts = bit.(t[1],1)
    laser = ttl_down(ts,250_000,1000);
    original = d[1];
    air_down,time_down = analog_down(original,250_000,1000);
    audio_data = d[2];
    @time sp = stft_params(audio=audio_data);
    @time spectrogram!(sp)
    @time sp.fout_denoised .= denoise_usv(sp.fout,sp.f)
    @time sp.usvON,sp.usvOFF = usv_whole(sp.fout_denoised,sp.dt_sng,sp.t_sng,sp.f)
    println("$(length(sp.usvON)) USVs are detected")

    return sp,laser,air_down,time_down
end
