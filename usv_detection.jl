
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

# use data itself for an argument
function spectrogram2(data::Array{Float64,1},nfft::Int64)
    npts = length(data)
    noverlap = Int(ceil(0.75*nfft)) # 384
    skip = nfft-noverlap # 128
    nblocks = Int(floor((npts-nfft)/skip) + 1)
    wa = hanning(nfft)
    nout = (nfft >> 1)+1

    # 
    fin = zeros(Float64,nfft)
    tmp = zeros(ComplexF64,nout)
    fout = zeros(Float64,nout,nblocks)
    plan = plan_rfft(fin)

    offset = 0
        @inbounds begin for k = 1:nblocks
        period = (1:nfft) .+ (k-1) *skip
        sig!(data,fin,period,wa)
        mul!(tmp,plan,fin) # calculate matrix product from LinearAlgebra
        fout!(fout,tmp,offset)
        offset+=nout
        end
    end

    fout
end

@with_kw mutable struct stft_params

    audio::Vector{Float64}
    npts::Int64 = length(audio)

    nfft::Int64 = 512
    noverlap::Int64 = Int(ceil(0.75*nfft))
    skip::Int64 = nfft-noverlap
    nblocks::Int64 = Int(floor((npts-nfft)/skip) + 1)
    wa::Vector{Float64} = hanning(nfft)
    nout::Int64 = (nfft >> 1)+1

    fs::Int64 = 250_000
    f = range(0,fs/2,length=nout)

    # for 1d audio signal
    dt::Float64 = 1.0/fs
    t::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = range(0,step=dt,length=npts)
    
    # for 2d spectrogram
    dt_sng::Float64 = skip/fs
    t_sng::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}} = range(0,step=dt_sng,length=nblocks)
    
    # 
    fin::Vector{Float64} = zeros(Float64,nfft)
    tmp = zeros(ComplexF64,nout)
    plan::FFTW.rFFTWPlan{Float64, -1, false, 1, UnitRange{Int64}} = plan_rfft(fin)

    fout::Matrix{Float64} = zeros(Float64,nout,nblocks)
    
    offset::Int64 = 0

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

## loading 
function loading(filename)
    @time d,t = loadNIDAQdataFromh5(filename);
    ## TTL downsampling
    laser = ttl_down(t[:,1],250_000,1000);
    original = d[:,1];
    air_down,time_down = analog_down(original,250_000,1000);
    audio_data = d[:,2];
    @time sp = stft_params(audio=audio_data);
    @time spectrogram!(sp)
    return sp,laser,air_down,time_down
end

