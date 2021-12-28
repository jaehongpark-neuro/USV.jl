# functions for wavesurfer reading

# voltage scales for the analog singals
function scaledDoubleAnalogData(data::Matrix{Int16},Scales::Matrix{Float64},sc::Matrix{Float64})

    inverseChannelScales=1.0./Scales;
    nScans,nChannels = size(data);
    nCoefficients = size(sc,1) ;
    scaledData=zeros(nScans,nChannels) ;

    if nScans>0 && nChannels>0 && nCoefficients>0
        for j = 1:nChannels
            @inbounds for i = 1:nScans
                datumAsADCCounts = Float64(data[i,j]) ;
                datumAsADCVoltage = sc[nCoefficients,j] ;
                @inbounds for k = (nCoefficients-1):-1:1
                    datumAsADCVoltage = sc[k,j] + datumAsADCCounts*datumAsADCVoltage ;
                end
                scaledData[i,j] = inverseChannelScales[j] * datumAsADCVoltage ;
            end
        end
    end
    return scaledData
end

# load .h5 file from wavesurfer
function loadNIDAQdataFromh5(name)
    fid = h5open(name,"r")
    sweep = name[findlast("_",name)[1]+1:findlast(".",name)[1]-1] # check the sweeps

    # matrix
    timestamps = read(fid[string("sweep_",sweep,"/digitalScans")]) # TTL inputs e.g. laer, camera, fiberphotometry
    analogs = read(fid[string("sweep_",sweep,"/analogScans")])
    channelScales = read(fid["header/AIChannelScales"])
    scalingCoefficients = read(fid["header/AIScalingCoefficients"])

    datas = scaledDoubleAnalogData(analogs,channelScales,scalingCoefficients);
    return datas,timestamps # Matrix for analog data and TTL inputs
end

# downsampling for analog signals and create a new timevector for it
function analog_down(anal::Vector{Float64},sample::Int64,down::Int64)
    # sample = sampling rate of the original signal e.g. 250khz = 250_000
    # down = the downsampled rate e.g. 1khz = 1000
    rate = down/sample
    anal_resample = resample(anal,rate)
    len = Int(length(anal)*rate)
    time_resample = range(0,step=1/down,length=len)

    return anal_resample,time_resample
end

# downsampling for ttl signals

function ttl_down(ttl::Vector{UInt8},sample::Int64,down::Int64)

    rate  = Int(sample/down)
    new_ttl_length = Int(length(ttl)/rate)
    new_ttl = zeros(UInt8,new_ttl_length,)

    for n = 1:new_ttl_length
        new_ttl[n] = ttl[rate*(n-1)+1]
    end
    return new_ttl
end