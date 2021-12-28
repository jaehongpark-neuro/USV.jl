# tutorial for plotting usv data with Makie.jl and Observables
# 12252021 updated

using GLMakie

include("usv_temp.jl")
# freely PAG-RAm
filename = "F:\\USV\\20211221\\freely\\PAG-RAm-ChR2(left)-freely-100hz-1_2021-12-21_0002.h5";
# headfix PAG-RAm
filename = "F:\\USV\\20211221\\PAG-RAm-ChR2(left)_2021-12-21_0002.h5";

# some USVs
# best; many short USV syllables
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-5ms-15ms-5s_2021-10-01_0002.h5"

filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-5ms-15ms-5s_2021-10-01_0003.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-5ms-15ms-5s_2021-10-01_0001.h5"

# sqeuak
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-cont-5s_2021-10-01_0003.h5"

filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-cont-5s_2021-10-01_0004.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-cont-10s_2021-10-01_0001.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-cont-10s_2021-10-01_0002.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-cont-10s_2021-10-01_0003.h5"

filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-5ms-10ms-5s_2021-10-01_0002.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-5ms-10ms-5s_2021-10-01_0004.h5"
filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4mm-10ms-20ms-5s_2021-10-01_0001.h5"

filename = "F:\\USV\\20211001\\CANE-RAm-ChR2-eyfp-(bilateral)-leftstim-4.5mm-cont-10s_2021-10-01_0001.h5"



filename = "F:\\USV\\20211014\\CANE-RAm-ChR2-mCh-L-11.40-burst-200-200_2021-10-14_0001.h5"


sp,laser,air_down,time_down = loading(filename);

#####
# within a function
# plot_usv(sp,laser,air_down,time_down)

# within a function
f = plot_usv(sp,laser,air_down,time_down,3)
save("F:\\labmeeting\\2022\\2022-1st\\20211227\\PAG-RAm-freely-example.png",f)

##

@time d,t = loadNIDAQdataFromh5(filename)
audio_data = d[:,1];
@time sp = stft_params(audio=audio_data);
@time spectrogram!(sp)

laser = ttl_down(t[:,1],250_000,1000);
original = d[:,1];
air_down,time_down = analog_down(original,250_000,1000);

