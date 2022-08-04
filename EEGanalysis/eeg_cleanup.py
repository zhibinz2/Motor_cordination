#%% Step 1: LOAD MODULES
# load modules
import sys
import os
import numpy as np
import pyarrow.feather as feather
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal
from sklearn.decomposition import FastICA

# from hdfstorage import savemat #UCI version
#%% Step 2: ADD TMSI python code to the path.  We should include this in some better way in the future.
# add path
p = os.getcwd()
print(p)
sys.path.insert(0, os.path.abspath(p + "/tmsi-python-interface-main/"))
from TMSiSDK.file_readers import Poly5Reader

# %% Step 3: Set the Data File Path
datapath = "/home/ramesh/Projects/flanker/PDC/p016/"
fnameList = os.listdir(datapath)

# %% Step 4: get the list of all EEG files in directory
# list eeg files
files = [f for f in fnameList if "Poly5" in f]
print("eeg files: ", files)

#%% Step 5: Setup parameters for the EEG extraction and variance and ICA cleanup.
plot = "y"  # make plots as you go along
diagnostics = "y"  # make diagnostic plots, in case something is wrong.
nEEGchan = 32  # number of EEG channels
eyechans = [0, 1, 2]  # list of eyechannels.  Must be at least 1 channel.
# I suggest Fp1,Fpz, and Fp2
stimchan = 40  # channel that has the stimulus information for segmentation.
# TMSi channel named Digi if you use digital inputs from StimTracker.
hipass = 0.5  # 0.25 is a good choice if you want readiness potentials or P300
lowpass = 50  # you could notch the 60 Hz and go to 100 Hz. this is simpler.
prestim = 0.75  # pre stimulus interval in seconds (always pad with extra 0.25 s)
poststim = 1.25  # post stimulus interval in seconds. (always pad with extra 0.25 s)
var_threshold = 2.5  # normalized variance threshold to reject trials.
chan_threshold = 2.5  #to identify 
# (smaller values are stricter)
corr_threshold = 0.4  # threshold for identifying whether an ICA component contains
# eye movement. (smaller values are stricter)
#%% Step 6: Load and merge datafiles
for j in range(len(files)):
    filename = datapath + files[j]
    file_obj = open(filename, "rb")
    polydata = Poly5Reader(filename)
    samples = polydata.samples
    sr = polydata.sample_rate
    channels = [c._Channel__name for c in polydata.channels]
    samples = samples.T
    if j == 0:
        data = samples
    else:
        data = np.concatenate((data, samples))

#%% Step 7: identify stimulus presentation times - this procedure varies with EEG
#   systems and whether you used analog or digital markers.
dval = 64
stimraw = np.array(np.where(data[:, stimchan] == dval)[0])
stimdiff = np.diff(stimraw)
stimind = np.where(stimdiff != 1)[0]
stimtime = np.append(stimraw[stimind], stimraw[-1])
ntrials = len(stimtime)  # this is how many trials we have.
nchan = data.shape[1]  # This is the total number of samples
#%% Step 8: filter EEG data.
#function that builds a filter 
def makefiltersos(sr,fp,fs,gp=3,gs=20):
	""" 	Wrapper function around scipy filter functions.  
	Makes it convenient by providing frequency parameters in terms of 
	frequencies in Hz.   
	INPUT: 	sr - sampling rate in Hz. 
		fp - pass frequency in Hz
		fs - stop frequency in Hz
		gp - pass band ripple in dB, default 3 dB
		gs - stop band attenuation in dB, default 20 dB

	OUTPUT: sos filter coefficients. 
			w,h for making bode plot 
	Automatically detects the type of filter.  if fp < fs the filter
	is low pass but if fp > fs the filter is highpass.  """
#
#set up filter parameters

	fn = sr/2
	wp = fp/fn
	ws = fs/fn
#get the filter order

	n,wn = signal.buttord(wp,ws,gp,gs);                                                            
#design the filter

#lowpass 
	if fp < fs:
		sos = signal.butter(n,wn,btype='lowpass',output='sos')
#highpass
	if fs < fp:
		sos = signal.butter(n,wn,btype='highpass',output='sos')
#get filter respons function	
	w,h = signal.sosfreqz(sos,fs=sr)
	return sos,w,h

# make a high pass digital filter
sos_hf, w, h = makefiltersos(sr, hipass, hipass / 2, gp=3, gs=12)
# make a low pass digital filter
sos_lf, w, h = makefiltersos(sr, lowpass, lowpass * 1.1, gp=3, gs=20)
data[:, 0:nEEGchan] = signal.sosfiltfilt(
    sos_lf, data[:, 0:nEEGchan], axis=0, padtype="odd"
)
data[:, 0:nEEGchan] = signal.sosfiltfilt(
    sos_hf, data[:, 0:nEEGchan], axis=0, padtype="odd"
)
#%% Step 9: segment filtered EEG by the stimulus onset.
# segment data
prestim = int(prestim * sr)  # convert prestim to samples
poststim = int(poststim * sr)  # convert poststim to samples
epochdata = np.zeros((ntrials, prestim + poststim + 1, nchan))
for j in range(len(stimtime)):
    epochdata[j, :, :] = data[stimtime[j] - prestim : stimtime[j] + poststim + 1, :]
if plot == "y":
    erp_noisy = np.mean(epochdata, axis=0)
    time = 1000 * np.arange(-prestim, poststim + 1) / sr
    plt.plot(time, erp_noisy[:, 0:nEEGchan])
    plt.title("ERP noisy")

#%% Step 10: Identify bad trials by standard deviation

eegstd = np.std(epochdata, axis=1)
chanstd = np.sum(eegstd[:, 0:nEEGchan], axis=0)
trialstd = np.sum(eegstd[:, 0:nEEGchan], axis=1)

if plot == "y":
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    ax[0].plot(chanstd / np.median(chanstd), "ro")
    ax[0].set_title("Channel Variability")
    ax[0].set_xlabel("Channels")
    ax[0].set_ylabel("Normalized Standard Deviation")
    ax[0].grid()
    ax[1].plot(trialstd / np.median(trialstd), "bo")
    ax[1].set_xlabel("Trials")
    ax[1].set_ylabel("Normalized Standard Deviation")
    ax[1].set_title("Trial Variability")
    ax[1].grid()
    plt.show(block=False)


# threshhold trials by standard deviation criteria and remove them.
# in principle you could threshold channels this way too.  But, I
# think with 32 channels you need to avoid that.  With 128 you could.

badtrials_eeg = np.where(trialstd / np.median(trialstd) > var_threshold)[0]
goodtrials = np.setdiff1d(range(ntrials), badtrials_eeg)
badchan_eeg = np.where(chanstd/np.median(chanstd) > chan_threshold)[0]
#%% Step 11: Detect no response trials and remove them.
responsesbytrial = np.sum(epochdata[:, prestim + 1 :, 40], axis=1)
badtrials_noresponse = np.where(responsesbytrial == 0)[0]

goodtrials = np.setdiff1d(goodtrials, badtrials_noresponse)

#%% Step 12: RUN ICA analysis to find and remove components correlated to eye blinks

# reshape gooddata for ICA

goodepochdata = epochdata[goodtrials, :, :]
goodepochdata = np.reshape(
    goodepochdata,
    (goodepochdata.shape[0] * goodepochdata.shape[1], goodepochdata.shape[2]),
)

# run ICA
ICA = FastICA(n_components=nEEGchan, whiten=True)
S = ICA.fit_transform(goodepochdata[:, 0:nEEGchan])
A = ICA.mixing_
W = ICA.components_
dubious_chans = np.unique(np.concatenate((np.array(eyechans),badchan_eeg)))
# compute correlations with eye channels
corrs = np.zeros((len(dubious_chans), nEEGchan))
for j in range(nEEGchan):
    for k in range(len(dubious_chans)):
        corrs[k, j] = np.corrcoef(S[:, j], goodepochdata[:, dubious_chans[k]])[0, 1]
corrsmax = np.max(corrs, axis=0)
if diagnostics == 'y':
    plt.plot(corrsmax)
    plt.title('Correlation Between ICA and Frontal Channels')
    plt.show()
# detect which components are not too correlated with eye channels
goodcomponents = abs(corrsmax) < corr_threshold
chancomponents = np.zeros(nEEGchan)
B = np.zeros(A.shape)
for j in range(nEEGchan):
    B[:,j] = A[:,j]**2/np.sum(A[:,j]**2)
    chancomponents[j] = np.max(B[:,j])
    if any(B[:,j] > 0.8):
        goodcomponents[j] = False

# recombine data without bad components.
AT = A.T
cleandata = S[:, goodcomponents] @ AT[goodcomponents, :]
# replace EEG data with clean data
goodepochdata[:, 0:nEEGchan] = cleandata
#%% Step 13: Save outputs to appropriate format.

# Do some data handling to restore trials to original order in file
goodepochdata = np.reshape(
    goodepochdata, (len(goodtrials), epochdata.shape[1], epochdata.shape[2])
)
finaldata = np.zeros(
    epochdata.shape
)  # An alternative would be to replace the data in epochdata
finaldata[goodtrials, :, :] = goodepochdata

if plot == "y":
    erp_clean = np.mean(goodepochdata, axis=0)
    time = 1000 * np.arange(-prestim, poststim + 1) / sr
    plt.plot(time, erp_clean[:, 0:nEEGchan])
    plt.title("ERP clean")

# fix channel names if needed
eeg_channels = [
    "Fp1",
    "Fpz",
    "Fp2",
    "F7",
    "F3",
    "Fz",
    "F4",
    "F8",
    "FC5",
    "FC1",
    "FC2",
    "FC6",
    "M1",
    "T7",
    "C3",
    "Cz",
    "C4",
    "T8",
    "M2",
    "CP5",
    "CP1",
    "CP2",
    "CP6",
    "P7",
    "P3",
    "Pz",
    "P4",
    "P8",
    "POz",
    "O1",
    "Oz",
    "O2",
]
channels[0:nEEGchan] = eeg_channels

# OUTPUT FOR UCI: create a dictionary for python use.

# data = dict()
# data['samplingrate'] = sr
# data['goodtrials'] = goodtrials
# data['badtrials_eeg'] = badtrials_eeg
# data['badtrials_noresponse'] = badtrials_noresponse
# data['prestim'] = prestim
# data['poststim'] = poststim
# data['data'] = finaldata
# data['stimulus'] = prestim
# data['channels'] = channels


# reshape for data frame for the feather file.
finaldata = np.reshape(finaldata, (ntrials * epochdata.shape[1], epochdata.shape[2]))

# create a stimulus channel for the feather file
stimulus = np.zeros(prestim + poststim + 1)
stimulus[prestim] = 1
stimulus = np.tile(stimulus, (1, ntrials))


#% make a new data frame for the feather file,

df_final = pd.DataFrame(finaldata, columns=channels)
df_final["stimulus"] = stimulus.T
feather.write_feather(df_final, datapath + files[0][:-6] + "_clean.ftr")

#%% Step 14: create a dataframe with trial level information.
goodtrialindex = np.ones(ntrials)
goodtrialindex[badtrials_eeg] = -1
goodtrialindex[badtrials_noresponse] = 0
df_trial = pd.DataFrame(goodtrialindex.T, columns=["goodtrialindex"])
feather.write_feather(df_trial, datapath + files[0][:-6] + "_trial.ftr")