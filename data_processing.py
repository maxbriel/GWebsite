#
# Author: Max Briel
#
# Several errors will be thrown, but these can be savely ignored.
# This code contains code from the LIGO collaborations
# QuickNotebook and Tutorial IPython notebooks.

# packages for plotting
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm, rc

# packages for reading files
import h5py
import healpy
import os

# Data analysis packages
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
from scipy.io import wavfile
from datetime import datetime, timedelta

# LIGO-specific readligo.py
import readligo as rl

# GW event times in GPS time
times = {
"GW150914":1126259462.44,
"LVT151012":1128678900.44,
"GW151226":1135136350.65,
"GW170104":1167559936.6,
"GW170608":1180922494.49,
"GW170814":1186741861.53,
"GW170817":1187008882.43
}

# CODE FROM LIGO Open Science Center QuickNotebook!
# -----------------------------------------------------------------------------
# function to writen data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back,
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

# function to keep the data within integer limits, and write to wavfile:
def write_wavfile(filename,fs,data):
    d = np.int16(data/np.max(np.abs(data)) * 32767 * 0.9)
    wavfile.write(filename,int(fs), d)

# function that shifts frequency of a band-passed signal
def reqshift(data,fshift=100,sample_rate=4096):
    """Frequency shift the signal by constant
    """
    x = np.fft.rfft(data)
    T = len(data)/float(sample_rate)
    df = 1.0/T
    nbins = int(fshift/df)
    # print T,df,nbins,x.real.shape
    y = np.roll(x.real,nbins) + 1j*np.roll(x.imag,nbins)
    y[0:nbins]=0.
    z = np.fft.irfft(y)
    return z

# -----------------------------------------------------------------------------

# Get the sky fits files
cwd = os.getcwd()
sky_fits = os.listdir(cwd+"/FITS")

# Remove Apple DS_Store from file list.
if ".DS_Store" in sky_fits:
    sky_fits.remove(".DS_Store")

# make a new folder for the images if it's not already present
if "events" not in os.listdir(cwd):
    os.mkdir("events")

# location for the output
storage = cwd+"/events/"

#setup color maps for the skymaps
color_map = cm.plasma
color_map.set_under("w", alpha=0)
color_map.set_bad("w", alpha=0)

# export the graticule if non-present
if "graticule.png" not in os.listdir(storage):
    # take first fit and empty it to plot a same size figure.
    map, header = healpy.read_map("FITS/"+sky_fits[0], h=True, verbose=False)
    for i in range(len(map)):
        map[i] = np.nan
    output = healpy.mollview(map, flip="geo", cbar=False, title="",  return_projected_map=True, cmap=color_map)
    healpy.graticule(dpar=20)
    plt.savefig(storage+"/graticule.png",figsize= (1000,500), dpi=300, transparent=True, bbox_inches='tight')
else:
    print("Graticule already present. Remove from file system to replot")

# Export the plots of the events
for i in sky_fits:
    filename = "FITS/"+i
    event_name = i[:-5]
    print("Working on event: "+event_name)
    # create an event folder in not already present
    if event_name not in os.listdir(storage):
        os.mkdir(storage+event_name)
    event_folder = storage+event_name+"/"

    if event_name+".png" not in os.listdir(event_folder):
        map, header = healpy.read_map(filename, h=True, verbose=False)
        for j in range(len(map)):   # remove low probability values
            if map[j] < 1e-5:
                map[j] = np.nan
        # Mark the maximum probability
        maximum = np.nanargmax(map)
        angle = healpy.pix2ang(healpy.get_nside(map), maximum)

        # plot the remaining probability and maximum
        output = healpy.mollview(map, flip="geo",cbar=False,title="",
                                 return_projected_map=True, cmap=color_map)

        healpy.projplot(angle, marker='*', color="#0091A5", markersize=18,
                        markeredgecolor="#003C45", markeredgewidth=0.5)

        plt.savefig(event_folder+event_name+".png",figsize= (1000,500),
                    dpi=300, transparent=True, bbox_inches='tight',  coord="G",notext=True)
        plt.close()
    else:
        print("Location information already present. Remove first to replot")


    t0 = times[event_name]
    # Get the current time from GPS time to UTC
    # in 1980: 19 leapseconds, in 2015: 36 leapseconds, 2017: 37
    if event_name[-6:-4] == "15":
        utc = datetime(1980, 1, 6) + timedelta(seconds=t0 - (36 - 19))
    else:
        utc = datetime(1980, 1, 6) + timedelta(seconds=t0 - (37 - 19))

    detector = "H1"
    data_file = "data/"+event_name+".hdf5"
    strain, time, chan_dict_H1 = rl.loaddata(data_file, 'H1')
    dt = time[1] - time[0]
    fs = int(np.round(1/dt))
    rel_time = time - t0
    #-- How much data to use for the ASD?
    deltat = 15  # Number of seconds on each side of data
    N_samp = deltat*fs

    # -- Center the ASD segment on the requested time
    indx = np.where(np.abs(rel_time) < dt)[0][0]

    strain_seg = strain[indx-N_samp : indx+N_samp]
    time_seg = rel_time[indx-N_samp : indx+N_samp]

    # number of sample for the fast fourier transform:
    NFFT = 1*fs
    fmin = 10
    fmax = 2000
    Pxx, freqs = mlab.psd(strain_seg, Fs = fs, NFFT=NFFT,
                          noverlap=NFFT/2, window=np.blackman(NFFT))

    # We will use interpolations of the ASDs computed above for whitening:
    psd = interp1d(freqs, Pxx)

    # Plot the ASD
    fig1 = plt.figure(figsize=(10,5))
    plt.loglog(freqs, np.sqrt(Pxx),'r',label='{0} strain'.format(detector),
                color="#D75749", linewidth=4)
    plt.axis([fmin, fmax, 1e-24, 1e-19])
    plt.xticks([10,1e2,1e3],[1,2,3], fontsize=15)
    plt.yticks([1e-19,1e-21, 1e-23], [-19,-21,-23], fontsize=15)
    plt.title("ASD", fontsize=20)

    if "ASD.png" not in os.listdir(event_folder):
        plt.savefig(event_folder+"ASD.png", dpi=300, transparent=True, bbox_inches='tight')
    else:
        print("ASD already present. Remove first to replot")
    plt.close(fig1)

    # now whiten the data
    strain_whiten = whiten(strain_seg,psd,dt)

    # We need to suppress the high frequencies with some bandpassing:
    high_freq = 600.
    low_freq  = 30.
    bb, ab = butter(4, [low_freq*2./fs, high_freq*2./fs], btype='band')
    strain_whitenbp = filtfilt(bb, ab, strain_whiten)

    # Plot the Waveform
    fig2 = plt.figure(figsize=(10,5))
    plt.plot(time_seg,strain_whitenbp,'r',label='H1 strain', linewidth=4,
                color="#D75749")
    plt.xlim([-0.15,0.05])
    plt.ylim([-5,5])
    plt.yticks([-5,0,5], fontsize=15)
    plt.xticks([-0.1,-0.05,0], fontsize=15)
    plt.plot([0,0], [-5,5], color='black')
    plt.title('Waveform at {:%H:%M:%S} UTC'.format(utc), fontsize=20)

    if "waveform.png" not in os.listdir(event_folder):
        plt.savefig(event_folder+"waveform.png", dpi=300, transparent=True, bbox_inches='tight')
    else:
        print("Waveform already present. Remove first to replot")
    plt.close(fig2)

    # pick a shorter FTT time interval, like 1/16 of a second:
    NFFT = int(fs/16.)
    # and with a lot of overlap, to resolve short-time features:
    NOVL = int(NFFT*15./16)
    # and choose a window that minimizes "spectral leakage"
    # (https://en.wikipedia.org/wiki/Spectral_leakage)
    window = np.blackman(NFFT)

    spec_cmap='viridis'

    # Calculate the whitened spectrogram
    plt.figure()
    spec_H1, freqs, bins, im = plt.specgram(strain_whiten, NFFT=NFFT, Fs=fs, window=window,
      noverlap=NOVL, xextent=[-deltat,deltat], cmap=spec_cmap, vmin=0, vmax=0.01, scale='linear')

    # Plot the Spectrogram
    fig3, ax = plt.subplots(figsize=(10,6))
    Z = np.flipud(spec_H1) / np.median(spec_H1)
    extent = -deltat, deltat, freqs[0], freqs[-1]
    im = ax.imshow(Z,  cmap=spec_cmap, extent=extent, vmin=0, vmax=25)

    ax.axis('auto')
    plt.axis([-0.5, 0.5, 0, 512])
    plt.xticks([-0.5,0,0.5],[0,0.5,1],fontsize=20)
    plt.yticks([256,512],[256,512], fontsize=20)

    if "spectrogram.png" not in os.listdir(event_folder):
        plt.savefig(event_folder+"spectrogram.png", dpi=300, transparent=True, bbox_inches='tight')
    else:
        print("Spectrogram already present. Remove first to replot")
    plt.close(fig3)
    deltat_sound = 2.                     # seconds around the event
    # index into the strain time series for this time interval:
    indxd = np.where((time_seg >= -deltat_sound) & (time_seg < deltat_sound))
    # parameters for frequency shift
    fshift = 400.
    speedup = 1.
    fss = int(float(fs)*float(speedup))

    # shift frequency of the data
    strain_H1_shifted = reqshift(strain_whitenbp,fshift=fshift,sample_rate=fs)
    # write the files:
    if "sound.wav" not in os.listdir(event_folder):
        write_wavfile(event_folder+"sound.wav",int(fs), strain_H1_shifted[indxd])
    else:
        print("Soundfile already present. Remove first to remake")
