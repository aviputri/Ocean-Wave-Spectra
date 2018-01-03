#!/bin/python2

import wave_util
import os
import sys
import numpy as np
import math

"""
our main process, which is called through the function at the bottom of this file. if you want to use the functions separately, ie just zerodown_data or synthetic, 'import main' in python interpreter and enter the functions manually.
"""
def main(argv):
    # if a data directory is specified from the command line, use it
    if len(argv) > 1:
        datpath = argv[1]
    # if not, look in the directory the program was launched from
    else:
        curpath = os.path.dirname(os.path.abspath(__file__))
        datpath = os.path.join(curpath, 'data')

    # if a sampling rate is specified use it, if not use a default value of 0.5
    if len(argv) > 2:
        sr = argv[2]
    else:
        sr = 2

    # zerodown all files in the data path
    dmean, dmax, dhts, dht10 = zerodown_data(datpath, sr)
    # write zerodown data to the specified files
    write('mean_data', 'max_data', 'hts_data', 'ht10_data', dmean, dmax, dhts, dht10)
    # generate synthetic data from specified measurement data file
    smean, smax, shts, sht10, data = synthetic('hts_data', 2, sr)
    # write zerodown data to the specified output files
    write('mean_synt', 'max_synt', 'hts_synt', 'ht10_synt', smean, smax, shts, sht10)


    #spectrm = perweekfft(datpath, 2.0)

    """
    sp = hugefft(datpath, sr)
    f = np.linspace(0,0.5,100)

    se = wave_util.genspec()

    import matplotlib.pyplot as plt
    plt.plot(f,sp,f,se)
    plt.show()
    """


"""
use this function to read measurement data from files in a folder, calculate various parameters of the individual waves, and write these parameters to file
"""
def zerodown_data(datpath, sr):
    # lists of relevant hourly data
    hmean = []
    hmax = []
    hhts = []
    hht10 = []

    # loop through all our data files and compute relevant hourly data
    """
    The weakness to this approach is that individual waves spanning files are not identified (they could but it would be a hassle). It also means we need to handle border cases (zero-dimension matrices) in our zerodown function.
    """
    for file in sorted(os.listdir(datpath)):
        # read raw data from file
        with open(os.path.join(datpath, file), 'r') as f:
            raw = f.read().split('\n')

        # remove first line (header) and any empty lines
        raw = [r for r in raw[1:] if r]
        # retrieve elevation data and convert into meters
        data = [r.split(',')[1].replace('-\.','-0\.') for r in raw]
        data = np.asarray(data, float)*0.01
        # zerodowncross
        mean, max, hts, ht10 = wave_util.zerodown(data, sr)
        # add to our hourly list
        print file, hts
        hmean.append(mean)
        hmax.append(max)
        hhts.append(hts)
        hht10.append(ht10)

    return (hmean, hmax, hhts, hht10)
    
"""
generates synthetic data from Hs.Ts pairs
"""
def synthetic(hsfile, spectrum, sr):
    curpath = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(curpath,hsfile),'r') as f:
        raw = f.read().split('\n')

    data = [[float(r.split(',')[0]),float(r.split(',')[1])] for r in raw if r]

    """
    constants for spectrum and sea water level generation. we define them outside the loop to save some time
    """
    max_t = 3600                # maximum duration of time series
    nyq_f = sr/2                # nyquist frequency of spectra: highest frequency bin
    del_f = 1.0/max_t           # interval between frequency bins
    f = np.linspace(del_f,nyq_f,max_t)    # array containing all frequency bins
    t = np.linspace(0,max_t,max_t*sr)
    #t = np.reshape(t, (t.size,1))

    # lists of relevant hourly data
    hmean = []
    hmax = []
    hhts = []
    hht10 = []

    index = 0

    # jonswap
    if spectrum == 0:
        for dat in data:
            index += 1
            # generate spectrum
            S = wave_util.jonswap(dat[0], dat[1], 3600, f)
            # generate elevation data
            E = wave_util.ema(S, 0.5, del_f)
            # zerodowncross
            mean, max, hts, ht10 = wave_util.zerodown(E, sr)
            print index, [dat[0] - hts[0], dat[1] - hts[1]]
            # add to our hourly list
            #print file, mean, max, hts, ht10
            hmean.append(mean)
            hmax.append(max)
            hhts.append(hts)
            hht10.append(ht10)

    # wallops
    elif spectrum == 1:
        for dat in data:
            index += 1
            # generate spectrum
            S = wave_util.wallops(dat[0], dat[1], 10.0, f)
            # generate elevation data
            E = wave_util.ema(S, 0.5, del_f)
            # zerodowncross
            mean, max, hts, ht10 = wave_util.zerodown(E, sr)
            print index, [dat[0] - hts[0], dat[1] - hts[1]]
            # add to our hourly list
            #print file, mean, max, hts, ht10
            hmean.append(mean)
            hmax.append(max)
            hhts.append(hts)
            hht10.append(ht10)

    # custom spectrum -- based on bert-schneider
    elif spectrum == 2:
        for dat in data:
            index += 1
            #generate spectrum
            S = wave_util.genspec(hs=dat[0], ts=dat[1]-1., f=f, foffset=0, A=0.072, C=0)

            E = wave_util.ema(S, 0.5, del_f)
            mean, max, hts, ht10 = wave_util.zerodown(E, sr)
            print index, [dat[0] - hts[0], dat[1] - hts[1]]
            # add to our hourly list
            #print file, mean, max, hts, ht10
            hmean.append(mean)
            hmax.append(max)
            hhts.append(hts)
            hht10.append(ht10)

    data[:] = np.asarray(data)
    hhts[:] = np.asarray(hhts)
    print np.mean(data), np.mean(hhts)

    return (hmean, hmax, hhts, hht10, data)

"""
generates synthetic data from Hs.Ts pairs
"""
def maximizem(sr):
    curpath = os.path.dirname(os.path.abspath(__file__))

    with open(os.path.join(curpath,'hts_data'),'r') as f:
        hts_raw = f.read().split('\n')
    dhts = [[float(r.split(',')[0]),float(r.split(',')[1])] for r in hts_raw if r]

    with open(os.path.join(curpath,'mean_data'),'r') as f:
        mean_raw = f.read().split('\n')
    dmean = [[float(r.split(',')[0]),float(r.split(',')[1])] for r in mean_raw if r]

    with open(os.path.join(curpath,'max_data'),'r') as f:
        max_raw = f.read().split('\n')
    dmax = [[float(r.split(',')[0]),float(r.split(',')[1])] for r in max_raw if r]

    with open(os.path.join(curpath,'ht10_data'),'r') as f:
        ht10_raw = f.read().split('\n')
    dht10 = [[float(r.split(',')[0]),float(r.split(',')[1])] for r in ht10_raw if r]

    """
    constants for spectrum and sea water level generation. we define them outside the loop to save some time
    """
    max_t = 3600                # maximum duration of time series
    nyq_f = sr/2                # nyquist frequency of spectra: highest frequency bin
    del_f = 1.0/max_t           # interval between frequency bins
    f = np.linspace(del_f,nyq_f,max_t)    # array containing all frequency bins
    t = np.linspace(0,max_t,max_t*sr)
    #t = np.reshape(t, (t.size,1))

    # list of hourly error rates
    emean = []
    emax = []
    ehts = []
    eht10 = []

    for m in xrange(2,10):
        # lists of relevant hourly data
        hmean = []
        hmax = []
        hhts = []
        hht10 = []

        index = 0
        for dat in dhts:
            index += 1
            # generate spectrum
            S = wave_util.wallops(dat[0], dat[1], m, f)
            # generate elevation data
            E = wave_util.ema(S, 0.5, del_f)
            # zerodowncross
            mean, max, hts, ht10 = wave_util.zerodown(E, sr)
            print index, [dat[0] - hts[0], dat[1] - hts[1]]
            # add to our hourly list
            #print file, mean, max, hts, ht10
            hmean.append(mean)
            hmax.append(max)
            hhts.append(hts)
            hht10.append(ht10)

        hhts[:] = np.asarray(hhts)
        hmean[:] = np.asarray(hmean)
        hmax[:] = np.asarray(hmax)
        hht10[:] = np.asarray(hht10)

        meandmean = np.mean(np.asarray(dmean))
        meandmax = np.mean(np.asarray(dmax))
        meandhts = np.mean(np.asarray(dhts))
        meandht10 = np.mean(np.asarray(dht10))

        emean.append((np.mean(hmean) - meandmean) / meandmean)
        emax.append((np.mean(hmax) - meandmax) / meandmax)
        ehts.append((np.mean(hhts) - meandhts) / meandhts)
        eht10.append((np.mean(hht10) - meandht10) / meandht10)
            
    return emean, emax, ehts, eht10

def perweekfft(datpath, sr):
    files = sorted(os.listdir(datpath))
    allfft = []
    for file in files[:2184]:
        with open(os.path.join(datpath, file), 'r') as f:
            raw = f.read().split('\n')

        # remove first line (header) and any empty lines
        raw = [r for r in raw[1:] if r]
        # retrieve elevation data and convert into meters
        data = [r.split(',')[1].replace('-\.','-0\.') for r in raw]
        data = np.asarray(data, float)*0.01
        data -= np.mean(data)
        allfft.append(np.abs(np.fft.rfft(data)))

    wfft = np.mean(np.vstack(allfft).reshape((13, 168, 3601)), axis=1)

    for n in xrange(wfft.shape[0]):
        with open('week_' + str(n) + '.sp', 'w') as f:
            [f.write(str(line) + '\n') for line in wfft[n]]

    return wfft

def hugefft(datpath, sr, smoothing_factor=100):
    import matplotlib.pyplot as plt

    files = sorted(os.listdir(datpath))
    allfft = []
    for file in files[:2184]:
        with open(os.path.join(datpath, file), 'r') as f:
            raw = f.read().split('\n')

        # remove first line (header) and any empty lines
        raw = [r for r in raw[1:] if r]
        # retrieve elevation data and convert into meters
        data = [r.split(',')[1].replace('-\.','-0\.') for r in raw]
        data = np.asarray(data, float)*0.01
        allfft.append(np.abs(np.fft.rfft(data)))
        
    avgsp = np.mean(np.vstack(allfft), axis=0)
    sp = avgsp[1:]
    sp = sp.reshape((smoothing_factor, avgsp.size/smoothing_factor))
    sp = np.mean(sp, axis=-1)

    #spectrum = np.fft.rfft(spec)
    #spectrum.tofile(os.join.path(datpath, 'full-spectrum'), format="text")

    return sp

def write(fmean, fmax, fhts, fht10, hmean, hmax, hhts, hht10):
    # write results to file
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fmean), 'w') as f:
        [f.write(str(mean[0]) + ',' + str(mean[1]) + '\n') for mean in hmean]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fmax), 'w') as f:
        [f.write(str(max[0]) + ',' + str(max[1]) + '\n') for max in hmax]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fhts), 'w') as f:
        [f.write(str(hts[0]) + ',' + str(hts[1]) + '\n') for hts in hhts]
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fht10), 'w') as f:
        [f.write(str(ht10[0]) + ',' + str(ht10[1]) + '\n') for ht10 in hht10]
    
"""
this is run first
"""
if __name__ == ('__main__'):
    main(sys.argv)
