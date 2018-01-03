#!/bin/python2
import numpy as np
import numpy.random as nprand
import numpy.fft as npf
import numexpr as ne
import math

"""
zerodown takes elevation data and its sampling ratio and computes values for individual wave height and period sorted descending by wave height. In addition it returns values for significant(1/3) and 1/10 individual wave height/period.
"""
def zerodown(elevation, sampling_ratio):
    # make no assumptions regarding the unit of time
    srd = 1.0/sampling_ratio
    time_length = len(elevation) * srd

    # 
    i_prev = 0

    #
    t = 0.0
    t_prev = 0.0

    ht = []

    # lagrange polynomial dimension
    n = np.arange(2)

    # normalize elevation
    elevation -= np.mean(elevation)

    # determine height and period of each individual wave
    for i in xrange(1,len(elevation)-1):
        # if elevation at current index is positive and elevation at the next index is negative
        if elevation.item(i) >= 0 and elevation.item(i+1) < 0:
            # calculate the interpolated time of zero crossing
            t = lagrpol(elevation[i:i+2], np.arange(i,i+2), n, 0.0)
            # calculate period of current individual wave
            T = srd * (t - t_prev)

            if T <= 30:
                # calculate height of current individual wave
                H = np.ptp(elevation[i_prev:i])

                # append to our ouput array
                ht.append([H, T])

            t_prev = t
            i_prev = i

    # sort descending and convert to numpy array
    ht.sort()
    ht[:] = ht[::-1]
    htsort = np.asarray(ht,float)

    # obtain mean, hs,ts and h10,t10 through running average pass
    n = len(htsort)
    n10 = n/10
    ns = n/3

    ht10 = [np.mean(htsort[:n10,0]), np.mean(htsort[n10:,1])]
    hts = [np.mean(htsort[:ns,0]), np.mean(htsort[:ns,1])]
    mean = [np.mean(htsort[:,0]), np.mean(htsort[:,1])]
    tmax = np.max(htsort[:,1])
    max = [htsort[0,0], tmax]

    return (mean, max, hts, ht10)


"""
lagrange interpolation
"""
def lagrpol(x, y, n, xx):
    fx = 0.0

    for i in n:
        l = 1.0
        for j in n:
            if j != i:
                l *= (xx - x.item(j)) / (x.item(i) - x.item(j))
        fx += l * y.item(i)

    return fx


"""
generate jonswap spectrum
"""
def jonswap(hs, ts, max_t, f):
    sp = []     # array for our output spectrum

    # spectrum generator constants
    gamma_1 = 3.3
    beta_i = ne.evaluate('0.0624*(1.094 - 0.01915*log(gamma_1))/(0.230 + 0.0336*gamma_1 - (0.185/(1.9 + gamma_1)))')
    divtstp = ne.evaluate('1 - 0.132*(gamma_1 + 0.2)**(-0.559)')
    divltfp = 2*0.07**2
    divgtfp = 2*0.09**2

    tp = ts/divtstp
    hssq = hs**2
    tppow4 = tp**(-4)

    # calculate fp (cutoff frequency) and figure out which frequency index it's closest to
    fp = int(max_t/tp)
    # split freqeuncy array into two separate subarrays at fp index
    flt, fgt = np.split(f,[fp])

    # calculate spectrum for each frequency subarray and combine them back: for frequencies less than fp (flt) use divltfp, otherwise (fgt) use divgtfp. in addition, use numexpr to speed things up
    sp.extend(ne.evaluate('beta_i*hssq*tppow4*flt**(-5)*exp(-1.25*(flt*tp)**(-4))*gamma_1**(exp(-((flt*tp - 1)**2)/divltfp))'))
    sp.extend(ne.evaluate('beta_i*hssq*tppow4*fgt**(-5)*exp(-1.25*(fgt*tp)**(-4))*gamma_1**(exp(-((fgt*tp - 1)**2)/divgtfp))'))

    return sp


"""
sea water elevation generation from spectrum
"""
def ema(sp, sr, del_f):
    pi2 = 2*math.pi
    srd = 1.0/sr

    # vectorize spectrum and weigh
    a = (np.asarray(sp,float)*2*del_f)**0.5

    # use inverse FFT to generate random phase time domain information from spectrum
    eta = npf.irfft(a * np.exp(1j * nprand.uniform(0, pi2, a.shape))) * a.size

    return eta

"""
generate wallops spectrum
"""
def wallops(hs, ts, m, f):
    sp = [] # initialize output array for spectrum

    # spectrum generator constants
    gamma = math.gamma(0.25*(m-1))
    tp = ne.evaluate('ts / (1 - 0.283*(m - 1.5)**(-0.684))')
    bw = ne.evaluate('0.0624*m**(0.25*(m-1)) * (1 - 0.748*(m + 2)**(-1.057)) / (4.0**(0.25*(m-5)) * gamma)')
    sp.extend(ne.evaluate('bw*hs**2*tp**(1-m)*f**(-m)*exp(-0.25*m*(tp*f)**(-4))'))

    return sp

def genspec(hs=1.5682, ts=3.8296, f=np.linspace(0,0.5,100), foffset=0.055, A=1.3, B=-0.25, C=14):
    sp = A*hs**2*ts**(-4)*(f+foffset)**(-5)*np.exp(B*(ts*(f+foffset))**(-4))+C
    return sp
