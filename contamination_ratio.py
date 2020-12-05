# -*- coding: utf-8 -*-

import numpy as np

    

def contamination_ratio(spike_train, sampling_rate, total_duration, 
                        refractory_time=1e-3, censored_time=.5e-3): 
    """
    Return ratio of contamination

    Parameters
    ----------
    spike_train : list,
        List of timing.
    sampling_rate : int,
        Sampling rate of the recording, if list_sample is is s set to 1.
    total_duration: float,
        Total lenght of the recording, same unit as spike_train.
    refractory_time : float,
        Duration of the refractory period, in second, Default is 1e-3.
    censored_time : float,
        Duration of the cesored period, in second. The min distance between two 
        spikes found with the spikesorter. Default is .5e-3


    Returns
    -------
    ratio_contamination : float
        Ratio of contamination.

    """
    
    assert refractory_time > censored_time, "censored period is higher than refractory period"
    
    refractory_sample = int(refractory_time * sampling_rate) # convert in sample
    censored_sample = int(censored_time * sampling_rate) # convert in sample
    Tf = refractory_sample - censored_sample # Duration on which contaminant spike are look after
    N = len(spike_train)*1. # size of the cluster

    isi = np.diff(spike_train)
    N_conta = np.count_nonzero((isi <= refractory_sample) *\
                               (isi >= censored_sample)) #count event between censored time and refractory time

    freq = N / (total_duration/sampling_rate) # freq of the unit
    
    ratio_contamination = 1 - (1 - N_conta / (Tf*N*freq))**0.5
    
    return ratio_contamination
