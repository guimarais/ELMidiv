import dd
import numpy as np
from getsig import getsig

def ddelmsync(shotnr, diag, signal, experiment='AUGD', edition=0,
              tBegin=0.0, tEnd=10.0, preft=0.001, suft=0.004,
              felm_min=0, felm_max=1000,
              elm_exper="AUGD", elm_edition=0):
    """Gets a selected 1-D signal and syncs it to the ELMs in the desired time interval
    
    Parameters
    ------------
    shotnr: int
        Number of the shot to analyse.
    diag: str
        Name of the shotfile containing the data.
    signal: str
        Name of the signal in 'diag' to analyse.
    experiment: str
        Naame of the experiment (username) containing the shotfile.
    edition: int
        Edition of the shotfile.
    tBegin: float
        Initial time of analysis.
    tEnd: float
        Final time for analysis.
    preft: float
        'Prefix time' to consider before ELM onset.
    suft: float
        'Suffix time' to consider after ELM ends.
    felm_min: float
        Minimum ELM frequency to include in analysis.
        Default value 0Hz considers all ELMs.        
    felm_max : float
        Maximum ELM frequency to include in analysis.
        Default value 1000Hz considers all ELMs.
    elm_exper: str
        User of the experiment. Default is public shotfile 'AUGD'.
    elm_edition: int
        Edition of the ELM shotfile. Default is latest edition, 0.
    
    Returns
    ------------
    synctime: np.array(float)
        Array of times resetted to the closest ELM.
    syncsig: np.array(float)    
        1D or 2D array of data ordered according to the sychronized times.
    """
    #################################
    ##Flag to check if data came from a signal group
    sgrp = False

    ###### Gets ELM data ############
    ELM = dd.shotfile("ELM", shotnr, experiment=elm_exper, edition=elm_edition)
    elmd = ELM("t_endELM", tBegin=tBegin, tEnd=tEnd)
    freq_ELM = ELM("f_ELM", tBegin=tBegin, tEnd=tEnd)
    t_endELM = elmd.data
    t_begELM = elmd.time
    ELM.close()
    ################################
    
    ################################
    ###### Get the signal in cause
    signal = getsig(shotnr, diag, signal, tBegin=tBegin, tEnd=tEnd, edition=edition)

    #################### Syncs the timebase to the ELM timebase     
    ###########################
    ###### Signal group
    ###########################
    syncsig = []#np.zeros_like(signal.data)
    synctime = []#np.zeros_like(signal.time)
        
    for elm in range(t_begELM.size):
        #Only accept ELMs at the chosen frequency window
        if (freq_ELM.data[elm]>=felm_min)&(freq_ELM.data[elm]<=felm_max):
            t1, t2 =t_begELM[elm]-preft, t_endELM[elm]+suft
            #Re-adjust ELM times so no overlap between consecutive ELMs occurs
            if (elm >=1 ) :
                tendprev = t_endELM[elm-1]                
                t1 = np.max([t1,tendprev])
            if  (elm<t_begELM.size-1):
                tstartnext =  t_begELM[elm+1]
                t2 = np.min([t2,tstartnext])

            elmind = np.where((signal.time >= t1) & (signal.time <=t2))
        
            synctime.append(signal.time[elmind]-t_begELM[elm])
        
            #Distinguish between 1D (signal) and 2D array (Signal group)
            if len(signal.data.shape)==1:
                syncsig.append(signal.data[elmind])
            elif len(signal.data.shape)==2:
                syncsig.append(signal.data[elmind,:])
            else:
                raise Exception('Array format not supported!')
        else:#Space left open for possible analysis
            a=0
            
    #Finally, return is again dependent on array dimensions
    if len(signal.data.shape)==1:
        synctime_return = np.concatenate(synctime)
        syncsig_return = np.concatenate(syncsig)
    if len(signal.data.shape)==2:
        synctime_return = np.concatenate(synctime)
        syncsig_return = np.concatenate(syncsig, axis=1)[0,:,:]
    
    return synctime_return, syncsig_return