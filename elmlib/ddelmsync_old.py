import dd
import numpy as np
from getsig import getsig

def ddelmsync_old(shotnr, diag, signal, edition=0,
              ti=0.0, tf=10.0, preft=0.001, suft=0.004,
              elm_exper="AUGD", elm_edition=0):
    """Gets a selected 1-D signal and syncs it to the ELMs in the desired time interval
    
    Parameters
    ------------
    shotnr: int
    diag: str
    signal: str
    edition: int
    ti: float
    tf: float
    preft: float
    suft: float
    elm_exper: str
    elm_edition: int
    
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
    elmd = ELM("t_endELM", tBegin=ti, tEnd=tf)
    t_endELM = elmd.data
    t_begELM = elmd.time
    ELM.close()
    ################################
    
    ################################
    ###### Get the signal in cause
    signal = getsig(shotnr, diag, signal, tBegin=ti, tEnd=tf, edition=edition)

    #################### Syncs the timebase to the ELM timebase     
    ###########################
    ###### Signal group
    ###########################
    syncsig = []#np.zeros_like(signal.data)
    synctime = []#np.zeros_like(signal.time)
        
    for elm in range(t_begELM.size):
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
        
    #Finally, return is again dependent on array dimensions
    if len(signal.data.shape)==1:
        synctime_return = np.concatenate(synctime)
        syncsig_return = np.concatenate(syncsig)
    if len(signal.data.shape)==2:
        synctime_return = np.concatenate(synctime)
        syncsig_return = np.concatenate(syncsig, axis=1)[0,:,:]
    
    return synctime_return, syncsig_return