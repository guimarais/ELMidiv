import dd
import numpy as np
from getsig import getsig

def ddelmsync(shotnr, diag, signal, edition=0,
              ti=0.0, tf=10.0, preft=0.001, suft=0.004,
              elm_exper="AUGD", elm_edition=0):
    #################################
    ### Time Base ELM Sync:
    ### Gets a selected 1-D signal and syncs it to the ELMs in the desired time interval
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
    #### Open with dd libraries ####
    ################################

    signal = getsig(shotnr, diag, signal, tBegin=ti, tEnd=tf, edition=edition)
    # Remove Offset (bolometers, etc)

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
        #synctime.append(signal.time[elmind]-signal.time[elm])
        #syncsig.append(signal.data[elmind])
        #plt.scatter(signal.time[elmind]-t_begELM[elm], signal.data[elmind], color='k')
        
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
        return np.concatenate(synctime), np.concatenate(syncsig)
    if len(signal.data.shape)==2:
        return np.concatenate(synctime), np.concatenate(syncsig, axis=1)[0,:,:]