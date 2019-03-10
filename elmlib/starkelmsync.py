#!/usr/bin/env python
#import dd
import dd
import numpy as np
from readStark import *

def starkelmsync(shotnr, tBegin=0.0, tEnd=10.0, item='nev', side='in',
                 preft=0.004, suft=0.008,
                 felm_min=0.0, felm_max=1000.0,
                 elm_exper='AUGD', elm_edition=0,
                 file=None):
    """Synchronizes DIVERTOR data to ELMs
    Data must be saved in a file that is output from DIVERTOR.
    
    Parameters
    ------------
    shotnr: int
        Number of the shot to analyse.       
    tBegin: float
        Initial time of analysis.
    tEnd: float
        Final time of analysis.
    item: str
    side: str
    preft: float
        'Prefix time' to consider before ELM onset.
    suft: float
        'Suffix time' to consider after ELM ends.
    felm_min: float
        Minimum ELM frequency to include in analysis.
        Default value 0Hz considers all ELMs.
    felm_max: float
        Maximum ELM frequency to include in analysis.
        Default value 1000Hz considers all ELMs.
    elm_exper: str
        User of the experiment. Default is public shotfile 'AUGD'.
    elm_edition: int
        Edition of the ELM shotfile. Default is latest edition, 0.
    file: str
        DIVERTOR data file. Default of 'None' will try to read a file stored in '$HOME/Divertor/#shotnr/'.
    
    Returns
    ------------
    synctime: np.array(float)
    deltas: np.array(float)
    syncmat: np.array(float, float)
    
    Example:
    import matplotlib.pyplot as plt
    t, s, m = starkelmsync(30554, ti=2.0, tf=3.0, item="jsat", side="in", suft=0.004, preft=0.002)
    plt.pcolormesh(t, s, m, shading='goraud')
    plt.show()
    """
    
    #Gets the filename where DIVERTOR data is stored
    if file is None:
        divertor_fname = getDivFname(shotnr, item, side)
    else:
        divertor_fname = file
    
    try:
        sdata = readDivData(divertor_fname)
        if (len(sdata.deltas) <= 2 | len(sdata.time) <= 2 ) :
            print "Returning dummy data from starkelmsync!"
            return sdata.time, sdata.deltas, sdata.data
    except:
        print "No such file!"
        raise

    analysis_mask = np.where((sdata.time >= tBegin) & (sdata.time <= tEnd))
    time_stark = sdata.time[analysis_mask]
    dummy = sdata.data[:, analysis_mask]
    ##MAS PORQUE CRLLLL!!!!!!???????
    mat = dummy[:, 0, :]

    ######## READ THE ELMS ########
    ELM = dd.shotfile("ELM", shotnr, experiment=elm_exper, edition=elm_edition)
    elmd = ELM("t_endELM", tBegin=tBegin, tEnd=tEnd)
    f_ELM = ELM('f_ELM')
    t_endELM = elmd.data
    t_begELM = elmd.time
    ELM.close()
    ###############################

    synctime = []
    dslen = len(sdata.deltas)
    syncmat = [[]]*dslen
    matT = mat.T

    for elm in range(t_begELM.size):
    #Only accept ELMs at the chosen frequency window
        if ((f_ELM.data[elm]>=felm_min) & (f_ELM.data[elm]<=felm_max)):
            t1, t2 =t_begELM[elm]-preft, t_endELM[elm]+suft
            #Re-adjust ELM times so no overlap between consecutive ELMs occurs
            if (elm >=1 ) :
                tendprev = t_endELM[elm-1]                
                t1 = np.max([t1,tendprev])
            if  (elm<t_begELM.size-1):
                tstartnext =  t_begELM[elm+1]
                t2 = np.min([t2,tstartnext])
            
            elmind = np.where((time_stark >= t_begELM[elm]-preft) & (time_stark <= t_endELM[elm]+suft))
            synctime.append(time_stark[elmind]-t_begELM[elm])

            for s in range(dslen):
                syncmat[s] = np.append(syncmat[s], [matT[elmind, s][0]])

        #dummy = np.append(mat[:, elmind], axis=0)
        #syncmat.append(dummy[:,0,:])

    synctime = np.concatenate(synctime)
    syncmat = np.array(syncmat)
    indt = np.argsort(synctime)

    return synctime[indt], sdata.deltas, syncmat[:, indt]


