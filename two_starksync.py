from sigelmsync import sigelmsync
from starkelmsync import starkelmsync
import matplotlib.cm as cm
import matplotlib.pylab as plt
import matplotlib.colors as colors
import numpy as np
from IPython import embed
from matplotlib import rc
rc('text', usetex=True)
import os

def smooth(data,N=5):
    return np.convolve(data,np.ones(N)/N,'same');

def smooth2d(data,N):
   for i in range(len(data[0,:])):
       data[:,i] = smooth(data[:,i],N)
   return data

############## PLOTS PAPER QUALITY IMAGE COMPARING HFS, JSAT AND SBD ELM-SYNCED DATA
def two_starksync(shotnr, ti, tf, preft=0.002, suft=0.008, refside = 'in', divside = 'in', divitem = 'jsat', elm_exper = "AUGD"):
    ######## Global settings ##############
    #shotnr = 30554
    #ti = 2.1
    #tf = 2.7
    #preft = 0.006
    #suft = 0.012

    mintime = -preft*1e3
    maxtime = suft*1e3

    ##Colorbar shrink
    cbshrink = 0.95
    ##Pad cb
    cbpad = 0.02

    #######################################
    fig = plt.figure(figsize=(5.0, 5.0), dpi=300)
    #fig = plt.figure(figsize=(3.0, 5.0), dpi=300)

    ticfont = 16
    font = {'family':'sans','size':ticfont}
    plt.rc('font',**font)
    plt.rc('legend',**{'fontsize':ticfont})

    ######################################
    ax1 = fig.add_subplot(211)

    ####Parameters for Refplot 
    if refside == 'in':
        refsgr = "HFSR"
        separatrix_signal = "Rin"
        minyval = 1.04
        maxyval = 1.17
    else:
        refsgr = "LFSR"
        separatrix_signal = "Raus"
        minyval = 2.08
        maxyval = 2.21

    th, rh = sigelmsync(shotnr, "RDL", refsgr, ti=ti, tf=tf, preft=preft, suft=suft, elm_exper=elm_exper)

    ##Don't use all data
    nchans = 11
    rh = rh[1:nchans]

    ## Check densities used (shortcut without dd)
    mindens = 0.5
    maxdens = 6.0
    dens = np.linspace(mindens, maxdens ,12)

    trin, rin = sigelmsync(shotnr, "FPG", separatrix_signal, ti=ti, tf=tf, preft=preft, suft=suft, elm_exper=elm_exper)

    indr = sorted(xrange(len(trin)), key=lambda ix: trin[ix])
    trinsort = trin[indr]
    rinsort = rin[indr]
    plt.fill_between(trinsort*1e3, rinsort, 1.5, color='black', alpha=0.25)

    ###Sep text
    #ax1.text(2.0, 1.14, r'$Separatrix$', color='#000000')

    #Plot Rin
    plt.scatter(trin*1e3, rin, color='k', edgecolors='k', s=1.0)

    #Plot the layers
    #colors = cm.rainbow(np.linspace(0, 1, len(rh)))

    fth = []
    frh = []
    fne = []

    for i in range(len(rh)):
        #    sc =plt.scatter(th*1e3, rh[i], c=colors[i], edgecolors=colors[i], s=1.0)
        #color = colors[i] * np.ones(len(th))
        fth.append(th*1e3)
        frh.append(rh[i])
        fne.append(dens[i] * np.ones(len(th)))

    ########
    fth = np.concatenate(fth)
    frh = np.concatenate(frh)
    fne = np.concatenate(fne)

    ########
    cmap = cm.get_cmap('jet')
    sc = plt.scatter(fth, frh, c=fne, s=6.0, lw=0, cmap=cmap)

    ###Colorbar
    cbticks = np.linspace(1, 5, 5)
    clb = plt.colorbar(sc, ticks = cbticks, shrink = cbshrink, pad = cbpad)
    clb.set_label("$\mathrm{n_{e}\,[10^{19}m^{-3}]}$", fontsize=ticfont)

    plt.ylim([minyval, maxyval])

    #Inner wall
    plt.hlines(1.045, mintime, maxtime, color='black', lw=3)
    
    plt.ylabel(r"$\mathrm{R\,[m]}$")

    plt.setp(ax1.get_xticklabels(), visible=False)
    ######################################
    
    ######################################
    #ax2 = fig.add_subplot(412, sharex=ax1)
    
    ##tp, sp = sigelmsync(shotnr, "XVS", "S2L1A10", ti=ti, tf=tf, preft=preft, suft=suft)#, file="/home/guimas/diags/30554/uic.lp")
    #tp, sp = sigelmsync(shotnr, "XVS", "S2L1A10", ti=ti, tf=tf, preft=preft, suft=suft, file="/home/guimas/diags/30554/uic.lp", offset=True)
    
    #plt.scatter(tp*1e3, sp, s=1.0, c="r", edgecolors="r")
    #plt.setp(ax2.get_xticklabels(), visible=False)
    ######################################

    ######################################
    ax4 = fig.add_subplot(212, sharex=ax1)

    t, s, m = starkelmsync(shotnr, ti=ti, tf=tf, item=divitem, side=divside, preft=preft, suft=suft, elm_exper=elm_exper)

    #Maximum vertical value for the emapping
    if divitem == 'nev':
        vmax = 80
        cbticks = np.linspace(20, vmax, 4)
        clblabelstr = "$\mathrm{n_{e,v}\,[10^{20}m^{-3}]}$"
    elif divitem == 'te':
        #vmax = np.floor(m.max())
        vmax = 20 #Limit it to 20eV
        cbticks = np.linspace(vmax/4, vmax, 4)
        clblabelstr = r"$\mathrm{t_{e}\,[eV]}$"
    elif divitem == 'net':
        vmax = np.min([24, np.floor(m.max())])
        cbticks = np.linspace(vmax/4, vmax, 4)
        clblabelstr = r"$\mathrm{n_{e,p}\,[10^{19}m^{-3}]}$"
    else: #Default to jsat
        vmax = 20
        cbticks = np.linspace(0, vmax, 5)
        clblabelstr = r"$\mathrm{\Gamma_{D^{+}}\,[10^{22}m^{-2}s^{-1}]}$"

    plt.pcolormesh(t*1e3, s, m, shading='goraud', vmax=vmax, cmap=cmap)
    ###Colorbar    
    clb = plt.colorbar(ticks = cbticks, shrink = cbshrink, pad = cbpad)
    clb.set_label(clblabelstr, fontsize=ticfont)

#    if divside=='in':
#        plt.ylim([min(s), max(s)])
#    else:
#        plt.ylim([-5, 20])

    plt.ylim([-5, 20])

    plt.ylabel(r"$\mathrm{\Delta S\,[cm]}$")
    plt.xlabel(r"$\mathrm{T-T_{ELM}\,[ms]}$")
    plt.hlines(0, mintime, maxtime, color='white')
    ######################################
    plt.subplots_adjust(left=0.16, bottom=0.12, right=0.96, top=0.98, wspace=0.10, hspace=0.05)
    #plt.subplots_adjust(left=0.28, bottom=0.12, right=0.84, top=0.98, wspace=0.10, hspace=0.05)

    plt.xlim([-preft*1e3, suft*1e3])
    #plt.xlim([mintime, maxtime])
    ######################################
    #plt.show()

    #Output file
    elmdir = os.path.expanduser('~') + "/ELM/"
    if not os.path.exists(elmdir):
        os.makedirs(elmdir)
    shotdir = elmdir + str(shotnr) + '/'
    if not os.path.exists(shotdir):
        os.makedirs(shotdir)

    filename = shotdir + "ss_" + str(shotnr) + "_" + str(ti) + "_" + str(tf) + "_" + refside + divside + ".png"
    plt.savefig(filename, format="png", dpi=300)
    print "Wrote: " + filename

if __name__=="__main__":
    sidestr = 'in'
    two_starksync(30554, 1.8, 2.8, preft=0.002, suft=0.008, refside=sidestr, divside=sidestr, elm_exper='guimas')
    sidestr = 'out'
    two_starksync(30554, 1.8, 2.8, preft=0.002, suft=0.008, refside=sidestr, divside=sidestr, elm_exper='guimas')

    sidestr = 'in'
    two_starksync(30554, 3.0, 3.5, preft=0.002, suft=0.008, refside=sidestr, divside=sidestr, elm_exper='guimas')
    sidestr = 'out'
    two_starksync(30554, 3.0, 3.5, preft=0.002, suft=0.008, refside=sidestr, divside=sidestr, elm_exper='guimas')

    #sidestr = 'in'
    #two_starksync(30554, 4.0, 4.5, preft=0.001, suft=0.003, refside=sidestr, divside=sidestr, elm_exper='guimas')
    #sidestr = 'out'
    #two_starksync(30554, 4.0, 4.5, preft=0.001, suft=0.003, refside=sidestr, divside=sidestr, elm_exper='guimas')
