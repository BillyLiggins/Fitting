import copy
import echidna
import echidna.output.plot as plot
import echidna.core.spectra as spectra
from echidna.output import store
import matplotlib.pyplot as plt

import argparse
import glob
import numpy as np
import os

def convertor(path):
    flist=np.array(glob.glob(path))
    for ntuple in flist:
        os.system("python ~/echidna/echidna/scripts/dump_spectra_ntuple.py -c ~/workspace/PhD/fitting/config.yml -f "+ str(ntuple)+" -s hdf5/")
    
def combinerNtuple(path,filename):
    flist=np.array(glob.glob(path))
    print flist
    first = True
    for hdf5 in flist:
        print hdf5
        if first:
            spectrum1 = store.fill_from_ntuple(hdf5)
            first = False
        else:
            spectrum2 = store.fill_from_ntuple(hdf5)
            spectrum1.add(spectrum2)
    store.dump(filename, spectrum1)



def combiner(path,filename):
    flist=np.array(glob.glob(path))
    print flist
    first = True
    for hdf5 in flist:
        print hdf5
        if first:
            spectrum1 = store.load(hdf5)
            first = False
        else:
            spectrum2 = store.load(hdf5)
            spectrum1.add(spectrum2)
    store.dump(filename, spectrum1)


"""The way you should do it is to define a lot of spectra and then plot them.
    You don't really know how to normlise the histrogram or indeed weather that is of any uses in the first
    place. 
"""

def slicer(spectrumPath,filler,nslice):
    for i in range(nslice):
       spectrum=store.load(spectrumPath)
       print spectrum.sum()
       shrink_dict = {"energy_reco_low": 0.,
                      "energy_reco_high": 0.6,
                      "radial_reco_low": i*6000.0/nslice,
                      "radial_reco_high": (i+1)*6000/nslice}
       spectrum.cut(**shrink_dict)
       spectrum.scale(1)
       spec2=copy.copy(spectrum)
       spec2._name=str(i*1000)+"mm to "+str((i+1)*1000)+"mm"
       print type(spec2)
       filler.append(spec2)

def slicerMC(spectrumPath,filler,nslice):
    for i in range(nslice):
       spectrum=store.load(spectrumPath)
       print spectrum.sum()
       shrink_dict = {"energy_mc_low": 0.,
                      "energy_mc_high": 1,
                      "radial_mc_low": i*6000.0/nslice,
                      "radial_mc_high": (i+1)*6000/nslice}
       spectrum.cut(**shrink_dict)
       spectrum.scale(1)
       spec2=copy.copy(spectrum)
       spec2._name="MC"
       print type(spec2)
       print "This gives the number os events in each window:"
       print "mc : "+str(i*6000.0/nslice)+"mm to "+str((i+1)*6000.0/nslice)+"mm : "+str(spec2.sum())
       filler.append(spec2)

def slicerReco(spectrumPath,filler,nslice):
    for i in range(nslice):
       spectrum=store.load(spectrumPath)
       print spectrum.sum()
       shrink_dict = {"energy_reco_low": 0.,
                      "energy_reco_high": 1.,
                      "radial_reco_low": i*6000.0/nslice,
                      "radial_reco_high": (i+1)*6000/nslice}
       spectrum.cut(**shrink_dict)
       spectrum.scale(1)
       spec2=copy.copy(spectrum)
       spec2._name="Reco"
       print type(spec2)
       print "This gives the number os events in each window:"
       print "reco : "+str(i*6000.0/nslice)+"mm to "+str((i+1)*6000.0/nslice)+"mm : "+str(spec2.sum())
       filler.append(spec2)


def signalPlotter(spectra,dim,name):
    i=0
    for spec in spectra:
        fig = plt.figure()
        ax= fig.add_subplot(1,1,1)

        par = spec.get_config().get_par(dim)
        width = par.get_width()
        bins = np.linspace(par._low,par._high, par._bins+1)
        x = bins[:-1] + 0.5*width

        plt.xlabel(str(dim)+ " [" + par.get_unit() + "]")
        plt.ylabel("Events per " + str(width) + " " + par.get_unit() + " bin")

        ax.set(title="Normalised energy spectrum in "+str(i*1000)+"mm to "+str((i+1)*1000)+"mm ",ylabel="Events per " + str(width) + " " + par.get_unit() + " bin", xlabel=str(dim)+" [" + par.get_unit() + "]")
        ax.hist(x,bins,weights=spec.project(dim),histtype="stepfilled", color="RoyalBlue",label=spec._name)
        fig.savefig("slice_"+str(name)+"_"+str(i*1000)+"_"+str((i+1)*1000)+".png")
        i=1+i

 
def combiPlotter(spectra,dim,name):
    i=0
    fig = plt.figure()
    ax= fig.add_subplot(1,1,1)
    for spec in spectra:
        par = spec.get_config().get_par(dim)
        width = par.get_width()
        bins = np.linspace(par._low,par._high, par._bins+1)
        x = bins[:-1] + 0.5*width

        plt.xlabel(str(dim)+ " [" + par.get_unit() + "]")
        plt.ylabel("Events per " + str(width) + " " + par.get_unit() + " bin")

        ax.set(title="Normalised energy spectrum in 1000mm slices",ylabel="Events per " + str(width) + " " + par.get_unit() + " bin", xlabel="energy_reco"+ " [" + par.get_unit() + "]")
        ax.hist(x,bins,weights=spec.project("energy_reco"),label=spec._name,histtype='step')
        ax.set_ylim([0,0.03])
        ax.set_xlim([0.2,0.7])
        ax.legend(loc="best")
    fig.savefig("combined_"+str(name)+".png")

def func(path,nslice,name):
    spectra=[]
    slicer(path,spectra,nslice)
    signalPlotter(spectra,"energy_reco",name)
    combiPlotter(spectra,"energy_reco",name)

def po210():

    convertor("po210_ntuple/*")
    combiner("hdf5/SolarPo**ntuple*","hdf5/SolarPo210_combined.hdf5")
     
    plotpath="plots/"
    
    func("hdf5/SolarPo210_combined.hdf5",6,"po210")

def bi210():
    convertor("bi210_ntuple/*")
    combiner("hdf5/SolarBi**ntuple*","hdf5/SolarBi210_combined.hdf5")
    plotpath="plots/"
    func("hdf5/SolarBi210_combined.hdf5",6,"bi210")

def compair(spectrumPathReco,spectrumPathMC,name):
    spectraReco=[]
    spectraMC=[]
    slicerReco(spectrumPathReco,spectraReco,6)
    slicerMC(spectrumPathMC,spectraMC,6)

    for i in range(0,len(spectraReco)):
        fig = plt.figure()
        ax= fig.add_subplot(1,1,1)

        par = spectraReco[i].get_config().get_par("energy_reco")
        width = par.get_width()
        bins = np.linspace(par._low,par._high, par._bins+1)
        x = bins[:-1] + 0.5*width

        ax.set(title="Normalised energy spectrum in "+str(i*1000)+"mm to "+str((i+1)*1000)+"mm ",ylabel="Events per " + str(width) + " " + par.get_unit() + " bin", xlabel="Energy [" + par.get_unit() + "]")
        ax.hist(x,bins,weights=spectraReco[i].project("energy_reco"),histtype="stepfilled",label=spectraReco[i]._name)

        par = spectraMC[i].get_config().get_par("energy_mc")
        width = par.get_width()
        bins = np.linspace(par._low,par._high, par._bins+1)
        x = bins[:-1] + 0.5*width

        ax.hist(x,bins,weights=spectraMC[i].project("energy_mc"),histtype="stepfilled",label=spectraMC[i]._name,alpha=0.75)
        ax.legend(loc=2)    
        fig.savefig("compare_"+str(name)+"_"+str(i*1000)+"_"+str((i+1)*1000)+".png")

        


if __name__=="__main__":

    print "You need to compare the recon against the mc"
    print "You should bin in bigger bins becuase you could then bin in 4d"
    """You need to plot the standard spectra"""
