# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#Author: Alex Sun
#Date: 11152015
#Purpose: generates pce model
#For simple leakage
# uncertain vars include
# perm [inzone, azmi]
# por [inzone, azmi]
# Aquitard thickness 
# Qinj
#Dependencies: chaospy
#The user needs to re-implement Problem and Geom for their problems
#==============================================================================
import numpy as np
import chaospy as cp
import os, sys
import pickle as pkl
from settings import Problem,Geom

def genLauncherFile(allSamples):
        #conver to lognormal
        Nvar=allSamples.shape[1]
        allSamples[:,0] = np.exp(allSamples[:,0])
        allSamples[:,1] = np.exp(allSamples[:,1])
        fo = open('paramlist', 'wb')
        counter=0
        for item in allSamples:
            runcmd='$WORK/anaconda2/bin/python elsamain.py --runnum {0} --vars '.format(counter)
            for i in range(Nvar):
                runcmd = runcmd + '{0} '.format(item[i])
            runcmd = runcmd + '\n'
            fo.write(runcmd)
            counter+=1
        fo.close()

class RandVar(object):
    '''
    This class describes a random variable
    Normal: (mean, std)
    Uniform: (lb, ub)
    '''
    def __init__(self, distType, param):
        '''
        @param distType, type of prob distribution
        @param param, parameters associated with the prob distribution
        '''
        self.distType = distType
        self.param = param


class VarSet(object):
    '''
    This class describes a random variable set
    for base case
    This is used in MC simulation and basecase
    '''
    def __init__(self):
        #inzone perm [m^2]
        var0 = RandVar('normal', [np.log(1e-13), 0.5])
        #above zone perm [m^2]
        var1 = RandVar('normal', [np.log(5e-14), 0.5])
        #inzone phi [-]
        var2 = RandVar('uniform', [0.1, 0.2])
        #above zone phi [-]
        var3 = RandVar('uniform', [0.05, 0.3])
        #aquitard thickness [m]
        var4 = RandVar('uniform', [10, 30])
        #injection rate [Mt/yr]
        var5 = RandVar('uniform', [0.5, 5])
        self.varset = [var0, var1, var2, var3, var4, var5]

        self.Nvar = len(self.varset)

        self.jointDist = cp.J(
            cp.Normal(mu=var0.param[0], sigma=var0.param[1]),
            cp.Normal(mu=var1.param[0], sigma=var1.param[1]),
            cp.Uniform(lo=var2.param[0], up=var2.param[1]),
            cp.Uniform(lo=var3.param[0], up=var3.param[1]),
            cp.Uniform(lo=var4.param[0], up=var4.param[1]),
            cp.Uniform(lo=var5.param[0], up=var5.param[1]),
        )


class PCEDriverSimple(object):
    def __init__(self, varsetInst, regen=True):
        self.varset = varsetInst.varset
        self.Nvar = varsetInst.Nvar
        self.jointDist = varsetInst.jointDist
        order = 3   #order of ortho polynomial    
        if regen:
            self.polynomials = cp.orth_ttr(order, self.jointDist)
            self.nodes, self.weights = cp.generate_quadrature(order+1, self.jointDist, rule="G", sparse=True)
            pkl.dump([self.polynomials, self.nodes,self.weights], open('pcecoeff.dump', 'wb'))
        else:
            [self.polynomials, self.nodes,self.weights]=pkl.load(open('pcecoeff.dump', 'rb')) 
        self.nsamples = self.nodes.shape[1]
        self.co2Model = None
        
    def genParamFile(self):
        '''
        This function generates the parameterlist file needed by Launcher
        The calculation is done in /wrangler/elsaproject/basecase 
        '''
        genLauncherFile(np.copy(self.nodes.T))

    def postProcess(self, Tmax=20, reLoad=True):
        #load the results files
        Tmax = Tmax*365.25
        leakCVol = None
        if reLoad: 
            print 'loading results from {0} files...'.format(self.nsamples)
            for runnum in range(self.nsamples):
                print runnum
                [monP,monh,cumMassCO2, cumMassBrine] = pkl.load(open('time%s_%s.dump'% (Tmax, runnum), 'rb'))
                if (leakCVol is None):
                    leakCVol = np.zeros((cumMassCO2.shape[0], self.nsamples))
                    leakBVol = np.zeros((cumMassBrine.shape[0], self.nsamples))
                leakCVol[:,runnum] = cumMassCO2
                leakBVol[:,runnum] = cumMassBrine
            print 'finished loading results'
            pkl.dump([leakCVol, leakBVol], open('bigmat.dump', 'wb'))
        else:
            [leakCVol,leakBVol] = pkl.load(open('bigmat.dump', 'rb')) 
        
        layer = 2 
        nlw = int(0.5*leakCVol.shape[0])
        ilw = (layer-1)*nlw+0    
        leak = leakCVol[ilw,:]

        #Note: the sample matrix needs to be nsamples * N_predictedpoints
        models = cp.fit_quadrature(self.polynomials, self.nodes, self.weights, leakCVol.T)
        #simple test
        asample = self.jointDist.sample(1,'H') 
        print 'leaked co2=', models(asample[0],asample[1],asample[2],asample[3],asample[4],asample[5])


    def getProblemGeom(self):
        '''        
        @return:  problem configuration for the current problem
        '''
        return  Problem()
        
    def getMetaModels(self):
        '''
        This function generates meta models using pre-run samples  
        '''
        [leakCVol,leakBVol] = pkl.load(open('bigmat.dump', 'rb')) 
        #this is the total numer of leak wells in each layer
        self.nlw = int(0.5*leakCVol.shape[0])
        #
        #Note: the sample matrix needs to be nsamples * N_leakywells
        #
        self.co2Model = cp.fit_quadrature(self.polynomials, self.nodes, self.weights, leakCVol.T)
        self.brineModel = cp.fit_quadrature(self.polynomials, self.nodes, self.weights, leakBVol.T)
        
    
    def getTotalLeak4Web(self, asample):
        '''
        Print out total leak mass from all leaky wells in the above zone for a given realization
        @param asample, a realization from the joint distribution
        '''
        if (self.co2Model is None):
            self.getMetaModels()
            
        metaLeakCO2Mass = self.co2Model(asample[0],asample[1],asample[2],asample[3],asample[4],asample[5])       
        metaLeakBrineMass = self.brineModel(asample[0],asample[1],asample[2],asample[3],asample[4],asample[5])

        return metaLeakCO2Mass[self.nlw:], metaLeakBrineMass[self.nlw:]
    
    def getMCResults4Web(self):
        '''
        Use surrogate model to perform Monte Carlo simulation
        @param nmc, number of Monte Carlo realizations 
        '''
        if (self.co2Model is None):
            self.getMetaModels()
        
        meanCO2 = cp.E(self.co2Model, self.jointDist)
        stdCO2 = cp.Std(self.co2Model, self.jointDist)
        
        meanBrine = cp.E(self.brineModel, self.jointDist)
        stdBrine = cp.Std(self.brineModel, self.jointDist)
 
        return [meanCO2[self.nlw:], stdCO2[self.nlw:], meanBrine[self.nlw:], stdBrine[self.nlw:]]

class MCSimple():
    def __init__(self, varsetInst, nsample):
        self.varset = varsetInst.varset
        self.Nvar = varsetInst.Nvar
        jointDist = varsetInst.jointDist

        self.randSamples = jointDist.sample(nsample, 'H')

    def genParamFile(self):
        '''
        The calculation is done in /wrangler/elsaproject/basemc        
        '''        
        print 'reached here'
        genLauncherFile(np.copy(self.randSamples.T))

if __name__ == "__main__":
    METHOD='MC'
    
    if METHOD=='PCE':
        myDriver = PCEDriverSimple(VarSet(),regen=False)
        #max simulation time in yr
        Tmax = 20
        #myDriver.postProcess(Tmax, reload=False)
        asample = myDriver.jointDist.sample(1, 'H')
        myDriver.getTotalLeak4Web(asample)
        myDriver.getMCResults4Web()

    elif METHOD=='MC':
        myDriver = MCSimple(VarSet(), 1000)
        myDriver.genParamFile()
