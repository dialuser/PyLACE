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
#Author: Alex sun
#Date: 11152016
#Date 3/27/2017
#Use this script to generate problem configuration
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
from core import Param, Assembler
import sklearn.preprocessing as skp
import pickle as pkl
import sys
from matplotlib import rcParams

#set up font for plot
params = {'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'xtick.major.size': 1.5,      # major tick size in points
          'xtick.minor.size': .5,      # minor tick size in points
          'ytick.major.size': 1.5,      # major tick size in points
          'ytick.minor.size': .5,      # minor tick size in points
          'xtick.major.pad': 1,      # distance to major tick label
          'xtick.minor.pad': 1,
          'ytick.major.pad': 2,      # distance to major tick label
          'ytick.minor.pad': 2,
          'axes.labelsize': 10,
          'axes.linewidth': .5,
          'font.size': 10,
          'lines.markersize': 4,            # markersize, in points
          'legend.fontsize': 11,
          'legend.numpoints': 2,
          'legend.handlelength': 1.
          }
rcParams.update(params)

class Problem(object):
    def __init__(self):
        self.uvec = [5.14025510793e-15, 2.57012755396e-14, 0.0921850137466, 0.142268126477, 20.0, 4.1546984073]
        self.param = Param(self.uvec)
        self.W = 4000 #[m], domain width
        self.L = 4000 #[m], domain length
        self.dx = 1000 #[m], x dist between leaky wells
        self.dy = 1000 #[m], y dist between leaky wells
        self.mx = 200
        self.my = 200 
        self.geom = Geom(width=self.W, length=self.L, dx=self.dx, dy=self.dy, mx=self.mx, my=self.my)
    
    def plotProblemConfig(self):
        '''
        Plots geometry for the web
        '''
        self.geom.plot()
        
    def setupAndRun(self):
        Tmax = 365.25*1
        
        self.param.setSimulationTime(Tmax) 
        
        self.param.setLeakWells(self.geom)
        self.param.setInjWells(Minj=self.uvec[-1], geom=self.geom)
        self.param.setMonWells(self.geom)

        RERUN = False
        if (RERUN):
            mytest = Assembler(self.param)
            mytest.iterate(isPlot=False)
            monP,monh,cumMassCO2, cumMassBrine = mytest.getResults()        
            pkl.dump([monP,monh,cumMassCO2, cumMassBrine], open('time%s.dump'%Tmax, 'wb'))
        else:
            [monP,monh,cumMassCO2, cumMassBrine] = pkl.load(open('time%s.dump'%Tmax, 'rb'))
            print len(monP)
            
        fig = plt.figure(1, figsize=(6,6))
        layer = 2 # 1 inzone, 2 above zone
        self.plotResult(monP, layer, fig, self.geom)
        fig = plt.figure(2, figsize=(6,6))
        self.plotResult(monh, layer, fig, self.geom)
        
    def plotResult(self, monres, layer, fig=None, geom=None):
        if (not fig):
            plt.figure()
        fig.hold(True)
        geom.plot(fig)                 
        msize = skp.minmax_scale(monres[(layer-1)*geom.Nm:layer*geom.Nm,-1],feature_range=(1,20))
        for i in range(geom.Nm):
            plt.plot(geom.monwellx[i], geom.monwelly[i], 'o', color='#51555b', markersize=msize[i])
                
class Geom(object):
    def __init__(self, width, length, dx, dy, mx, my):
        #form leak locations
        X,Y = np.meshgrid(np.arange(-0.5*width, 0.5*width+dx, dx), np.arange(-0.5*length, 0.5*length+dy, dy))
        X = np.reshape(X, (X.shape[0]*X.shape[1]))
        Y = np.reshape(Y, (Y.shape[0]*Y.shape[1]))
        
        #set inj locs
        self.injx = [0.0]
        self.injy = [0.0]

        #remove inj locs from leak locs
        clocs=[]        
        for ix in range(len(self.injx)):
            indx = np.where(X==self.injx[ix]) 
            indy = np.where(Y==self.injy[ix])
            ind = np.intersect1d(indx, indy)
            if (ind.size):
                ind=ind[0]
                clocs.append(ind)
        self.wellx=np.delete(X, clocs)
        self.welly=np.delete(Y, clocs)
        self.Nleaks = len(self.wellx)
        
        # set monitoring locations
        X,Y = np.meshgrid(np.arange(-0.5*width-mx, 0.5*width+2*mx, mx), np.arange(-0.5*length-my, 0.5*length+2*my, my))
        X = np.reshape(X, (X.shape[0]*X.shape[1]))
        Y = np.reshape(Y, (Y.shape[0]*Y.shape[1]))
        # check coincidence with leak and injector locations and do adjustment
        rw = 0.15
        toler = 0.01
        for ix in range(len(self.injx)):
            indx = np.where(abs(X-self.injx[ix])<toler) 
            indy = np.where(abs(Y-self.injy[ix])<toler)
            ind = np.intersect1d(indx, indy)
            if (ind.size):
                ind=ind[0]
                X[ind] = X[ind] + rw
                Y[ind] = Y[ind] + rw

        for ix in range(len(self.wellx)):            
            indx = np.where(abs(X-self.wellx[ix])<toler) 
            indy = np.where(abs(Y-self.welly[ix])<toler)
            ind = np.intersect1d(indx, indy)

            if (ind.size):
                ind=ind[0]
                X[ind] = X[ind] + rw
                Y[ind] = Y[ind] + rw

        self.monwellx = X 
        self.monwelly = Y
        self.Nm = len(self.monwellx)
        
    def plot(self, fig=None):
        if (not fig):
            plt.figure(dpi=200, figsize=(8,8))
        plt.axes().set_aspect('equal')
        plt.hold(True)
        plt.plot(self.injx, self.injy, 'bs', markersize=8, label='Injector')
        plt.plot(self.monwellx, self.monwelly, 'g.', markersize=5, label='Monitoring Well')
        plt.plot(self.wellx, self.welly, 'ro', label='Leaky Well')

        #plt.legend(bbox_to_anchor=(0.1, 1.02, 0.9, .052), loc=3,
        #   ncol=3, mode="expand", borderaxespad=0.)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True, numpoints=1)
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.xlim((-3000,3000))
        plt.ylim((-3000,3000))
        plt.savefig('problemSimple.eps', pad_inches=0.0, bbox_inches='tight')
        
if __name__ == "__main__":
    problem = Problem()
    problem.plotProblemConfig()
    plt.show()
    