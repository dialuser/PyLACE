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
#@author: Alex Sun
#@date: 11/1/2016
#@desc: implements multilayer solution according to Bau paper, 2015 comput geoscie
#http://co2interface.princeton.edu/index.html
#chose the simple model
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt 
import sys

class Well(object):
    def __init__(self, x,y):
        '''
        @param x, y location of the well
        @param kw, permeability of the well
        '''
        self.x = x
        self.y = y
        self.rw = 0.15
        
    def getDist(self, x0, y0):
        '''
        @param x0,y0 injector location
        '''
        return np.sqrt((self.x-x0)*(self.x-x0)+(self.y-y0)*(self.y-y0))

class LeakWell(Well):
    def __init__(self,x,y,kw):
        '''
        @param, kw
        '''
        super(LeakWell, self).__init__(x,y)
        self.kw= kw
        
class MonWell(Well):
    def __init__(self, x, y):
        super(MonWell, self).__init__(x,y)


class InjWell(Well):
    def __init__(self, x, y, Q):
        '''
        @param Q: injection rate [m3/s]
        ''' 
        super(InjWell, self).__init__(x,y)
        self.Q = Q

    def getRate(self, t):
        return self.Q
    
    def getCumVol(self,t):
        return abs(self.Q*t)
 
 
class InjWellComplex(Well):
    '''
    11152016: for variable injection rate
    not fully implemented
    '''
    def __init__(self, x, y, Qrate, Tarray):
        super(InjWellComplex, self).__init__(x, y)
        self.Qrate = Qrate
        self.Tarray = Tarray 
 
    def setRate(self, Tarray, Qrate):
        '''
        @param Tarray: discretized time axis        
        '''      
        self.Qrate = Qrate
        self.Tarray = Tarray
        
    def getRate(self, t):
        return np.interp(t, self.Tarray, self.Qrate)

    def getCumVol(self, t):
        return np.trapz(self.Qrate, self.Tarray)
               
        
class Layer(object):
    def __init__(self, H, D):
        '''
        @param H, thickness of the layer
        @param D, depth from the surface to the bottom of the current layer      
        '''
        self.H = H
        self.D = D

class Aquifer(Layer):
    def __init__(self, H, D, k, phi):
        super(Aquifer, self).__init__(H, D)
        self.k = k
        self.phi = phi    
        self.P0 = None

            
class Param:
    '''
    Parameters for the problem
        L: number of aquifers, meaning there are L-1 aquitards
        M: number of injectors
        N: number of leak wells
        Nm: number of monitoring locations
        
    '''
    def __init__(self, uvec=None):
        self.uvec = uvec
        self.dict ={
            'rho_c': 479., #kg/m3
            'rho_b': 1045.,
            'mu_b': 0.0002535, #Pa*S 
            'mu_c': 3.95e-5,
            'cf': 4.2e-10,  #assumed to be the same as brine???          
            'sr_b': 0.2,
            'sr_c': 0.0, 
            'h1min': 0.45,
            'krcmax': 0.55, #max co2 relative perm, i.e, kc = krcmax*f(S)
            'g': 9.81,
            'rw': 0.15, 
            'kw':1e-12, #m^2
        }
        if not uvec:
            perm = [5e-14, 2e-14]
            H = [10., 10.] #aquifer thickness
            Ha = [10.]    #aquitard thickness
            phi = [0.15, 0.15]
        else:
            perm = [uvec[0], uvec[1]]
            H = [10., 10.] #aquifer thickness
            Ha = [uvec[4]]    #aquitard thickness
            phi = [uvec[2], uvec[3]]
                
        #this is the default   
        self.setSimulationTime(Tmax=365)
        #
        # set the number of layers
        #
        self.L = 2
        #depth to the bottom of the bottom most aquifer
        DL = 3010.0         
        # define aquifers
        self.aquifers=[]
        D0 = DL        
        for i in range(self.L):
            aquifer = Aquifer(H[i], D0, perm[i], phi[i])
            aquifer.P0 = self.dict['rho_b']*self.dict['g']*D0
            self.aquifers.append(aquifer)
            #move to next aquifer bottom
            if (i<self.L-1):
                D0 = D0-Ha[i]-H[i]
        # define aquitard
        self.aquitards=[]
        D0 = DL - H[0]
        for i in range(self.L-1):
            self.aquitards.append(Layer(Ha[i], D0))
            D0 = D0 - Ha[i]-H[i+1]                        
        #
        # define injectors
        #
        self.setInjWells()
        #
        # define leaky wells
        #
        self.setLeakWells()
        #
        # define monitoring locations
        #
        self.setMonWells()
        

    def setSimulationTime(self, Tmax=100):
        '''
        form total simulation time [d]
        @param Tmax, maximum simulation time in days
        '''
        self.Tmax = Tmax #[d]
        if (Tmax<=100.0):
            #self.T = np.array(np.arange(0.1,self.Tmax, 0.5)*86400) #in [s]    
            self.T = np.array(np.arange(0.1,self.Tmax, 0.5)) #in [d]
        else:
            self.T = np.concatenate((np.arange(1.0,100, 1.0)*86400, np.arange(100, Tmax, 5.)*86400))
            
    def setLeakWells(self, geom=None):
        '''
        Set leak locations
        @param wx, coords of well x
        @param wy, coords of well y
        '''
        self.leaks = []

        if (not geom):
#             self.N = 24
#             wx = [-2000., -1000.,     0.,  1000.,  2000., -2000., -1000.,     0.,  1000.,  2000.,
#                     -2000., -1000.,  1000.,  2000., -2000., -1000.,     0.,  1000.,  2000., -2000.,
#                         -1000.,     0.,  1000.,  2000.]
#             wy = [-2000., -2000., -2000., -2000., -2000., -1000., -1000., -1000., -1000., -1000.,
#                   0.,     0.,     0.,     0.,  1000.,  1000.,  1000.,  1000.,  1000.,  2000.,
#                   2000.,  2000.,  2000.,  2000.]
            self.N=1
            wx = [-500.]
            wy = [-500.]
            for i in range(self.N):
                self.leaks.append(LeakWell(wx[i],wy[i], self.dict['kw']))            
        else:
            #set defaults
            self.N = geom.Nleaks
            for i in range(self.N):
                self.leaks.append(LeakWell(geom.wellx[i], geom.welly[i], self.dict['kw']))

    def setInjWells(self, Minj=1, geom=None):
        self.injw = []

        Qinj = Minj*1e9/(86400.0*365.25) #kg/s, Qinj = 31.68808781 #kg/s = 1Mt per year
        Qinj = Qinj/self.dict['rho_c'] #m3/s, volumetric rate                

        if not geom:
            self.M = 1
    
            ComplexWell = False
            if ComplexWell:
                r = np.random.RandomState(1234)
                rvec = r.randn(len(self.T))
                sigm = 0.10
                Qrate = Qinj*(1+sigm*rvec)
                for i in range(self.M):
                    self.injw.append(InjWellComplex(0, 0, Qrate, self.T))                    
            else:
                for i in range(self.M):
                    self.injw.append(InjWell(0, 0, Qinj))        
        else:
            self.M= len(geom.injx)
            for i in range(self.M):
                self.injw.append(InjWell(geom.injx[i], geom.injy[i], Qinj))
    
    def setMonWells(self, geom=None):
        self.monwells = []

        if not geom:
            self.Nm = 1
            wx = [-1000.15, ]
            wy = [-1000.15, ]
            for i in range(self.Nm):
                self.monwells.append(MonWell(wx[i],wy[i]))        
        else:
            self.Nm= geom.Nm
            for i in range(self.Nm):
                self.monwells.append(MonWell(geom.monwellx[i], geom.monwelly[i]))
                
class Assembler:
    '''
    #sign convention: upward Q is positive, downward negative
    Main class for solving P and S iteratively
    
    '''
    def __init__(self, param):
        self.param = param
        self.injw = self.param.injw
        self.leakw = self.param.leaks
        self.monw = self.param.monwells
        
    def iterate(self, isPlot=False):
        L = self.param.L
        N = self.param.N
        Nm = self.param.Nm
        
        Pvec = np.zeros((L*N))
        hvec = np.zeros((L*N))
        Qvec0 = np.zeros(((L-1)*N))
        Pmon0 = np.zeros((L*Nm))
        #
        # cumulative mass of co2, brine and total
        # note: vol is a misnormer
        #
        self.cumVolC = np.zeros((L*N))
        self.cumVolB = np.zeros((L*N))
        self.cumVol = np.zeros((L*N))         
        
        Qold = np.zeros(((L-1)*N))
        cumVol0 = np.zeros((L*N))
        cumVolC0 = np.zeros((L*N))
        cumVolB0 = np.zeros((L*N))

        #Time discretization
        T = self.param.T
        #record time history 
        self.Phist = np.zeros((L*N, len(T)))
        self.Pmonhist = np.zeros((L*Nm, len(T)))        
        self.hmonhist = np.zeros((L*Nm, len(T)))

        #maximum number of iterations
        Iter = 50
        #convergence criterion
        epsil =0.001
        #get initial conditions
        for l in range(L):
            for iw in range(N):
                j = l*N+iw
                Pvec[j] = self.param.aquifers[l].P0
            for iw in range(Nm):
                j = l*Nm+iw
                Pmon0[j] = self.param.aquifers[l].P0
        P0 = np.copy(Pvec)
        
        
        dampfact = 0.4
                
        counter=0
        deltaT = 0.05*86400 #[s]
        t0 = 0.0
        deltaTGood = 0.0
        t = 0.05*86400 #[s]
        tfail = 0.001*86400 #[s]]
        Tmax = max(T)
        while t<Tmax:         
            it = 1
            print '****** t=%s ********' % (t/86400.0)
            while True:
                #print 'iter=%s' % it
                Pvec, hvec = self.getPS(hvec, t)
                #check for pressure buildup                    
                Qvec = self.getQ(Pvec, hvec, Qvec0, cumVol0, cumVolC0, cumVolB0, t, (t-t0))
                Qvec = dampfact*Qvec + (1-dampfact)*Qold 
                #deal with zero flow
                validQ = np.where(Qvec>1e-9)[0]

                if validQ.size:                    
                    err = abs(np.divide((Qold[validQ]-Qvec[validQ]), Qvec[validQ]))
                    print max(err)
                    if (max(err) < epsil):
                        break
                else:
                    break
                it +=1
                if it >Iter:
                    print 'bad iteration'
                    #raise Exception('did not converge')
                    #cut time step
                    print 'warning: cutting time step ' 
                    deltaT = deltaT*0.6   
                    t= t0+deltaT    
                    it = 1    
                    if (deltaT<tfail):
                        print 'failed values ', self.param.uvec
                        raise Exception('failed to converge')    
                Qold = Qvec
                
            cumVol0 = np.copy(self.cumVol)
            cumVolC0 = np.copy(self.cumVolC)
            cumVolB0 = np.copy(self.cumVolB)
            Qvec0 = np.copy(Qvec)
            if ((T[counter]-t)<(1e-3)*86400):
                #output
                self.Phist[:,counter] = np.divide(Pvec - P0, 1e6)
                Pvec, h1mon = self.updateMonWells( hvec, t)
                self.Pmonhist[:, counter] = np.divide(Pvec-Pmon0, 1e6)
                self.hmonhist[:, counter] = h1mon
                counter+=1
            t0 = t
            deltaTGood = max([1.2*deltaT, deltaTGood])
            print 'dt=', deltaT
            deltaT = min([T[counter]-t, deltaTGood])
            t += deltaT
        
        if isPlot:
            plt.figure()
            #show Delta P
            plt.hold(True)
            for l in range(L):
                for iw in range(N):
                    plt.plot(T/86400, self.Phist[l*N+iw,:], label='layer%s, iw=%s' % (l, iw))
            plt.legend(loc=0)
                    
            plt.figure()
            plt.hold(True)
            for l in range(L):
                for iw in range(Nm):
                    plt.plot(T/86400, self.Pmonhist[l*Nm+iw,:], label='layer%s, iw=%s' % (l, iw))
            
            plt.legend(loc=0)
    
    def getResults(self):
        return self.Pmonhist, self.hmonhist, self.cumVolC*self.param.dict['rho_c'], self.cumVolB*self.param.dict['rho_b']
    
    def updateMonWells(self, hvec, t):
        '''
        Get pressure for each monitoring well        
        '''
        L = self.param.L
        Nm = self.param.Nm #number of monitoring well
        N = self.param.N
        M = self.param.M
     
        dgamma = (self.param.dict['rho_b']-self.param.dict['rho_c'])*self.param.dict['g']
        
        P = np.zeros((L*Nm))
        h1 = np.zeros((L*Nm))
        for l in range(L):
            aquifer = self.param.aquifers[l]
            for ilw in range(Nm):                
                currentMon = self.monw[ilw]
                j = l*Nm+ilw                
                P[j]=0.0
                hAll=[]       
                if (l==0):
                    #sources are the injectors
                    #This is the boundary condition and should stay fixed
                    for iw in range(M):
                        currentInj = self.injw[iw]
                        r = currentMon.getDist(currentInj.x, currentInj.y)
                        dP, htemp = self.__getP1(r, t, aquifer, currentInj.getCumVol(t), abs(currentInj.getRate(t)), self.param.dict)                        
                        P[j]+=dP
                        hAll.append(htemp)
                    #sinks
                    for jlw in range(N):
                        if (abs(self.cumVol[jlw])>0):
                            rhoeff = hvec[jlw]*self.param.dict['rho_c']+(1-hvec[jlw])*self.param.dict['rho_b']
                            cumVol  = self.cumVol[jlw]
                            Qavg =  cumVol/t                            
                            aLeak = self.leakw[jlw]
                            r = currentMon.getDist(aLeak.x, aLeak.y)
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, Qavg, self.param.dict, hvec[jlw])
                            P[j]+= dP
                else:
                    #sources are the leak wells        
                    for jlw in range(N):
                        if (abs(self.cumVol[l*N+jlw])>0):  
                            cumVol  = self.cumVol[l*N+jlw]
                            Qavg =  cumVol/t                            
                            aLeak = self.leakw[jlw]
                            r = currentMon.getDist(aLeak.x, aLeak.y)                                
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, Qavg, self.param.dict, hvec[(l-1)*N+jlw])
                            #the co2 saturaiton should only be possible when plume reaches here
                            htemp = min([htemp, hvec[(l-1)*N+jlw]]) #by this time, hvec should have been updated

                            hAll.append(htemp)
                            P[j]+=dP
                    #sinks
                    if (not l==L-1):
                        if (abs(self.cumVol[(l+1)*N+jlw])>0):
                            #average leak out rate
                            rhoeff = hvec[l*N+jlw]*self.param.dict['rho_c']+(1-hvec[l*N+jlw])*self.param.dict['rho_b']
                            cumVol  = self.cumVol[(l+1)*N+jlw]
                            Qavg =  cumVol/t
                            
                            aLeak = self.leakw[jlw]
                            r = currentMon.getDist(aLeak.x, aLeak.y)
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, -Qavg, self.param.dict, hvec[l*N+jlw])
                            P[j]+=dP
                            
                P[j] = dgamma*aquifer.H*P[j] + aquifer.P0
                if not hAll:
                    h1[j] = 0.0
                else:
                    h1[j] = max(hAll)  
        return P, h1
            
    
    def getQ(self, Pvec, hvec, Qvec0, cumVol0, cumVolC0, cumVolB0, t, dt):
        '''
        calculate flux across each leaky well
        @param: Pvec, pressure for the current iteration
        @param: Qvec0, flux for the previous step
        @param: cumVol0, cumulative total volume up to the previous step
        @param: cumVolC0, cumulative co2 leak volume up to the previous step
        @param: t, current time
        @param: dt,current time step  
        '''
        A = np.pi*self.param.dict['rw']*self.param.dict['rw']
        
        L = self.param.L
        N = self.param.N
        g_const = self.param.dict['g']
        rhoc = self.param.dict['rho_c']
        rhob = self.param.dict['rho_b']
        gammac = rhoc*g_const
        gammab = rhob*g_const        
                     
        Qvec=np.zeros((L-1)*N)
        h1min = self.param.dict['h1min']
        Qmax = 0.1*self.injw[0].Q
        
        for l in range(L-1):
            aquitard = self.param.aquitards[l]
            aquifer = self.param.aquifers[l]
            for ilw in range(N):
                mc,mb = self.getLambda(hvec[l*N+ilw])

                Pdn = Pvec[l*N+ilw] 
                Pup = Pvec[(l+1)*N+ilw]
                
                #no upconing
                if hvec[l*N+ilw]>h1min:
                    Qc = A*self.leakw[ilw].kw*mc*(-Pup+Pdn-gammac*(aquitard.H+aquifer.H))/aquitard.H
                    Q = Qc
                elif hvec[l*N+ilw]<=h1min and hvec[l*N+ilw]>0.01:
                    Qb = A*self.leakw[ilw].kw*mb*(-Pup+Pdn-gammab*(aquitard.H+aquifer.H))/aquitard.H
                    Qc = A*self.leakw[ilw].kw*mc*(-Pup+Pdn-gammac*(aquitard.H+aquifer.H))/aquitard.H 
                    Q = Qb+Qc     
                else:
                    Qb = A*self.leakw[ilw].kw*mb*(-Pup+Pdn-gammab*(aquitard.H+aquifer.H))/aquitard.H
                    Qc = 0.0
                    Q = Qb                                  
                Q = min(Q, Qmax)
                #average between the current and the previous time step           
                dM = dt*0.5*(Qvec0[l*N+ilw]+Q)
                self.cumVolC[l*N+ilw] = cumVolC0[l*N+ilw] - Qc*dt
                self.cumVolC[(l+1)*N+ilw] = cumVolC0[(l+1)*N+ilw] + Qc*dt
                self.cumVolB[l*N+ilw] = cumVolB0[l*N+ilw] - Qb*dt
                self.cumVolB[(l+1)*N+ilw] = cumVolB0[(l+1)*N+ilw] + Qb*dt
                
                self.cumVol[l*N+ilw] =  cumVol0[l*N+ilw] - dM
                self.cumVol[(l+1)*N+ilw] = cumVol0[(l+1)*N+ilw] + dM
                
                Qvec[l*N+ilw] = Q  
        return Qvec           

    def getLambda(self, hval=None):
        '''
        hval (normalized by aquifer thickness) is equivalent to actual saturation 
        '''
        krcmax = self.param.dict['krcmax']
        krbmax = 1.0
        muc = self.param.dict['mu_c']
        mub = self.param.dict['mu_b']
        # co2 saturation at brine residual
        Srb = self.param.dict['sr_b']    
        Src = self.param.dict['sr_c']
        h1min = self.param.dict['h1min']
        
        if (hval is None):
            #for injector (h always == 1) 
            krc = krcmax
            krb = krbmax
        else:
            #these should be for leaky wells            
            if hval > h1min:
                #single phase co2  
                krc = krcmax
                krb = 1.0
            elif hval>0.05 and hval<=h1min:
                #two phase, linear relative perm       
                #these should be for leaky wells
                Sb = 1-hval
                S = (Sb-Srb)/(1-Srb-Src)            
                if S >= 0.05 and S < 1.0:
                    #two phase, linear relative perm                            
                    #krc = krcmax*np.power(1-S,2)*(1-np.power(S,2))
                    #krb = krbmax*np.power(S, 2)
                    krc = krcmax*(1-S)
                    krb = krbmax*S
                elif S < 0.05: 
                    #single phase brine
                    krc = krcmax
                    krb = 0.5
                elif S>= 1.0:
                    krc = 0.01
                    krb = krbmax                                     
            else:
                #single phase brine
                krc = 0.005
                krb = 1.0
                
        mc = krc/muc
        mb = krb/mub
        return mc,mb        
    
    def getPS(self, hvec0, t):
        '''
        Get pressure and co2 column height for each leaky well
        @param hvec0: co2 column height from the previous step
        '''
        L = self.param.L
        N = self.param.N
        M = self.param.M
     
        dgamma = (self.param.dict['rho_b']-self.param.dict['rho_c'])*self.param.dict['g']
        
        P = np.zeros((L*N))
        hvec = np.zeros((L*N))
        for l in range(L):
            aquifer = self.param.aquifers[l]
            for ilw in range(N):
                currentLeak = self.leakw[ilw]
                j = l*N+ilw                
                P[j]=0.0
                hAll=[]
               
                if (l==0):
                    #sources are the injectors
                    #This is the boundary condition and should stay fixed
                    for iw in range(M):
                        currentInj = self.injw[iw]
                        r = currentLeak.getDist(currentInj.x, currentInj.y)
                        dP, htemp = self.__getP1(r, t, aquifer, currentInj.getCumVol(t), abs(currentInj.getRate(t)), self.param.dict)                        
                        P[j]+=dP
                        hAll.append(htemp)
                    #sinks
                    for jlw in range(N):
                        if (abs(self.cumVol[jlw])>0):
                            cumVol  = self.cumVol[jlw]
                            Qavg =  cumVol/t
                            if (not jlw==ilw):
                                anotherLeak = self.leakw[jlw]
                                r = currentLeak.getDist(anotherLeak.x, anotherLeak.y)
                            else:
                                r = currentLeak.rw
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, Qavg, self.param.dict, hvec0[jlw])
                            P[j]+= dP 
                else:
                    #sources are the leak wells themselves        
                    for jlw in range(N):
                        if (abs(self.cumVol[l*N+jlw])>0):                    
                            cumVol  = self.cumVol[l*N+jlw]
                            Qavg =  cumVol/t
                            
                            #see Nordbotten 2009, p.746                             
                            if (not jlw==ilw):
                                anotherLeak = self.leakw[jlw]
                                r = currentLeak.getDist(anotherLeak.x, anotherLeak.y)                                
                            else:
                                r = currentLeak.rw
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, Qavg, self.param.dict, hvec0[l*N+jlw])
                            htemp = min([htemp, hvec[(l-1)*N+jlw]]) #by this time, hvec should have been updated
                            P[j]+=dP
                            hAll.append(htemp)
                    #sinks
                    if (not l==L-1):
                        if (abs(self.cumVol[(l+1)*N+jlw])>0):
                            cumVol  = self.cumVol[(l+1)*N+jlw]
                            Qavg =  cumVol/t

                            if (not ilw==jlw):
                                anotherLeak = self.leakw[jlw]
                                r = currentLeak.getDist(anotherLeak.x, anotherLeak.y)
                            else:
                                r = currentLeak.rw
                            dP, htemp = self.__getP1(r, t, aquifer, cumVol, -Qavg, self.param.dict, hvec0[l*N+jlw])
                            P[j]+=dP
                        
                P[j] = dgamma*aquifer.H*P[j] + aquifer.P0
                if not hAll:
                    hvec[j] = 0.0
                else:
                    hvec[j] = max(hAll)     
        return P, hvec        
                
    def __getP1(self, r, t, aquifer, cumVol, Q, param, hval=None):
        '''
        Note: hval is zero in single phase region
        @param r, distance from injector, in [m]
        @param t, time in [day]
        @param aquifer, current aquifer
        @param cumVol, cumulative mass
        @param Q, injection rate for inj, or average leak flux for leak well 
        @return P, absolute pressure in Pa
        '''
        H = aquifer.H
        phi = aquifer.phi
        k = aquifer.k
        g_const = param['g']
        Swr  = param['sr_b']
        rhow = param['rho_b']
        rhoc = param['rho_c']
        cf = param['cf']
        
        mc,mw = self.getLambda(hval)
            
        #Q:in is positive  
        chi = 2.*np.pi*H*phi*(1.0-Swr)*r*r/abs(cumVol) #eq. (3)
        
        gamma = 2*np.pi*(rhow-rhoc)*g_const*k*mw*H*H/Q #eq. 4

        psi = 4.5*np.pi*H*phi*k*mw*(1.0-Swr)/(abs(Q)*cf) #eq. 5 

        # mobility ratio
        mr=mc/mw        
        P1 = 0.0

        def __getRPlume():
            return np.sqrt((mr*cumVol)/(np.pi*(1-Swr)*phi*H))
        
        def P2(X):
            return -np.log(X/psi)/(2.0*gamma)
        
        def P3(X):
            return (1.0-np.sqrt(X/(2.0*mr)))/gamma + P2(2.*mr) + F(X)
        
        def P4(X):
            return -np.log(X*mr*0.5)/(2.0*mr*gamma) + P3(2.0/mr)
        
        def F(X):                 
            return -mr/(mr-1.0)*(getH(X) - np.log((mr-1.0)*getH(X)+1.0)/(mr-1.0)) # Eq.24 of Nordbotten 2009

        def getH(X):
            if (mr>0):
                if X<= 2/mr:
                    h1 = 1.0
                elif X<2*mr:                
                    h1 = (np.sqrt(2.0*mr/X)-1.0)/(mr-1.0) #eq. 6    
                    h1 = max([min([h1, 1.0]), 0.0])    
                else:
                    h1 = 0.0
            else:
                h1 = 0.0            
            
            return h1


        h1 = getH(chi)       
        if (chi>=psi):
            P = P1
        elif (chi>=2*mr):
            P = P2(chi)
        elif (chi>2/mr):            
            P = P3(chi)     
        else:
            P = P4(chi)
        return P, h1

    def getRoutP(self,t):
        '''
        Get the outer radius of plume buildup for inzone
        '''
        return np.sqrt((2.25*self.param.aquifers[0].k*t)/(self.param.dict['mu_b']*self.param.dict['cf']))
        
    def inzone(self):
        '''
        11/6/2016: validated inzone solution against
        http://co2interface.princeton.edu/index.html, chose the simple model
        '''
        dgamma = (self.param.dict['rho_b']-self.param.dict['rho_c'])*self.param.dict['g']
        t =365.25*86400  #[d]
        Rmax = self.getRoutP(t)
        rarray =np.arange(1, Rmax, 10)
        rarray =np.arange(0.1, 2000, 10)
        aquifer = self.param.aquifers[0]
        P = []
        H = []
        currentInj = self.injw[0]

        for r in rarray:
            Ptemp, htemp=self.__getP1(r, t, aquifer, currentInj.getCumVol(t), (currentInj.Q), self.param.dict)    
            P.append(dgamma*aquifer.H*Ptemp) #this is the pressure build-up only
            H.append((1-htemp)*aquifer.H)
        P = np.divide(P, 1e6) #convert o MPa
        plt.figure()
        plt.subplot(211)
        plt.plot(rarray, P)
        plt.xlabel('r')
        plt.subplot(212)
        plt.plot(rarray, H)
        plt.xlabel('r')
        
        r = 2000.
        tarray = np.arange(0.1, 365.25, 0.1)*86400 
        P = []
        H = []
        for t in tarray:
            Ptemp, htemp=self.__getP1(r, t, aquifer, currentInj.getCumVol(t), (currentInj.Q), self.param.dict)    
            P.append(dgamma*aquifer.H*Ptemp) #this is the pressure build-up only
            H.append((htemp)*aquifer.H)
        P = np.divide(P, 1e6) #convert o MPa
        plt.figure()
        plt.subplot(211)
        plt.plot(tarray/86400, P)
        plt.xlabel('t')
        plt.subplot(212)
        plt.plot(tarray/86400, H)
        plt.xlabel('t')
    
    def inzoneFunc(self, t, rarray):
        dgamma = (self.param.dict['rho_b']-self.param.dict['rho_c'])*self.param.dict['g']
        aquifer = self.param.aquifers[0]
        currentInj = self.injw[0]
        P = []
        
        for r in rarray:
            Ptemp, _ =self.__getP1(r, t, aquifer, currentInj.getCumVol(t), (currentInj.Q), self.param.dict)    
            P.append(dgamma*aquifer.H*Ptemp) #this is the pressure build-up only

        P = np.divide(P, 1e6) #convert o MPa
        
         
if __name__ == "__main__":
    
    
    mytest = Assembler(Param())
    mytest.inzone()
    plt.show()
    sys.exit()
    mytest.iterate(isPlot=True)
    #mytest.inzone()
    a,b,c,d = mytest.getResults()
    print b
    plt.show()
    