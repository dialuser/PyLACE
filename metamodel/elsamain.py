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
#Author: Alex Sun
#Date: 12/9/2017
#Main interface for Elsa
#==============================================================================
from core import Param, Assembler
from settings import Geom
import pickle as pkl
import argparse
import sys
import matplotlib.pyplot as plt
import sklearn.preprocessing as skp

def main(argv):

    try:
        parser = argparse.ArgumentParser(description="run elsa model") 
        
        #use '+' to gather all commands in one list
        parser.add_argument('--runnum', type=int, help="run number")
        parser.add_argument('--vars', nargs='+', help='list of uncertain variables')
        
        args = parser.parse_args()
        runnum=args.runnum
        ustring = args.vars
        #convert to float
        print ustring
        uval=[float(item) for item in ustring]
    except Exception:
        parser.print_usage()
        sys.exit(2)

    if (uval):
        print runnum
        print uval
        param = Param(uval)
        W = 4000 #[m], domain width
        L = 4000 #[m], domain length
        dx = 1000 #[m], x dist between leaky wells
        dy = 1000 #[m], y dist between leaky wells
        mx = 200
        my = 200 
        geom = Geom(width=W, length=L, dx=dx, dy=dy, mx=mx, my=my)
        
        Tmax = 365.25*1
        
        param.setSimulationTime(Tmax)     
        param.setLeakWells(geom)
        #set the injection rate here
        param.setInjWells(Minj=uval[-1], geom=geom)
        param.setMonWells(geom)
            
        mytest = Assembler(param)
        mytest.iterate(isPlot=False)
        monP,monh,cumMassCO2, cumMassBrine = mytest.getResults()        
        pkl.dump([monP,monh,cumMassCO2, cumMassBrine], open('time%s_%s.dump'% (Tmax, runnum), 'wb'))
                          
if __name__ == "__main__":
   main(sys.argv[1:])