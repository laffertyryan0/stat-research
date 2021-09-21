from scipy import special
import numpy as np
import os
import shutil
import pandas as pd
import json

#I_j is integral from -inf to t of x^j e^(-1/2 x^2)
def gen_I(t,max_j = 10):
    if(t == np.inf):
        Ijm2 = np.sqrt(2*np.pi)
        Ijm1 = 0       
    else:
        Ijm2 = np.sqrt(np.pi/2)*(special.erf(t/np.sqrt(2))+1)
        Ijm1 = -np.exp(-t**2/2)
    yield Ijm2
    yield Ijm1
    j = 2
    while(j<=max_j):
        Ij = (0 if t == np.inf else -t**(j-1)) * np.exp(-t**2/2) + (j-1)*Ijm2
        yield Ij
        j = j + 1
        Ijm2 = Ijm1
        Ijm1 = Ij


#need to integrate from 0 to t of (x+a)^k e^(-x^2/2) where k is even

def F(a,k):
    norm = sum(special.binom(k,j)*a**(k-j)*I_j
                   for j,I_j in enumerate(gen_I(np.inf,max_j = k)))
    def outFunction(t):
        return sum(special.binom(k,j)*a**(k-j)*I_j
                   for j,I_j in enumerate(gen_I(t,max_j = k)))/norm
    return outFunction

def tabulateF(kMax = 20,
              aStep = .01,
              aMin = -10,
              aMax = 10,
              xMin = -5,
              xMax = 5,
              xStep = .01):
    args = locals()
    with open("./json/configTable.json",mode = "w") as file:
        file.write(str(args).replace("'","\""))
    dirname = "table"
    shutil.rmtree("table") if os.path.isdir(dirname) else 0
    os.mkdir(dirname)
    for k in range(2,kMax+2,2):
        filepath = dirname + "/tablekequals" + str(k) + ".csv"
        table = []
        a_list = []
        for a in np.arange(aMin,aMax + aStep,aStep):
            a_list = a_list + ["a="+str(np.round(a,5))]
            table.append([np.round(F(a,k)(x),5) for x in np.arange(xMin,xMax + xStep,xStep)])
            x_list = ["x="+str(np.round(x,5)) for x in np.arange(xMin,xMax + xStep,xStep)]
        pd.DataFrame(data = table,
                     index = a_list,
                     columns = x_list).to_csv(filepath,header = True,index = True)

tabulateF(kMax = 6,aStep = .1,xStep = .1)
