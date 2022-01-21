import xlrd
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import integrate
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
import pandas as pd
import os
import string
from matplotlib.ticker import FormatStrFormatter

#Do you want graphs on or off? 0 = off, 1 = on, 2 = only recent plots
graphToggle = 0
#Open raw data files
fileLoc = ("/Users/abbielger/Desktop/Research/rawdata.csv")
CF88Loc = ("/Users/abbielger/Desktop/Research/CF88C12C12.csv")
alphaLoc = ("/Users/abbielger/Desktop/Research/alpha.csv")
dataframe = read_csv(fileLoc,header = 0)
dataframeCF88 = read_csv(CF88Loc,header = 0)
dataframeAlpha = read_csv(alphaLoc, header=0)
data= dataframe.values
CF88Data = dataframeCF88.values
alphaData = dataframeAlpha.values

#Create CF88 data frames of temperature and reaction rates
CF88temp, CF88RR = CF88Data[:,0], CF88Data[:,1]
#Create data frames for the alpha channel
E_a0, modSa0, E_a1, modSa1 = alphaData[:,0], alphaData[:,1], alphaData[:,2], alphaData[:,3]
Ecm, modSp0, modSp1, modSp0p1 = data[:,0], data[:,1], data[:,2], data[:,3]
#Ecm= data[:,0]

#modSp0 = np.full(len(Ecm),3e16)
#modSp1 = np.full(len(Ecm),3e16)
#modSp0p1 = np.full(len(Ecm),3e16)
#Find the cross-section of p0, p1 and p0+p1 by converting from modified S-factors
#to cross-section in barns

import functions
modSp0p1Extend = [3E+16]*37
4.56E+17
#modSp0p1Extend2 = [5.5E+16]*10
#modSp0p1Extend3 = [5.5E+16]*10
#modSp0p1Extend = [9.54361111111111E+15,	9.49277777777778E+15,	9.44194444444444E+15,	9.39111111111111E+15,	9.34027777777778E+15,	9.28944444444444E+15,	9.23861111111111E+15,	9.18777777777778E+15,	9.13694444444444E+15,	9.08611111111111E+15,	9.03527777777778E+15,	8.98444444444444E+15,	8.93361111111111E+15,	8.88277777777778E+15,	8.83194444444444E+15,	8.78111111111111E+15,	8.73027777777778E+15,	8.67944444444444E+15,	8.62861111111111E+15,	8.57777777777778E+15,	8.52694444444444E+15,	8.47611111111111E+15,	8.42527777777778E+15,	8.37444444444444E+15,	8.32361111111111E+15,	8.27277777777778E+15,	8.22194444444444E+15,	8.17111111111111E+15,	8.12027777777778E+15,	8.06944444444444E+15]
#modSp0p1Extend = np.append(modSp0p1Extend1,modSp0p1Extend2)
#modSp0p1Extend = np.append(modSp0p1Extend,modSp0p1Extend3)
modSp0p1 = np.append(modSp0p1,modSp0p1Extend)
modSp0Extend = [4.17E+15]*37
modSp0 = np.append(modSp0,modSp0Extend)
modSp1Extend = [5.38E+15]*37
modSa0Extend = [15E+15]*37
modSa1Extend = [18E+15]*37
modSp1 = np.append(modSp1,modSp1Extend)
modSa0 = np.append(modSa0,modSa0Extend)
modSa1 = np.append(modSa1,modSa1Extend)

extendE = [2.8, 	2.9, 	3, 	3.1, 	3.2, 	3.3, 	3.4, 	3.5, 	3.6, 	3.7, 	3.8, 	3.9, 	4, 	4.1, 	4.2, 	4.3, 	4.4, 	4.5, 	4.6, 	4.7, 	4.8, 	4.9, 	5, 	5.1, 	5.2, 	5.3, 	5.4, 	5.5, 	5.6, 	5.7, 5.8 , 5.9, 6.1, 6.2,6.3,6.4,6.5]
Ecm = np.append(Ecm,extendE)
E_a0 = np.append(E_a0,extendE)
E_a1 = np.append(E_a1,extendE)

#print(len(extendE))
#Ecm from the dataframe is in MeV
EcmkeV,Ea0keV,Ea1keV = functions.convertMevTokeV(Ecm),functions.convertMevTokeV(E_a0), functions.convertMevTokeV(E_a1)
twopieta,twopieta_a0,twopieta_a1  = functions.twoPiEta(EcmkeV),functions.twoPiEta(Ea0keV),functions.twoPiEta(Ea1keV)
sigmap0,sigmap1,sigmap0p1 = functions.convertModStoCrossSection(modSp0,Ecm,twopieta),functions.convertModStoCrossSection(modSp1,Ecm,twopieta),functions.convertModStoCrossSection(modSp0p1,Ecm,twopieta)
sigmaa0,sigmaa1 = functions.convertModStoCrossSection(modSa0,E_a0,twopieta_a0),functions.convertModStoCrossSection(modSa1,E_a1,twopieta_a1)

print(E_a0)

#Fitting alpha verses energy so I can use the energy values of the proton channel
"""fig, ax = plt.subplots()
ax.plot(E_a1, sigmaa0,color= '#26532B',linestyle = "solid", label = '$\mathregular{\sigma_{p0}}$')
ax.set_yscale('log')
plt.show()
"""




#Sigma is in barns, so we convert to centimeters for stopping power calculations
sigmap0cm = functions.convertBarnsToCentimeters(sigmap0)
sigmap1cm = functions.convertBarnsToCentimeters(sigmap1)
sigmap0p1cm = functions.convertBarnsToCentimeters(sigmap0p1)
sigmaa0cm = functions.convertBarnsToCentimeters(sigmaa0)
sigmaa1cm = functions.convertBarnsToCentimeters(sigmaa1)

#3.342509998425694e-16
#2.4754786900865917e-17
#1.2349379884113887e-18
#6.180403073623777e-19

fittedStoppingPower = functions.fitStoppingPower(Ecm)
fittedStoppingPowera0 = functions.fitStoppingPower(E_a0)
fittedStoppingPowera1 = functions.fitStoppingPower(E_a1)


if graphToggle == 1:
#Plot cross section verses the energy
    fig, ax = plt.subplots()
    ax.plot(Ecm, sigmap0,color= '#26532B',linestyle = "solid", label = '$\mathregular{\sigma_{p0}}$')
    ax.plot(Ecm, sigmap1, color= '#A3D0AA', linestyle = "solid",label = '$\mathregular{\sigma_{p1}}$')
    ax.plot(Ecm, sigmap0p1,color=  '#2C8C99',linestyle = "solid", label = '$\mathregular{\sigma_{p0+p1}}$')
    ax.plot(E_a0, sigmaa0, color= '#B95F89', linestyle = "solid",label = r'$\mathregular{\sigma_ {\alpha 0}}$')
    ax.plot(E_a1, sigmaa1, color= '#DD6E42', linestyle = "solid",label = r'$\mathregular{\sigma_{\alpha 1}}$')

    #plt.xlim([0.75,2])
    #plt.ylim([0,2e19])
    plt.xlabel("$\mathregular{E_{cm}}$ [MeV]")
    plt.ylabel("Cross-Section [barns]")
    #Comment out the log if you'd prefer to ignore Coulomb Barrier
    ax.set_yscale('log')
    legend = ax.legend(loc='lower right')
    plt.title("Center of Mass Energy vs. Cross-Section")
    plt.show()

#Plot astrophysical s-factor verses the energy
    fig, ax = plt.subplots()
    ax.plot(Ecm, modSp0,color= '#26532B',linestyle = "solid", label = 'S$^*_{p0}$')
    ax.plot(Ecm, modSp1, color= '#A3D0AA', linestyle = "solid",label = 'S$^*_{p1}$')
    ax.plot(Ecm, modSp0p1,color=  '#2C8C99',linestyle = "solid", label = 'S$^*_{p0p1}$')
    ax.plot(E_a0, modSa0, color= '#B95F89', linestyle = "solid",label = 'S$^*_{a0}$')
    ax.plot(E_a1, modSa1, color= '#DD6E42', linestyle = "solid",label = 'S$^*_{a1}$')

    #plt.xlim([0.75,2])
    #plt.ylim([0,2e19])
    plt.xlabel("$\mathregular{E_{cm}}$ [MeV]")
    plt.ylabel("Astrophysical S-Factor")
    #Comment out the log if you'd prefer to ignore Coulomb Barrier
    ax.set_yscale('log')
    legend = ax.legend(loc='upper right')
    plt.title("Center of Mass Energy vs. Cross-Section")
    plt.show()
maxEcm = 2.7
index = np.where(Ecm == 2.7)


yp0,yieldp0 = functions.calculateYield(sigmap0cm,fittedStoppingPower,Ecm)
yp1,yieldp1 = functions.calculateYield(sigmap1cm,fittedStoppingPower,Ecm)
yp0p1,yieldp0p1 = functions.calculateYield(sigmap0p1cm,fittedStoppingPower,Ecm)
ya0,yielda0 = functions.calculateYield(sigmaa0cm,fittedStoppingPowera0,E_a0)
ya1,yielda1 = functions.calculateYield(sigmaa1cm,fittedStoppingPowera1,E_a1)

#yieldp0,yieldarrayp0 = file.calculateYield(sigmap0cm,fittedStoppingPower,Ecm)

#Find the yields for all energies

yieldarrayp0 = [0]*(len(Ecm))
x = len(Ecm)
for i in range(x):
    if i == 0:
        yieldarrayp0[i] = 0
    else:
        yieldarrayp0[i] = scipy.integrate.simps(yp0[0:i],Ecm[0:i])

#sigma,e = sigmap0cm,fittedStoppingPower

yieldarrayp0p1 = [0]*(len(Ecm))
for i in range(x):
    if i == 0:
        yieldarrayp0p1[i] = 0
    else:
        yieldarrayp0p1[i] = scipy.integrate.simps(yp0p1[0:i],Ecm[0:i])
#7.504515876328869e-12
############### Getting experimental data and other data sources
EZickesfoose, yieldZickesfoose, Efriend, yieldFriend = functions.getData()

EZickesfoose = np.flip(EZickesfoose)
yieldZickesfoose = np.flip(yieldZickesfoose)
yieldFriend = np.divide(yieldFriend,8)
z = np.divide(yieldarrayp0p1,8)

if graphToggle == 1:
    fig, ax = plt.subplots()
    ax.plot(Ecm, z, 'bo', markersize = 4,label = '$\mathregular{\sigma_{p0+p1}}$')
    #ax.plot(x, yieldarrayp0, 'r+', markersize = 3,label = '$\mathregular{\sigma_{p0}}$')
    ax.plot(Efriend,yieldFriend, 'r+', markersize=3,  label = 'Friends Integral')
    ax.plot(EZickesfoose,yieldZickesfoose, 'y*',  label = 'Zickefoose Data')
    plt.xlim([2,2.7])
    #plt.ylim([1e-31,1e-28])
    plt.xlabel("Maximum $\mathregular{E_{cm}}$ [MeV]")
    plt.ylabel("Yield")
    #Comment out the log if you'd prefer to ignore Coulomb Barrier
    ax.set_yscale('log')
    legend = ax.legend(loc='upper right')
    plt.title("Yield vs Maximum Energy")
    plt.show()

#constants
NA = 6.02214e23 #avogadros constant [1/mol]
mu = np.sqrt(24/144) #reduced mass for 2 carbon 12 atoms
Tlist = (0.11,.12, .13, .14, .15, .16, .18, .2, .25, .3, .35, .4, .45, .5,.6,.7,.8,.9,1,1.25,1.5,1.75,2,2.5, 3,3.5, 4,5,6,7,8,9,10)
T= np.array(Tlist)

y = [0]*len(T)
ya0 = [0]*len(T)
ya1 = [0]*len(T)
yLow = [0]*len(T)
RR=[0]*len(T)
RRa0=[0]*len(T)
RRa1=[0]*len(T)
RRLow=[0]*len(T)
#sigmap0 = np.divide(sigmap0,8)
#sigmap1 = np.divide(sigmap1,8)


for i in range(len(T)):
    y = Ecm*sigmap0p1*np.exp(-11.605*Ecm/T[i])
    ya0 = E_a0[0:178]*sigmaa0[0:178]*np.exp(-11.605*E_a0[0:178]/T[i])
    ya1 = E_a1*sigmaa1*np.exp(-11.605*E_a1/T[i])
    #RR[i] = np.add((1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(y,Ecm),(1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(h,Ecm))
    RR[i] = (1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(y,Ecm)
    RRa0[i] = (1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(ya0[0:178],E_a0[0:178])
    RRa1[i] = (1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(ya1,E_a1)


#print(test)
#Find ratios of RR from Ecm=0 to Ecm=some value over RR from Ecm=0 to Ecm = 2.7
#Currently just for T=0.6
#Already have total RR -> RR[14] is RR at T=0.6
RRLowBound=[0] * len(Ecm)
RRLowBound2=[0] * len(Ecm)
yLowBound =[0] * len(Ecm)

if graphToggle == 1:
    fig, ax = plt.subplots()
    for k in range(len(T)):
        for i in range(len(Ecm)):
            yLowBound = Ecm[0:i+1]*sigmap0p1[0:i+1]*np.exp(-11.605*Ecm[0:i+1]/T[k])
            #RR[i] = np.add((1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(y,Ecm),(1/2)*(3.7318*(10**10)/(T[i])**(3/2))*mu*scipy.integrate.simps(h,Ecm))
            RRLowBound[i] = (1/2)*(3.7318*(10**10)/(T[k])**(3/2))*mu*scipy.integrate.simps(yLowBound,Ecm[0:i+1])
        frac = RRLowBound/RR[k]

        legend = ax.legend(loc='lower right', prop={'size': 10})
        ax.plot(Ecm, frac,label=T[k])
    plt.xlabel('Energy')
    plt.ylabel('Integral from 0 to E/Total Integral')
    plt.title('Ratio')
    plt.show()
"""
frac = np.divide(RRLow,RR)
fig, ax = plt.subplots()
ax.plot(Ecm,frac)
plt.xlim(2,2.7)
plt.ylim(.8,1.1)
plt.show()
"""
someTemps = np.linspace(0.1,0.6,20)
EcmLow = Ecm[0:120]
sigmap0p1Low = sigmap0p1[0:120]



y = [0]*len(someTemps)

for i in range(len(someTemps)):
    y[i] = Ecm*sigmap0p1*np.exp(-11.605*Ecm/T[i])

frac = np.divide(z,y)


if graphToggle == 1:
    fig, ax = plt.subplots()
    for i in range(len(someTemps)):
        ax.plot(Ecm, y[i], label = round(T[i],2))
    legend = ax.legend(loc='upper left')
    plt.xlabel("Ecm")
    plt.ylabel("Integrand")
    #plt.ylim([0,0.3e-9])
    #plt.xlim([2,2.75])
    ax.set_yscale('log')
    plt.show()
#print(Ecm[0:120])

"""
f = np.linspace(0.5,1,30)
for i in someTemps:
    m = Ecm*sigmap0p1*np.exp(-Ecm/(11.605*f[i]))
    fig, ax = plt.subplots()
    ax.plot(Ecm, m)
    plt.show()
"""
RR = np.add(RR,RRa0)
RR = np.add(RR,RRa1)

RR3 = np.multiply(RR,3.333)
RRa = np.multiply(RR,2)
RRcorrected = np.multiply(RR,3.75)
RR8 = np.divide(RR,8)
#ydiff = np.subtract(CF88RR,RR[1:])
#ydiff=np.abs(ydiff)
#print(sum(ydiff))
#print(RR)
#print(RRa0)

ratioRR = np.divide(RR,CF88RR)
ratioRR3 = np.divide(RR3,CF88RR)
ratioRRa = np.divide(RRa,CF88RR)
ratioRR8 = np.divide(RR8,CF88RR)
ratioRRcorrected = np.divide(RRcorrected,CF88RR)

if graphToggle == 1 or 2:
    #Plotting both our data vs temperature and CF88 vs temperature
    fig, ax = plt.subplots()
    ax.plot(T, RRa, 'ys', markersize = 5,label = 'Reaction Rate Scaled by 2')
    ax.plot(T, RR3, 'b+', markersize = 4,label = 'Reaction Rate Scaled by 3.333')
    ax.plot(T, RR, 'ro', markersize = 3,label = 'Reaction Rate')
    ax.plot(T, RR8, 'g*', markersize = 2,label = 'Reaction Rate Scaled by 1/8')
    ax.plot(CF88temp, CF88RR, 'g*', markersize = 3,label = 'CF88 Data')
    ax.set_yscale('log')
    legend = ax.legend(loc='lower right')
    plt.xlabel("Temperature [GK]")
    plt.ylabel("Reaction Rate [$\sigma v$]")
    plt.title("Reaction Rate vs Temperature")
    plt.show()

    #Plot ratio of CF88/Ours vs temperature
    fig, ax = plt.subplots()
    ax.plot(CF88temp[0:25], ratioRR[0:25], 'k-', markersize = 4,label = 'Ratio')
    ax.plot(CF88temp[0:25], ratioRR3[0:25], 'g-', markersize = 4,label = 'Ratio scale of 3.333')
    ax.plot(CF88temp[0:25], ratioRRcorrected[0:25], 'm-', markersize = 4,label = 'Ratio scale of 3.75')
    ax.plot(CF88temp[0:25], ratioRRa[0:25], 'r-', markersize = 4,label = 'Ratio scale of 2')
    ax.plot(CF88temp[0:25], ratioRR8[0:25], 'y-', markersize = 4,label = 'Ratio scale of 1/8')
    plt.axvspan(0.4, 0.5, color='purple', alpha=0.3,label='Superbursts')
    plt.axvspan(0.6, 1.2, color='red', alpha=0.3,label='Hydrostatic Burning')
    plt.axvspan(1.8, 2.6, color='green', alpha=0.3,label='Explosive Burning')
    ax.axhline(1, color='red', linestyle='--')
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.tick_params(axis='x', which='minor')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
    legend = ax.legend(loc='upper right',prop={'size': 15})
    plt.ylim(.1,1000)
    plt.xlabel("Temperature [GK]")
    plt.ylabel("R$_{THM}$/R$_{CF88}$")
    plt.title("$^{12}$C + $^{12}$C Reaction Rate Ratio")
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(CF88temp[14:23], RR[14:23], 'k-', markersize = 4,label = 'THM')
    ax.plot(CF88temp[14:23], CF88RR[14:23], 'b-', markersize = 4,label = 'CF88')
    #ax.plot(CF88temp[0:25], ratioRR3[0:25], 'g+', markersize = 4,label = 'Ratio scale of 3.333')
    #ax.plot(CF88temp[0:25], ratioRRa[0:25], 'r-', markersize = 4,label = 'Ratio scale of 2')
    #ax.plot(CF88temp[0:25], ratioRR8[0:25], 'y-', markersize = 4,label = 'Ratio scale of 1/8')
    ax.set_yscale('log')
    #plt.xlim(10E-15,10E-3)
    legend = ax.legend(loc='upper left')
    #plt.ylim(.1,1000)
    plt.xlabel("Temperature [GK]")
    plt.ylabel("Reaction Rate (cm$^3$mole$^{-1}$sec$^{-1}$)")
    plt.title("Reaction Rates")
    plt.show()

#Save note file


#6.079037556095241e-07
#1.1926743325931574e-07
RRvsT = np.vstack((T,CF88RR))

#print(RRvsT)
header = 'T          RR'
np.savetxt('/Users/abbielger/Desktop/Research/RRvsT',RRvsT, fmt="%e", header=header)
