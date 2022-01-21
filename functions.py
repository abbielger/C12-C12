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

def convertMevTokeV(E):
    a = E*1000
    return a

def twoPiEta(EkeV):
    a = 31.29 * 6 * 6 *np.sqrt((144/24)/EkeV)
    return a

def convertModStoCrossSection(modS,E,pieta):
    a = (modS/E)*np.exp(-pieta-0.46*E)
    return a

def convertBarnsToCentimeters(sigma):
    a = sigma*10**(-24)
    return a

def fitStoppingPower(Ecm):
    srimLoc = ("/Users/abbielger/Desktop/Research/srimdata.csv")

    srimdataframe = read_csv(srimLoc,header = 0)
    srimdata = srimdataframe.values
    srimEcm, stoppingPower = srimdata[:,0], srimdata[:,3]
    srimEcmMeV = srimEcm * 10**(-3)
    stoppingPowerMeV = stoppingPower*10**(-21)

    def poly(x, a, b, c, d, e, f, g):
        return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + (f * x**6)  + g

    popt, _ = curve_fit(poly, srimEcmMeV, stoppingPowerMeV)
    a,b,c,d,e,f,g = popt
    #Fit the data for the stopping power
    fittedStoppingPower = (a * Ecm) + (b * Ecm**2) + (c * Ecm**3) + (d * Ecm**4) + (e * Ecm**5) + (f * Ecm**6)  + g

    return fittedStoppingPower

def calculateYield(sigma, e, Ecm):
    I = np.divide(sigma,e)
    y = scipy.integrate.simps(I,Ecm)

    return I,y

def getData():
    locZickesfoose = ("/Users/abbielger/Desktop/Research/zickefoose.xls")
    wbZickesfoose = xlrd.open_workbook(locZickesfoose)
    sheetZickesfoose = wbZickesfoose.sheet_by_index(0)

    EZickesfoose = [0]*(sheetZickesfoose.nrows-1)
    yieldZickesfoose = [0]*(sheetZickesfoose.nrows-1)


    for i in range(sheetZickesfoose.nrows-1):
        EZickesfoose[i] = (sheetZickesfoose.cell_value(i+1,0))

    for j in range(sheetZickesfoose.nrows-1):
        yieldZickesfoose[j] = (sheetZickesfoose.cell_value(j+1,1))


    locFriend = ("/Users/abbielger/Desktop/Research/friendTTYdata.xls")
    wbFriend = xlrd.open_workbook(locFriend)
    sheetFriend = wbFriend.sheet_by_index(0)

    Efriend = [0]*(sheetFriend.nrows-1)
    yieldFriend = [0]*(sheetFriend.nrows-1)


    for i in range(sheetFriend.nrows-1):
        Efriend[i] = (sheetFriend.cell_value(i,0))

    for j in range(sheetFriend.nrows-1):
        yieldFriend[j] = (sheetFriend.cell_value(j,1))

    return EZickesfoose, yieldZickesfoose, Efriend, yieldFriend
