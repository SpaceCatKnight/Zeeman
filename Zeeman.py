# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:22:09 2017

@author: HP
"""


import numpy as np
import matplotlib.pyplot as plt
 

# Konstanten


h = 6.626*10**(-34)         #Planck Konstante
c = 2.9979*10**8            #Lichtgeschw.
mb = 9.2740*10**(-24)       #Bohrsches Magneton
sa = np.pi*0.01002**2       #Querschnittsfläche Spule
sn = 127                    #Windungen Spule 
n = 1.5172                  #Brechungsindex Lummerplatte
dn = -0.0025/(40*10**(-9))  #Ableitung Brechungsindex bei lambda0
lambda0 = 540*10**(-9)      #Unverschobene Wellenlänge
d = 3.213*10**-3            #Dicke der Lummerplatte
V = 4.372                   #Durchschnittlich gemessene Voltzahl
S = 230.42                  #Durchschnittliche gemmessener Abstand Drehpunkt-Mikrometerschraube
#m=0                        #Ordnung von Maximum von innen nach aussen 
                      
                      
                      
                      
# Funktionen 

def read_from_file(filename):                   #Liest Daten ein
    return np.loadtxt(filename, skiprows=1)

def Mittelwerte(dat):                           #Gibt Liste mit gemittelten Daten aus
    dat0 = dat-np.mean(dat)
    dat0m = []
    for i in range(0,18):
        dat0m.append(np.abs(np.mean(dat0[:,i])))
    return dat0m
    

def Abstand(dat, i):                            #Gibt gemittelten Abstand zwischen Maxima auf beiden Seiten aus
    A = 0
    dat0m = Mittelwerte(dat)
    A= dat0m[i]+dat0m[-(i+1)]  
    return A 
    
def Theta(dat,S,i):                             #Gibt Winkel des Maximums aus
    A = Abstand(dat, i)
    theta = np.arctan(A/(2*S)) 
    return theta    
    
def Maximumsordnung(dat, d, n, S, m, lambda0): #Gibt tatsächliche Ordnung des Maximums aus
    i = 3*m+1
    theta = Theta(dat,S,i)
    M = (2*d/lambda0)*np.sqrt(n**2-1+np.sin(theta)**2)
    return M 
        
def Frequenzunterschied(dat, d, lambda0, n, S, m, c, dn): #gibt Frequenzunterschied von unverschobenem zu verschobenen Übergängen
    
    ip = 3*m
    im = 3*m+2
    
    M = Maximumsordnung(dat, d, n, S, m, lambda0)
    thetap = Theta(dat,S,ip)
    thetam = Theta(dat,S,im)
    deltaf = (-c/lambda0**2)*(np.sin(thetam)**2-np.sin(thetap)**2)/(lambda0*M**2/d**2-4*n*dn)
    return deltaf
    
    

def Magnetfeld(V, sa, sn):                      #Berechnet Magnetfeldstärke
    B = V/(314.16*sa*sn)
    return B 

def Landefaktor(dat, d, lambda0, n, S, m, c, dn, V, sa, sn, mb, h): #Berechnet Landéfaktor
    deltaf = Frequenzunterschied(dat, d, lambda0, n, S, m, c, dn)
    B = Magnetfeld(V,sa,sn)
    gj = h*deltaf/(mb*B)
    return gj

   
#Test

daty = read_from_file('Gelb.txt')
datbg = read_from_file('BlauGruen.txt')

print (Landefaktor(daty, d, lambda0, n, S, 0, c, dn, V, sa, sn, mb, h)+Landefaktor(daty, d, lambda0, n, S, 1, c, dn, V, sa, sn, mb, h)+Landefaktor(daty, d, lambda0, n, S, 2, c, dn, V, sa, sn, mb, h))/3

print (Landefaktor(datbg, d, lambda0, n, S, 0, c, dn, V, sa, sn, mb, h)+Landefaktor(datbg, d, lambda0, n, S, 1, c, dn, V, sa, sn, mb, h)+Landefaktor(datbg, d, lambda0, n, S, 2, c, dn, V, sa, sn, mb, h))/3
