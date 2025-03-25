#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:44:57 2024

@author: jachym
"""

import scipy.io
import os
import pickle
import numpy as np

# vstup se nacita z input() pri spousteni

def get_id(path):
    # urceni id pro soubor podle simulovanych parametru v nazvu souboru
    name = path.split("/")[-1]
    name_array = list(name)
    integer_id = name_array[0] + name_array[2] + name_array[3] + name_array[5] + name_array[6] + name_array[9] + name_array[10] + name_array[12] + name_array[13]
    return integer_id


def gammaFromFile(path1, path2):
    
    # nacte soubor do pythonu
    mat = scipy.io.loadmat(path1)
    mat.update(scipy.io.loadmat(path2))
        
    # urceni id pro soubor podle simulovanych parametru v nazvu souboru
    name = path1.split("/")[-1]
    type = name[0]
    integer_id = get_id(path1)
    
    # nacteni dulezitych promenych do dict
    wanted = []
    if type == "1":
        wanted = ["dens01","dens02","dt","edgeenergyflux01","edgeenergyflux02","edgeflux01","edgeflux02",
                  "vxav01","vxav02","vyav01","vyav02","vzav01","vzav02","Epar","Eperp","Potav",
                  "mksTe","mksB","tau","mksn0","mksmainionm","mksmainionq","Npc","alphaxz","alphayz"]
    if type == "2":
        wanted = ["dens01","dens02","dens03","dt","edgeenergyflux01","edgeenergyflux02","edgeenergyflux03","edgeflux01",
                  "edgeflux02","edgeflux03","vxav01","vxav02","vxav03","vyav01","vyav02","vyav03","vzav01","vzav02","vzav03",
                  "Epar","Eperp","Potav","mksTe","mksB","tau","mksn0","mksmainionm","mksmainionq","Npc","alphaxz","alphayz"]
    dic_arr = [["id", integer_id]]
    for key in wanted:
        dic_arr.append([key, mat[key]])
    PlasmaDict = dict(dic_arr)
    
    
    
    
    
    
    # nacitani parametru do promennych
    Te = PlasmaDict["mksTe"]
    B = PlasmaDict["mksB"]
    tau = PlasmaDict["tau"]
    n0 = PlasmaDict["mksn0"]
    mi = PlasmaDict["mksmainionm"]
    qi = PlasmaDict["mksmainionq"]
    # konstanty
    eps0 = 8.8541e-12
    kB = 1.3806e-23
    e = 1.6021e-19
    amu = 1.6e-27
    me = 9.109e-31

    Te = Te * e / kB
    qi = qi*e
    mi = mi*amu
    
    # promenne plazmatu
    lamD = np.sqrt((eps0*kB*Te)/(n0*(e**2)))

    omeC = (abs(qi)*B)/mi

    # koncentrace plynu z id souboru
    HeC = int(PlasmaDict["id"][5]+PlasmaDict["id"][6])/100
    HC = int(PlasmaDict["id"][3]+PlasmaDict["id"][4])/100
    SPot = -int(PlasmaDict["id"][7] + PlasmaDict["id"][8])
    
    # denormalizace podle PDF navodu
    Q01 = PlasmaDict["edgeenergyflux01"] * n0 * mi * (lamD**3) * (omeC**3)
    Q02 = PlasmaDict["edgeenergyflux02"] * n0 * mi * (lamD**3) * (omeC**3)
    Gamma01 = PlasmaDict["edgeflux01"] * n0 * lamD * omeC
    Gamma02 = PlasmaDict["edgeflux02"] * n0 * lamD * omeC
    
    # pro soubory s dvemi casticemi
    if type == '2':
        Q03 = PlasmaDict["edgeenergyflux03"] * n0 * mi * (lamD**3) * (omeC**3)
        Gamma03 = PlasmaDict["edgeflux03"] * n0 * lamD * omeC
        mi = (mi*HeC)+(amu*HC)
        
    # pocitani gamma z denormalizovanych velicin
    def calculate_gamma(border, plot=False):

        Gamma = Gamma01+Gamma02+1e-20
        if type == '2': 
            Gamma = Gamma01+Gamma02+Gamma03+1e-20

        gamma = (Q01+Q02)/(Gamma * kB * Te)
        if type == '2':
            gamma = (Q01+Q02+Q03)/(Gamma * kB * Te)

        return np.median(gamma[border])
    
    # predikce gammy podle Stangebyho
    def gamma(tau, m_e, m_i, delta_e, ZEff):
        a = (2.5*tau)/ZEff + 2/(1-delta_e)
        b = 0.5*np.log((2*np.pi*(m_e/m_i))*(ZEff+tau)*(1-delta_e)**(-2))
        return a-b
    
    return [type, tau, HC, HeC, SPot, calculate_gamma(30), 
            gamma(tau, me, mi, 0, (HC + HeC*2))]
    
    

def gammaFromDir(inp):
    
    # prochazi vsechny soubory ve slozce a vola funkci gammaFromFile
    result_list = []
    ids_done = []   # pole ke kontrole t a o souboru
    for filename1 in os.listdir(inp):
        print(filename1)
        if get_id(inp+'/'+filename1) not in ids_done:
            for filename2 in os.listdir(inp):
                if filename1 != filename2 and get_id(inp+'/'+filename1) == get_id(inp+'/'+filename2):
                    result_list.append(gammaFromFile(inp+'/'+filename1, inp+'/'+filename2))
                    ids_done.append(get_id(inp+'/'+filename1))
    result_list.sort(key=lambda x: x[4])
    result_list.sort(key=lambda x: x[2])
    result_list.sort(key=lambda x: x[0])
    return result_list
    


if __name__ == "__main__":
    print('gamma calculation for THEH simulations')
    print('--------------------------------------')
    result = gammaFromDir(input('Path to matlab files (directory):  '))
    with open("gamma.pickle", "wb") as file:
        pickle.dump(result, file)
    print(result)

# unused definitions
    # potav = PlasmaDict["Potav"] * Te
    # # denormalizace rychlosti
    # normalizeV = lambda v, rho: (v*lamD*omeC)/(Npc*rho)
    # vx01norm = normalizeV(PlasmaDict["vxav01"], PlasmaDict["dens01"]+1e-20)   # proti deleni nulou
    # vy01norm = normalizeV(PlasmaDict["vyav01"], PlasmaDict["dens01"]+1e-20)
    # vz01norm = normalizeV(PlasmaDict["vzav01"], PlasmaDict["dens01"]+1e-20)
    # vx02norm = normalizeV(PlasmaDict["vxav02"], PlasmaDict["dens02"]+1e-20)
    # vy02norm = normalizeV(PlasmaDict["vyav02"], PlasmaDict["dens02"]+1e-20)
    # vz02norm = normalizeV(PlasmaDict["vzav02"], PlasmaDict["dens02"]+1e-20)density03 = PlasmaDict["dens03"] * n0
    # vx03norm = normalizeV(PlasmaDict["vxav03"], PlasmaDict["dens03"]+1e-20)
    # vy03norm = normalizeV(PlasmaDict["vyav03"], PlasmaDict["dens03"]+1e-20)
    # vz03norm = normalizeV(PlasmaDict["vzav03"], PlasmaDict["dens03"]+1e-20)
    # density01 = PlasmaDict["dens01"] * n0
    # density02 = PlasmaDict["dens02"] * n0
    # dt = PlasmaDict["dt"]/omeC    
    # Epar = PlasmaDict["Epar"] * Te / lamD
    #     Eperp = PlasmaDict["Eperp"] * Te / lamD
    #     phi = PlasmaDict["Potav"] * Te
    # cs0 = np.sqrt((kB*Te*(1+tau))/mi)
    # rL = (mi*cs)/(abs(qi)*B)
        # # proc tady mam if?
        # if PlasmaDict["alphaxz"][0][0] == 90 or PlasmaDict["alphaxz"][0][0] == -90:
        #     alphaxz = np.radians(PlasmaDict["alphaxz"])   # to radians
        # else:
        #     alphaxz = PlasmaDict["alphaxz"]

        # if PlasmaDict["alphayz"][0][0] == 90 or PlasmaDict["alphayz"][0][0] == -90:
        #     alphayz = np.radians(PlasmaDict["alphayz"])   # to radians
        # else:
    #     alphayz = PlasmaDict["alphayz"]
        
    # Npc = PlasmaDict["Npc"]
    
    # cs = np.sqrt((kB*Te)/mi)
