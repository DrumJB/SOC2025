#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 13:47:51 2023

@author: jachym
"""

import pickle
import numpy as np

class DenormalizedSimulation:
    
    def __init__(self, output, type):
        
        import numpy as np
        # loading the plasma parameters from given file
        Te = output["mksTe"]
        B = output["mksB"]
        tau = output["tau"]
        self.tau = tau
        n0 = output["mksn0"]
        self.n0=n0
        mi = output["mksmainionm"]
        qi = output["mksmainionq"]
        Npc = output["Npc"]
        # loading natural constants
        eps0 = 8.8541e-12
        kB = 1.3806e-23
        self.kB=kB
        e = 1.6021e-19
        amu = 1.6e-27
    
        Te = Te * e / kB
        qi = qi*e
        mi = mi*amu
        self.Te = Te
        self.mi = mi
        
        lamD = np.sqrt((eps0*kB*Te)/(n0*(e**2)))
        cs = np.sqrt((kB*Te)/mi)
        self.cs0 = np.sqrt((kB*Te*(1+tau))/mi)
        
        omeC = (abs(qi)*B)/mi
        rL = (mi*cs)/(abs(qi)*B)

        if output["alphaxz"][0][0] == 90 or output["alphaxz"][0][0] == -90:
            self.alphaxz = np.radians(output["alphaxz"])   # to radians
        else:
            self.alphaxz = output["alphaxz"]
      
        if output["alphayz"][0][0] == 90 or output["alphayz"][0][0] == -90:
            self.alphayz = np.radians(output["alphayz"])   # to radians
        else:
            self.alphayz = output["alphayz"]

        # retriving helium and hydrogen concentration from the id:
        self.HeC = int(output["id"][5]+output["id"][6])/100
        self.HC = int(output["id"][3]+output["id"][4])/100
        self.SPot = -int(output["id"][7] + output["id"][8])
        
        self.density01 = output["dens01"] * n0
        self.density02 = output["dens02"] * n0
        dt = output["dt"]/omeC
        self.Q01 = output["edgeenergyflux01"] * n0 * mi * (lamD**3) * (omeC**3)
        self.Q02 = output["edgeenergyflux02"] * n0 * mi * (lamD**3) * (omeC**3)
        self.Gamma01 = output["edgeflux01"] * n0 * lamD * omeC
        self.Gamma02 = output["edgeflux02"] * n0 * lamD * omeC
        self.potav = output["Potav"] * Te
        
        normalizeV = lambda v, rho: (v*lamD*omeC)/(Npc*rho)
        self.vx01norm = normalizeV(output["vxav01"], output["dens01"]+1e-20)   # I don't want to divide with zero
        self.vy01norm = normalizeV(output["vyav01"], output["dens01"]+1e-20)
        self.vz01norm = normalizeV(output["vzav01"], output["dens01"]+1e-20)
        self.vx02norm = normalizeV(output["vxav02"], output["dens02"]+1e-20)
        self.vy02norm = normalizeV(output["vyav02"], output["dens02"]+1e-20)
        self.vz02norm = normalizeV(output["vzav02"], output["dens02"]+1e-20)
        
        if type == 2:
            self.density03 = output["dens03"] * n0
            self.vx03norm = normalizeV(output["vxav03"], output["dens03"]+1e-20)
            self.vy03norm = normalizeV(output["vyav03"], output["dens03"]+1e-20)
            self.vz03norm = normalizeV(output["vzav03"], output["dens03"]+1e-20)
            self.Q03 = output["edgeenergyflux03"] * n0 * mi * (lamD**3) * (omeC**3)
            self.Gamma03 = output["edgeflux03"] * n0 * lamD * omeC
            self.mi = (mi*self.HeC)+(amu*self.HC)

        self.Epar = output["Epar"] * Te / lamD
        self.Eperp = output["Eperp"] * Te / lamD
        self.phi = output["Potav"] * Te

def return_GTPC(pc_file, calculate_Gamma=False, save_pc=True):
    with open(f"{pc_file}", "rb") as file:
        save_class = pickle.load(file)
    data = save_class
    data01 = data[0]
    data02 = data[1]
    
    def get_angle(alpha, beta):
        # zvolím si výšku trojúhelníku rovnou jedné, dopočítávám druhou odvěsnu
        x = 1/np.tan(alpha)
        y = 1/np.tan(beta)
        return np.arcsin(1/np.sqrt(1+np.sqrt(x**2+y**2)**2))[0][0]
    
    def calculate_flux(object, border, plot=False, ratio=True, type=1):     # border should be 30
        alpha = get_angle(object.alphaxz, object.alphayz)
        Gamma0 = (object.n0*object.cs0*np.sin(alpha))[0]
        # Zde místo abych použil vzoreček výše, budu počítat s nasimulovanou hustotou iontů a elektronů - density 01 a density 02.
        Gamma0 = ((object.density01.mean()-object.density02.mean())*object.cs0*np.sin(alpha))
        if type==2: Gamma0 = ((object.density01.mean()+object.density02.mean()-object.density03.mean())*object.cs0*np.sin(alpha))
        
        Gamma = object.Gamma01[border]
        if type==2: Gamma = object.Gamma01[border] + object.Gamma02[border]
        
        if ratio: return np.array(Gamma/Gamma0)[0]
        return Gamma0
    
    def calculate_gamma(object, border, type, plot=False, calculate_p_flux=False):

        Gamma = object.Gamma01+1e-20
        if type == 2: Gamma = object.Gamma01+object.Gamma02+1e-20
        if calculate_p_flux:
            Gamma = calculate_flux(object, border, plot=False, ratio=False)+1e-20
            if type == 2: Gamma = calculate_flux(object, border, plot=False, ratio=False, type=2)+1e-20
        
        if type == 1:
            gamma = object.Q01/(Gamma * object.kB * object.Te)
        elif type == 2:
            gamma = (object.Q01+object.Q02)/(Gamma * object.kB * object.Te)

        return np.median(gamma[border])
    
    def gamma(tau, m_e, m_i, delta_e, ZEff):
        a = (2.5*tau)/ZEff + 2/(1-delta_e)
        b = 0.5*np.log((2*np.pi*(m_e/m_i))*(ZEff+tau)*(1-delta_e)**(-2))
        return a-b
    
    def gamma_pot(tau, Te, m_e, m_i, ZEff, pot):
        a = -(1.6021e-19*pot)/(1.3806e-23*Te)+(2.5*tau)/ZEff
        b = 2*(((ZEff+tau)*(2*np.pi*(m_e/m_i)))**(-0.5))*(np.exp((-1.6021e-19*pot)/(1.3806e-23*Te)))
        return a+b
    
    def gamma_pot_de(tau, Te, m_e, m_i, ZEff, pot, delta_e):
        a = (-1.6021e-19*pot)/(1.3806e-23*Te)+(2.5*tau)/ZEff + 2/(1-delta_e)
        b = 2*(((ZEff+tau)*(1-delta_e)*(2*np.pi*(m_e/m_i)))**-0.5)*(np.e**((-1.6021e-19*pot)/(1.3806e-23*Te)))
        return a+b
    
    gamma_TPC_01 = []
    gamma_TPC_02 = []
    me = 9.1e-31
    delta_e = 0.4
    data01.sort(key=lambda x: x.SPot)
    data01.sort(key=lambda x: x.HC)
    data01.sort(key=lambda x: x.tau[0][0])
    data02.sort(key=lambda x: x.SPot)
    data02.sort(key=lambda x: x.HC)
    data02.sort(key=lambda x: x.tau[0][0])
    for ddata01 in data01:
        #ZEff = (ddata01.HC + ddata01.HeC*2)
        ZEff = 1
        gamma_TPC_01.append([ddata01.tau[0][0], ddata01.HC, ddata01.HeC, ddata01.SPot,
                             calculate_gamma(ddata01, 30, 1, calculate_p_flux=calculate_Gamma), 
                             gamma(ddata01.tau[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], delta_e, ZEff),
                             gamma_pot(ddata01.tau[0][0], ddata01.Te[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], ZEff, ddata01.SPot), 
                             gamma_pot_de(ddata01.tau[0][0], ddata01.Te[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], ZEff, ddata01.SPot, delta_e)])
        print([ddata01.tau[0][0], ddata01.HC, ddata01.HeC, ddata01.SPot,
                             calculate_gamma(ddata01, 30, 1, calculate_p_flux=calculate_Gamma), 
                             gamma(ddata01.tau[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], delta_e, ZEff),
                             gamma_pot(ddata01.tau[0][0], ddata01.Te[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], ZEff, ddata01.SPot), 
                             gamma_pot_de(ddata01.tau[0][0], ddata01.Te[0][0], ddata01.mi[0][0]/200, ddata01.mi[0][0], ZEff, ddata01.SPot, delta_e)])
    for i in range(len(data02)):
        ddata02 = data02[i]
        ddata01 = data01[i]
        #ZEff = (ddata02.HC + ddata02.HeC*2)
        ZEff = 1
        gamma_TPC_02.append([ddata02.tau[0][0], ddata02.HC, ddata02.HeC, ddata02.SPot,
                             calculate_gamma(ddata02, 30, 2, calculate_p_flux=calculate_Gamma), 
                             gamma(ddata02.tau[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], delta_e, ZEff), 
                             gamma_pot(ddata02.tau[0][0], ddata02.Te[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], ZEff, ddata02.SPot), 
                             gamma_pot_de(ddata02.tau[0][0], ddata02.Te[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], ZEff, ddata02.SPot, delta_e)])
        print([ddata02.tau[0][0], ddata02.HC, ddata02.HeC, ddata02.SPot,
                             calculate_gamma(ddata02, 30, 2, calculate_p_flux=calculate_Gamma), 
                             gamma(ddata02.tau[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], delta_e, ZEff), 
                             gamma_pot(ddata02.tau[0][0], ddata02.Te[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], ZEff, ddata02.SPot), 
                             gamma_pot_de(ddata02.tau[0][0], ddata02.Te[0][0], ddata02.mi[0][0]/200, ddata02.mi[0][0], ZEff, ddata02.SPot, delta_e)])
    if save_pc:
        name = "gamma_TPC_012-Pot-de"
        with open(f"{name}.pickle", "wb") as file:
            pickle.dump([gamma_TPC_01, gamma_TPC_02], file)
    return gamma_TPC_01, gamma_TPC_02

return_GTPC("gamma-TPC1223.pickle", calculate_Gamma=False)