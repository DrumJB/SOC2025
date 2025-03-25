#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 21:28:53 2023

@author: jachym
"""
import scipy.io
import os
import pickle

class save_class:
    
    def __init__(self):
        self.r = []

    def load_mat(self, mat_path):
        mat = scipy.io.loadmat(mat_path)

        # find out the type of file:
        keys = mat.keys()
        type = ""
        for key in keys:
            if key == "Eperp":
                type = "o"
            elif key == "irel":
                type = "t"
            elif key == "stype":
                type = "0"
        if type == "":
            print("Error: Can't identify the file type!")
        # identify angles:
        name = mat_path.split("/")[-1]
        name_array = list(name)
        type += name[0]
        integer_id = name_array[0] + name_array[2] + name_array[3] + name_array[5] + name_array[6] + name_array[9] + name_array[10] + name_array[12]
        # load the wanted variables:
        wanted = []
        if type == "o1":
            wanted = ["dens01","dens02","dt","edgeenergyflux01","edgeenergyflux02","edgeflux01","edgeflux02",
                      "vxav01","vxav02","vyav01","vyav02","vzav01","vzav02","Epar","Eperp","Potav"]
        if type == "t1":
            wanted = ["mksTe","mksB","tau","mksn0","mksmainionm","mksmainionq","Npc","alphaxz","alphayz"]
        if type == "o2":
            wanted = ["dens01","dens02","dens03","dt","edgeenergyflux01","edgeenergyflux02","edgeenergyflux03","edgeflux01",
                      "edgeflux02","edgeflux03","vxav01","vxav02","vxav03","vyav01","vyav02","vyav03","vzav01","vzav02","vzav03",
                      "Epar","Eperp","Potav"]
        if type == "t2":
            wanted = ["mksTe","mksB","tau","mksn0","mksmainionm","mksmainionq","Npc","alphaxz","alphayz"]
        dic_arr = [["id", integer_id]]
        for key in wanted:
            dic_arr.append([key, mat[key]])
        appended = False
        for dic in self.r:
            if dic["id"] == integer_id:
                dic.update(dict(dic_arr))
                appended = True
        if not appended:
            self.r.append(dict(dic_arr))

def return_data(save_pc=True, name="saved_pickle"):
    save_pickle = save_class()
    for filename in os.listdir("d"):
        print(filename)
        save_pickle.load_mat("d/"+filename)
    if save_pc:
        with open(f"{name}.pickle", "wb") as file:
            pickle.dump(save_pickle.r, file)
    return save_pickle.r

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

class DenormalizedFile:

    def __init__(self, data, name="datafile"):
        
        self.data01 = []
        self.data02 = []
        for d in data:
            print(d["id"])
            if "dens03" in list(d.keys()):
                self.data02.append(DenormalizedSimulation(d, 2))
            else:
                self.data01.append(DenormalizedSimulation(d, 1))
        with open(f"{name}.pickle", "wb") as file:
            pickle.dump([self.data01, self.data02], file)

if __name__ == "__main__":
    name = "gamma-TPC1223"
    data = return_data(save_pc=False, name=name)
    result = DenormalizedFile(data, name=name)