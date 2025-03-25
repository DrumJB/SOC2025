#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 12:07:56 2023

@author: jachym
"""

step = 20    # number of simulations in one file
name = "gamma-TPC1223"

taus = ["05", "10", "20", "35", "50"]
Hs = ["02", "05", "10", "40", "60", "90", "95", "98"]
Hes = ["98", "95", "90", "60", "40", "10", "05", "02"]
Pots = ["01", "03", "06", "10"]
namesA = []
namesB = []
for tau in taus:
    for Pot in Pots:
        for i in range(len(Hs)):
            namesA.append(f"1T{tau}H{Hs[i]}He{Hes[i]}P{Pot}.sh")
            namesB.append(f"2T{tau}H{Hs[i]}He{Hes[i]}P{Pot}.sh")

number_of_sim = len(namesA)
for i in range(int(number_of_sim/step)):
    with open(f"QS-{name}{i}.sh", "w+") as output:
        output.write('CHECK_MARK="✅"\n')
        output.write('echo -e "\n\e[4mRequest process\e[0m"\n')
        for j in range(step):
            output.write(f"qsub {name}/{namesA[j+i*step]}\n")
            output.write(f"qsub {name}/{namesB[j+i*step]}\n")
            output.write('echo -e "${CHECK_MARK} '+str(namesB[j+i*step][1:])+'"\n')
        output.write("All calculations has been requested.")