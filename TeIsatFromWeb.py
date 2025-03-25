import pandas as pd
import numpy as np
import requests
import io

def load_shot(shot_n, reference_shot=None, Te=False):
        # returns isat and plasma potential with time step 100 ns

    if reference_shot is None: reference_shot = shot_n

    mode_url = f"http://golem.fjfi.cvut.cz/shots/{shot_n}/Diagnostics/LangBallPenProbe/Parameters/s_mode_lp"
    IpDf = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{reference_shot}/Diagnostics/BasicDiagnostics/Results/Ip.csv', names=['time', 'Ip'])

    if Te: df = pd.read_csv(f'/home/jachym/SPICEdata/golem/{shot_n}.csv')
    else:

        oscilloscope_url = f"http://golem.fjfi.cvut.cz/shots/{shot_n}/Devices/Oscilloscopes/TektrMSO64-a/TektrMSO64_ALL.csv"

        response = requests.get(oscilloscope_url).content.decode("utf-8")
        header, found = 0, False
        while not found:
            if response[header:header+4] == 'TIME':
                found = True
            else:
                header += 1
            if header>1000:
                print('Error: Header not found')
                found=True
        df = pd.read_csv(io.StringIO(response[header:]), sep=",")

    mode = requests.get(mode_url).content.decode("utf-8")[:-1]
    

    start_index = IpDf[['Ip']].idxmax().values[0]+500
    end_index = IpDf[['Ip']].idxmax().values[0]+2500
    
    # scaling constants
    R_Langmuir = 100
    R_BallPen = 100
    Langmuir_data = R_Langmuir*df['CH2']-np.mean(df['CH2'][0:500])
    BallPen_data = R_BallPen*df['CH3']-np.mean(df['CH3'][0:500])

    return Langmuir_data[start_index:end_index], BallPen_data[start_index:end_index]


def get_Isat(shot_n, reference_shot=None, Te=False):
    # returns isat and plasma potential with time step 100 ns

    LP, BPP = load_shot(shot_n, reference_shot=reference_shot, Te=Te)
    R_Isat = 47

    Isat = LP/R_Isat * 1e3
    
    return Isat, BPP

def get_Te(shot_n, reference_shot=None):
    # returns Te with time step 100 ns

    LP, BPP = load_shot(shot_n, reference_shot=reference_shot, Te=True)
    alpha_golem = 2.2

    Te = (BPP-LP)/alpha_golem
    
    return Te

def get_r(n_shot):
    return float(requests.get(f"http://golem.fjfi.cvut.cz/shots/{n_shot}/Diagnostics/LangBallPenProbe/Parameters/r_lp_tip").content.decode("utf-8")[:-1])

def get_Ulp(n_shot):
    return float(requests.get(f"http://golem.fjfi.cvut.cz/shots/{n_shot}/Diagnostics/LangBallPenProbe/Parameters/u_lp").content.decode("utf-8")[:-1])

def get_LP(shot_n, reference_shot=None, Te=True):
    LP, BPP = load_shot(shot_n, reference_shot=reference_shot, Te=Te)
    return LP