def ProbeAnalysisTe(shotlist, rs, time, alpha=2.0, BPP_alpha=1.0, alpha_He=2.0, Te_shift=1.0, R_LP=100.0, R_BPP=100.0):
    ''' 
    function for analyzing the probe data when configured for Te measurement.
    Loads the probe data from golem shot homepage, creating xarray dataset.
    Website of the latest discharge is http://golem.fjfi.cvut.cz/shots/0

    Then uses the input calibration parameters to get the Te profile and plasma
    potential. Then estimates the separatrix position from plasma potential and
    corrects the profiles. 

    Parameters
    ----------
    shotlist(List-like): python list with discharge numbers
    rs(List-like): python list with probe radial position
    time: Uniformly spaced array of times to average on
    alpha(float): calibration constant for the combined ball-pen and Langmuir probe
    BPP_alpha(float): BPP calibration coefficient
    BPP_He: BPP correction for He
    Te_shifg(float): radial difference between Langmuir and ball-pen probe
    R_LP(float): Resistance of voltage divider for Langmuir probe
    R_BPP(float): Resistance of voltage divider for Ballpen probe
    
    
    Returns
    ----------
    ds(xr.Dataset): xarray Dataset with all the data and coressponding coordinates (t, r)
    '''

    import numpy as np
    import xarray as xr
    import pandas as pd

    ds_container = []  # list prepared to store data in the for cycle
    default_t = None  # variable to store time axis for the data (can be different for each discharge :( )

    # loop over the shotlis
    for shot in shotlist:
        # open data from the shot homepage using pandas read_csv function
        # ship first 10 rows, because these are not data, just osciloscope setup
        data = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Devices/Oscilloscopes/TektrMSO64-a/TektrMSO64_ALL.csv',
                           skiprows=10)
        # set the default time as the time coordinate of the very first discharge
        if default_t is None:
            default_t = data['TIME'] * 1e3
        
        # load probes and interpolate on common time
        LP = xr.DataArray(data['CH2'] * R_LP, dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        BPP = xr.DataArray(data['CH3'] * R_BPP, dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        
        # remove offset (the initiale phase of the signal from 0 to 1 ms)
        LP -= LP.sel(t=slice(0,1)).mean('t')
        BPP -= BPP.sel(t=slice(0,1)).mean('t')
        
        # create the dataset and append it to list (creating list of xarray datasets)
        ds = xr.Dataset({'Ubpp': BPP, 'Ulp': LP})
        ds_container.append(ds)

    # after going through all the discharges, we concat (connect them) along radial axis
    # the radial axis is created using pd.Index(rs, name='r')
    # this is done like this because simple puting rs would not create the axis name
    ds = xr.concat(ds_container, pd.Index(rs, name='r'))
    ds['Te'] = (ds['Ubpp'] - ds['Ulp']) / alpha
    ds['phi'] = ds['Ubpp'] + BPP_alpha * ds['Te']
    
    # LP shift
    # because probes are on a bit different radial coordinate, we sometimes need to manually adjust
    # the shift. For this, we have to interpolate one of the probes
    r_LP = ds.r - Te_shift
    LP_da = xr.DataArray(ds['Ulp'].data,dims=['r','t'],coords={'r':r_LP,'t':ds.t})
    LP_intp = LP_da.interp(r=ds.r,method='linear')
    Te_shifted = xr.DataArray((ds['Ubpp']-LP_intp)/alpha_He,dims=['r','t'], coords={'r':ds.r,'t':ds.t})
    ds['Te_shifted'] = Te_shifted
    ds['Ulp_intp'] = LP_intp
    
    # --- Normalize 
    # first of all we create an array of all the times in the discharge we might be interested in
    step = time[1] - time [0]
    t_bins = np.arange(0,20, step)
    # now, using xarray, we groupbu_bins using our t_bins array
    gb_ds = ds.groupby_bins('t', t_bins).mean('t')
    # gb_ds_Isat = ds_Isat.groupby_bins('t', t_bins).mean('t')
    # gb_n = n.groupby_bins('t', t_bins).mean('t')

    # doing the same for the plasma potential directly (because we want to normalize by maximum of plasma potential)
    gb = ds['phi'].groupby_bins('t', t_bins).mean('t')
    # after groupingby is done, lets find positions of maximums in each time interval
    VSL_bins = gb[gb.argmax('r')].r
    # now the magix, use this line to save positions of maximums
    VSL = VSL_bins.sel(t_bins=pd.IntervalIndex.from_tuples([(i,i+step) for i in time])).data

    
    
    data_time = []
    for i, t in enumerate(time):
        r_norm = rs - VSL[i]
        phi_norm = ds['phi'].sel(t=slice(t,t+step)).mean('t', skipna=True)
        phi_std = ds['phi'].sel(t=slice(t,t+step)).std('t', skipna=True)
        Te_norm = ds['Te_shifted'].sel(t=slice(t,t+step)).mean('t', skipna=True)
        Te_std = ds['Te_shifted'].sel(t=slice(t,t+step)).std('t', skipna=True)
        Ulp_intp_norm = ds['Ulp_intp'].sel(t=slice(t,t+step)).mean('t', skipna=True)
        Ulp_intp_std = ds['Ulp_intp'].sel(t=slice(t,t+step)).std('t', skipna=True)
        
        
        
        data ={
            'r_norm': r_norm,
            'phi' : phi_norm,
            'phi_std' : phi_std,
            'Te': Te_norm,
            'Te_std' : Te_std,
            'Ulp_intp': Ulp_intp_norm,
            'Ulp_intp_std': Ulp_intp_std
            }
        data_time.append(data)
    data_time = np.array(data_time)
    
    Norm_data ={
        't': time,
        'data': data_time
        }
    
    return Norm_data
    
    
def get_Isat(shotlist, rs, R_Isat=47):
    import numpy as np
    import xarray as xr
    import pandas as pd
    # prepare an empty list to store data in
    ds_container = []
    default_t = None

    # loop over discharges
    for shot in shotlist:
        data = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Devices/Oscilloscopes/TektrMSO64-a/TektrMSO64_ALL.csv',
                           skiprows=10)
        # set default time axis
        if default_t is None:
            default_t = data['TIME'] * 1e3
        
        # load signals (Isat and plasma potential)
        Isat = xr.DataArray(data['CH2'], dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        BPP = xr.DataArray(data['CH3'], dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        Isat = Isat / R_Isat * 1e3
        
        # remove offset
        Isat -= Isat.sel(t=slice(0,1)).mean('t')
        BPP -= BPP.sel(t=slice(0,1)).mean('t')
        
        # create the dataset and append to the list
        ds = xr.Dataset({'Ubpp': BPP, 'Isat': Isat})
        ds_container.append(ds)
    
    # concat the list of dataset into one dataset along r axis
    ds = xr.concat(ds_container, pd.Index(rs, name='r'))   
    
    return ds

def get_basic_probe_dataset(shotlist, rs, alpha=2.0, BPP_alpha=1.0, alpha_He=2.0, Te_shift=1.0, R_LP=100.0, R_BPP=100.0):
    ''' 
    function for loading data the probe datafrom golem shot homepage, creating xarray dataset.
    Website of the latest discharge is http://golem.fjfi.cvut.cz/shots/0

    Parameters
    ----------
    shotlist(List): python list with discharge numbers
    rs(List): python list with probe radial position
    alpha(float): calibration constant for the combined ball-pen and Langmuir probe
    Te_shifg(float): radial difference between Langmuir and ball-pen probe
    
    Returns
    ----------
    ds(xr.Dataset): xarray Dataset with all the data and coressponding coordinates (t, r)
    '''
    
    import numpy as np
    import xarray as xr
    import pandas as pd
    
    ds_container = []  # list prepared to store data in the for cycle
    default_t = None  # variable to store time axis for the data (can be different for each discharge :( )

    # loop over the shotlis
    for shot in shotlist:
        # open data from the shot homepage using pandas read_csv function
        # ship first 10 rows, because these are not data, just osciloscope setup
        data = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Devices/Oscilloscopes/TektrMSO64-a/TektrMSO64_ALL.csv',
                           skiprows=10)
        # set the default time as the time coordinate of the very first discharge
        if default_t is None:
            default_t = data['TIME'] * 1e3
        
        # load probes and interpolate on common time
        LP = xr.DataArray(data['CH2'] * R_LP, dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        BPP = xr.DataArray(data['CH3'] * R_BPP, dims=['t'], coords={'t': data['TIME'] * 1e3}).interp(t=default_t)
        
        # remove offset (the initiale phase of the signal from 0 to 1 ms)
        LP -= LP.sel(t=slice(0,1)).mean('t')
        BPP -= BPP.sel(t=slice(0,1)).mean('t')
        
        # create the dataset and append it to list (creating list of xarray datasets)
        ds = xr.Dataset({'Ubpp': BPP, 'Ulp': LP})
        ds_container.append(ds)

    # after going through all the discharges, we concat (connect them) along radial axis
    # the radial axis is created using pd.Index(rs, name='r')
    # this is done like this because simple puting rs would not create the axis name
    ds = xr.concat(ds_container, pd.Index(rs, name='r'))
    ds['Te'] = (ds['Ubpp'] - ds['Ulp']) / alpha
    ds['phi'] = ds['Ubpp'] + BPP_alpha * ds['Te']
    
    # LP shift
    # because probes are on a bit different radial coordinate, we sometimes need to manually adjust
    # the shift. For this, we have to interpolate one of the probes
    r_LP = ds.r - Te_shift
    LP_da = xr.DataArray(ds['Ulp'].data,dims=['r','t'],coords={'r':r_LP,'t':ds.t})
    LP_intp = LP_da.interp(r=ds.r,method='linear')
    Te_shifted = xr.DataArray((ds['Ubpp']-LP_intp)/alpha_He,dims=['r','t'], coords={'r':ds.r,'t':ds.t})
    ds['Te_shifted'] = Te_shifted
    ds['Ulp_intp'] = LP_intp

    
    return ds
    
def ImportShotData(shot_number):
    """
    Imports basic shot data (Ip, Uloop, Bt) from the input shot number.
    INPUTS
    shot_number: (int or iterable) shot number or list of shot numbers to fetch
    RETURNS
    ds: (xr.Dataset) data set with the data with axis ('t', 'shot')
    """
    
    # Import needed libraries
    import numpy as np
    import xarray as xr
    import pandas as pd
    
    if isinstance(shot_number, int) or isinstance(shot_number, np.int64):
        shot_number = [shot_number] # If we pass only one shot as int convert into iterable
    
    ds_container = []
    t_start = []
    t_end = []
    # Iterate over the shots we want to import
    for i, shot in enumerate(shot_number):
        # Fetch data from server
        Ip = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/BasicDiagnostics/Results/Ip.csv',
                        skiprows=0, names=['time[s]','Ip[kA]'])
        Uloop = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/BasicDiagnostics/Results/U_loop.csv',
                        skiprows=0, names=['time[s]','Uloop[V]'])
        Bt = pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/BasicDiagnostics/Results/Bt.csv',
                        skiprows=0, names=['time[s]','Bt[T]'])
        
        # This is spaguetti code for importint the shot times, should clean it up.
        t_start.append(pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/PlasmaDetection/Results/t_plasma_start',
                                   names=['t_start']).to_numpy()[0,0])
        t_end.append(pd.read_csv(f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/PlasmaDetection/Results/t_plasma_end',
                                 names=['t_end']).to_numpy()[0,0])
        # Set one uniform time for all the data
        if i == 0:
            default_t = Ip['time[s]']
            
        # Convert data into xarrays
        IpX = xr.DataArray(Ip['Ip[kA]'], dims=['t'], coords={'t': Ip['time[s]']}).interp(t=default_t)
        UloopX = xr.DataArray(Uloop['Uloop[V]'], dims=['t'], coords={'t': Uloop['time[s]']}).interp(t=default_t)
        BtX =xr.DataArray(Bt['Bt[T]'], dims=['t'], coords={'t': Bt['time[s]']}).interp(t=default_t)
        WpX = IpX * UloopX
        
        # Group xarrays into Dataset
        ds_i = xr.Dataset({'Ip': IpX, 'Uloop': UloopX, 'Bt': BtX, 'Wp': WpX})
        ds_container.append(ds_i)
    
    ds = xr.concat(ds_container, pd.Index(shot_number, name='shot'))
    t_startX = xr.DataArray(np.array(t_start),dims=['shot'], coords={'shot':shot_number})
    t_endX = xr.DataArray(np.array(t_end),dims=['shot'], coords={'shot':shot_number})
    
    ds['t_start'] = t_startX
    ds['t_end'] = t_endX
    
    # Compute mean values in the discharge for each shot and save to dataset
    Ip_mean = []
    Ip_std = []
    Bt_mean = []
    Bt_std = []
    Uloop_mean = []
    Uloop_std = []
    Wp_mean = []
    Wp_std = []
    
    for i, shot in enumerate(shot_number):
        time_slice = slice(ds['t_start'].sel(shot=shot).data, ds['t_end'].sel(shot=shot).data)

        Ip_mean.append(ds['Ip'].sel(t=time_slice,shot=shot).mean('t'))
        Ip_std.append(ds['Ip'].sel(t=time_slice,shot=shot).std('t'))
        Bt_mean.append(ds['Bt'].sel(t=time_slice,shot=shot).mean('t'))
        Bt_std.append(ds['Bt'].sel(t=time_slice,shot=shot).std('t'))
        Uloop_mean.append(ds['Uloop'].sel(t=time_slice,shot=shot).mean('t'))
        Uloop_std.append(ds['Uloop'].sel(t=time_slice,shot=shot).std('t'))
        Wp_mean.append(ds['Wp'].sel(t=time_slice,shot=shot).mean('t'))
        Wp_std.append(ds['Wp'].sel(t=time_slice,shot=shot).std('t'))
        
    ds['Ip_mean'] = xr.DataArray(np.array(Ip_mean),dims=['shot'], coords={'shot':shot_number})
    ds['Ip_std'] = xr.DataArray(np.array(Ip_std),dims=['shot'], coords={'shot':shot_number})
    ds['Bt_mean'] = xr.DataArray(np.array(Bt_mean),dims=['shot'], coords={'shot':shot_number})
    ds['Bt_std'] = xr.DataArray(np.array(Bt_std),dims=['shot'], coords={'shot':shot_number})
    ds['Uloop_mean'] = xr.DataArray(np.array(Uloop_mean),dims=['shot'], coords={'shot':shot_number})
    ds['Uloop_std'] = xr.DataArray(np.array(Uloop_std),dims=['shot'], coords={'shot':shot_number})
    ds['Wp_mean'] = xr.DataArray(np.array(Wp_mean),dims=['shot'], coords={'shot':shot_number})
    ds['Wp_std'] = xr.DataArray(np.array(Wp_std),dims=['shot'], coords={'shot':shot_number})

        
    return ds


def plasma_detection(shot_no=0):
    import requests
    import numpy as np
    import matplotlib.pyplot as plt
    import urllib
    
    plasma_start=float(requests.get(f'http://golem.fjfi.cvut.cz/shots/{shot_no}/Diagnostics/PlasmaDetection/Results/t_plasma_start').text)
    plasma_end=float(requests.get(f'http://golem.fjfi.cvut.cz/shots/{shot_no}/Diagnostics/PlasmaDetection/Results/t_plasma_end').text)

    return plasma_start, plasma_end


def icon_fig(shotlist=[0], save=False, axvline_x = 10):
    import requests
    import numpy as np
    import matplotlib.pyplot as plt
    import urllib
    
    
    figs, axs = plt.subplots(len(shotlist), 1, figsize=(10, 6 * len(shotlist)), sharex=True)

    for i, shot in enumerate(shotlist):
        try:
            url_rad = f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/FastCameras/plasma-film_Vertical2.png'
            url_vert = f'http://golem.fjfi.cvut.cz/shots/{shot}/Diagnostics/FastCameras/plasma-film_Radial2.png'
            urllib.request.urlretrieve(url_rad, 'plasma-film_Vertical2.png')
            urllib.request.urlretrieve(url_vert, 'plasma-film_Radial2.png')
        except:
            print('data not downloaded')

        vert = plt.imread('plasma-film_Vertical2.png')
        rad = plt.imread('plasma-film_Radial2.png')
        maxlen = min(vert.shape[1], rad.shape[1])
        stacked = np.vstack((rad[:, :maxlen, :], vert[:, :maxlen, :])) 

        # Dynamic time axis based on plasma_detection function
        t_start, t_end = plasma_detection(shot)
        time_axis = np.linspace(t_start, t_end, maxlen)  # Adjusting the time axis range


        # Show the stacked image
        axs[i].imshow(stacked, aspect='auto')

        # Set x-axis to represent time in milliseconds
        axs[i].set_xticks(np.linspace(0, maxlen, num=11))  # Adjust number of ticks based on time range
        axs[i].set_xticklabels([f'{t:.1f} ms' for t in np.linspace(t_start, t_end, 11)],
                          fontweight='bold')  # Set labels dynamically

        # Set labels and title
        axs[i].set_xlabel('t [ms]', fontweight='bold')
        axs[i].set_ylabel('Pixel Position', fontweight='bold')
        axs[i].set_title('Plasma image from fast cameras', fontweight='bold')

        # Save the figure
        if(save is True):
            plt.savefig(f'fast_cameras_{shot_no}.png', dpi=300)
            
        t_index = np.argmin(np.abs(t-8))
        axs[i].axvline(x=t_index, color='white', ls='--')
    return figs, axs
















    
    