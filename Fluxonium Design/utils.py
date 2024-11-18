"""
This file contains fitting functions which are available in Andersen Lab Analysis,
but are copied into this repository to ensure it is standalone. 
"""
import os
import numpy as np
import h5py
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Ensure font sizes are large enough to read
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)

# Ensure matplotlib NEVER, EVER plots an ugly
# figure with invisible axes.
plt.rcParams['axes.facecolor']=(1,1,1,1)
plt.rcParams['savefig.facecolor']=(1,1,1,1)


def fig_prepare(xlabel, ylabel, xscale="linear", yscale="linear"):
    """Prepare a matplotlib figure

    Parameters
    ----------
    xlabel : str
        Label for x-axis
    ylabel : str 
        Label for y-axis
    xscale : str
        Valid xscale label for matplotlib x-axis
    yscale : str
        Valid xscale label for matplotlib y-axis

    Returns
    -------
    Matplotlib figure and axis
    """
    fig,ax = plt.subplots()
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid(color="#A2A2A2", linestyle="--", linewidth=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig,ax


def get_data_folders(base_dir: str, 
                     start: int, 
                     stop: int):
    """
    Finds a range of data folders based on their timestamp.

    Parameters
    ----------
    base_dir : str
        Directory to search in
    
    start : int
        Start time stamp

    stop : int
        Stop time stamp

    Returns
    -------
    List of data folders
    """
    all_data_folders = next(os.walk(base_dir))[1]
    folder_numbers = np.array([int(x.split("-")[1]) for x in all_data_folders])
    idx = np.argwhere((start<=folder_numbers)&(folder_numbers<=stop))
    data_folders = np.squeeze(np.array(all_data_folders)[idx], axis=1)
    data_folders = [base_dir+x for x in data_folders]

    return data_folders


def interp_fr(base_dir, date_dir, exp_number):
    data_folders = get_data_folders(base_dir+date_dir, start=exp_number, stop=exp_number)

    with h5py.File(data_folders[0]+"/analysis_ResonatorFluxScanAnalysis/dataset_processed.hdf5", "r") as f:
        data_bias_current = np.array(list(f["bias_current"]))
        data_frequency = np.array(list(f["frequency"]))
        data_mag = np.array(list(f["magnitude"]))

    return interp1d(data_bias_current, data_frequency[np.argmin(data_mag, axis=0)])


def asymmetric_lorentzian(f, a0, a1, fr, kappa, gamma, alpha):
    """
    Fitting function typically used for resonators.

    Parameters
    ----------
    a0, a1 : float
        Parameters of linear background. Background is of the form
        a0 + a1*f
    fr : float
        Resonance frequency
    kappa : float
        Coupling loss. kappa = 2pi*fr/Qc
    gamma : float
        Internal loss. gamma = 2pi*fr/Qi
    alpha : float
        Resonator assymetry
    """
    S21 = 1 - kappa*(1+alpha*1j) / (2j*2*np.pi*(f-fr)+kappa+gamma)
    p_abs = np.poly1d([a0, a1])
    return p_abs(f)*np.abs(S21)


def double_asymmetric_lorentzian(f, a0, a1, fr1, fr2, kappa1, kappa2, gamma1, gamma2, alpha, A):
    """
    Fitting function typically used for resonators. We assume that a0, a1 and the assymmetry
    alpha are the same for both resonators.

    Parameters
    ----------
    a0, a1 : float
        Parameters of linear background. Background is of the form
        a0 + a1*f
    fr1 : float
        Resonance frequency 1
    fr2 : float
        Resonance frequency 1
    kappa : float
        Coupling loss. kappa = 2pi*fr/Qc
    gamma : float
        Internal loss. gamma = 2pi*fr/Qi
    alpha : float
        Resonator assymetry
    """
    S21_1 = 1 - kappa1*(1+alpha*1j) / (2j*2*np.pi*(f-fr1)+kappa1+gamma1)
    S21_2 = 1 - kappa1*(1+alpha*1j) / (2j*2*np.pi*(f-fr2)+kappa2+gamma2)
    p_abs = np.poly1d([a0, a1])
    return p_abs(f)*np.abs(A*S21_1 + (1-A)*S21_2)
