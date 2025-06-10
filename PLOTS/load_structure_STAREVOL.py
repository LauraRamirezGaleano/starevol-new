#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
*******************************************************************************
  load_structure_STAREVOL : Load structure data from STAREVOL results.
  Version 1.01, 14/06/2024
  VOJE Thomas <thomas.voje@umontpellier.fr>
*******************************************************************************
"""

__author__  = "VOJE Thomas <thomas.voje@umontpellier.fr>"
__date__    = "14/06/2024"
__version__ = "1.01"

# =============================================================================
#   Imports
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import warnings

from matplotlib.collections import LineCollection
from time import time

# =============================================================================
#   Main parameters 
# =============================================================================

# Path where the structure files are loaded (RESULTS/DATA/)
loading_path = r"/home/thomasvoje/Codes/starevol/RESULTS/DATA/"

# Path where the figures are saved
saving_path = r"/home/thomasvoje/Figures/Reunion/"

# If True, skip the first line of the files when loading to skip the header
header = True

# If True, ignore warnings
ignore_warnings = True

# Colors used by default in the code (by default matplotlib's basic colors)
base_colors = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple",
               "tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]
base_colormap = "viridis"

# =============================================================================
#   Initialization
# =============================================================================

# Apply header and ignore_warnings
skiprows = 0
if header:
    skiprows = 1
if ignore_warnings:
    warnings.filterwarnings("ignore")
    
# Define tables
profile_names = []
colors = []
labels = []
linestyles = []

# Initializing quantity tables
# -----------------------------------------------------------------------------
# All the quantities loaded by this program are defined here.

mass, L, R, age, seq, dt = [r"m$_{\rm r}$ $[M_{\odot}]$"],[],[],[],[],[]
nb_shell, yzi = [r"shell number"],[r"Mixing type"]
dm, lum = [r"dm$_{\rm r}$ $[M_{\odot}]$"],[r"L$_{\rm r}$ $[L_{\odot}]$"]
lum_rad = [r"L$_{\rm rad}$ $[L_{\odot}]$"]
r, T, rho = [r"R $[R_{\odot}]$"],[r"T $[K]$"],[r"$\rho$ $[cgs]$"]
P, beta, eta, lnf = [r"P $[cgs]$"],[r"$\beta$"],[r"$\eta$"],[r"ln(f)"]
s, accel = [r"Entropy"],[r"acceleration $[cgs]$"]
u, xmr = [r"velocity $[cgs]$"],[r"mass coordinate"]
tau, kappa, mu, mue = [r"$\tau$"],[r"$\kappa$ $[cgs]$"],[r"$\mu$"],[r"$\mu_e$"]
kappa_e = [r"$\kappa_e$ $[cgs]$"]
abad, abmu = [r"$\nabla_{\rm ad}$"],[r"$\nabla_{\rm \mu}$"]
abrad, abla = [r"$\nabla_{\rm rad}$"],[r"$\nabla$"]
cs, cp = [r"sound speed $[cgs]$"],[r"c$_{\rm P}$ $[cgs]$"]
cg = [r"isothermal sound speed $[cgs]$"]
gamma1, pvisc = [r"$\Gamma_{1}$"],[r"P$_{\rm visc}$ $[cgs]$"]
e_nucl, e_int = [r"e$_{\rm nucl}$ $[cgs]$"],[r"e$_{\rm int}$ $[cgs]$"]
e_nupla, e_nunucl = [r"e$_{\rm \nu}$ $[cgs]$"],[r"e$_{\rm \nu,nucl}$ $[cgs]$"]
e_grav = [r"e$_{\rm grav}$ $[cgs]$"]
D_conv, D_turb = [r"D$_{\rm conv}$ $[cgs]$"],[r"D$_{\rm turb}$ $[cgs]$"]
v_conv, t_conv = [r"$v_{conv}$ $[cgs]$"],[r"$\tau_{conv}$ $[cgs]$"]
F_conv = [r"F$_{\rm conv}$/F$_{\rm tot}$"]
F_rad = [r"F$_{\rm rad}$/F$_{\rm tot}$"]
compression_rate, heat_rate = [r"compression rate"],[r"heat rate"]
work_e_grav = [r"contribution of work in e$_{\rm grav}$ $[cgs]$"]
heat_e_grav = [r"contribution of heat in e$_{\rm grav}$ $[cgs]$"]
dlum_dt = [r"dln(L)/dt"]
H1, H2, He3, He4 = [r"$^{1}$H"],[r"$^{2}$H"],[r"$^{3}$He"],[r"$^{4}$He"]
Li6, Li7, Be7, B8 = [r"$^{6}$Li"],[r"$^{7}$Li"],[r"$^{7}$Be"],[r"$^{8}$B"]
Be9, B10, B11, C12 = [r"$^{9}$Be"],[r"$^{10}$B"],[r"$^{11}$B"],[r"$^{12}$C"]
C13, N13, C14, N14 = [r"$^{13}$C"],[r"$^{13}$N"],[r"$^{14}$C"],[r"$^{14}$N"]
N15, O15, O16, O17 = [r"$^{15}$N"],[r"$^{15}$O"],[r"$^{16}$O"],[r"$^{17}$O"]
O18, F18, F19, F20 = [r"$^{18}$O"],[r"$^{18}$F"],[r"$^{19}$F"],[r"$^{20}$F"]
Ne20, Ne21, Ne22 = [r"$^{20}$Ne"],[r"$^{21}$Ne"],[r"$^{22}$Ne"]
Na22, Ne23, Na23 = [r"$^{22}$Na"],[r"$^{23}$Ne"],[r"$^{23}$Na"]
Na24, Mg24, Na25 = [r"$^{24}$Na"],[r"$^{24}$Mg"],[r"$^{25}$Na"]
Mg25, Mg26, Al26m = [r"$^{25}$Mg"],[r"$^{26}$Mg"],[r"$^{26m}$Al"]
Al26g, Mg27, Al27 = [r"$^{26g}$Al"],[r"$^{27}$Mg"],[r"$^{27}$Al"]
Si28, Si29, Si30 = [r"$^{28}$Si"],[r"$^{29}$Si"],[r"$^{30}$Si"]
P31, S32, S33 = [r"$^{31}$P"],[r"$^{32}$S"],[r"$^{33}$S"]
S34, S35, Cl35 = [r"$^{34}$S"],[r"$^{35}$S"],[r"$^{35}$Cl"]
S36, Cl36, Cl37 = [r"$^{36}$S"],[r"$^{36}$Cl"],[r"$^{37}$Cl"]
neutron, heavy, sumX = [r"$^{1}$n"],[r"heavy"],[r"Sum mass fractions"]
metallicity, T_eff, n_eff = [r"metallicity"],[r"T$_{\rm eff}$ $[K]$"],["neff"]
Hp, gmr = [r"H$_{\rm P}$ $[cgs]$"],[r"g$_{\rm r}$ $[cgs]$"]
phi_eos, delta_eos, N2 = [r"$\phi$"],[r"$\delta$"],[r"N$^2$"]
N2_thermal, N2_chemical = [r"N$^2$ (thermal)"],[r"N$^2$ (chemical)"]
L_edd = [r"L$_{\rm Edd}$ $[L_{\odot}]$"]
# See Jiang et al. (2015,ApJ,813,74)
tau_crit = [r"$\tau_{\rm crit}$"]
tau_0 = [r"$\tau_{\rm 0}$"]

# -----------------------------------------------------------------------------

# =============================================================================
#   Console display and input
# =============================================================================

print(__doc__)
print(" loading results in ..", loading_path)
print(" saving figures in ...", saving_path)
N = int(input("\n Number of STAREVOL profiles to load : "))
for n in range(N):
    profile_name_tmp = input(f"\n profile name {n+1} : ")
    profile_names.append(profile_name_tmp)

# YOU CAN ALSO LOAD FILES BY ADDING A profile_names LIST BELOW THEN RUNNING 
# THE PROGRAM WITH N = 0.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

N = len(profile_names)

# =============================================================================
#   Define plotting functions
# =============================================================================

def plot(x, y, fig=None, ax=None, output=False, domain=None, figsize=None,
         invert_xaxis=False, invert_yaxis=False, xlabel=None, ylabel=None,
         log_xaxis=False, log_yaxis=False, save_fig=False, save_name=None,
         mapping=None, cmap=base_colormap, log_mapping=False, cbarlabel=None,
         cbar_min=None, cbar_max=None, title=None):

    ''' 
    Main plotting function. Plot y as a function of x. Allows mapping (ie. 
    coloring the plot according to the value of a specified quantity).
    
    Parameters 
    ----------
    
    x, y : lists of ndarray
        Quantities to be plotted. 

    fig, ax : Figure and Axes [None]
        If specified, modify the figure given as an entry.

    output : bool [False]
        If True, returns the figure and axes.

    domain : list, ndarray or tuple [None]
        If specified, form needed is (xmin,xmax,ymin,ymax). Plot using the 
        values in domain as axis limits.
        
    figsize : tuple [None]
        If specified, size of the figure.
        
    invert_xaxis : bool [False]
        If True, invert the x-axis.
        
    invert_yaxis : bool [False]
        If True, invert the y-axis.
        
    xlabel : str [None]
        If specified, label of the x-axis.
        
    ylabel : str [None]
        If specified, label of the y-axis.
        
    log_xaxis : bool [False]
        If True, plot with log(x) instead of x.
        
    log_yaxis : bool [False]
        If True, plot with log(y) instead of y.
        
    save_fig : bool [False]
        If True, save the figure in the saving_path folder.
        
    save_name : str [None]
        If specified, name under which the figure is saved.
        
    mapping : list of ndarray [None]
        If specified, quantity according to which the plot is colored.
        
    cmap : str ["viridis"]
        Name of the colormap.
        
    log_mapping : bool [False]
        If True, plot with log(mapping) instead of mapping.
        
    cbarlabel : str [None]
        If specified, label of the colormap.
        
    cbar_min : float [None]
        If specified, minimum value for the colorbar.
        
    cbar_max : float [None]
        If specified, maximum value for the colorbar.
        
    title : str [None]
        If specified, title of the figure.
        
    Returns
    -------
    
    fig, ax : Figure and Axes (returned only is output = True).
    '''

    # Initializing the figure and axes
    if fig is None or ax is None:
        if figsize is None:
            fig, ax = plt.subplots(figsize=(10,7))
        else:
            fig, ax = plt.subplots(figsize=figsize)

    # Main loop for plotting the data
    for n in range(N):

        X = x[n]
        Y = y[n]

        # Mapping : Plotting an other quantity as the colormap of the plot
        # (if specified)
        if mapping is None:
            linestyle = linestyles[n]
            if (linestyle=='solid' or linestyle=='dashed' or linestyle=='-' 
                or linestyle=='--' or linestyle=='-.'  or linestyle=='dotted'
                or linestyle=='dashdot' or linestyle=='--' or linestyle==':'
                or linestyle==' '):
                ax.plot(X, Y, ls=linestyle, color=colors[n], label=labels[n])
            else:
                ax.plot(X, Y, marker=linestyle, color=colors[n], 
                        label=labels[n])
        else:
            # See :
            # https://matplotlib.org/stable/gallery/lines_bars_and_
            # markers/multicolored_line.html
            points = np.array([X, Y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            mapp = mapping[n]
            if log_mapping:
                mapp = np.log10(mapping[n])
            if cbar_min is None or cbar_max is None:
                mappmin, mappmax = np.nanmin(mapp), np.nanmax(mapp)
            else:
                mappmin, mappmax = cbar_min, cbar_max
            norm = plt.Normalize(mappmin, mappmax)
            lc = LineCollection(segments, cmap=cmap, norm=norm)
            lc.set_array(mapp)
            line = ax.add_collection(lc)
            ax.plot(X[-1], Y[-1], color="k", 
                    label=labels[n])

    # Colorbar of the mapping (if specified)
    if mapping is not None:
        cbar = fig.colorbar(line, ax=ax)
        if cbarlabel is None:
            cbar.set_label(label=mapping[-1], size=14)
        else:
            cbar.set_label(label=cbarlabel, size=14)

    # Computing the axis limits (if unspecified)
    if domain is None:
        xmin = x[0][0]
        if log_xaxis:
            xmin = np.log10(x[0][0])
        for n in range(N):
            X = np.nanmin(x[n])
            if log_xaxis:
                X = np.log10(np.nanmin(x[n]))
            if xmin>X:
                xmin = X
        xmax = x[0][0]
        if log_xaxis:
            xmax = np.log10(x[0][0])
        for n in range(N):
            X = np.nanmax(x[n])
            if log_xaxis:
                X = np.log10(np.nanmax(x[n]))
            if xmax<X:
                xmax = X
        ymin = y[0][0]
        if log_yaxis:
            ymin = np.log10(y[0][0])
        for n in range(N):
            Y = np.nanmin(y[n])
            if log_yaxis:
                Y = np.log10(np.nanmin(y[n]))
            if ymin>Y:
                ymin = Y
        ymax = y[0][0]
        if log_yaxis:
            ymax = np.log10(y[0][0])
        for n in range(N):
            Y = np.nanmax(y[n])
            if log_yaxis:
                Y = np.log10(np.nanmax(y[n]))
            if ymax<Y:
                ymax = Y
        width_factor = 0.1*abs(xmax-xmin)
        height_factor = 0.1*abs(ymax-ymin)
        xmin = xmin - width_factor
        xmax = xmax + width_factor
        ymin = ymin - height_factor
        ymax = ymax + height_factor
        if log_xaxis:
            xmin, xmax = 10**xmin, 10**xmax
        if log_yaxis:
            ymin, ymax = 10**ymin, 10**ymax
        domain = (xmin,xmax,ymin,ymax)

    # Setting the axis limits
    ax.axis(domain)

    # Inverting axis (if specified)
    if invert_xaxis:
        ax.invert_xaxis()
    if invert_yaxis:
        ax.invert_yaxis()

    # Logarithmic scale (if specified)
    if log_xaxis:
        ax.set_xscale("log")
    if log_yaxis:
        ax.set_yscale("log")

    # Adding labels to axes (if specified)
    if xlabel is None:
        ax.set_xlabel(x[-1], fontsize=14)
    else:
        ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is None:
        ax.set_ylabel(y[-1], fontsize=14)
    else:
        ax.set_ylabel(ylabel, fontsize=14)
    
    if title is not None:
        ax.set_title(title, y=1.02, fontsize=15)

    # Adding legend and ticks and displaying the figure
    i = 0
    for label in labels:
        if label is None:
            i+=1
    if i!=N:
        ax.legend()
    ax.tick_params(axis='both', labelsize=12)
    fig.tight_layout()
    fig.show()
    
    if x is yzi or y is yzi or mapping is yzi:
        print("\n Mixing type")
        print(" -3 : semi-convective shell (Ledoux off)")
        print(" -2 : convective shell")
        print(" -1 : atmospheric shell")
        print("  1 : overshoot shell")
        print("  2 : thermohaline shell")
        print("  3 : semi-convective shell (Ledoux on)")
        print("  4 : radiative shell")

    # Saving figure (if specified)
    if save_fig:
        print("\n Saving plot.")
        if save_name is None:
            fig.savefig(saving_path+"plot_structure_starevol.png", dpi=200)
        else:
            fig.savefig(saving_path+save_name, dpi=200)
            
    # Returning the figure and axis (if specified)
    if output:
        return fig, ax

# =============================================================================

def multiplot(x, y, n=1, fig=None, ax=None, figsize=None, domain=None,
              invert_xaxis=False, invert_yaxis=False, xlabel=None, title=None,
              log_xaxis=False, log_yaxis=False, colors=None, output=False,
              save_fig=False, save_name=None):
    
    '''
    Multiple plotting function. Plot each quantity in the list y as a function
    of x for the selected loaded STAREVOL result.
    
    Parameters
    ----------
    
    n : int [1]
        Number of the result plotted.
    
    The rest is identical to plot() except that y is a list of quantities
    instead of a quantity.
    
    Returns
    -------
    
    Identical to plot().
    '''
    
    n = n-1
    
    # Initializing the figure and axes
    if fig is None or ax is None:
        if figsize is None:
            fig, ax = plt.subplots(figsize=(10,7))
        else:
            fig, ax = plt.subplots(figsize=figsize)
            
    if colors is None:
        colors = base_colors
    X = x[n]
    for i in range(len(y)):
        Y = y[i][n]
        linestyle = linestyles[n]
        if (linestyle=='solid' or linestyle=='dashed' or linestyle=='-' 
            or linestyle=='--' or linestyle=='-.'  or linestyle=='dotted'
            or linestyle=='dashdot' or linestyle=='--' or linestyle==':'
            or linestyle==' '):
            ax.plot(X, Y, ls=linestyle, color=colors[i%len(colors)],
                    label=y[i][-1])
        else:
            ax.plot(X, Y, marker=linestyle,
                    color=colors[i%len(colors)])
            
    # Computing the axis limits (if unspecified)
    if domain is None:
        xmin = np.nanmin(x[n])
        if log_xaxis:
            xmin = np.log10(np.nanmin(x[n]))
        xmax = np.nanmax(x[n])
        if log_xaxis:
            xmax = np.log10(np.nanmax(x[n]))
        ymin = y[0][n][0]
        if log_yaxis:
            ymin = np.log10(y[0][n][0])
        for i in range(len(y)):
            Y = np.nanmin(y[i][n])
            if log_yaxis:
                Y = np.log10(np.nanmin(y[i][n]))
            if ymin>Y:
                ymin = Y
        ymax = y[0][n][0]
        if log_yaxis:
            ymax = np.log10(y[0][n][0])
        for i in range(len(y)):
            Y = np.nanmax(y[i][n])
            if log_yaxis:
                Y = np.log10(np.nanmax(y[i][n]))
            if ymax<Y:
                ymax = Y
        width_factor = 0.1*abs(xmax-xmin)
        height_factor = 0.1*abs(ymax-ymin)
        xmin = xmin - width_factor
        xmax = xmax + width_factor
        ymin = ymin - height_factor
        ymax = ymax + height_factor
        if log_xaxis:
            xmin, xmax = 10**xmin, 10**xmax
        if log_yaxis:
            ymin, ymax = 10**ymin, 10**ymax
        domain = (xmin,xmax,ymin,ymax)

    # Setting the axis limits
    ax.axis(domain)
    
    # Inverting axis (if specified)
    if invert_xaxis:
        ax.invert_xaxis()
    if invert_yaxis:
        ax.invert_yaxis()

    # Logarithmic scale (if specified)
    if log_xaxis:
        ax.set_xscale("log")
    if log_yaxis:
        ax.set_yscale("log")

    # Adding labels to axes (if specified)
    if xlabel is None:
        ax.set_xlabel(x[-1], fontsize=14)
    else:
        ax.set_xlabel(xlabel, fontsize=14)
        
    if title is not None:
        ax.set_title(title, y=1.02, fontsize=15)

    # Adding legend and ticks and displaying the figure
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    fig.tight_layout()
    fig.show()
    
    # Saving figure (if specified)
    if save_fig:
        print("\n Saving plot.")
        if save_name is None:
            fig.savefig(saving_path+"plot_structure_starevol.png", dpi=200)
        else:
            fig.savefig(saving_path+save_name, dpi=200)
            
    # Returning the figure and axis (if specified)
    if output:
        return fig, ax

# =============================================================================

def rho_T(output=False, domain=None, cmap=base_colormap, mapping=None, 
          cbarlabel=None, figsize=None, log_mapping=False, 
          save_fig=False, save_name=None):
    
    fig, ax = plot(log10(rho), log10(T), output=True, 
                   domain=domain, figsize=figsize,
                   cmap=cmap, mapping=mapping, cbarlabel=cbarlabel,
                   log_mapping=log_mapping, save_fig=save_fig,
                   save_name=save_name)
        
    # Returning the figure and axis (if specified)
    if output:
        return fig, ax

# =============================================================================
#   Define data management functions 
# =============================================================================

def _load(n):
    
    ''' Private function. Load model number n. '''
    
    profile_name = profile_names[n]

    # Load .p0 file (general informations)
    L.append(np.loadtxt(fname=loading_path+profile_name+".p0", 
                        skiprows=skiprows, usecols=(2)))
    R.append(np.loadtxt(fname=loading_path+profile_name+".p0", 
                   skiprows=skiprows, usecols=(3)))
    age.append(np.loadtxt(fname=loading_path+profile_name+".p0", 
                     skiprows=skiprows, usecols=(4)) / (365.2425*24*3600))
    dt.append(np.loadtxt(fname=loading_path+profile_name+".p0", 
                         skiprows=skiprows, usecols=(11)))
    seq.append(int(np.loadtxt(fname=loading_path+profile_name+".p0", 
                         skiprows=skiprows, usecols=(12))))
    
    # Load .p1 file (main profiles)
    nb_shell.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                     skiprows=skiprows, usecols=(0)))
    yzi.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                     skiprows=skiprows, usecols=(1)))
    r.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                   skiprows=skiprows, usecols=(2)))
    T.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                   skiprows=skiprows, usecols=(3)))
    rho.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                     skiprows=skiprows, usecols=(4)))
    P.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                   skiprows=skiprows, usecols=(5)))
    beta.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                           skiprows=skiprows, usecols=(6)))
    eta.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                         skiprows=skiprows, usecols=(7)))
    lnf.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                          skiprows=skiprows, usecols=(8)))
    s.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                          skiprows=skiprows, usecols=(9)))
    u.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                          skiprows=skiprows, usecols=(11)))
    xmr.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                     skiprows=skiprows, usecols=(12)))
    accel.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p1", 
                            skiprows=skiprows, usecols=(13)))
    
    # Load .p2 file (other thermodynamic quantities)
    tau.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                     skiprows=skiprows, usecols=(1)))
    kappa.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                     skiprows=skiprows, usecols=(2)))
    mu.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                    skiprows=skiprows, usecols=(3)))
    mue.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                   skiprows=skiprows, usecols=(4)))
    abad.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                      skiprows=skiprows, usecols=(5)))
    abmu.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                        skiprows=skiprows, usecols=(6)))
    abrad.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                        skiprows=skiprows, usecols=(8)))
    abla.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                          skiprows=skiprows, usecols=(9)))
    cs.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                    skiprows=skiprows, usecols=(10)))
    cp.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                skiprows=skiprows, usecols=(11)))
    gamma1.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                             skiprows=skiprows, usecols=(12)))
    e_int.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                       skiprows=skiprows, usecols=(13)))
    pvisc.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p2", 
                 skiprows=skiprows, usecols=(14)))
    
    # Load .p3 file (other thermodynamic quantities)
    t_conv.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(21)))
    v_conv.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(22)))
    F_conv.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(23)))
    F_rad.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(24)))
    e_nucl.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(25)))
    e_nupla.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(26)))
    e_grav.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(27)))
    e_nunucl.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(28)))
    D_conv.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(29)))
    D_turb.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p3", 
                       skiprows=skiprows, usecols=(30)))

    # Load .p4-9 file (mass fractions)
    neutron.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                              skiprows=skiprows, usecols=(1)))
    H1.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                         skiprows=skiprows, usecols=(2)))
    H2.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                         skiprows=skiprows, usecols=(3)))
    He3.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(4)))
    He4.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(5)))
    Li6.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(6)))
    Li7.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(7)))
    Be7.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(8)))
    B8.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                         skiprows=skiprows, usecols=(9)))
    Be9.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p4", 
                          skiprows=skiprows, usecols=(10)))
    B10.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(1)))
    B11.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(2)))
    C12.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(3)))
    C13.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(4)))
    N13.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(5)))
    C14.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(6)))
    N14.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(7)))
    N15.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(8)))
    O15.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(9)))
    O16.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p5", 
                          skiprows=skiprows, usecols=(10)))
    O17.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                          skiprows=skiprows, usecols=(1)))
    O18.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                          skiprows=skiprows, usecols=(2)))
    F18.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                          skiprows=skiprows, usecols=(3)))
    F19.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                          skiprows=skiprows, usecols=(4)))
    F20.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                          skiprows=skiprows, usecols=(5)))
    Ne20.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                           skiprows=skiprows, usecols=(6)))
    Ne21.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                           skiprows=skiprows, usecols=(7)))
    Ne22.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                           skiprows=skiprows, usecols=(8)))
    Na22.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                           skiprows=skiprows, usecols=(9)))
    Ne23.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p6", 
                           skiprows=skiprows, usecols=(10)))
    Na23.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(1)))
    Na24.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(2)))
    Mg24.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(3)))
    Na25.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(4)))
    Mg25.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(5)))
    Mg26.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(6)))
    Al26m.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                            skiprows=skiprows, usecols=(7)))
    Al26g.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                            skiprows=skiprows, usecols=(8)))
    Mg27.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(9)))
    Al27.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p7", 
                           skiprows=skiprows, usecols=(10)))
    Si28.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                           skiprows=skiprows, usecols=(1)))
    Si29.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                           skiprows=skiprows, usecols=(2)))
    Si30.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                           skiprows=skiprows, usecols=(3)))
    P31.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(4)))
    S32.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(5)))
    S33.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(6)))
    S34.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(7)))
    S35.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(8)))
    Cl35.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                           skiprows=skiprows, usecols=(9)))
    S36.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p8", 
                          skiprows=skiprows, usecols=(10)))
    Cl36.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p9", 
                           skiprows=skiprows, usecols=(1)))
    Cl37.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p9", 
                           skiprows=skiprows, usecols=(2)))
    heavy.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p9", 
                            skiprows=skiprows, usecols=(3)))
    sumX.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p9", 
                           skiprows=skiprows, usecols=(4)))
    
    compression_rate.insert(-1, 
                            np.loadtxt(fname=loading_path+profile_name+".p11", 
                                       skiprows=skiprows, usecols=(1)))
    heat_rate.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p11", 
                                    skiprows=skiprows, usecols=(2)))
    work_e_grav.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p11", 
                                      skiprows=skiprows, usecols=(3)))
    heat_e_grav.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p11", 
                                      skiprows=skiprows, usecols=(4)))
    dlum_dt.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p11", 
                                  skiprows=skiprows, usecols=(5)))
    
    Hp.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(1)))
    gmr.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(2)))
    phi_eos.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(5)))
    delta_eos.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(6)))
    N2_thermal.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(7)))
    N2_chemical.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(8)))
    N2.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p13", 
                           skiprows=skiprows, usecols=(9)))
    
    # Compute other quantities
    tau_index = 0
    for i in range(len(tau[n])):
        if abs(tau[n][i]-2/3)<abs(tau[n][tau_index]-2/3):
            tau_index = i
    T_eff.insert(-1, T[n][tau_index])
    n_eff.insert(-1, tau_index)
    kappa_e.insert(-1, 0.2005*(1.0+H1[n]+H2[n]))
    metallicity.insert(-1, 1-(H1[n]+H2[n]+He3[n]+He4[n]))
    mass.insert(-1, np.loadtxt(fname=loading_path+profile_name+".p0", 
                      skiprows=skiprows, usecols=(1))*xmr[n])
    dm_tmp = [0.]
    for i in range(1,len(mass[n])):
        dm_tmp.append(mass[n][i]-mass[n][i-1])
    dm.insert(-1, dm_tmp)
    lum.insert(-1, 16*5.6704e-5*6.67259e-8*T[n]**4*mass[n]*1.989e33*abrad[n]/ \
               (3*kappa[n]*P[n])/3.846e33*4*np.pi+1e-5)
    lum_rad.insert(-1, 16*5.6704e-5*6.67259e-8*T[n]**4*mass[n]*1.989e33* \
                   abla[n]/(3*kappa[n]*P[n])/3.846e33*4*np.pi+1e-5)
    L_edd.insert(-1, 4*np.pi*2.99792458e10*6.67259e-8*mass[n]/kappa[n]* \
                 1.989e33/3.846e33+1e-5)
    cg.insert(-1, np.sqrt(P[n]/rho[n]))
    tau_crit.insert(-1, 2.99792458e10/cg[n])
    tau_0.insert(-1, kappa[n]*rho[n]*Hp[n])
    

# =============================================================================

def change_color():
    
    ''' Change the color of input specified model. '''
    
    global colors
    
    n = int(input("\n profile number : "))
    color = input(" new color : ")
    if color=="":
        color = base_colors[(n-1)%len(base_colors)]
    colors[n-1] = color
        
# =============================================================================

def change_label():
    
    ''' Change the label of input specified model. '''
    
    global labels
    
    n = int(input("\n profile number : "))
    label = input(" new label : ")
    if label=="":
        label = profile_names[n-1]
    elif label=="None" or label=="none":
        label = None
    labels[n-1] = label
        
# =============================================================================

def change_linestyle():
    
    ''' Change the linestyle of an input specified model. '''
    
    global linestyles
    
    n = int(input("\n model number : "))
    linestyle = input(" new linestyle : ")
    if linestyle=="":
        linestyle = "solid"
    linestyles[n-1] = linestyle
    
# =============================================================================

def change_all():
    
    ''' Change all the features of one model. '''
    
    n = int(input("\n model number : "))
    color = input(" new color : ")
    if color=="":
        color = base_colors[(n-1)%len(base_colors)]
    colors[n-1] = color
    label = input(" new label : ")
    if label=="":
        label = profile_names[n-1]
    elif label=="None" or label=="none":
        label = None
    labels[n-1] = label
    linestyle = input(" new linestyle : ")
    if linestyle=="":
        linestyle = "solid"
    linestyles[n-1] = linestyle
 
# =============================================================================
        
def add_model():
    
    '''  Load an input specified supplementary model. '''
    
    global N
    
    N = N+1
    profile_name = input(f"\n profile name {N} : ")
    color_tmp = input(f" color {N} : ")
    if color_tmp=="":
        color_tmp = base_colors[(N-1)%len(base_colors)]
    label_tmp = input(f" label {N} : ")
    if label_tmp=="":
        label_tmp = profile_name
    profile_names.append(profile_name)
    colors.append(color_tmp)
    labels.append(label_tmp)
    
    try:
        _load(N-1)
        print("\n Model added succesfully.")
    except:
        N = N-1
        print("\n Error when loading the file, model not added.")
       
# =============================================================================

def infos():
    
    ''' Print informations on all loaded models. '''

    for n in range(N):
        print(f"\n model number : {n+1}")
        print(" profile name :", profile_names[n])
        print(" color        :", colors[n])
        print(" label        :", labels[n])
        print(" linestyle    :", linestyles[n])
        
# =============================================================================
#   Define mathematical operations on quantities functions
# =============================================================================

def add(x,y):
    
    '''
    Return a quantity corresponding to the quantity x added to the quantity
    y.
    '''
    
    add_x_y = []
    for n in range(N):
        add_x_y.append(x[n]+y[n])
    string_x = x[-1]
    string_y = y[-1]
    add_x_y.append(string_x+"+"+string_y)
    return add_x_y

# =============================================================================

def minus(x):
    
    '''
    Return a quantity corresponding to the quantity x multiplied by -1.
    '''
    
    minus_x = []
    for n in range(N):
        minus_x.append(-1*x[n])
    string = x[-1]
    minus_x.append("-"+string)
    return minus_x

# =============================================================================

def multiply(x,y):
    
    '''
    Return a quantity corresponding to the quantity x multiplied by the
    quantity y.
    '''
    
    multiply_x_y = []
    for n in range(N):
        multiply_x_y.append(x[n]*y[n])
    multiply_x_y.append(x[-1]+"*"+y[-1])
    return multiply_x_y

# =============================================================================

def divide(x,y):
    
    '''
    Return a quantity corresponding to the quantity x divided by the quantity
    y.
    '''
    
    divide_x_y = []
    for n in range(N):
        divide_x_y.append(x[n]/y[n])
    string_x = x[-1]
    string_y = y[-1]
    divide_x_y.append(string_x+"/"+string_y)
    return divide_x_y

# =============================================================================

def log10(x):
    
    '''
    Return a quantity corresponding to the decimal logarithm of the quantity x.
    '''

    log10_x = []
    for n in range(N):
        log10_x.append(np.log10(x[n]))
    string = x[-1]
    log10_x.append(r"log("+string+r")")
    return log10_x

# =============================================================================

def absol(x):
    
    '''
    Return a quantity corresponding to the absolute value of the quantity x.
    '''

    abs_x = []
    for n in range(N):
        abs_x.append(np.abs(x[n]))
    string = x[-1]
    for i in range(len(string)):
        if string[i-1:i+1] == " $":
            break
    if i != len(string)-1:
        abs_x.append(r"|"+string[:i-1]+r"|"+string[i-1:])
    else:
        abs_x.append(r"|"+string+r"|")
    return abs_x

# =============================================================================

def derivative(x,y):
    
    '''
    Return a quantity corresponding to the derivative of x with respect to y.
    '''
    
    dx_dy = []
    for n in range(N):
        tmp_dx_dy = np.zeros_like(x[n])
        for i in range(1,len(x[n])):
            if (y[n][i]-y[n][i-1])!=0:
                tmp_dx_dy[i] = (x[n][i]-x[n][i-1])/(y[n][i]-y[n][i-1])
        tmp_dx_dy[0] = tmp_dx_dy[1]
        dx_dy.append(tmp_dx_dy)
    dx_dy.append("d("+x[-1]+")/d("+y[-1]+")")
    return dx_dy

# =============================================================================

def cumulative(x):
    
    '''
    Return a quantity corresponding to the cumulative value of x.
    '''
    
    cumulative_x = []
    for n in range(N):
        tmp_cumulative_x = []
        tmp_cumulative_x_value = 0
        for i in range(len(x[n])):
            tmp_cumulative_x_value += x[n][i]
            tmp_cumulative_x.append(tmp_cumulative_x_value)
        cumulative_x.append(tmp_cumulative_x)
    cumulative_x.append("cumulative "+x[-1])
    return cumulative_x

# =============================================================================
#   Load files
# =============================================================================

# Get start time to compute total runtime
start_time = time()

# Main loop
try:
    for n in range(N):
        _load(n)
        color_tmp = base_colors[n%len(base_colors)]
        label_tmp = profile_names[n]
        colors.append(color_tmp)
        labels.append(label_tmp)
        linestyles.append("solid")
except:
    N = 0
    print("\n Error when loading the files.")

end_time = time()
print("\n Total loading time : {:.1f} s".format(end_time-start_time))

# End of loading, console display
if N>0:
    print("\n STAREVOL results loaded succesfully.")
else :
    print("\n Exiting load_structure_STAREVOL without loading STAREVOL"+\
          " results.")

# =============================================================================
#   Acknowledgements
# =============================================================================

# This software uses the following Python packages : 

# NumPy : Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., 
#         Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., 
#         Smith, N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, 
#         M. H., Brett, M., Haldane, A., del Rio, J. F., Wiebe, M., 
#         Peterson, P., Gérard-Marchant, P., Sheppard, K., Reddy, T.,
#         Weckesser, W., Abbasi, H., Gohlke, C., and Oliphant, T. E. 
#         (2020). Array programming with NumPy. Nature, 585(7825):357–362

# matplotlib : Hunter, J. D. (2007). Matplotlib: A 2D Graphics Environment. 
#              Computing in Science and Engineering, 9(3):90–95.
        
# =============================================================================
#   End
# =============================================================================