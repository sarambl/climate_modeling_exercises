# Code up of North and Coakley's seasonal EBM model
# Simplified to eliminate the ocean domain
# Designed to run without a seasonal cycle and hence
# to stop once an equilibrium solution is reached.
# The model uses an implicit trapezoidal method
# so the timestep can be long.

import numpy as np
from ebm_models.albedo import albedo
import matplotlib.pyplot as plt
import pandas as pd
from ebm_models.setupfastM import setupfastM


# %%
# size of domain.

# Choose parameters.
# scaleQ = dQ/Q0
# A
def run_1d_ebm(D = 0.44,
               scaleQ=1,
               A=203.3,
               B=2.09,
               coldstartflag=False,
               hadleyflag=False,
               noalbedoflag=False,
               jmx=151,
               ):
    """

    :param D: heat diffusion coefficient
    :param scaleQ: dQ/Q0
    :param A:  OLR constant
    :param B:  OLR coef.
    :param coldstartflag: If start cold
    :param hadleyflag: If hadley cell simulation
    :param noalbedoflag: If no albedo feedback
    :param jmx: number of points along latitude
    :return:
    """
    """
    jmx = 151
    Dmag = 0.44;
    scaleQ = 1;
    A = 203.3;
    B = 2.09;
    coldstartflag = False

    hadleyflag = False
    albedoflag = False
    """
    print(f'Running 1-D model with settings: \n'
          f'A={A}, B={B}, Dmag={D} \n'
          f'scaleQ={scaleQ}\n')
    # heat diffusion coefficient.
    Toffset = 0.
    # if (exist('coldstartflag')==1);
    if coldstartflag:
        Toffset = -40
        # Uncomment for warm start
        # Toffset = +40

    # heat capacity over land.
    Cl = 0.2  # something small to make it equilibriate quickly

    # time step in fraction of year
    delt = 1. / 50
    NMAX = 1000

    # set up x array.
    delx = 2.0 / jmx

    x = np.arange(-1.0 + delx / 2, 1.0, delx)
    #  [-1.0+delx/2:delx:1.0-delx/2];
    phi = np.arcsin(x) * 180 / np.pi
    # obtain annual array of daily averaged-insolation.
    # [insol] = sun(x);
    # Legendre polynomial realizatin of mean annual insol.
    Q = 338.5
    S = Q * (1 - 0.241 * (3 * x ** 2 - 1))
    S = scaleQ * S
    # S=S[:]

    # set up inital T profile
    T = 20 * (1 - 2 * x ** 2)
    # T=T(:);
    T = T + Toffset
    # load T_final
    Tinit = T

    # setup D(x) if simulating the Hadley Cell
    # and calculate the matrix Mh and invM.
    if hadleyflag:
        # Simulate Hadley Cell with Lindzen and Farrell plan

        xmp = np.arange(-1, 1 + delx, delx)
        len(xmp)
        Diff = D * (1 + 9 * np.exp(-(xmp / np.sin(25 * np.pi / 180)) ** 6))
        [invM, Mh] = setupfastM(delx, jmx, Diff, B, Cl, delt)
    else:
        Diff = D * np.ones([jmx + 1])
        [invM, Mh] = setupfastM(delx, jmx, Diff, B, Cl, delt)
    print('')

    # Boundary conditions
    # Set up initial value for h.
    alb = albedo(T, jmx, x, noalbedoflag)
    len(alb)
    len(S)
    src = (1 - alb) * S / Cl - A / Cl
    h = np.matmul(Mh, T) + src
    # Global mean temperature
    # Note that the grid is set up so that the grid cells cover equal area
    Tglob = np.mean(T)
    #, weights=np.cos(x/180*np.pi)) # /np.mean(np.cos(x/180))

    # Timestepping loop
    for n in range(0, NMAX):

        Tglob_prev = Tglob

        # Calculate src for this loop.
        alb = albedo(T, jmx, x, noalbedoflag)
        src = ((1 - alb) * S - A) / Cl
        # src=src(:)

        # Calculate new T.
        T = -np.matmul(invM, (0.5 * (h + src) + T / delt))

        # Calculate h for next loop.
        h = np.matmul(Mh, T) + src

        # Check to see if global mean temperature has converged
        Tglob = np.mean(T)
        Tchange = Tglob - Tglob_prev
        if abs(Tchange) < 1.0e-12:
            break

    # save T_final.mat T

    # compute meridional heat flux and its convergence
    a = 6.37e+6  # earth radius in meters
    [invM, Mh] = setupfastM(delx, jmx, Diff, 0., 1.0, delt)
    Dmp = 0.5 * (Diff[1:jmx + 1] + Diff[0:jmx])
    divF = np.matmul(Mh, T)

    F = -2 * np.pi * a ** 2 * np.sqrt(1 - x ** 2) * Dmp * np.gradient(T, delx, edge_order=2)

    plot_results(A, B, D, F, S, T, Tglob, Toffset, alb, divF, phi, scaleQ)
    # %%
    output_df = pd.DataFrame(index=phi)
    output_df.index.name = 'latitude'
    output_df['T_C'] = T
    output_df[r'Delta_F'] = divF
    output_df['albedo'] = alb
    output_df['F'] = F

    return output_df
    # %%


def plot_results(A, B, Dmag, F, S, T, Tglob, Toffset, alb, divF, phi, scaleQ):
    fig, axs = plt.subplots(3, 1, figsize=[8, 8], sharex='col')
    ax = axs[0]
    ax.plot(phi, T, '-',marker='.', linewidth=1.5)
    # plot(phi,T,'.-','linewidth',1.5)
    ax.set_ylabel(r'Temperature [$^\circ$C]')
    # set(gca,'position',[0.1300    0.71    0.7750    0.21]);
    ax.set_title(f'Global mean temperature is {Tglob:.2f}')
    ax.grid()
    ax = axs[1]
    ax.plot(phi, F * 1e-15, '-', linewidth=1.5)
    ax.set_ylabel('Poleward Heat Flux [10$^{15}$ W]')
    ax.grid()

    ax = axs[2]
    ax.plot(phi, divF, '-', label=r'$\Delta$ F', linewidth=1.5)
    ax.plot(phi, (1 - alb) * S, '--', label='SWd', linewidth=1.5)
    ax.plot(phi, A + B * T, '-.', label='LWu', linewidth=1.5)
    #            pphi,A+B*T,'.','linewidth',1.5)
    ax.set_ylabel('Energy Balance Terms [W $^{-2}$]')
    ax.set_xlabel('latitude')
    ax.grid()

    plt.legend()
    fig.suptitle(f'D = {Dmag}, '
                 f' Q/Qo ={scaleQ},'
                 f'A = {A},'
                 f'B= {B}, '
                 f'Toffset = {Toffset}'
                 )
    plt.show()

# fast()
# %%
