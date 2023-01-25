import numpy as np


def setupfastM(delx, jmx, D, B, Cl, delt):
    """
    Setup inverse for solving
    :param delx: delta latitude
    :param jmx: number of points along latitude
    :param D: Diffusion?
    :param B: ORL coef.
    :param Cl:heat capacity over land.
    :param delt: time step in fraction of year
    :return:
    """
    # set up lambda array.

    # lam=(1-[-1:delx:1]'.^2)/delx^2;
    lam = (1 - np.arange(-1, 1 + delx, delx) ** 2) / delx ** 2
    lam[-1] = 0
    lam = D[:] * lam[:]
    M = np.zeros([jmx, jmx])
    M[0, 0] = -B - lam[1]
    M[0, 1] = lam[1]

    M[-1, -2] = lam[-2]
    M[-1, -1] = -B - lam[-2]
    for j in range(1, jmx - 1):
        M[j, j - 1] = lam[j]
        M[j, j] = -B - (lam[j + 1] + lam[j])
        M[j, j + 1] = lam[j + 1]
    # M[-2,-1] = lam[-1]

    # add in heat capacities
    M = M / Cl

    # calculate the inverse of M, the matrix operator.
    Mh = M.copy()
    M = 0.5 * M
    for j in range(0, jmx):
        M[j, j] = M[j, j] - 1. / delt

    invM = np.linalg.inv(M, )
    return [invM, Mh]
