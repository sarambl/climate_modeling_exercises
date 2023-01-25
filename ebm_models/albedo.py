import numpy as np


def albedo(T, jmx, x, NoAlbFB_Flag):
    """

 :param T: Temperature
 :param jmx: number of points along latitude
 :param x: latitude
 :param NoAlbFB_Flag: if True, ignores albedo feedback
 :return: albedo
 """
    alb = np.ones([jmx]) * 0.3

    # alternative albedo that depends on latitude (zenith angle)
    # alb=0.31+0.08*(3*x^2-1)/2
    if NoAlbFB_Flag:
        k = np.argwhere(np.abs(x) >= 0.95)
    else:
        k = np.argwhere(T <= -10)

    alb[k] = 0.6

    return alb

# %%
