# Stefan-Boltzmann constant (sigma = 5.67×10^−8 W m^−2 K^−4)
sigma = 5.67e-8


def ebm0(S=1361, alpha=0.3, epsilon=1):
    """
    EBM0 - Calculates the equilibrium temperature T
    by solving the planetary radiation balance equation:

    S/4 * (1-alpha) = epsilon * sigma * T^4

    The user is allowed to change:
    total solar irradiance (S), albedo (alpha) and emissivity (epsilon)

    sigma is the Stefan-Boltzmann constant

    Created by Anders Moberg, June 2013
    anders.moberg@natgeo.su.se
    Translated and edited by Sara Blichner, January 2023
    sara.blichner@gmail.com

    :param S: Total solar irradiance. 1361 Wm^2 is the new TSI value from Kopp & Lean, GRL, 2011, doi:10.1029/GL045777
    :param alpha: 0.3 standard value of Earth's albedo
    :param epsilon: 1 blackbody emissivity = 1
    :return:
  """
    print('Calculating temperature for equilibrium conditions: \n'
          'Incoming solar radiation = Outgoing terrestrial radiation \n'
          'With equation: \n'
          'S/4 * (1-alpha) = epsilon * sigma * T^4 \n'
          f'where sigma = 5.67x10^(-8), S={S} Wm^2, alpha={alpha}\n'
          f'and epsilon = {epsilon}')
    T_K = (1 / sigma * S / 4 * (1 - alpha) / epsilon) ** .25
    # Convert to Celcius:
    T_C = T_K - 273.15
    print(f'The calculated temperature is {T_C} degrees C')
    return T_C
