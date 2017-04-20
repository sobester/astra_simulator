# coding=utf-8

"""
A collection of small flight-specific tools used throughout the simulation.

These tools generally perform basic performance calculations, such as balloon
diameter and drag.

University of Southampton
"""
from math import pi


# Constants
MOLECULAR_MASS = {'air': 0.02896, # kg/mol
                  'helium': 0.004002602, # kg/mol
                  'hydrogen': 0.0020159 # kg/mol
                 }

# Not the same as the pure gas molecular mass, due to mixing fractions
MIXEDGAS_MOLECULAR_MASS = {
    "Helium" : (0.945 * MOLECULAR_MASS['helium'] +
                0.055 * MOLECULAR_MASS['air']),
    "Hydrogen": (0.985 * MOLECULAR_MASS['hydrogen'] +
                0.015 * MOLECULAR_MASS['air'])
    }

def liftingGasMass(nozzleLift, balloonMass, ambientTempC,
    ambientPressMbar, gasMolecularMass, excessPressureCoefficient):
    """
    Calculates the gas mass, balloon volume and balloon diameter.

    Parameters
    ----------
    nozzleLift : scalar (kg)
    balloonMass : scalar (kg)
    ambientTempC : scalar (Celsius)
    ambientPressMbar : scalar (millibar)
    gasMolecularMass : scalar (kg/mol)
    excessPressCoefficient : scalar (dimension-less)

    Returns
    -------
    gas_mass, balloon_volume, balloon_diameter : tuple (float, float, float
    """
    R = 8.31447 # m^3Pa/mol K
    gasDensity = excessPressureCoefficient * ambientPressMbar * 100 * gasMolecularMass / (
        R * (ambientTempC + 273.15))
    airDensity = ambientPressMbar * 100 * MOLECULAR_MASS['air'] / (R * (ambientTempC + 273.15))

    balloonVolume = (nozzleLift + balloonMass) / (airDensity - gasDensity)
    gasMass = balloonVolume * gasDensity
    balloonDiameter = (6. * balloonVolume / pi) ** (1. / 3)

    return [gasMass, balloonVolume, balloonDiameter]


def gasMassForFloat(currentAltitude, floatingAltitude,
    gasMassAtInflation, gasMassAtFloatingAltitude, ventStart=500):
    """
    Returns the gas mass profile to simulate the valves venting air to reach
    floating conditions.

    This method should only be used for floating flights.

    Parameters
    ----------
    currentAltitude : scalar
        the altitude for which the gas mass is required [m]
    floatingAltitude : scalar
        the target altitude for float [m]
    gasMassAtInflation : scalar
        the initial gas mass, before venting [kg]
    gasMassAtFloatingAltitude : scalar
        the final gas mass required to reach floating conditions [kg]
    ventStart : scalar
        how many meters below the target floating altitude should the venting
        start. Default is 500m

    Returns
    -------
    gasMassAtFloatingAltitude : float
    """
    if currentAltitude < (floatingAltitude - ventStart):
        # Below the threshold to start venting. The gas mass is unchanged
        return gasMassAtInflation
    elif currentAltitude < floatingAltitude:
        # The balloon is approaching the target altitude, start venting

        # The following is the equation for linear interpolation between the
        # start and the end of the venting.
        return (gasMassAtInflation - gasMassAtFloatingAltitude) / float(ventStart) * (
            floatingAltitude - currentAltitude) + gasMassAtFloatingAltitude
    else:
        # The balloon has reached the floating altitude
        return gasMassAtFloatingAltitude


def nozzleLiftForFloat(nozzleLiftAtInflation, airDensity, gasDensity,
    balloonVolume, balloonMass, currentAltitude, floatingAltitude,
    ventStart=500):
    """Returns the nozzle lift of a balloon in floating configuration.

    This method should only be used for floating flights. Below the venting
    threshold (floatingAltitude-ventStart), the nozzle lift is unchanged.
        Otherwise, it's recalculated.

    Parameters
    ----------
    nozzleLiftAtInflation : scalar
        nozzle lift of the balloon at inflation [kg]
    airDensity : scalar
        density of air at current altitude [kg/m3]
    gasDensity : scalar
        density of the gas at current altitude [kg/m3]
    balloonVolume : scalar
        volume of the balloon at current altitude [m3]
    balloonMass : scalar
        total balloon mass at current altitude [kg]
    currentAltitude : scalar
        altitude at which the nozzle lift is required [m]
    floatingAltitude : scalar
        target altitude for floating conditions [m]
    ventStart : scalar (default 500)
        how many meters below the target floating altitude should the
        venting start. [m]

    Returns
    -------
    CurrentNozzleLift : float [kg].
    """
    if currentAltitude < (floatingAltitude - ventStart):
        # Below the threshold to start venting. The nozzle lift is unchanged.
        return nozzleLiftAtInflation
    else:
        # The balloon has started venting gas, nozzle lift needs to be recalculated.
        currentNozzleLift = (airDensity - gasDensity) * balloonVolume - balloonMass

        return currentNozzleLift


def balloonDrag(diameter, ascentRate, density, viscosity,
    lowCD, highCD, ReBand, transition):
    """Calculates and returns the balloon drag, using model from [1].

    Parameters
    ----------
    diameter : scalar. balloon diameter [m]
    ascent rate : scalar. balloon vertical velocity [m/s]
    density : scalar. current air density [kg/m3]
    viscosity : scalar. current air viscosity [kg/m3]
    lowCd : the lower 

    Returns
    -------
    balloon_drag : float [Newtons].

    References
    ----------
    [1] [1] Sobester, A., Czerski, H.,  Zapponi, N., and Castro, I. P., "High
    Altitude Gas Balloon Trajectory Prediction - a Monte Carlo Model". AIAA
    Journal, 52, (4), 2014, pp. 832-842.
    (doi:10.2514/1.J052900 <http://dx.doi.org/10.2514/1.J052900>).

    """
    reynoldsNumber = density * abs(ascentRate) * diameter / float(viscosity)

    # Balloon CD
    if reynoldsNumber < transition * 1e5:
        balloonCD = highCD
    elif reynoldsNumber < (transition + ReBand) * 1e5:
        balloonCD = highCD - (highCD - lowCD) * (
            (reynoldsNumber - transition * 1e5) / float(ReBand * 1e5) )
    else:
        balloonCD = lowCD

    return 0.5 * density * (abs(ascentRate) ** 2) * (pi * (diameter ** 2) / 4.) * balloonCD


def parachuteDrag(descentRate, density, parachuteAref, parachuteCD):
    """
    Calculates and returns the drag generated by the parachute on descent.

    Parameters
    ----------
    descentRate : float
        balloon vertical velocity [m/s]
    density : float
        current air density [kg/m3]

    Returns
    -------
        Returns a float with the parachute drag in [N].
    """
    return 0.5 * density * descentRate ** 2 * parachuteAref * parachuteCD
