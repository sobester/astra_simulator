# coding=utf-8

"""
available_balloons.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Available Balloons
This file contains key data of all the balloons currently supported by the simulator.
The dictionary returns a tuple of four items: (weight, nominal burst diameter, weibull lambda, weibull k).

Warning: this module is used privately and should not be used per se.



University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 23/02/2013

Data courtesy of Dr Andras Sobester
"""
__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

# Which balloon weights are supported
# Inside the tuples, the first value represents the nominal burst diameter, the second value the Weibull lambda
# and the third value the Weibull k.
balloons = {
    "TA10": (0.01, 0.45, 0.5051, 14.3577),
    "TA20": (0.02, 0.7, 0.7857, 14.3577),
    "TA30": (0.03, 0.88, 0.9877, 14.3577),
    "TA45": (0.045, 1.1, 1.2346, 14.3577),
    "TA100": (0.1, 1.96, 2.1999, 14.3577),
    "TA200": (0.2, 3, 3.3672, 14.3577),
    "TA300": (0.3, 3.78, 4.2427, 14.3577),
    "TA350": (0.35, 4.12, 4.6243, 14.3577),
    "TA450": (0.45, 4.72, 5.2977, 14.3577),
    "TA500": (0.5, 4.99, 5.6008, 14.3577),
    "TA600": (0.6, 6.02, 6.7569, 14.3577),
    "TA700": (0.7, 6.53, 7.3293, 14.3577),
    "TA800": (0.8, 7, 7.8568, 14.3577),
    "TA1000": (1.0, 7.86, 8.8221, 14.3577),
    "TA1200": (1.2, 8.63, 9.6863, 14.3577),
    "TA1500": (1.5, 9.44, 10.5955, 14.3577),
    "TA2000": (2, 10.54, 11.8301, 14.3577),
    "TA3000": (3, 13, 14.5913, 14.3577)
}

meanToNominalBurstRatio = 1.08116