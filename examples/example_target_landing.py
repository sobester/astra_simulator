# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-04 18:21:52
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-06-19 11:33:41
"""
example_forecast.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Example: Forecast based Simulation

University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""
import logging
import numpy as np
import os
from astra.target_landing import targetFlight
from datetime import datetime, timedelta

logging.basicConfig(level=logging.DEBUG)


if __name__ == "__main__":
    from astra.simulator import *

    np.random.seed(62)

    # Environment parameters
    launch_datetime = datetime.now() + timedelta(days=1)
    # launch site in Stoney Cross, New Forest
    launchLat, launchLon, launchElev = (50.903824, -1.63697, 114.0)
    # Lat lon for somewhere near winchester.
    targetLat = 51.077214
    targetLon = -1.1866875398775423
    targetElev = 137

    # Launch setup
    simulator = targetFlight(start_dateTime=launch_datetime,
                             targetLat=targetLat,
                             targetLon=targetLon,
                             targetElev=targetElev,
                             launchSites=[(launchLat, launchLon, launchElev)],
                             balloonModel='TA100',
                             balloonGasType="Helium",
                             parachuteModel=None,
                             nozzleLift=1,
                             trainEquivSphereDiam=0.1,
                             inflationTemperature=0.0,
                             payloadTrainWeight=0.38,
                             windowDuration=72,
                             HD=False,
                             maxFlightTime=18000,
                             debugging=True,
                             log_to_file=False,
                             progress_to_file=False,
                             outputFile=os.path.join(''))


    # Run the simulation
    bestProfile_bf, X, Y, distances = simulator.bruteForce()
