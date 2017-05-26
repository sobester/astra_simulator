# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-04 18:21:52
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-05-09 19:31:10
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
logging.basicConfig(level=logging.DEBUG)

if __name__ == "__main__":
    from datetime import datetime, timedelta
    from astra.simulator import *

    np.random.seed(62)

    # Environment parameters
    # Launch site: Daytona Beach, FL
    #        time: tomorrow, this time
    # Note that the forecast duration is set to 24 hours
    launch_datetime = datetime.now() + timedelta(days=1)
    simEnvironment = forecastEnvironment(launchSiteLat=29.2108,      # deg
                                         launchSiteLon=-81.0228,     # deg
                                         launchSiteElev=4,           # m
                                         dateAndTime=launch_datetime,
                                         forceNonHD=True,
                                         debugging=True,
                                         forecastDuration=24)


    # Launch setup
    # simFlight = flight(environment=simEnvironment,
    #                    balloonGasType='Helium',
    #                    balloonModel='TA800',
    #                    nozzleLift=1,                                # kg
    #                    payloadTrainWeight=0.433,                    # kg
    #                    parachuteModel='SPH36',
    #                    numberOfSimRuns=10,
    #                    trainEquivSphereDiam=0.1,                    # m
    #                    floatingFlight=False,
    #                    excessPressureCoeff=1,
    #                    outputFile=os.path.join('.', 'astra_output'),
    #                    debugging=True,
    #                    log_to_file=True)

    # simFlight.maxFlightTime = 5*60*60

    # Run the simulation
    # simFlight.run()
