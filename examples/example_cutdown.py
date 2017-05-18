# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-18 17:55:28
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-05-18 18:15:48
if __name__ == "__main__":
    from datetime import datetime, timedelta
    import sys
    import os
    sys.path.append(os.path.abspath('../astra'))
    from simulator import *


    simEnvironment = forecastEnvironment(debugging=True, log_to_file=False)
    simFlight = flight(debugging=True, log_to_file=False)

    # Environment parameters
    # Launch site: Daytona Beach, FL
    #        time: tomorrow, this time
    simEnvironment.launchSiteLat = 29.2108                  # deg
    simEnvironment.launchSiteLon = -81.0228                 # deg
    simEnvironment.launchSiteElev = 4                       # m
    simEnvironment.dateAndTime = datetime.now() + timedelta(days=1)
    simEnvironment.forceNonHD = True

    # The forecast download is automatically called by the flight object below.
    # However, if you'd like to download it in advance, uncomment the following line.
    # The flight object will automatically detect it and won't download the forecast twice.
    # simEnvironment.loadForecast()

    # Launch parameters
    simFlight.environment = simEnvironment
    simFlight.balloonGasType = 'Helium'
    simFlight.balloonModel = 'TA100'
    simFlight.nozzleLift = 0.6                                # kg
    simFlight.payloadTrainWeight = 0.38                    # kg
    simFlight.parachuteModel = 'SPH36'
    simFlight.numberOfSimRuns = 10
    simFlight.trainEquivSphereDiam = 0.1                    # m
    simFlight.cutdown = True
    simFlight.cutdownAltitude = 14000
    simFlight.excessPressureCoeff = 1
    #simFlight.maxFlightTime = 5*60*60
    simFlight.outputFile = os.path.join('.', 'astra_output_cutdown')


    # Run the simulation
    simFlight.run()
