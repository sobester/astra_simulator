"""
example_forecast.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Example: Forecast based Simulation

University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from datetime import datetime, timedelta

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
simFlight.balloonModel = 'TA800'
simFlight.nozzleLift = 1                                # kg
simFlight.payloadTrainWeight = 0.433                    # kg
simFlight.parachuteModel = 'SPH36'
simFlight.numberOfSimRuns = 10
simFlight.trainEquivSphereDiam = 0.1                    # m
simFlight.floatingFlight = False
simFlight.floatingAltitude = 30000                      # m
simFlight.excessPressureCoeff = 1
#simFlight.maxFlightTime = 5*60*60
simFlight.outputFile = '/Users/nick/Desktop/Simulation_Results'


# Run the simulation
simFlight.run()