"""
example_sounding.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Example: Sounding based Simulation

University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, FINAL DATE GOES HERE
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from datetime import datetime

from simulator import *


simEnvironment = soundingEnvironment(debugging=False,log_to_file=False)
simFlight = flight(debugging=False,log_to_file=False)

# Environment parameters
# Launch site: University of Southampton
#        time: 3 hours from sounding's time
simEnvironment.launchSiteLat = 50.2245                  # deg
simEnvironment.launchSiteLon = -5.3069                  # deg
simEnvironment.launchSiteElev = 60                      # m
simEnvironment.distanceFromSounding = 0                 # km
simEnvironment.timeFromSounding = 3                     # hours
simEnvironment.inflationTemperature = 10.5              # degC
simEnvironment.dateAndTime = datetime.now()
simEnvironment.UTC_offset = 0
simEnvironment.loadSounding('/path/to/myFile.sounding')

# Launch parameters
simFlight.environment = simEnvironment
simFlight.balloonGasType = 'Helium'
simFlight.balloonWeight = 0.8                           # kg
simFlight.nozzleLift = 1                                # kg
simFlight.payloadTrainWeight = 0.433                    # kg
simFlight.parachuteType = 3
simFlight.numberOfSimRuns = 10
simFlight.trainEquivSphereDiam = 0.1                    # m
simFlight.floatingFlight = False
simFlight.floatingAltitude = 30000                      # m
simFlight.excessPressureCoeff = 1
#simFlight.maxFlightTime = 5*60*60
simFlight.outputFile = '/path/to/Simulation_Results'


# Run the simulation
simFlight.run()