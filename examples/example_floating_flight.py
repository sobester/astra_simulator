import logging
logging.basicConfig(level=logging.DEBUG)
from datetime import datetime, timedelta
from astra.simulator import *


if __name__ == "__main__":
    # Environment parameters
    # Launch site: Daytona Beach, FL
    #        time: tomorrow, this time
    launch_datetime = datetime.now() + timedelta(days=1)
    simEnvironment = forecastEnvironment(launchSiteLat=29.2108,      # deg
                                         launchSiteLon=-81.0228,     # deg
                                         launchSiteElev=4,           # m
                                         dateAndTime=launch_datetime,
                                         forceNonHD=True,
                                         debugging=True)

    # The forecast download is automatically called by the flight object below.
    # However, if you'd like to download it in advance, uncomment the following
    # line. The flight object will automatically detect it and won't download
    # the forecast twice.
    # simEnvironment.loadForecast()

    # Launch setup
    simFlight = flight(environment=simEnvironment,
                       balloonGasType='Helium',
                       balloonModel='TA800',
                       nozzleLift=1,                                # kg
                       payloadTrainWeight=0.433,                    # kg
                       parachuteModel='SPH36',
                       numberOfSimRuns=10,
                       trainEquivSphereDiam=0.1,                    # m
                       floatingFlight=True,
                       floatingAltitude=30000,                      # m
                       excessPressureCoeff=1,
                       outputFile=os.path.join('.', 'astra_output_floating'),
                       debugging=True,
                       log_to_file=True,
                       floatDuration=30*60,  #seconds
                       cutdownTimeout=3*3600)


    #simFlight.maxFlightTime = 5*60*60


    # Run the simulation
    simFlight.run()
