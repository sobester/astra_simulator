# coding=utf-8

"""
sim_manager.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Simulator Manager
This is a command line interface for the simulator system.
The job of the simulator manager is to collect all the input parameters, set up the necessary objects for the
simulation according to the flight type (standard or floating), initialize all the objects correctly and run the
simulation.


USAGE
--------------

python sim_manager.py [arguments]

Arguments are:
    -x [longitude], required
    -y [latitude], required
    -z [elevation], required
    -t [date and time], required.  NOTE: Format should be exactly 'dd/mm/yyyy HH:mm +0000', where +0000 indicated the
        time zone offset from UTC (for example, Central European Summer time is +0200). This is important in order to
        download the correct weather forecast
    -b [balloon weight in kg], required
    -g ['Helium' or 'Hydrogen'], required
    -n [nozzle lift in kg], required
    -p [payload weight in kg], required
    -c [parachute type (3 or 4)], required
    -h [train equivalent sphere diameter in m], required
    -r [number of sim runs], required
    -o [output filename], required
    --noHD forces the forecast to be downloaded at standard definition, optional
    --debug activates detailed output, optional

    Additionally, the following can be used for a sounding-based atmosphere:
        -s|--sounding [file path] -e [distance from sounding] -i [time from sounding] -k [temperature]

    Or, for a floating flight:
        -f|--float -a [altitude in m] -l [maximum flight time limit in hours]


Type python sim_manager.py --help for documentation.

University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 12/02/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

import sys
import getopt
from datetime import datetime

from simulator import *


def main(argv):
    global debug

    soundingFile = None

    try:
        opts, args = getopt.getopt(argv, "x:y:z:t:b:g:n:p:c:o:r:s:e:i:k:fa:h:l:",
                                   ["sounding=", "float", "help", "debug", "noHD"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    if '--debug' in argv:
        debug = True
    else:
        debug = False

    # Initialize objects
    soundingWxr = soundingEnvironment(debugging=debug, log_to_file=True)
    forecastWxr = forecastEnvironment(debugging=debug, log_to_file=True)
    simFlight = flight(debugging=debug, log_to_file=True, progress_to_file=True)

    # Read arguments and assign them to the appropriate object parameters
    for opt, arg in opts:
        if opt in '-x':
            forecastWxr.launchSiteLon = float(arg)
            soundingWxr.launchSiteLon = float(arg)
        elif opt in '-y':
            forecastWxr.launchSiteLat = float(arg)
            soundingWxr.launchSiteLat = float(arg)
        elif opt in '-z':
            forecastWxr.launchSiteElev = float(arg)
            soundingWxr.launchSiteElev = float(arg)
        elif opt in '-t':
            forecastWxr.UTC_offset = int(arg[-5:-2]) + float(arg[-2:]) / 60
            forecastWxr.dateAndTime = datetime.strptime(arg[:-6], '%d/%m/%Y %H:%M')
            soundingWxr.dateAndTime = datetime.strptime(arg[:-6], '%d/%m/%Y %H:%M')
        elif opt in '-b':
            simFlight.balloonWeight = float(arg)
        elif opt in '-g':
            simFlight.balloonGasType = arg
        elif opt in '-n':
            simFlight.nozzleLift = float(arg)
        elif opt in '-p':
            simFlight.payloadTrainWeight = float(arg)
        elif opt in '-c':
            simFlight.parachuteType = int(arg)
        elif opt in '-h':
            simFlight.trainEquivSphereDiam = float(arg)
        elif opt in '-r':
            simFlight.numberOfSimRuns = int(arg)
        elif opt in ('-s', '--sounding'):
            soundingFile = arg
        elif opt in '-k':
            soundingWxr.inflationTemperature = float(arg)
        elif opt in '-e':
            soundingWxr.distanceFromSounding = float(arg)
        elif opt in '-i':
            soundingWxr.timeFromSounding = float(arg)
        elif opt in ('-f', '--float'):
            simFlight.floatingFlight = True
        elif opt in '-a':
            simFlight.floatingAltitude = float(arg)
        elif opt in '-o':
            simFlight.outputFile = arg
        elif opt in '-l':
            simFlight.maxFlightTime = int(arg) * 3600.
            forecastWxr.maxFlightTime = simFlight.maxFlightTime
        elif opt in '--noHD':
            forecastWxr.forceNonHD = True
        elif opt in '--help':
            usage()
            sys.exit(0)
        elif opt in '--debug':
            continue
        else:
            usage()
            sys.exit(2)

    # Select if simulating a flight based on forecast or sounding
    if soundingFile is None:
        simFlight.environment = forecastWxr
    else:
        soundingWxr.loadSounding(soundingFile)
        simFlight.environment = soundingWxr

    # Run the simulation
    simFlight.run()


def usage():
    print __doc__


if __name__ == "__main__":
    main(sys.argv[1:])