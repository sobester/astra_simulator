# coding=utf-8

"""
simulator.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Simulator Module
This module is responsible for the whole flight simulation and the processing of results.
With all the flight parameters, the flight class defines the framework of operations necessary to accomplish
the balloon flight simulation.

Note: This should be the first and only module manually imported by the user. All the other necessary modules are
automatically imported by the simulator as needed.


USAGE
--------------

Instantiate a new flight object to deal with the whole simulation.
When a new flight object is instantiated, it is empty and needs to be configured before running the simulation.
This includes providing the object with all the flight parameters and environmental data required for the flight
simulation to run.

Upon instantiation, the following three optional parameters can be passed to the object:
    debugging           (type: bool)
    log_to_file         (type: bool)
    progress_to_file    (type: bool)
These parameters are used to configure the debugging and error logging system as well as the progress updater.
By default, they are all False, which means that errors and progress information are displayed on the terminal or
command line. See the NOTE below for more information.

The standard syntax to instantiate and initialize a flight object is as follows:

    import simulator
    my_flight_simulation = simulator.flight()

Once the object is initialized, the flight can be configured.
The following parameters of the flight object need to be defined (parameters in [brackets] are optional):
    environment         an environment object (either soundingEnvironment or forecastEnvironment) already configured
        with the environmental model. See the Weather Module for more information.
    balloonGasType      'Helium'|'Hydrogen' (type: string)
    balloonModel        model of the balloon. Accepted values can be found in available_balloons_parachutes.py. (type: string)
    nozzleLift          nozzle lift in kg (type: float)
    payloadTrainWeight  weight of the payload train in kg (type: float)
    parachuteModel      model of the parachute. Accepted values can be found in available_balloons_parachutes.py. (type: string)
    trainEquivSphereDiam  the diameter (in meters) of sphere that generates the same drag as the train (type: float)
    numberOfSimRuns     total number of simulations to be run. If 1, a deterministic simulation is run. If >1,
        Monte Carlo perturbations are applied to each simulation (type: int)
    [excessPressureCoeff] the coefficient of excess pressure of the balloon on the gas inside. Default is 1. Currently
        not implemented (type: float)
    [floatingFlight]    TRUE if the balloon is going to vent air and float rather than burst. Default is FALSE (type: bool)
    [floatingAltitude]  the target altitude for the balloon to float, in meters. Ignored if floatingFlight is FALSE
        (type: float)
    [ventingStart]      how many meters below the target floating altitude should the venting start. Ignored if
        floatingFlight is FALSE. Default is 1000m (type: float)
    [maxFlightTime]     what is the maximum duration, in seconds, of a flight to be simulated. Warning: setting a high
        maxFlightTime during a floatingFlight could take a long time to simulate! Default is 14400 seconds (4 hours)
        (type: int)
    outputFile          the path of the output file containing all the simulation data. See NOTE for available formats.
        (type: string)

After all the flight parameters have been entered, the simulation is ready to be performed.

To simply run the simulation, the run() method can be used.
The standard syntax is as follows:

    my_flight_simulation.run()

The run() method takes all the necessary steps to run the simulation and wraps together all the simulation-specific
methods found below. For standard users, the use of the run() method is RECOMMENDED.
This method returns 0 if the simulation succeeded. If there were any errors, they will be displayed depending on how the
logging system has been set up. By default, they will be displayed on the terminal or command line.
See the NOTE for debugging and progress for more information.

For advanced users, the simulation can be performed using more specific methods.

These are the available methods to control the simulation:
    preflight()
        This runs a series of pre-flight checks and calculations to verify the consistency of the flight
        parameters entered and to prepare all the data required for simulation. It's MANDATORY to execute this method
        before the simulation (otherwise the simulation will throw an error)
        If successful, no errors are thrown. Enable debugging for detailed information.
    fly(flightNumber)
        This executes a single simulation. It should be run N times, where N is the number of simulation runs specified
        upon configuration in the numberOfSimRuns variable. flightNumber should have values 0...N-1
        If successful, no errors are thrown. Enable debugging for detailed information.
    postflight(outputFormat)
        After all the simulations have been executed, this method puts the results together, processes them and stores
        them in the outputFormat required. See NOTE for output formats.
    reset(keepParameters=False)
        This method resets the simulation. If keepParameters is FALSE, all the parameters are deleted and the flight
        object is ready to be reconfigured. If keepParameters is TRUE, all the results are deleted but the initial flight
        configuration is preserved. preflight() should still be run again before performing anymore simulations.
        Default is FALSE.
    updateProgress(value,action)
        This updates the progress file with both the value entered and the action being performed. See NOTE for debugging
        and progress (value types: value:float,action:int)


NOTE: AVAILABLE OUTPUT FORMATS
--------------

The format of the file to be generated depends on the extension of the file passed in the outputFile.
Accepted extensions are:
    'json'    JavaScript data structure, used to provide data to the web interface;
    'kml'     Standard format for geographical data. It can be opened by Google Maps, Google Earth, etc;
    'kmz'     Zipped kml, good for online use (eg Google Maps). The file size is significantly reduced;
    'csv'     Comma separated values, the preferred format for further data processing with any other software;
    'csv.zip' Zipped csv file;
    'web'     This stores json, kml and csv.zip in the same folder requested. This is used by the web interface to
              prepare files for export.
    none      If no extension is found, a new folder is created and ALL output formats are stored.


NOTE: DEBUGGING SYSTEM AND PROGRESS UPDATER
--------------

Both the debugging system and the progress updater are very easy to setup and use.
For the debugging system, here is how to do it:
When initializing the flight object, two optional boolean parameters can be passed:
    debugging,
    log_to_file.

If debugging is set to TRUE, all the information regarding the simulation will be logged. If it's set to FALSE, only
errors will be logged.
log_to_file determines the location where the errors and debug messages will be logged. If TRUE, an error.log file will
be created in your current folder and everything will go in there (make sure you have permissions to write in the
current folder, or the simulator will not work!)
If it's FALSE, all the errors and debug messages will be displayed on the terminal or command line.

The progress updater works pretty much in the same way: there is an optional boolean parameter than can be passed upon
initialization of the flight object:
    progress_to_file.

If TRUE, a progress .json file will be created in the current folder and will be constantly updated.
If FALSE, progress information about the simulation will be displayed on the terminal or command line.

By default, all the three optional parameters are set to FALSE.


EXAMPLES
--------------

A typical usage of the Simulator would be as follows (see the Weather Module documentation to setup the
environment object):

    import simulator
    my_flight_simulation = simulator.flight()

    my_flight_simulation.environment = my_environment_object
    my_flight_simulation.balloonGasType = 'Helium'
    my_flight_simulation.balloonModel = 'TA800'
    my_flight_simulation.nozzleLift = 1.1
    my_flight_simulation.payloadTrainWeight = 0.43
    my_flight_simulation.parachuteModel = 'SPH36'
    my_flight_simulation.trainEquivSphereDiam = 0.1
    my_flight_simulation.numberOfSimRuns = 10
    my_flight_simulation.floatingFlight = True
    my_flight_simulation.floatingAltitude = 28000
    my_flight_simulation.ventingStart = 500
    my_flight_simulation.outputFile = 'output.csv'

    my_flight_simulation.run()


An advanced usage of the Simulator could be as follows:

    import simulator
    my_flight_simulation = simulator.flight()

    my_flight_simulation.environment = my_environment_object
    my_flight_simulation.balloonGasType = 'Helium'
    my_flight_simulation.balloonModel = 'TA800'
    my_flight_simulation.nozzleLift = 1.1
    my_flight_simulation.payloadTrainWeight = 0.43
    my_flight_simulation.parachuteModel = 'SPH36'
    my_flight_simulation.trainEquivSphereDiam = 0.1
    my_flight_simulation.numberOfSimRuns = 10
    my_flight_simulation.floatingFlight = True
    my_flight_simulation.floatingAltitude = 28000
    my_flight_simulation.ventingStart = 500
    my_flight_simulation.outputFile = 'output.csv'

    my_flight_simulation.preflight()
    for flightSimulation in range(10):
        my_flight_simulation.fly(flightSimulation)
        # extract any results or data from the simulation at this point
    my_flight_simulation.postflight()


University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""
from math import pi
from datetime import timedelta
from sys import stdout
import os
import logging
from six.moves import range
import numpy
from scipy.integrate import odeint
from .flight_tools import flight_tools
from .weather import *
from . import drag_helium
from . import available_balloons_parachutes

# Error and warning logger
logger = logging.getLogger(__name__)


class flight(object):
    def __init__(self,
                 debugging=False,
                 log_to_file=False,
                 progress_to_file=False):
        """
        Initialize all the parameters of the object and setup the debugging if
        required.
        """
        # User defined variables
        self.environment = None             # weather object
        self.balloonGasType = ''
        self.balloonModel = None            # See available_balloons
        self.nozzleLift = 0.0               # kg
        self.payloadTrainWeight = 0.0       # kg
        self.parachuteModel = None
        self.numberOfSimRuns = 0
        self.trainEquivSphereDiam = 0.0     # m
        self.floatingFlight = False
        self.floatingAltitude = 0.0         # m
        self.ventingStart = 1000            # meters below the target altitude
        self.maxFlightTime = 18000
        self.excessPressureCoeff = 1
        self.outputFile = ''
        # self.launchSiteLat = 0.0   --- These will be automatically fetched
        # self.launchSiteLon = 0.0   --- during preflight from the environment object.
        # self.launchSiteElev = 0.0  --- DO NOT enter them here, they will be ignored.

        # Output
        self.results = []

        # Simulation precision - not user defined!
        self.samplingTime = 3               # seconds

        # Private variables - do not edit!
        self._preflightCompleted = False
        self._hasRun = False
        self._usingGFS = False
        self._debugging = debugging
        self._progressToFile = progress_to_file
        self._progressFile = ''
        self._balloonWeight = 0.0
        self._meanBurstDia = 0.0
        self._lowCD = []
        self._highCD = []
        self._transition = []
        self._ReBand = []
        self._burstDiameter = []
        self._parachuteCD = []
        self._balloonReturnFraction = []
        self._gasMassAtInflation = 0.0
        self._gasMolecularMass = 0.0
        self._balloonVolumeAtInflation = 0.0
        self._balloonDiaAtInflation = 0.0
        self._gasMassAtFloat = 0.0
        self._balloonVolumeAfFloat = 0.0
        self._balloonDiaAtFloat = 0.0
        self._totalAscendingMass = 0.0
        self._totalDescendingMass = 0.0
        self._lastFlightBurstAlt = 0.0
        self.flightTools = flight_tools()
        self.cutdown = False
        self.cutdownAltitude = None

        self._totalStepsForProgress = 0

        if debugging:
            log_lev = logging.DEBUG
        else:
            log_lev = logging.WARNING

        if log_to_file:
            logging.basicConfig(filename='error.log',
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=log_lev)
        else:
            logger.setLevel(log_lev)

    @profile
    def run(self):
        """
        This method takes all the necessary steps to run the simulation and
        wraps together all the simulation-specific methods found below. This
        method is RECOMMENDED for standard users.

        Returns 0 if the simulation succeeded, or the error number if there
        were errors. Check the error.log and the simulator's documentation for
        a full description of the error(s).
        """
        # Prepare progress file: try and create the file
        self._progressFile = os.path.splitext(self.outputFile)[0] +\
            '_progress.json'
        self.updateProgress(1.0, 2)
        self._totalStepsForProgress = self.numberOfSimRuns + 1

        # _________________________________________________________________ #
        # PREFLIGHT SEQUENCE

        self.preflight()

        if not self._preflightCompleted:
            # Check whether the preflight was completed successfully
            logger.error("""There was an error while performing the preflight
                validations and calculations.""")
            logger.error('The simulation was interrupted.')
            return
        self.updateProgress(0.0, 0)

        # _________________________________________________________________ #
        # RUN THE FLIGHT SIMULATION
        for flightNumber in range(self.numberOfSimRuns):
            logger.debug('SIMULATING FLIGHT %d' % (flightNumber + 1))
            self.fly(flightNumber)
            self.updateProgress(
                float(flightNumber + 1) / self._totalStepsForProgress, 0)

        # _________________________________________________________________ #
        # POSTFLIGHT HOUSEKEEPING AND RESULT PROCESSING
        self._hasRun = True
        self.postflight()
        self.updateProgress(1.0, 0)
        return 0

    @profile
    def preflight(self):
        """
        Run a series of pre-flight checks and calculations to verify the
        consistency of the flight parameters entered and to prepare all the
        data required for simulation.

        It is MANDATORY to execute this method before the simulation (otherwise
        the simulation will throw an error)

        If successful, no errors are thrown. Enable debugging for detailed
        information.
        """

        logger.debug('Preflight sequence starting...')

        # _________________________________________________________________ #
        # Fetch launch site information from the environment data
        self.launchSiteLat = self.environment.launchSiteLat
        self.launchSiteLon = self.environment.launchSiteLon
        self.launchSiteElev = self.environment.launchSiteElev

        # _________________________________________________________________ #
        # Variable validation
        toBreak = False
        if self.environment is None:
            logger.error('Flight environment not defined!')
            toBreak = True

        if not self.environment._weatherLoaded\
           and isinstance(self.environment, soundingEnvironment):
            logger.error('Flight environment has not loaded weather data!')
            toBreak = True

        if self.balloonGasType not in ['Helium', 'Hydrogen']:
            logger.warning('{} is an invalid gas type'.format(
                self.balloonGasType))
            logger.warning("The calculation will carry on with Helium.")
            self.balloonGasType = 'Helium'

        if self.balloonModel not in available_balloons_parachutes.balloons:
            logger.error(
                '{} is an invalid balloon model'.format(self.balloonModel))
            logger.error('Supported models are %s. Please correct it.',
                         available_balloons_parachutes.balloons.keys())
            toBreak = True

        if self.parachuteModel is not None and self.parachuteModel not in\
           available_balloons_parachutes.parachutes:
            logger.error('An invalid parachute model was found.')
            logger.error('Supported models are None or %s. Please correct it.',
                         available_balloons_parachutes.parachutes.keys())
            toBreak = True

        if self.nozzleLift == 0.0:
            logger.error('Nozzle lift cannot be zero!')
            toBreak = True

        if self.payloadTrainWeight == 0.0:
            logger.error('Payload train weight cannot be zero!')
            toBreak = True

        if self.numberOfSimRuns == 0:
            logger.error('The number of sim runs cannot be zero!')
            toBreak = True

        if self.launchSiteLat == 0.0 or self.launchSiteLon == 0.0:
            logger.warning('Warning: Are you sure you set the launch site location correctly?')
        # if self.trainEquivSphereDiam == 0.0:
        #     logger.error('The train equivalent sphere diameter cannot be zero!')
        #     toBreak = True
        if self.nozzleLift <= self.payloadTrainWeight:
            logger.error(
                'The nozzle lift is too low for the balloon to climb! Adjust the nozzle lift before continuing')
            toBreak = True

        # If the forecast hasn't been downloaded yet, do it now.
        if not self.environment._weatherLoaded and isinstance(self.environment, forecastEnvironment):
            self.environment.maxFlightTime = self.maxFlightTime
            self.environment.loadForecast(self.updateProgress)
            if not self.environment._weatherLoaded:
                logger.error('Cannot load forecast!')
                toBreak = True

        # Check if the output file can be written and create it. If not, stop the preflight.
        try:
            logger.debug("Creating output file {}".format(self.outputFile))
            out = open(self.outputFile, 'w')
            out.close()
            os.remove(self.outputFile)
        except IOError:
            if not os.path.isdir(self.outputFile):
                logger.error('The output file cannot be created.\n')
                toBreak = True

        if toBreak:
            logger.error('The simulation was interrupted.')
            return

        logger.debug('Flight parameters validation succeeded.')
        logger.debug('Beginning Preflight...')

        # ____________________________________________________________ #
        # Preflight calculations

        # Balloon performance estimation
        # According to the balloon weight entered, select the related mean
        # burst diameter and its Weibull coefficients.
        self._balloonWeight = \
            available_balloons_parachutes.balloons[self.balloonModel][0]
        self._meanBurstDia = (
            available_balloons_parachutes.meanToNominalBurstRatio *
            available_balloons_parachutes.balloons[self.balloonModel][1]
        )
        self._weibull_lambda = \
            available_balloons_parachutes.balloons[self.balloonModel][2]
        self._weibull_k = \
            available_balloons_parachutes.balloons[self.balloonModel][3]

        logger.debug("""Balloon performance: Mean burst diameter: %.4f,
            Lambda: %.4f, k: %4f""" % (
            self._meanBurstDia, self._weibull_lambda, self._weibull_k))

        # ____________________________________________________________ #
        # Variable initialization

        # Load the flight tools and prep them for simulation

        self.results = []
        self._lowCD = []
        self._highCD = []
        self._transition = []
        self._ReBand = []
        self._burstDiameter = []
        self._parachuteCD = []
        self._balloonReturnFraction = []

        # Check if GFS is being used
        if isinstance(self.environment, forecastEnvironment):
            self._usingGFS = True
        else:
            self._usingGFS = False

        # Initialize variables subject to perturbation (Monte Carlo simulation)
        if self.numberOfSimRuns == 1:
            # Deterministic simulation: use standard values, no perturbations
            self._lowCD.append(0.225)
            self._highCD.append(0.425)
            self._transition.append(3.296)
            self._ReBand.append(0.363)
            self._burstDiameter.append(self._meanBurstDia)
            self._parachuteCD.append(0.9)
            self._balloonReturnFraction.append(0.03)
        else:
            # Monte Carlo simulation: perturb values
            for _ in range(self.numberOfSimRuns):
                mcIndex = numpy.random.random_integers(0, (
                    numpy.size(drag_helium.transitions[:, 0])) - 1)
                self._lowCD.append(drag_helium.transitions[mcIndex, 0])
                self._highCD.append(drag_helium.transitions[mcIndex, 1])
                self._transition.append(drag_helium.transitions[mcIndex, 2])
                self._ReBand.append(drag_helium.transitions[mcIndex, 3])
                self._balloonReturnFraction.append(
                    0.03 + numpy.random.random() * (1 - 0.03))
                self._parachuteCD.append(0.9 + 0.2 * numpy.random.random() *
                                         (-1) ** round(numpy.random.random()))
                self._burstDiameter.append(
                    self._weibull_lambda *
                    numpy.random.weibull(self._weibull_k))
                # Perturb the wind for Monte Carlo.
            self.environment.perturbWind(self.numberOfSimRuns)

        # _________________________________________________________________ #
        # Lifting gas mass calculations
        if self.balloonGasType == 'Helium':
            self._gasMolecularMass = (
                0.945 * self.flightTools.HeMolecularMass +
                0.055 * self.flightTools.airMolecularMass
            )
        else:
            self._gasMolecularMass = (
                0.985 * self.flightTools.HydrogenMolecularMass +
                0.015 * self.flightTools.airMolecularMass
            )

        # Use flight tools to calculate all preliminary balloon calculations
        # (gas mass, balloon volume and diameter at inflation)
        (self._gasMassAtInflation, self._balloonVolumeAtInflation,
            self._balloonDiaAtInflation) = self.flightTools.liftingGasMass(
                self.nozzleLift,
                self._balloonWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                             self.launchSiteLon,
                                             self.launchSiteElev,
                                             self.environment.dateAndTime),
                self._gasMolecularMass,
                self.excessPressureCoeff
        )

        logger.debug('Lifting gas mass calcs completed!')
        logger.debug('Gas molecular mass: %.4f' % self._gasMolecularMass)
        logger.debug('Gas mass: %.4f' % self._gasMassAtInflation)
        logger.debug('Balloon volume at inflation: %.4f' %
                     self._balloonVolumeAtInflation)
        logger.debug('Balloon diameter ad inflation: %.4f' %
                     self._balloonDiaAtInflation)

        # In case it's a floating flight, calculate the gas mass and the
        # balloon volume needed at the target altitude for the balloon to
        # float.

        # WARNING: There is an approximation here, in case GFS is being used:
        # this data is calculated as if the balloon were to float at the launch
        # site, rather than at its true location. It can probably be fixed by
        # continuously calculating these during flight simulation, but this
        # would slow the whole simulation down.

        if self.floatingFlight:
            # Calculate the balloon characteristics at the floating altitude
            (self._gasMassAtFloat, self._balloonVolumeAtFloat,
                self._balloonDiaAtFloat) = self.flightTools.liftingGasMass(
                    # This is the Nozzle Lift, which has to be equal to the
                    # payload train weight for the balloon to float. If they
                    # are, the sum of the (vertical) forces is 0.
                    self.payloadTrainWeight,
                    self._balloonWeight,
                    self.environment.getTemperature(
                        self.launchSiteLat,
                        self.launchSiteLon,
                        self.floatingAltitude,
                        self.environment.dateAndTime),
                    self.environment.getPressure(self.launchSiteLat,
                                                 self.launchSiteLon,
                                                 self.floatingAltitude,
                                                 self.environment.dateAndTime),
                    self._gasMolecularMass,
                    self.excessPressureCoeff
            )

            logger.debug('Lifting gas mass calcs for floating:')
            logger.debug('Gas molecular mass: %.4f' % self._gasMolecularMass)
            logger.debug('Gas mass: %.4f' % self._gasMassAtFloat)
            logger.debug('Balloon volume at floating alt: %.4f' %
                         self._balloonVolumeAtFloat)
            logger.debug('Balloon diameter ad floating alt: %.4f' %
                         self._balloonDiaAtFloat)

        # Configure the flight tools for the correct parachute
        if self.parachuteModel is not None:
            self.flightTools.parachuteAref = \
                available_balloons_parachutes.parachutes[self.parachuteModel]
        else:
            self.flightTools.parachuteAref = 0.0

        self._preflightCompleted = True

        logger.debug('Preflight completed!')

    @profile
    def fly(self, flightNumber):
        """
        Execute a single simulation.
        It should be run N times, where N is the number of simulation runs
        specified upon configuration in the numberOfSimRuns variable.
        flightNumber should have values 0...N-1

        If successful, no errors are thrown. Enable debugging for detailed
        information.
        """

        # Check whether the preflight sequence was performed. If not, stop the
        # simulation.
        if not self._preflightCompleted:
            logger.error("""Preflight sequence needs to be performed before the
                actual flight! Simulation interrupted.""")
            return

        # Flight-specific variables initialization
        # Prepare the flight tools to deliver results specific to this flight.
        self.flightTools.highCD = self._highCD[flightNumber]
        self.flightTools.lowCD = self._lowCD[flightNumber]
        self.flightTools.ReBand = self._ReBand[flightNumber]
        self.flightTools.transition = self._transition[flightNumber]
        self.flightTools.parachuteCD = self._parachuteCD[flightNumber]
        self._totalAscendingMass = (self.payloadTrainWeight +
                                    self._gasMassAtInflation +
                                    self._balloonWeight)
        self._totalDescendingMass = (self.payloadTrainWeight +
            self._balloonWeight * self._balloonReturnFraction[flightNumber])
        # Use the correct wind profile according to the simulation type: a
        # standard one for deterministic simulations or a perturbed one for
        # Monte Carlo simulations.
        if self.numberOfSimRuns == 1:
            currentFlightWindDirection = self.environment.getWindDirection
            currentFlightWindSpeed = self.environment.getWindSpeed
        else:
            currentFlightWindDirection =\
                self.environment.getMCWindDirection[flightNumber]
            currentFlightWindSpeed = \
                self.environment.getMCWindSpeed[flightNumber]

        logger.debug('Flight variables initialized.')
        logger.debug('Low CD: %.4f' % self.flightTools.lowCD)
        logger.debug('High CD: %.4f' % self.flightTools.highCD)
        logger.debug('Transition: %.4f' % self.flightTools.transition)
        logger.debug('Re Band: %.4f' % self.flightTools.ReBand)
        logger.debug('Burst Diameter: %.4f' %
                     self._burstDiameter[flightNumber])
        logger.debug('Balloon Return Fraction: %.4f' %
                     self._balloonReturnFraction[flightNumber])
        logger.debug('Parachute CD: %.4f' % self.flightTools.parachuteCD)

        logger.debug('Total ascending mass: %.4f' % self._totalAscendingMass)
        logger.debug('Total descending mass: %.4f' % self._totalDescendingMass)

        logger.debug('Beginning simulation of flight %d' % (flightNumber + 1))

        self._lastFlightBurst = False
        self._lastFlightBurstAlt = 0.0
        self._lastBurstIndex = 0
        self._currentLatPosition = self.launchSiteLat
        self._currentLonPosition = self.launchSiteLon
        self._currentTime = 0
        self._totalRuns = 0

        print(self.cutdown, self.cutdownAltitude)

        def ode(y, t):
            """
            This is the right-hand side of the 2nd order ODE defining the
            vertical motion of the balloon.
            """
            self._totalRuns += 1
            # Extract inputs from array y
            altitude = y[0]
            ascentRate = y[1]

            if self.cutdown and not self._lastFlightBurst:
                if altitude > self.cutdownAltitude:
                    # Burst the balloon
                    logger.debug('Bursting the balloon at {}m altitude'.format(
                        altitude))
                    self._lastFlightBurst = True
                    self._lastFlightBurstAlt = altitude
                    return numpy.array([0.0, 0.0])

            # Calculate current position in time and space
            # This is used to fetch the correct atmospheric data if the GFS is
            # being used.
            currentTime = self.environment.dateAndTime + timedelta(seconds=t)

            # Convert wind direction and speed to u- and v-coordinates
            windLon, windLat = tools.dirspeed2uv(
                currentFlightWindDirection(self._currentLatPosition,
                                           self._currentLonPosition,
                                           altitude,
                                           currentTime),
                currentFlightWindSpeed(self._currentLatPosition,
                                       self._currentLonPosition,
                                       altitude,
                                       currentTime) * 0.514444
            )

            # Calculate how much the drift has been from the previous to the
            # current iteration, convert it to degrees and add it up to the
            # previous lat,lon position to find the current one.
            dLat, dLon = tools.m2deg(windLat * (t - self._currentTime),
                                     windLon * (t - self._currentTime),
                                     self._currentLatPosition)
            self._currentLatPosition += dLat
            self._currentLonPosition += dLon
            # Check if within bounds, correct otherwise
            if self._currentLatPosition > 90:
                self._currentLatPosition = 180 - self._currentLatPosition
                self._currentLonPosition += 180
            elif self._currentLatPosition < -90:
                self._currentLatPosition = -180 - self._currentLatPosition
                self._currentLonPosition += 180
            if self._currentLonPosition > 180:
                self._currentLonPosition -= 360
            elif self._currentLonPosition <= -180:
                self._currentLonPosition += 360

            # Update the current time, to be used in the next iteration to work
            # out the time difference.
            self._currentTime = t

            if not self._lastFlightBurst:
                # THE BALLOON IS CURRENTLY ASCENDING

                # If floating flight, check whether the balloon is approaching
                # the target altitude to start venting gas out. The altitude at
                # which venting starts is floatingAltitude - ventingStart.
                if self.floatingFlight:
                    gasMass = self.flightTools.gasMassForFloat(
                        altitude,
                        self.floatingAltitude,
                        self._gasMassAtInflation,
                        self._gasMassAtFloat,
                        ventStart=self.ventingStart
                    )
                else:
                    gasMass = self._gasMassAtInflation


                # Calculate current balloon diameter to check for burst
                gasDensity = self.excessPressureCoeff * self.environment.getPressure(self._currentLatPosition,
                                                                                     self._currentLonPosition, altitude,
                                                                                     currentTime) * 100 * self._gasMolecularMass / (
                                 8.31447 * tools.c2kel(
                                     self.environment.getTemperature(self._currentLatPosition, self._currentLonPosition,
                                                                     altitude, currentTime)))
                balloonVolume = gasMass / gasDensity
                balloonDiameter = (6 * balloonVolume / pi) ** (1. / 3)

                # If floating flight, calculate the nozzle lift if the gas is being vented.
                if self.floatingFlight:
                    nozzleLift = self.flightTools.nozzleLiftForFloat(
                        self.nozzleLift,
                        self.environment.getDensity(self._currentLatPosition, self._currentLonPosition, altitude,
                                                    currentTime),
                        gasDensity,
                        balloonVolume,
                        self._balloonWeight,
                        altitude,
                        self.floatingAltitude,
                        ventStart=self.ventingStart
                    )
                else:
                    nozzleLift = self.nozzleLift

                if balloonDiameter < self._burstDiameter[flightNumber]:
                    # THE BALLOON DIDN'T BURST

                    # Calculate balloon and train drag
                    currentDensity = self.environment.getDensity(self._currentLatPosition, self._currentLonPosition,
                                                                 altitude, currentTime)
                    currentViscosity = self.environment.getViscosity(self._currentLatPosition, self._currentLonPosition,
                                                                     altitude, currentTime)
                    balloonDrag = self.flightTools.balloonDrag(balloonDiameter, ascentRate, currentDensity,
                                                               currentViscosity)
                    trainDrag = self.flightTools.balloonDrag(self.trainEquivSphereDiam, ascentRate, currentDensity,
                                                             currentViscosity)

                    # External Forces
                    externalForces = (nozzleLift - self.payloadTrainWeight) * 9.81 - balloonDrag - trainDrag

                    # Derivatives
                    dvdt = externalForces / self._totalAscendingMass
                    dhdt = ascentRate

                    return numpy.array([dhdt, dvdt])

                else:
                    # THE BALLOON HAS BURST

                    # Floating flight is set to false because if the balloon has burst, the flight is now standard.
                    self._lastFlightBurst = True
                    self._lastFlightBurstAlt = altitude
                    return numpy.array([0.0, 0.0])

            else:
                # THE BALLOON IS CURRENTLY DESCENDING


                # Calculate parachute Drag
                currentDensity = self.environment.getDensity(self._currentLatPosition, self._currentLonPosition,
                                                             altitude, currentTime)
                currentViscosity = self.environment.getViscosity(self._currentLatPosition, self._currentLonPosition,
                                                                 altitude, currentTime)
                if self.parachuteModel == None:
                    parachuteDrag = 0
                else:
                    parachuteDrag = self.flightTools.parachuteDrag(abs(ascentRate), currentDensity)

                # Train Drag
                trainDrag = abs(self.flightTools.balloonDrag(self.trainEquivSphereDiam, abs(ascentRate), currentDensity,
                                                             currentViscosity))

                # External Forces
                externalForces = -self.payloadTrainWeight * 9.81 * (
                    1 + self._balloonReturnFraction[flightNumber]) + parachuteDrag + trainDrag

                # Derivatives
                dvdt = externalForces / self._totalDescendingMass
                dhdt = ascentRate

                return numpy.array([dhdt, dvdt])

        # Define the initial conditions, the time vector at which we want simulation data to be stored, and run the
        # integration.
        # Note: the simulation carries on all the way to the maxFlightTime, even if the altitude becomes negative.
        # Negative values of altitude will be trimmed later on.
        initialConditions = numpy.array([self.launchSiteElev, 0.0])
        timeVector = numpy.arange(0, self.maxFlightTime + self.samplingTime, self.samplingTime)

        logger.debug('Beginning integration.')

        ### INTEGRATION ###
        if self._usingGFS:
            solution = odeint(ode, initialConditions, timeVector, rtol=1e-3, atol=1e-3)
        else:
            solution = odeint(ode, initialConditions, timeVector)
        ###################

        logger.debug('Integration completed. Post-processing...')

        # Extract altitude and ascent rate data from the solution array
        solution_altitude = numpy.array(solution[:, 0])
        #solution_ascrate = numpy.array(solution[:,1]) ### Currently not used.

        # Trim negative altitude values from results and then trim time to the same length.
        if self._lastFlightBurst:
            solution_altitude = solution_altitude[solution_altitude > 0]
            solution_altitude[-1] = 0.0

        #solution_ascrate = solution_ascrate[:len(solution_altitude)] ### Currently not used.
        timeVector = timeVector[:len(solution_altitude)]

        # Calculate drift
        lastDriftLat = 0.0
        lastDriftLon = 0.0
        latitudeProfile = [self.launchSiteLat]
        longitudeProfile = [self.launchSiteLon]
        for eachAlt, eachTime in zip(solution_altitude[1:], timeVector[1:]):
            # Calculate the position of the balloon at each point and use it to work out its location
            currentTime = self.environment.dateAndTime + timedelta(seconds=float(eachTime))
            # Gather wind speed and convert to m/s
            windSpeed = currentFlightWindSpeed(latitudeProfile[-1], longitudeProfile[-1], eachAlt,
                                               currentTime) * 0.514444
            # Convert the wind to u- and v-coordinates
            windLon, windLat = tools.dirspeed2uv(
                currentFlightWindDirection(latitudeProfile[-1], longitudeProfile[-1], eachAlt, currentTime), windSpeed)
            # Store the drift in meters (this is the distance between the LAUNCH SITE and the current location)
            lastDriftLat += windLat * self.samplingTime
            lastDriftLon += windLon * self.samplingTime
            # Convert it to degrees
            dLat, dLon = tools.m2deg(lastDriftLat, lastDriftLon, latitudeProfile[-1])
            # Store the new latitude and longitude
            latitudeProfile.append(self.launchSiteLat + dLat)
            longitudeProfile.append(self.launchSiteLon + dLon)


        # Check that latitude and longitude are within bounds and correct if they are not (for example, if the balloon
        # flew over the North or the South Pole).
        for i in range(len(latitudeProfile)):
            if latitudeProfile[i] > 90:
                latitudeProfile[i] = 180 - latitudeProfile[i]
                longitudeProfile[i] += 180
            elif latitudeProfile[i] < -90:
                latitudeProfile[i] = -180 - latitudeProfile[i]
                longitudeProfile[i] += 180

            if longitudeProfile[i] > 180:
                n = int(longitudeProfile[i]-180)/360 + 1
                longitudeProfile[i] -= n * 360
            elif longitudeProfile[i] <= -180:
                n = int(abs(longitudeProfile[i])-180)/360 + 1
                longitudeProfile[i] += n * 360


        # Find burst point index or "target altitude reached" index
        if self._lastFlightBurst:
            # Store burst point index
            index = tools.find_nearest_index(solution_altitude, self._lastFlightBurstAlt)
        else:
            if solution_altitude[-1] < self.floatingAltitude - 100:
                # In this case, the balloon hasn't reached the target altitude. This is probably because the
                # maxFlightTime is too low. Show an error.
                index = -1
            else:
                # Store target altitude reached index
                index = tools.find_nearest_index(solution_altitude, self.floatingAltitude)

        # STORE RESULTS OF CURRENT SIMULATION
        if self._lastFlightBurst:
            # The results are:   flight number  time vector latitude profile longitude profile   altitude    burst index   burst altitude   has burst
            self.results.append(
                [flightNumber + 1, timeVector, latitudeProfile, longitudeProfile, solution_altitude, index,
                 self._lastFlightBurstAlt, True])
        else:
            # The results are:   flight number  time vector latitude profile longitude profile   altitude   burst index  target altitude  has burst
            self.results.append(
                [flightNumber + 1, timeVector, latitudeProfile, longitudeProfile, solution_altitude, index,
                 self.floatingAltitude, False])

        logger.debug('Simulation completed.')

    @profile
    def postflight(self):
        """
        After all the simulations have been executed, this method puts the results together, processes them and stores
        them in the outputFormat required.
        The format of the file to be generated depends on the extension given in the outputFile.
        Accepted extensions are:
            'json'  JavaScript data structure, used to pass data to the web interface;
            'kml'   Standard format for geographic data. It can be opened by Google Maps, Google Earth, etc;
            'kmz'   Zipped kml, good for online use (eg Google Maps). The file size is significantly reduced;
            'csv'   Comma separated values, the preferred format for further data processing with any other software;
            'web'   This stores json, kml and csv in the same folder requested. This is used by the web interface to
                    prepare files for export.
            none    If no extension is found, a new folder is created and ALL output formats are stored.

        The results are stored in the outputFile specified upon configuration.
        """

        def store_data(outputPath, data_format):
            """
            Function responsible for storing the data in the given format.
            """

            logger.debug('Storing data... Requested format: %s' % data_format)

            if data_format in ('json', 'JSON'):
                # GENERATE JSON FILE OUT OF RESULTS

                # Every how many points should we store one (this is used to reduce the size of the json file)
                shrinkFactor = 5


                # Initialize the data lists
                jsonBurstMarkers = []
                jsonLandingMarkers = []
                jsonFloatMarkers = []
                jsonPaths = []

                # Go through the results of each flight simulation
                for flightResult in self.results:
                    thisFlightHasBurst = flightResult[7]

                    # Build up the flight path data
                    jsonPath = ['{ "points" : [\n']
                    jsonPoints = []
                    pointNumber = 0
                    for eachPoint in numpy.transpose(flightResult[2:5]):
                        if pointNumber % shrinkFactor == 0:
                            jsonPoints.append('{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f }\n' % (
                                numpy.clip(eachPoint[0], -85.0511, 85.0511), eachPoint[1], eachPoint[2]))
                        pointNumber += 1

                    jsonPath.append(','.join(jsonPoints))
                    jsonPath.append('\n] }')
                    jsonPaths.append(''.join(jsonPath))

                    if thisFlightHasBurst:
                        # BALLOON HAS BURST DURING THIS SIMULATION
                        # Build up the landing markers data
                        flightHrs, flightMins, flightSecs = tools.prettySeconds(flightResult[1][-1])
                        landTime = self.environment.dateAndTime + timedelta(seconds=float(flightResult[1][-1]))
                        jsonLandingMarkers.append(
                            '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Payload Landing Site", "simNumber" : "%d", "otherData" : "Flight time: %d hours, %d mins and %d secs <br /> ETA: %s" }\n' % (
                                numpy.transpose(flightResult[2:5])[-1][0], numpy.transpose(flightResult[2:5])[-1][1],
                                numpy.transpose(flightResult[2:5])[-1][1], flightResult[0], flightHrs, flightMins,
                                flightSecs, landTime.strftime("%d/%m/%Y %H:%M:%S")))

                        # Build up the burst markers data
                        jsonBurstMarkers.append(
                            '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Balloon Burst", "simNumber" : "%d", "otherData" : "Burst altitude: %.0f m"  }\n' % (
                                numpy.transpose(flightResult[2:5])[flightResult[5]][0],
                                numpy.transpose(flightResult[2:5])[flightResult[5]][1], flightResult[6],
                                flightResult[0],
                                flightResult[6]))
                    else:
                        # BALLOON HASN'T BURST DURING THIS SIMULATION

                        # Build up the float markers data
                        if flightResult[5] == -1:
                            # This is the case when the target altitude was not reached. Show an error to the user.
                            jsonFloatMarkers.append(
                                '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Target Altitude NOT Reached", "simNumber" : "%d", "otherData" : "Max altitude: %.0f m <br />Try to increase the flight duration." }\n' % (
                                    numpy.transpose(flightResult[2:5])[flightResult[5]][0],
                                    numpy.transpose(flightResult[2:5])[flightResult[5]][1], flightResult[6],
                                    flightResult[0], numpy.transpose(flightResult[2:5])[flightResult[5]][2]))
                        else:
                            jsonFloatMarkers.append(
                                '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Target Altitude Reached", "simNumber" : "%d", "otherData" : "Target altitude: %.0f m" }\n' % (
                                    numpy.transpose(flightResult[2:5])[flightResult[5]][0],
                                    numpy.transpose(flightResult[2:5])[flightResult[5]][1], flightResult[6],
                                    flightResult[0], flightResult[6]))

                # Put all the non-empty lists above together
                if len(jsonBurstMarkers) != 0:
                    jsonBurstMarkers = ['"burstMarkers" : [\n', ','.join(jsonBurstMarkers), '],']

                if len(jsonLandingMarkers) != 0:
                    jsonLandingMarkers = ['"landingMarkers" : [\n', ','.join(jsonLandingMarkers), '],']

                if len(jsonFloatMarkers) != 0:
                    jsonFloatMarkers = ['"floatMarkers" : [\n', ','.join(jsonFloatMarkers), '],']

                jsonPaths = ['"flightPaths" : [\n', ','.join(jsonPaths), ']']

                # Put them all together in one single array
                jsonToAdd = [
                    '{',
                    ''.join(jsonBurstMarkers),
                    ''.join(jsonFloatMarkers),
                    ''.join(jsonLandingMarkers),
                    ''.join(jsonPaths),
                    '}'
                ]

                # Try to gain write permission on the output file specified upon configuration
                try:
                    jsonFile = open(outputPath, 'w')
                except IOError:
                    logger.error('Cannot create output file.')
                    return

                # Write all the data to the file
                jsonFile.write(''.join(jsonToAdd))
                jsonFile.close()


            elif data_format in ('kml', 'kmz', 'KML', 'KMZ'):
                # GENERATE KML OUT OF RESULTS

                # Author of the KML file
                kmlAuthor = 'ASTRA High Altitude Balloon Flight Planner, University of Southampton'
                launchPinURL = 'http://maps.google.com/mapfiles/ms/micons/red-dot.png'
                burstPinURL = 'http://maps.google.com/mapfiles/ms/micons/yellow-dot.png'
                landingPinURL = 'http://maps.google.com/mapfiles/ms/micons/red-dot.png'

                kmlPaths = []
                kmlMarkers = []

                for flightResult in self.results:
                    flightHrs, flightMins, flightSecs = tools.prettySeconds(flightResult[1][-1])
                    thisFlightHasBurst = flightResult[7]

                    kmlPath = [
                        '<Placemark>\n<name>Simulation %d</name>\n<styleUrl>#stratoLine</styleUrl>\n<LineString>\n<coordinates>\n' % (
                            flightResult[0])]

                    pointNumber = 0
                    for eachPoint in numpy.transpose(flightResult[2:5]):
                        if thisFlightHasBurst:
                            if pointNumber % 10 == 0:
                                kmlPath.append('%.5f,%.5f,%.5f\n' % (eachPoint[1], eachPoint[0], eachPoint[2]))
                        else:
                            kmlPath.append('%.5f,%.5f,%.5f\n' % (eachPoint[1], eachPoint[0], eachPoint[2]))
                        pointNumber += 1

                    kmlPath.append(
                        '</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</LineString>\n</Placemark>\n')

                    kmlPaths.append(''.join(kmlPath))

                    # Add balloon launch point
                    kmlMarkers.append(
                        '<Placemark>\n<name>Balloon Launch</name>\n<styleUrl>#launchPin</styleUrl>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                            numpy.transpose(flightResult[2:5])[0][1], numpy.transpose(flightResult[2:5])[0][0],
                            numpy.transpose(flightResult[2:5])[0][2]))

                    if thisFlightHasBurst:
                        # Add balloon burst point
                        kmlMarkers.append(
                            '<Placemark>\n<name>Balloon burst (Simulation %d)</name>\n<styleUrl>#burstPin</styleUrl>\n<description>Burst altitude: %.0f m.' % (
                                flightResult[0], flightResult[6]))
                        kmlMarkers.append(
                            '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                                numpy.transpose(flightResult[2:5])[flightResult[5]][1],
                                numpy.transpose(flightResult[2:5])[flightResult[5]][0], flightResult[6]))

                        # Add balloon landing point
                        kmlMarkers.append(
                            '<Placemark>\n<name>Payload landing (Simulation %d)</name>\n<styleUrl>#landingPin</styleUrl>\n<description>Flight time: %d hours, %d minutes and %d seconds.' % (
                                flightResult[0], flightHrs, flightMins, flightSecs))
                        kmlMarkers.append(
                            '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                                numpy.transpose(flightResult[2:5])[-1][1], numpy.transpose(flightResult[2:5])[-1][0],
                                numpy.transpose(flightResult[2:5])[-1][2]))
                    else:
                        # Add target altitude reached point
                        kmlMarkers.append(
                            '<Placemark>\n<name>Target altitude reached (Simulation %d)</name>\n<styleUrl>#burstPin</styleUrl>\n<description>Target altitude: %.0f m.' % (
                                flightResult[0], flightResult[6]))
                        kmlMarkers.append(
                            '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                                numpy.transpose(flightResult[2:5])[flightResult[5]][1],
                                numpy.transpose(flightResult[2:5])[flightResult[5]][0], flightResult[6]))

                kmlToAdd = [
                    '<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n<Document>\n<name>Flight Simulations</name>\n<atom:author>\n<atom:name>%s</atom:name>\n</atom:author>\n' % kmlAuthor,
                    '<Style id="landingPin">\n<IconStyle>\n<scale>1.0</scale>\n<Icon>\n<href>%s</href>\n</Icon>\n</IconStyle>\n</Style>\n' % landingPinURL,
                    '<Style id="launchPin">\n<IconStyle>\n<scale>1.0</scale>\n<Icon>\n<href>%s</href>\n</Icon>\n</IconStyle>\n</Style>\n' % launchPinURL,
                    '<Style id="burstPin">\n<IconStyle>\n<scale>1.0</scale>\n<Icon>\n<href>%s</href>\n</Icon>\n</IconStyle>\n</Style>\n' % burstPinURL,
                    '<Style id="stratoLine">\n<LineStyle>\n<width>1.0</width>\n</LineStyle>\n</Style>\n',
                    ''.join(kmlPaths),
                    ''.join(kmlMarkers),
                    '</Document>\n</kml>'
                ]

                if data_format in ('kmz', 'KMZ'):
                    import zipfile, tempfile

                    try:
                        outputKml = tempfile.NamedTemporaryFile(mode='w')
                    except IOError:
                        logger.error('Error: cannot create a new temporary file')
                        return

                    outputKml.write(''.join(kmlToAdd))
                    outputKml.flush()


                    # Convert KML to KMZ by zipping it.
                    try:
                        import zlib

                        zipCompression = zipfile.ZIP_DEFLATED
                    except:
                        zipCompression = zipfile.ZIP_STORED

                    with zipfile.ZipFile(outputPath, 'w') as zipKmz:
                        zipKmz.write(outputKml.name, arcname='output.kml', compress_type=zipCompression)
                    zipKmz.close()
                    outputKml.close()

                    logger.debug(('KMZ file generated! ', self.outputFile))

                else:
                    # Write the KML file.

                    try:
                        outputKml = open(outputPath, 'w')
                    except IOError:
                        logger.error('Error: cannot create the output file in the given directory')
                        return

                    outputKml.write(''.join(kmlToAdd))
                    outputKml.close()


            elif data_format in ('csv', 'CSV', 'csv.zip', 'CSV.ZIP'):
                # GENERATE CSV FILE OUT OF RESULTS

                # Calculate how many rows we need in total
                totalRows = 0
                for flightResult in self.results:
                    totalRows += len(flightResult[1])

                # Columns are Flight #, Time from launch, Lat, Lon, Alt, Remarks
                totalColumns = 6

                # Add one column for the header and generate an empty matrix
                csvMatrix = numpy.zeros((totalRows + 2) * totalColumns).reshape(totalRows + 2, totalColumns).astype(
                    '|S100')
                # Clear the remarks column (it was full of zeros)
                csvMatrix[..., 5] = ''

                currentPosition = 2
                # Populate matrix with data
                for flightResult in self.results:
                    numberOfPoints = len(flightResult[1])
                    thisFlightHasBurst = flightResult[7]

                    # Flight number, Time, Lat, Lon, Alt
                    for i in range(5):
                        csvMatrix[currentPosition:currentPosition + numberOfPoints, i] = flightResult[i]

                    # Remarks: Launch
                    csvMatrix[currentPosition, 5] = 'Balloon Launch'

                    if thisFlightHasBurst:
                        # Remarks: Burst
                        csvMatrix[currentPosition + flightResult[5], 5] = 'Balloon Burst'

                        # Remarks: Landing
                        csvMatrix[currentPosition + numberOfPoints - 1, 5] = 'Balloon Landed'
                    else:
                        # Remarks: Target Altitude Reached
                        csvMatrix[currentPosition + flightResult[5], 5] = 'Balloon has reached target altitude'

                    currentPosition += numberOfPoints

                # Add header
                csvMatrix[0] = [
                    'Data generated by the ASTRA High Altitude Balloon Flight Planner - University of Southampton', '',
                    '', '', '', '']

                csvMatrix[1] = ['Simulation #', 'Time from launch [s]', 'Latitude', 'Longitude', 'Altitude [m]',
                                'Remarks']

                # Save file
                if data_format in ('csv', 'CSV'):

                    numpy.savetxt(outputPath, csvMatrix, delimiter=',', fmt='%s')

                else:
                    import zipfile, tempfile

                    try:
                        outputCsv = tempfile.NamedTemporaryFile()
                    except IOError:
                        logger.error('Error: cannot create a new temporary file')
                        return

                    # Create a zipped version of the file, since it's probably big (approx 150KB/simulation not compressed)
                    numpy.savetxt(outputCsv, csvMatrix, delimiter=',', fmt='%s')
                    outputCsv.flush()

                    # Zip CSV file
                    try:
                        import zlib

                        zipCompression = zipfile.ZIP_DEFLATED
                    except:
                        zipCompression = zipfile.ZIP_STORED

                    with zipfile.ZipFile(outputPath, 'w') as zipCsv:
                        zipCsv.write(outputCsv.name, arcname='ASTRA Simulation Results.csv', compress_type=zipCompression)
                    zipCsv.close()
                    outputCsv.close()

                    logger.debug(('CSV-ZIP file generated! ', self.outputFile))



            else:
                logger.error('Output data format not understood! Is the specified extension correct?')

            # Done
            logger.debug('Output file generated.')

        ######################################################
        # Figure out which formats to output and export them #

        # Extract the output format from the outputFile path.
        fileName, fileExtension = os.path.splitext(self.outputFile)

        # Remove the dot from the extension
        fileExtension = fileExtension[1:]

        # Check a folder path has not been given and, if it has, remove slash
        if fileName[-1] == '/':
            fileName = fileName[:-1]

        if fileExtension == '':
            # If file has no extension, store json, kml, kmz, and csv

            try:
                os.mkdir(self.outputFile)
            except OSError:

                if not os.path.isdir(self.outputFile):
                    logger.error('The specified output path already exists. Change it or add an extension.')

            for data_format in ['json', 'kml', 'kmz', 'csv', 'csv.zip']:
                path = fileName + '/' + fileName.split('/')[-1] + '.' + data_format
                store_data(path, data_format)

        elif fileExtension == 'web':
            # If file has extension web, store json, kml and csv.gz

            for data_format in ['json', 'kml', 'csv.zip']:
                path = fileName + '.' + data_format
                store_data(path, data_format)

        else:
            # Only store the required one

            store_data(self.outputFile, fileExtension)


    def reset(self, keepParameters=False):
        """
        Reset the simulation.

        If keepParameters is FALSE, all the parameters are deleted and the flight object is ready to be reconfigured.
        If keepParameters is TRUE, all the results are deleted but the initial flight configuration is preserved.
        preflight() should still be run again before performing anymore simulations.
        Default is FALSE.
        """

        if keepParameters:
            self.results = []
            self._lastFlightBurst = False
            self._lastBurstIndex = 0
            self._lastFlightBurstAlt = 0.0
            self._hasRun = False
            self._preflightCompleted = False
        else:
            self.__init__()

    @profile
    def updateProgress(self, value, action):
        """
        Update the progress file with the ratio of value and the total steps of the simulation calculated when
        the method run() is executed.

        Note: if run() is not being used, the object's parameter _totalStepsForProgress should be defined before
        executing this method for it to work properly.
        """

        if self._progressToFile or self._debugging:
            # Try to gain write permissions to the progress file
            try:
                progFile = open(self._progressFile, 'w')
            except IOError:
                logger.error('Cannot create output file.')
                return

            # Write to file
            progFile.write('{ "progress" : %f, "action" : %d }' % (value, action))
            progFile.close()
        else:
            if action == 0:
                if value < 1:
                    stdout.write('\rThe simulation is running. Current progress: %d%%' % int(value * 100))
                else:
                    stdout.write('\rThe simulation is running. Current progress: 100%')
                    stdout.write('\nSimulation completed.\n')
            elif action == 1:
                if value < 1:
                    stdout.write('\rDownloading weather forecast: %d%%' % int(value * 100))
                else:
                    stdout.write('\rDownloading weather forecast. 100%')
                    stdout.write('\nWeather downloaded.\n')
            elif action == 2:
                stdout.write('Preparing simulation\n')

            stdout.flush()