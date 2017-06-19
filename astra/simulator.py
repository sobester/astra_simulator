# coding=utf-8

"""
This module contains classes for primary flight simulation and the processing
of results.

With all the flight parameters, the flight class defines the framework of
operations necessary to accomplish the balloon flight simulation. See examples
for usage.

University of Southampton
"""
from math import pi
from datetime import timedelta
from sys import stdout
import os
import logging
import logging.handlers
from six.moves import range, builtins
import numpy
from scipy.integrate import odeint
from . import flight_tools as ft
from .weather import forecastEnvironment, soundingEnvironment
from . import global_tools as tools
from . import drag_helium
from . import available_balloons_parachutes
import json

# Pass through the @profile decorator if line profiler (kernprof) is not in use
try:
    builtins.profile
except AttributeError:
    def profile(func):
        return func

# Error and warning logger
logger = logging.getLogger(__name__)


class flightProfile(object):
    """
    Parameters
    ----------
    highestAltitude : scalar
        The highest altitude reached by the flight. This could be the burst
        altitude, if bursting occurred, or the floating altitude if it didn't
    """
    def __init__(self,
                 launchDateTime,
                 nozzleLift,
                 flightNumber,
                 timeVector,
                 latitudeProfile,
                 longitudeProfile,
                 altitudeProfile,
                 highestAltIndex,
                 highestAltitude,
                 hasBurst):
        self.launchDateTime = launchDateTime
        self.nozzleLift = nozzleLift
        self.flightNumber = flightNumber
        self.timeVector = timeVector
        self.latitudeProfile = latitudeProfile
        self.longitudeProfile = longitudeProfile
        self.altitudeProfile = altitudeProfile
        self.highestAltIndex = highestAltIndex
        self.highestAltitude = highestAltitude
        self.hasBurst = hasBurst

    def getJsonPath(self):
        # Every how many points should we store one (this is used to
        # reduce the size of the json file)
        shrinkFactor = 5

        # Build up the flight path data
        jsonPath = ['{ "points" : [\n']
        jsonPoints = []
        for pointNumber, eachPoint in enumerate(numpy.transpose([self.latitudeProfile, self.longitudeProfile, self.altitudeProfile])):
            if pointNumber % shrinkFactor == 0:
                jsonPoints.append('{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f }\n' % (
                    numpy.clip(eachPoint[0], -85.0511, 85.0511), eachPoint[1], eachPoint[2]))

        jsonPath.append(','.join(jsonPoints))
        jsonPath.append('\n] }')

        return jsonPath

    def getKMLPath(self):
        """
        """
        kmlPath = []
        flightHrs, flightMins, flightSecs = tools.prettySeconds(self.timeVector[-1])
        thisFlightHasBurst = self.hasBurst

        kmlPath = [
            '<Placemark>\n<name>Simulation %d</name>\n<styleUrl>#stratoLine</styleUrl>\n<LineString>\n<coordinates>\n' % (
                self.flightNumber)]

        pointNumber = 0
        for eachPoint in numpy.transpose([self.latitudeProfile, self.longitudeProfile, self.altitudeProfile]):
            if thisFlightHasBurst:
                if pointNumber % 10 == 0:
                    kmlPath.append('%.5f,%.5f,%.5f\n' % (eachPoint[1], eachPoint[0], eachPoint[2]))
            else:
                kmlPath.append('%.5f,%.5f,%.5f\n' % (eachPoint[1], eachPoint[0], eachPoint[2]))
            pointNumber += 1

        kmlPath.append(
            '</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</LineString>\n</Placemark>\n')

        return kmlPath

    # def getCSVMarkers(self):

    def getKMLMarkers(self):
        """
        """
        kmlMarkers = []
        # Add balloon launch point
        kmlMarkers.append(
            '<Placemark>\n<name>Balloon Launch</name>\n<styleUrl>#launchPin</styleUrl>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                self.longitudeProfile[0], self.latitudeProfile[0],
                self.altitudeProfile[0]))

        if self.hasBurst:
            # Add balloon burst point
            flightHrs, flightMins, flightSecs = tools.prettySeconds(self.timeVector[-1])

            kmlMarkers.append(
                '<Placemark>\n<name>Balloon burst (Simulation %d)</name>\n<styleUrl>#burstPin</styleUrl>\n<description>Burst altitude: %.0f m.' % (
                    self.flightNumber, self.highestAltitude))
            kmlMarkers.append(
                '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                    self.longitudeProfile[self.highestAltIndex],
                    self.latitudeProfile[self.highestAltIndex],
                    self.highestAltitude))

            # Add balloon landing point
            kmlMarkers.append(
                '<Placemark>\n<name>Payload landing (Simulation %d)</name>\n<styleUrl>#landingPin</styleUrl>\n<description>Flight time: %d hours, %d minutes and %d seconds.' % (
                    self.flightNumber, flightHrs, flightMins, flightSecs))
            kmlMarkers.append(
                '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                    self.longitudeProfile[-1], self.latitudeProfile[-1],
                    self.altitudeProfile[-1]))
        else:
            # Add target altitude reached point
            kmlMarkers.append(
                '<Placemark>\n<name>Target altitude reached (Simulation %d)</name>\n<styleUrl>#burstPin</styleUrl>\n<description>Target altitude: %.0f m.' % (
                    self.flightNumber, self.highestAltitude))
            kmlMarkers.append(
                '</description>\n<Point>\n<coordinates>\n%.5f,%.5f,%.5f\n</coordinates>\n<altitudeMode>absolute</altitudeMode>\n</Point>\n</Placemark>\n' % (
                    self.longitudeProfile[self.highestAltIndex],
                    self.latitudeProfile[self.highestAltIndex], self.highestAlt))
        return kmlMarkers


    def getJsonMarkers(self):
        jsonFloatMarkers = []
        jsonBurstMarkers = []
        jsonLandingMarkers = []
        if self.hasBurst:
            # BALLOON HAS BURST DURING THIS SIMULATION
            # Build up the landing markers data
            flightHrs, flightMins, flightSecs = tools.prettySeconds(self.timeVector[-1])
            landTime = self.launchDateTime + timedelta(seconds=float(self.timeVector[-1]))
            jsonLandingMarkers.append(
                '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Payload Landing Site", "simNumber" : "%d", "otherData" : "Flight time: %d hours, %d mins and %d secs <br /> ETA: %s" }\n' % (
                    self.latitudeProfile[-1], self.longitudeProfile[-1],
                    self.altitudeProfile[-1], self.flightNumber, flightHrs, flightMins,
                    flightSecs, landTime.strftime("%d/%m/%Y %H:%M:%S")))

            # Build up the burst markers data
            jsonBurstMarkers.append(
                '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Balloon Burst", "simNumber" : "%d", "otherData" : "Burst altitude: %.0f m"  }\n' % (
                    self.latitudeProfile[self.highestAltIndex],
                    self.longitudeProfile[self.highestAltIndex],
                    self.highestAltitude,
                    self.flightNumber,
                    self.highestAltitude))
        else:
            # BALLOON HASN'T BURST DURING THIS SIMULATION

            # Build up the float markers data
            if self.highestAltIndex == -1:
                # This is the case when the target altitude was not reached. Show an error to the user.
                jsonFloatMarkers.append(
                    '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Target Altitude NOT Reached", "simNumber" : "%d", "otherData" : "Max altitude: %.0f m <br />Try to increase the flight duration." }\n' % (
                        self.latitudeProfile[self.highestAltIndex],
                        self.longitudeProfile[self.highestAltIndex],
                        self.highestAltitude,
                        self.flightNumber,
                        self.highestAltitude))
            else:
                jsonFloatMarkers.append(
                    '{ "lat" : %.5f, "lng" : %.5f, "alt" : %.5f, "label" : "Target Altitude Reached", "simNumber" : "%d", "otherData" : "Target altitude: %.0f m" }\n' % (
                        self.latitudeProfile[highestAltIndex],
                        self.longitudeProfile[highestAltIndex],
                        self.highestAltitude,
                        self.flightNumber,
                        self.highestAltitude))

        return jsonBurstMarkers, jsonFloatMarkers, jsonLandingMarkers

class flight(object):
    """Primary Balloon flight simulation class.

    Provides methods for solving the ascent rate equation [ref needed], 

    Parameters
    ----------
    environment : :obj:`environment`
        either soundingEnvironment or forecastEnvironment) already configured
        with the environmental model. See the Weather Module for more
        information.
    balloonGasType : string
        'Helium'|'Hydrogen'
    balloonModel : string
        model of the balloon. Accepted values can be found in
        available_balloons_parachutes.py.
    nozzleLift : scalar
        nozzle lift in kg
    payloadTrainWeight : scalar
        weight of the payload train in kg
    parachuteModel : string
        model of the parachute. Accepted values can be found in
        available_balloons_parachutes.py.
    trainEquivSphereDiam : scalar
        the diameter (in meters) of sphere that generates the same drag as the
        train
    numberOfSimRuns : int
        total number of simulations to be run. If 1, a deterministic simulation
        is run. If >1, Monte Carlo perturbations are applied to each simulation
    [excessPressureCoeff] : scalar (default 1)
        TODO: the coefficient of excess pressure of the balloon on the gas
        inside. Currently unused
    [floatingFlight] : bool (default False)
        TRUE if the balloon is going to vent air and float rather than burst.
    [floatingAltitude] : scalar
        the target altitude for the balloon to float, in meters. Ignored if
        floatingFlight is FALSE
    [ventingStart] : scalar (default 1000m)
        how many meters below the target floating altitude should the venting
        start. Ignored if floatingFlight is FALSE.
    [maxFlightTime] : scalar (default 18000 seconds)
        Maximum duration, in seconds, of a flight to be simulated. Warning:
        setting a high maxFlightTime during a floatingFlight could take a long
        time to simulate!
    [outputFile] : string
        the path of the output file containing all the simulation data. The
        format of the file to be generated depends on the extension. Available
        formats are,
        'json' : JavaScript data structure, used to provide data to the
            web interface
        'kml' : Standard format for geographical data. It can be opened by
            Google Maps, Google Earth, etc
        'kmz' : Zipped kml, good for online use (eg Google Maps). The
            file size is significantly reduced;
        'csv' : Comma separated values, the preferred format for further
            data processing with any other software;
        'csv.zip' : Zipped csv file;
        'web' : This stores json, kml and csv.zip in the same folder
            requested. This is used by the web interface to prepare files
            for export. If no extension is found, a new folder is created
            and ALL output formats are stored.
    [debugging] : bool (default False)
        if TRUE, all the information regarding the simulation will be logged.
        If it's set to FALSE, only errors will be logged.
    [log_to_file] : bool
        determines the location where the errors and debug messages will be
        logged. If TRUE, an error.log file will be created in your current
        folder and everything will go in there (make sure you have permissions
        to write in the current folder, or the simulator will not work!)
        If FALSE, all the errors and debug messages will be displayed on the
        terminal or command line.
    [progress_to_file] : bool
        If TRUE, a progress .json file will be created in the current folder
        and will be updated during preflight and then once per 'flight'.
        If FALSE, progress information about the simulation will be displayed
        on the terminal or command line.

    Notes
    -----
    * The run() method takes all the necessary steps to run the simulation.
      For standard users, the use of the run() method is RECOMMENDED.

    * To simply run the simulation, the run() method can be used. :Example:

        >>> my_flight_simulation.run()
    """
    def __init__(self,
                 balloonGasType,
                 balloonModel,
                 nozzleLift,
                 payloadTrainWeight,
                 environment=None,
                 maxFlightTime=18000,
                 parachuteModel=None,
                 numberOfSimRuns=10,
                 trainEquivSphereDiam=0.1,
                 floatingFlight=False,
                 floatingAltitude=None,
                 ventingStart=1000,
                 excessPressureCoeff=1.,
                 cutdown=False,
                 cutdownAltitude=None,
                 outputFile='',
                 debugging=False,
                 log_to_file=False,
                 progress_to_file=False):
        """
        Initialize all the parameters of the object and setup the debugging if
        required.
        """
        self.reset()
        if debugging:
            log_lev = logging.DEBUG
        else:
            log_lev = logging.WARNING

        if log_to_file:
            # Reset the app logger handlers and reset basic config
            # file handler (this may not be the best way to do this for all
            # situations, but it should allow decent logging on the web app
            # while not interrupting with general use of astra as a standalone
            # library
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)

            # Set a maximum log file size of 5MB:
            handler = logging.handlers.RotatingFileHandler('astra_py_error.log',
                mode='a',
                maxBytes=10*1024*1024)
            formatter = logging.Formatter('%(asctime)s.%(msecs)d %(levelname)s:%(name)s %(message)s')
            handler.setFormatter(formatter)

            stream = logging.StreamHandler()
            stream.setFormatter(formatter)

            logging.basicConfig(handlers=[handler, stream], level=log_lev)

            # rootlogger.addHandler(handler)
            # rootlogger.setLevel(log_lev)
        else:
            logger.setLevel(log_lev)

        self._progressToFile = progress_to_file
        self._debugging = debugging
        self.outputFile = outputFile

        # User defined variables
        # Note: As setters are used in some of these variables, there is a
        # small amount of order dependency here: e.g., self.payloadtrainWeight
        # must be defined before nozzleLift 
        self.environment = environment                # weather object
        self.balloonGasType = balloonGasType
        self.balloonModel = balloonModel              # See available_balloons
        self.payloadTrainWeight = payloadTrainWeight  # kg
        self.nozzleLift = nozzleLift                  # kg
        self.parachuteModel = parachuteModel
        self.numberOfSimRuns = numberOfSimRuns
        self.trainEquivSphereDiam = trainEquivSphereDiam     # m
        self.floatingFlight = floatingFlight
        self.floatingAltitude = floatingAltitude             # m
        self.ventingStart = ventingStart     # m (below the target altitude)
        self.excessPressureCoeff = excessPressureCoeff
        self.maxFlightTime = maxFlightTime
        self.cutdown = cutdown
        self.cutdownAltitude = cutdownAltitude

        # Simulation precision - not user defined!
        self._samplingTime = 3               # seconds

        # Note: the launch site latitude, longitude and elevation will be
        # fetched and stored as attributes from the environment object at
        # runtime. Do not enter values here, as they will be ignored

    # ----------------------------------------------------------------------
    # Properties
    # ----------------------------------------------------------------------
    @property
    def samplingTime(self):
        """Note: No setter exists for this, as it should not be set by the user
        """
        return self._samplingTime

    @property
    def launchSiteLat(self):
        return self._launchSiteLat

    @launchSiteLat.setter
    def launchSiteLat(self, new_launchSiteLat):
        if new_launchSiteLat == 0.0:
            logger.warning('Launch Site Latitude set to zero.')
        self._launchSiteLat = new_launchSiteLat

    @property
    def launchSiteLon(self):
        return self._launchSiteLon

    @launchSiteLon.setter
    def launchSiteLon(self, new_launchSiteLon):
        if new_launchSiteLon == 0.0:
            logger.warning('Launch Site Longitude set to zero.')
        self._launchSiteLon = new_launchSiteLon

    @property
    def environment(self):
        return self._environment
    
    @environment.setter
    def environment(self, new_environment):
        if new_environment:

            if not new_environment._weatherLoaded:
                try:
                    new_environment.load(self.updateProgress)
                except:
                    logger.exception("Environment loading failed: see traceback")
                    raise

            self.launchSiteLat = new_environment.launchSiteLat
            self.launchSiteLon = new_environment.launchSiteLon
            self.launchSiteElev = new_environment.launchSiteElev

            # Check if GFS is being used
            if isinstance(new_environment, forecastEnvironment):
                self._usingGFS = True
            else:
                self._usingGFS = False
        else:
            logger.warning("Empty environment loaded. Please load a valid environment")
        self._environment = new_environment


    @property
    def balloonGasType(self):
        return self._balloonGasType

    @balloonGasType.setter
    def balloonGasType(self, new_balloonGasType):
        try:
            self._gasMolecularMass = ft.MIXEDGAS_MOLECULAR_MASS[new_balloonGasType]
        except:
            logger.exception('{} is an invalid gas type'.format(
                new_balloonGasType))
            raise

        self._balloonGasType = new_balloonGasType

    @property
    def balloonModel(self):
        return self._balloonModel

    @balloonModel.setter
    def balloonModel(self, new_balloonModel):

        try:
            balloon_params = available_balloons_parachutes.balloons[new_balloonModel]
        except KeyError:
            logger.exception("""{} is an invalid balloon model.
                See astra.available_balloons_parachutes""".format(
                new_balloonModel))
            raise
        self._balloonWeight = balloon_params[0]

        self._meanBurstDia = (
            available_balloons_parachutes.meanToNominalBurstRatio *
            balloon_params[1]
        )
        self._weibull_lambda = balloon_params[2]
        self._weibull_k = balloon_params[3]

        logger.debug("""Balloon performance: Mean burst diameter: %.4f,
            Lambda: %.4f, k: %4f""" % (
            self._meanBurstDia, self._weibull_lambda, self._weibull_k))

    @property
    def parachuteModel(self):
        return self._parachuteModel

    @parachuteModel.setter
    def parachuteModel(self, new_parachuteModel):
        # Configure the flight tools for the correct parachute
        try:
            self.parachuteAref = \
                available_balloons_parachutes.parachutes[new_parachuteModel]
        except KeyError:
            logger.exception("{} is not a valid parachute name. See astra.available_balloons_parachutes.parachutes")
        self._parachuteModel = new_parachuteModel

    @property
    def nozzleLift(self):
        return self._nozzleLift

    @nozzleLift.setter
    def nozzleLift(self, new_nozzleLift):
        if new_nozzleLift <= 0.0:
            logger.error('Nozzle lift cannot be negative or zero.')
            raise ValueError('Nozzle lift cannot be negative or zero.')
        elif new_nozzleLift <= self.payloadTrainWeight:
            logger.error('The nozzle lift is too low for the balloon to climb! Adjust the nozzle lift.')
            raise ValueError('The nozzle lift is too low for the balloon to climb! Adjust the nozzle lift.')
        self._nozzleLift = new_nozzleLift

    @property
    def payloadTrainWeight(self):
        return self._payloadTrainWeight

    @payloadTrainWeight.setter
    def payloadTrainWeight(self, new_payloadTrainWeight):
        if new_payloadTrainWeight <= 0:
            logger.error('Payload train weight cannot be zero!')
            raise ValueError('Payload train weight cannot be zero!')
        self._payloadTrainWeight = new_payloadTrainWeight

    @property
    def numberOfSimRuns(self):
        return self._numberOfSimRuns

    @numberOfSimRuns.setter
    def numberOfSimRuns(self, new_numberOfSimRuns):
        if new_numberOfSimRuns <= 0:
            logger.error('Number of sim runs cannot be negative or zero.')
            raise ValueError('Number of sim runs cannot be negative or zero.')
        self.initMonteCarloParams(new_numberOfSimRuns)
        self._numberOfSimRuns = new_numberOfSimRuns

    @property
    def outputFile(self):
        return self._outputFile
    
    @outputFile.setter
    def outputFile(self, new_outputFile):
        if new_outputFile:
            try:
                out = open(new_outputFile, 'w')
                out.close()
                os.remove(new_outputFile)
            except IOError:
                if not os.path.isdir(new_outputFile):
                    logger.exception('The output file cannot be created.\n')
                    raise
        # Prepare progress file: try and create the file
        self._progressFile = os.path.splitext(new_outputFile)[0] +\
            '_progress.json'
        self._outputFile = new_outputFile
    
    # ----------------------------------------------------------------------
    def reset(self, keepParameters=False):
        """
        Reset the simulation.

        Parameters
        ----------
        keepParameters : bool (default False)
            If FALSE, all the parameters are deleted and the flight object is
            ready to be reconfigured. If keepParameters is TRUE, all the
            results are deleted but the initial flight configuration is
            preserved. preflight() should still be run again before performing
            anymore simulations.
        """
        if keepParameters:
            self.results = []
            self._lastFlightBurst = False
            self._lastBurstIndex = 0
            self._lastFlightBurstAlt = 0.0
            self._hasRun = False
            self._preflightCompleted = False
        else:
            self.results = []
            self._lastFlightBurst = False
            self._lastBurstIndex = 0
            self._lastFlightBurstAlt = 0.0
            self._hasRun = False
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
            self._balloonVolumeAtInflation = 0.0
            self._balloonDiaAtInflation = 0.0
            self._gasMassAtFloat = 0.0
            self._balloonVolumeAfFloat = 0.0
            self._balloonDiaAtFloat = 0.0
            self._totalAscendingMass = 0.0
            self._totalDescendingMass = 0.0
            self._lastFlightBurstAlt = 0.0
            self._totalStepsForProgress = 0

    @profile
    def run(self):
        """
        This method takes all the necessary steps to run the simulation and
        wraps together all the simulation-specific methods found below. This
        method is RECOMMENDED for standard users.

        Returns 0 if the simulation succeeded, or the error number if there
        were errors. Check the error.log and the simulator's documentation for
        a full description of the error(s)

        Errors occurring during preflight will be raised, stopping execution
        and delivering the full stack trace to both the stdout and the module
        logger (named astra.simulator).
        """
        self.updateProgress(1.0, 2)
        self._totalStepsForProgress = self.numberOfSimRuns + 1

        # _________________________________________________________________ #
        # PREFLIGHT SEQUENCE
        self.results = []
        # ____________________________________________________________ #
        # Variable initialization

        try:
            self._preflight(self.environment.dateAndTime)
        except:
            logger.exception(
                "Error during preflight validations and calculations:")
            raise

        self.updateProgress(0.0, 0)

        # _________________________________________________________________ #
        # RUN THE FLIGHT SIMULATION
        for flightNumber in range(self.numberOfSimRuns):
            logger.debug('SIMULATING FLIGHT %d' % (flightNumber + 1))
            result, _ = self.fly(flightNumber, self.environment.dateAndTime, runPreflight=False)
            self.results.append(result)
            self.updateProgress(
                float(flightNumber + 1) / self._totalStepsForProgress, 0)

        # _________________________________________________________________ #
        # POSTFLIGHT HOUSEKEEPING AND RESULT PROCESSING
        self._hasRun = True
        self.postflight()
        self.updateProgress(1.0, 0)
        return None

    def initMonteCarloParams(self, numberOfSimRuns):
        # Monte Carlo params
        self._lowCD = []
        self._highCD = []
        self._transition = []
        self._ReBand = []
        self._burstDiameter = []
        self._parachuteCD = []
        self._balloonReturnFraction = []

        # Initialize variables subject to perturbation (Monte Carlo simulation)
        if numberOfSimRuns == 1:
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
            for _ in range(numberOfSimRuns):
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
            self.environment.perturbWind(numberOfSimRuns)

    @profile
    def _preflight(self, launchDateTime):
        """
        Run a series of pre-flight checks and calculations to verify the
        consistency of the flight parameters entered and to prepare all the
        data required for simulation.

        It is MANDATORY to execute this method before the simulation (otherwise
        the simulation will throw an error)

        If successful, no errors are thrown. Enable debugging for detailed
        information.

        Parameters
        ----------
        launchDateTime : :obj:`datetime.datetime`
        """
        logger.debug('Preflight sequence starting...')

        # _________________________________________________________________ #
        # Lifting gas
        # Use flight tools to calculate all preliminary balloon calculations
        # (gas mass, balloon volume and diameter at inflation)
        (self._gasMassAtInflation, self._balloonVolumeAtInflation,
            self._balloonDiaAtInflation) = ft.liftingGasMass(
                self.nozzleLift,
                self._balloonWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                             self.launchSiteLon,
                                             self.launchSiteElev,
                                             launchDateTime),
                self._gasMolecularMass,
                self.excessPressureCoeff
        )

        logger.debug('Lifting gas mass calcs completed!')
        logger.debug('Gas molecular mass: %.4f' % self._gasMolecularMass)
        logger.debug('Gas mass: %.4f' % self._gasMassAtInflation)
        logger.debug('Balloon volume at inflation: %.4f' %
                     self._balloonVolumeAtInflation)
        logger.debug('Balloon diameter at inflation: %.4f' %
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
                self._balloonDiaAtFloat) = ft.liftingGasMass(
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
                                                 launchDateTime),
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

        self._preflightCompleted = True

        logger.debug('Preflight completed!')

    @profile
    def fly(self, flightNumber, launchDateTime, runPreflight=True):
        """
        Execute a single simulation.
        It should be run N times, where N is the number of simulation runs
        specified upon configuration in the numberOfSimRuns variable.
        flightNumber should have values 0...N-1

        If successful, no errors are thrown. Enable debugging for detailed
        information.

        Parameters
        ----------
        flightNumber : int
            the number of the flight, to be used as the index for getting
            Monte Carlo parameters initialised in preflight
        launchDateTime : :obj:`datetime.datetime
            The date and time of the launch (note that this may be different
            from the start window of the self.environment object)
        [runPreflight] : bool, default True
            Runs the preflight calculations if True. Useful to switch this off
            if a uniform launchDateTime is used across all flights (as is the
            case with the Monte Carlo simulation in self.run)
        Returns
        -------
        result : list of numpy array
            arrays for flight number, time, 
        """
        # Check whether the preflight sequence was performed. If not, stop the
        # simulation.
        if runPreflight:
            try:
                self._preflight(self.environment.dateAndTime)
            except:
                logger.exception(
                    "Error during preflight validations and calculations:")
                raise
        else:
            assert(self._preflightCompleted),\
                "'preflight' must be run at least once before attempting a flight"
        # Flight-specific variables initialization
        # Prepare the flight tools to deliver results specific to this flight.
        highCD = self._highCD[flightNumber]
        lowCD = self._lowCD[flightNumber]
        ReBand = self._ReBand[flightNumber]
        transition = self._transition[flightNumber]
        parachuteCD = self._parachuteCD[flightNumber]
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
        logger.debug('Low CD: %.4f' % lowCD)
        logger.debug('High CD: %.4f' % highCD)
        logger.debug('Transition: %.4f' % transition)
        logger.debug('Re Band: %.4f' % ReBand)
        logger.debug('Burst Diameter: %.4f' %
                     self._burstDiameter[flightNumber])
        logger.debug('Balloon Return Fraction: %.4f' %
                     self._balloonReturnFraction[flightNumber])
        logger.debug('Parachute CD: %.4f' % parachuteCD)

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
                if altitude >= self.cutdownAltitude:
                    # Burst the balloon
                    logger.debug('Bursting the balloon at {}m altitude'.format(
                        altitude))
                    self._lastFlightBurst = True
                    self._lastFlightBurstAlt = altitude
                    return numpy.array([0.0, 0.0])

            # Calculate current position in time and space
            # This is used to fetch the correct atmospheric data if the GFS is
            # being used.
            currentTime = launchDateTime + timedelta(seconds=t)

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
                    gasMass = ft.gasMassForFloat(
                        altitude,
                        self.floatingAltitude,
                        self._gasMassAtInflation,
                        self._gasMassAtFloat,
                        ventStart=self.ventingStart
                    )
                else:
                    gasMass = self._gasMassAtInflation


                # Calculate current balloon diameter to check for burst
                gasTemp = tools.c2kel(self.environment.getTemperature(
                    self._currentLatPosition, self._currentLonPosition,
                    altitude, currentTime))
                gasDensity = (self.excessPressureCoeff *
                    self.environment.getPressure(self._currentLatPosition,
                        self._currentLonPosition, altitude,currentTime) * 100
                        * self._gasMolecularMass / (8.31447 * gasTemp))
                balloonVolume = gasMass / gasDensity
                balloonDiameter = (6 * balloonVolume / pi) ** (1. / 3)

                # If floating flight, calculate the nozzle lift if the gas is
                # being vented.
                if self.floatingFlight:
                    nozzleLift = ft.nozzleLiftForFloat(
                        self.nozzleLift,
                        self.environment.getDensity(self._currentLatPosition,
                                                    self._currentLonPosition,
                                                    altitude,
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
                    currentDensity = self.environment.getDensity(
                        self._currentLatPosition, self._currentLonPosition,
                        altitude, currentTime)
                    currentViscosity = self.environment.getViscosity(
                        self._currentLatPosition, self._currentLonPosition,
                        altitude, currentTime)
                    balloonDrag = ft.balloonDrag(balloonDiameter,
                                                 ascentRate,
                                                 currentDensity,
                                                 currentViscosity,
                                                 lowCD,
                                                 highCD,
                                                 ReBand,
                                                 transition
                                                 )
                    trainDrag = ft.balloonDrag(self.trainEquivSphereDiam,
                                                   ascentRate,
                                                   currentDensity,
                                                   currentViscosity,
                                                   lowCD,
                                                   highCD,
                                                   ReBand,
                                                   transition)

                    # External Forces
                    externalForces = ((nozzleLift - self.payloadTrainWeight) *
                        9.81) - balloonDrag - trainDrag

                    # Derivatives
                    dvdt = externalForces / self._totalAscendingMass
                    dhdt = ascentRate

                    return numpy.array([dhdt, dvdt])

                else:
                    # THE BALLOON HAS BURST

                    # Floating flight is set to false because if the balloon
                    # has burst, the flight is now standard.
                    self._lastFlightBurst = True
                    self._lastFlightBurstAlt = altitude
                    return numpy.array([0.0, 0.0])

            else:
                # THE BALLOON IS CURRENTLY DESCENDING


                # Calculate parachute Drag
                currentDensity = self.environment.getDensity(
                    self._currentLatPosition, self._currentLonPosition,
                    altitude, currentTime)
                currentViscosity = self.environment.getViscosity(
                    self._currentLatPosition, self._currentLonPosition,
                    altitude, currentTime)
                if self.parachuteModel == None:
                    parachuteDrag = 0
                else:
                    parachuteDrag = ft.parachuteDrag(abs(ascentRate),
                        currentDensity, self.parachuteAref, parachuteCD)
                # Train Drag
                trainDrag = abs(ft.balloonDrag(self.trainEquivSphereDiam,
                                                   abs(ascentRate),
                                                   currentDensity,
                                                   currentViscosity,
                                                   lowCD, highCD,
                                                   ReBand,
                                                   transition))

                # External Forces
                externalForces = -self.payloadTrainWeight * 9.81 * (
                    1 + self._balloonReturnFraction[flightNumber])\
                    + parachuteDrag + trainDrag

                # Derivatives
                dvdt = externalForces / self._totalDescendingMass
                dhdt = ascentRate

                return numpy.array([dhdt, dvdt])

        # Define the initial conditions, the time vector at which we want
        # simulation data to be stored, and run the integration.
        # Note: the simulation carries on all the way to the maxFlightTime,
        # even if the altitude becomes negative. Negative values of altitude
        # will be trimmed later on.
        initialConditions = numpy.array([self.launchSiteElev, 0.0])
        timeVector = numpy.arange(0, self.maxFlightTime + self.samplingTime,
            self.samplingTime)

        logger.debug('Beginning integration.')

        ### INTEGRATION ###
        if self._usingGFS:
            solution = odeint(ode, initialConditions, timeVector, rtol=1e-3,
                atol=1e-3)
        else:
            solution = odeint(ode, initialConditions, timeVector)
        ###################

        logger.debug('Integration completed. Post-processing...')

        # Extract altitude and ascent rate data from the solution array
        solution_altitude = numpy.array(solution[:, 0])
        # solution_ascrate = numpy.array(solution[:,1]) ### Currently not used.

        # Trim negative altitude values from results and then trim time to the
        # same length.
        if self._lastFlightBurst:
            solution_altitude = solution_altitude[solution_altitude > 0]
            solution_altitude[-1] = 0.0

        ### Currently not used.
        #solution_ascrate = solution_ascrate[:len(solution_altitude)]

        timeVector = timeVector[:len(solution_altitude)]

        # Calculate drift
        lastDriftLat = 0.0
        lastDriftLon = 0.0
        latitudeProfile = [self.launchSiteLat]
        longitudeProfile = [self.launchSiteLon]
        for eachAlt, eachTime in zip(solution_altitude[1:], timeVector[1:]):
            # Calculate the position of the balloon at each point and use it to
            # work out its location
            currentTime = launchDateTime + timedelta(
                seconds=float(eachTime))
            # Gather wind speed and convert to m/s
            windSpeed = currentFlightWindSpeed(latitudeProfile[-1],
                longitudeProfile[-1], eachAlt, currentTime) * 0.514444
            # Convert the wind to u- and v-coordinates
            windLon, windLat = tools.dirspeed2uv(
                currentFlightWindDirection(latitudeProfile[-1],
                    longitudeProfile[-1], eachAlt, currentTime),
                windSpeed)
            # Store the drift in meters (this is the distance between the
            # LAUNCH SITE and the current location)
            lastDriftLat += windLat * self.samplingTime
            lastDriftLon += windLon * self.samplingTime
            # Convert it to degrees
            dLat, dLon = tools.m2deg(lastDriftLat, lastDriftLon,
                latitudeProfile[-1])
            # Store the new latitude and longitude
            latitudeProfile.append(self.launchSiteLat + dLat)
            longitudeProfile.append(self.launchSiteLon + dLon)


        # Check that latitude and longitude are within bounds and correct if
        # they are not (for example, if the balloon flew over the North or the
        # South Pole).
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
            index = tools.find_nearest_index(solution_altitude,
                self._lastFlightBurstAlt)
        else:
            if solution_altitude[-1] < self.floatingAltitude - 100:
                # In this case, the balloon hasn't reached the target altitude.
                # This is probably because the maxFlightTime is too low:
                #  Show an error.
                index = -1
            else:
                # Store target altitude reached index
                index = tools.find_nearest_index(solution_altitude,
                    self.floatingAltitude)

        # STORE RESULTS OF CURRENT SIMULATION
        if self._lastFlightBurst:
            highestAltitude = self._lastFlightBurstAlt
        else:
            highestAltitude = self.floatingAltitude


        resultProfile = flightProfile(launchDateTime, self.nozzleLift,
            flightNumber + 1, timeVector, latitudeProfile,
                  longitudeProfile, solution_altitude, index,
                  highestAltitude, self._lastFlightBurst)

        logger.debug('Simulation completed.')

        return resultProfile, solution

    def write_JSON(self, filename):
        """
        """
        # GENERATE JSON FILE OUT OF RESULTS

        # Initialize the data lists
        jsonBurstMarkers = []
        jsonLandingMarkers = []
        jsonFloatMarkers = []
        jsonPaths = []

        # Go through the results of each flight simulation
        for profile in self.results:
            thisFlightHasBurst = profile.hasBurst

            jsonPath = profile.getJsonPath()
            jsonPaths.append(''.join(jsonPath))

            burstMarkers, floatMarkers, landingMarkers = profile.getJsonMarkers()
            jsonBurstMarkers.extend(burstMarkers)
            jsonLandingMarkers.extend(burstMarkers)
            jsonFloatMarkers.extend(burstMarkers)


        jsonPaths = ['"flightPaths" : [\n', ','.join(jsonPaths), ']']

        # Put all the non-empty lists above together
        if len(jsonBurstMarkers) != 0:
            jsonBurstMarkers = ['"burstMarkers" : [\n', ','.join(jsonBurstMarkers), '],']

        if len(jsonLandingMarkers) != 0:
            jsonLandingMarkers = ['"landingMarkers" : [\n', ','.join(jsonLandingMarkers), '],']

        if len(jsonFloatMarkers) != 0:
            jsonFloatMarkers = ['"floatMarkers" : [\n', ','.join(jsonFloatMarkers), '],']

        # Put them all together in one single array
        jsonToAdd = [
            '{',
            ''.join(jsonBurstMarkers),
            ''.join(jsonLandingMarkers),
            ''.join(jsonFloatMarkers),
            ''.join(jsonPaths),
            '}'
        ]

        # Try to gain write permission on the output file specified upon configuration
        try:
            jsonFile = open(filename, 'w')
        except IOError:
            logger.error('Cannot create output file.')
            return

        # Write all the data to the file
        jsonFile.write(''.join(jsonToAdd))
        jsonFile.close()

    def write_KML(self, filename, zipped=False):
        """GENERATE KML OUT OF RESULTS

        Also allows zipped kml (kmz) if .kmz is provided in the extension
        """

        # Author of the KML file
        kmlAuthor = 'ASTRA High Altitude Balloon Flight Planner, University of Southampton'
        launchPinURL = 'http://maps.google.com/mapfiles/ms/micons/red-dot.png'
        burstPinURL = 'http://maps.google.com/mapfiles/ms/micons/yellow-dot.png'
        landingPinURL = 'http://maps.google.com/mapfiles/ms/micons/red-dot.png'

        kmlPaths = []
        kmlMarkers = []

        for profile in self.results:
            kmlPaths.append(''.join(profile.getKMLPath()))
            kmlMarkers.append(''.join(profile.getKMLMarkers()))


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

        if zipped:
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

            with zipfile.ZipFile(filename, 'w') as zipKmz:
                zipKmz.write(outputKml.name, arcname='output.kml',
                    compress_type=zipCompression)
            zipKmz.close()
            outputKml.close()

            logger.debug(('KMZ file generated! ', filename))

        else:
            # Write the KML file.

            try:
                outputKml = open(filename, 'w')
            except IOError:
                logger.error('Error: cannot create the output file in the given directory')
                return

            outputKml.write(''.join(kmlToAdd))
            outputKml.close()

    def write_CSV(self, filename, zipped=False):
        """GENERATE CSV FILE OUT OF RESULTS
        """
        # Calculate how many rows we need in total
        totalRows = 0
        for profile in self.results:
            totalRows += len(profile.timeVector)

        # Columns are Flight #, Time from launch, Lat, Lon, Alt, Remarks
        totalColumns = 6

        # Add one column for the header and generate an empty matrix
        csvMatrix = numpy.zeros((totalRows + 2, totalColumns), dtype='|S100')
        # Clear the remarks column (it was full of zeros)
        # csvMatrix[:, 5] = ''

        currentPosition = 2
        # Populate matrix with data
        for profile in self.results:
            numberOfPoints = len(profile.latitudeProfile)
            thisFlightHasBurst = profile.hasBurst

            # Flight number, Time, Lat, Lon, Alt
            # for i in range(5):
            csvMatrix[currentPosition:currentPosition + numberOfPoints, :5] = \
            numpy.column_stack([numpy.ones(numberOfPoints)*profile.flightNumber, profile.timeVector,
            profile.latitudeProfile, profile.longitudeProfile,
            profile.altitudeProfile])

            # Remarks: Launch
            csvMatrix[currentPosition, 5] = 'Balloon Launch'

            if thisFlightHasBurst:
                # Remarks: Burst
                csvMatrix[currentPosition + profile.highestAltIndex, 5] = 'Balloon Burst'

                # Remarks: Landing
                csvMatrix[currentPosition + numberOfPoints - 1, 5] =\
                    'Balloon Landed'
            else:
                # Remarks: Target Altitude Reached
                csvMatrix[currentPosition + self.highestAltIndex, 5] =\
                    'Balloon has reached target altitude'

            currentPosition += numberOfPoints

        # Add header
        csvMatrix[0] = [
            'Data generated by the ASTRA High Altitude Balloon Flight Planner - University of Southampton', '',
            '', '', '', '']

        csvMatrix[1] = ['Simulation #', 'Time from launch [s]', 'Latitude',
            'Longitude', 'Altitude [m]', 'Remarks']

        # Save file
        if not zipped:

            numpy.savetxt(filename, csvMatrix, delimiter=',', fmt='%s')

        else:
            import zipfile, tempfile

            try:
                outputCsv = tempfile.NamedTemporaryFile()
            except IOError:
                logger.error('Error: cannot create a new temporary file')
                return

            # Create a zipped version of the file, since it's probably big
            # (approx 150KB/simulation not compressed)
            numpy.savetxt(outputCsv, csvMatrix, delimiter=',', fmt='%s')
            outputCsv.flush()

            # Zip CSV file
            try:
                import zlib

                zipCompression = zipfile.ZIP_DEFLATED
            except:
                zipCompression = zipfile.ZIP_STORED

            with zipfile.ZipFile(filename, 'w') as zipCsv:
                zipCsv.write(outputCsv.name,
                    arcname='ASTRA Simulation Results.csv',
                    compress_type=zipCompression)
            zipCsv.close()
            outputCsv.close()

            logger.debug(('CSV-ZIP file generated! ', self.outputFile))

    def write(self, filename):
        """
        Function responsible for storing the data in the given format.

        The format of the file to be generated depends on the extension given
        in the outputFile.

        Parameters
        ----------
        filename : string
            the name of the output file. Accepted extensions are:
            'json'  JavaScript data structure, used to pass data to the web
                interface;
            'kml'   Standard format for geographic data. It can be opened by
                Google Maps, Google Earth, etc;
            'kmz'   Zipped kml, good for online use (eg Google Maps). The file
                size is significantly reduced;
            'csv'   Comma separated values, the preferred format for further
                data processing with any other software;
            'web'   This stores json, kml and csv in the same folder requested.
                This is used by the web interface to prepare files for export.
            none    If no extension is found, a new folder is created and ALL
                output formats are stored.
        """
        # Extract the output format from the outputFile path.
        fileExtension = os.path.splitext(filename)[1]

        # Remove the dot from the extension
        data_format = fileExtension[1:]

        logger.debug('Attempting to store data... Requested format: %s' % data_format)

        # {extension: (write_function, kwargs)}
        # formats = {'json':(write_JSON, {}),
        #            'kml': (write_KML, {'zipped': False}),
        #            'kmz': (write_KML, {'zipped': True})
        #            'csv': (write_CSV, {'zipped': False})
        #            'zip': (write_CSV, {'zipped': True})
        #            }
        # write_func, kwargs = formats[data_format.lower()]

        if data_format in ('json', 'JSON'):
            self.write_JSON(filename)

        elif data_format in ('kml', 'KML'):
            self.write_KML(filename, zipped=False)

        elif data_format in ('kmz', 'KMZ'):
            self.write_KML(filename, zipped=True)

        elif data_format in ('csv', 'CSV'):
            self.write_CSV(filename, zipped=False)

        elif data_format in ('zip', 'ZIP'):
            # This is used for .csv.zip double extension:
            self.write_CSV(filename, zipped=True)

        else:
            logger.error('Output data format not understood! Is the specified extension correct?')

        # Removing this temporarily: needs further success checks, but
        # os.path.isfile doesn't check if the file has been written recently.
        # logger.debug('Output file {} generated.')

    def postflight(self):
        """
        After all the simulations have been executed, this method puts the
        results together, processes them and stores them in the outputFormat
        required.

        The results are stored in the outputFile specified upon configuration.

        See Also
        --------
        astra.simulator.flight, self.outputFile
        """
        baseName, data_format = os.path.splitext(self.outputFile)

        if data_format == '':
            # In the case that no extension is provided, use the base name as
            # the folder name as write all types to this file
            try:
                os.mkdir(self.outputFile)
            except OSError:
                if not os.path.isdir(self.outputFile):
                    logger.error('The specified output path already exists. Change it or add an extension.')

            for data_format in ['json', 'kml', 'kmz', 'csv', 'csv.zip']:
                path = os.path.join(baseName, 'out' + '.' + data_format)
                try:
                    self.write(path)
                except:
                    logger.exception("flight simulator failed to write output {} file".format(data_format))
                    raise
        elif data_format == '.web':
            # If file has extension web, store json, kml and csv.gz

            for data_format in ['json', 'kml', 'csv.zip']:
                path = baseName + '.' + data_format
                try:
                    self.write(path)
                except:
                    logger.exception("flight simulator failed to write output {} file".format(data_format))
                    raise

        else:
            # Only store the required one
            self.write(self.outputFile)

    def updateProgress(self, value, action):
        """
        Update the progress file with the ratio of value and the total steps of
        the simulation calculated when the method run() is executed.

        Notes
        -----
        * if run() is not being used, the object's parameter
        _totalStepsForProgress should be defined before executing this method
        for it to work properly.
        """
        if self._progressToFile or self._debugging:
            # Try to gain write permissions to the progress file
            try:
                progFile = open(self._progressFile, 'w')
            except IOError:
                # In this case, opt for error rather than exception as failing
                # to open the progress log is not a fatal error
                logger.error('Cannot open progress file.')
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
