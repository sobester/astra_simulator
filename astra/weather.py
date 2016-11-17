# coding=utf-8

"""
weather.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Weather Module
This module is responsible for dealing with weather and the modelling of the atmosphere. It works closely together
with the Simulator Module, providing it with all the atmospheric data throughout the entire simulation.

This module contains the environment class, providing a common interface to all its subclasses. This allows the
Simulator Module not to care of what kind of environment is being used. If the interface with the other objects
is common, the Simulator Module can deal with them in exactly the same way. This approach also provides great
flexibility for extending the atmospheric model in the future: a new environment subclass can be defined and be
implemented very quickly in the current program, without needing to modify any code in any other module.

There are two environment subclasses:
    forecastEnvironment: responsible for a GFS-based atmospheric model
    soundingEnvironment: responsible for an atmospheric model generated from a sounding file


USAGE
--------------

The environment class defines a common interface for the Simulator Module to deal with.
Such interface has the following variables and methods:
    VARIABLES
        launchSiteLat: latitude of the launch site [deg] (type: float)
        launchSiteLon: longitude of the launch site [deg] (type: float)
        launchSiteElev: elevation of the launch site above Mean Sea Level [m] (type: float)
        dateAndTime: the launch local time (type: datetime.datetime)
        UTC_offset: the local time zone's offset from UTC, in hours. See NOTE below about automatic UTC offset
            (type: float)
        inflationTemperature: the ambient temperature during the balloon inflation [degC] (type: float)
    METHODS
        getTemperature(lat,lon,alt,time): request the temperature for the point at the given location at the given time.
            Returns a float, temperature in [degrees Celsius]
        getPressure(lat,lon,alt,time): request the pressure for the point at the given location at the given time.
            Returns a float, pressure in [millibar]
        getWindDirection(lat,lon,alt,time): request the wind direction for the point at the given location at the given time.
            Returns a float, direction in [degrees] clockwise from north
        getWindSpeed(lat,lon,alt,time): request the wind speed for the point at the given location at the given time.
            Returns a float, speed in [knots]
        getDensity(lat,lon,alt,time): request the density for the point at the given location at the given time.
            Returns a float, density in [kg/m3]
        getViscosity(lat,lon,alt,time): request the viscosity for the point at the given location at the given time.
            Returns a float, viscosity in [Pa s]

        Note: for all the methods above, lat [deg], lon [deg], alt [m], time [datetime.datetime].

Refer to the individual classes and subclasses for details on how to use them.

NOTE: AUTOMATIC UTC OFFSET
--------------

The UTC offset can be automatically retrieved for any LAND location by using Google Maps API's Time Zone service.
By default, the UTC offset is automatically retrieved using the launch site GPS location if the UTC_offset is set to 0.

See documentation for the global_tools.getUTCOffset() function for more information.



University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from datetime import timedelta
import logging

import numpy
from scipy.interpolate import UnivariateSpline

import global_tools as tools
import GFS

# Error and warning logger
logger = None


class environment:
    """
    This class defines a common interface for the Simulator Module to deal with.
    Such interface has the following variables and methods:
    VARIABLES
        launchSiteLat: latitude of the launch site [deg] (type: float)
        launchSiteLon: longitude of the launch site [deg] (type: float)
        launchSiteElev: elevation of the launch site above Mean Sea Level [m] (type: float)
        dateAndTime: the launch time (type: datetime.datetime)
        UTC_offset: the offset in hours between the current time zone and UTC (for example, Florida in winter has a
            UTC_offset = -5) (type: float)
    METHODS
        getTemperature(lat,lon,alt,time): request the temperature for the point at the given location at the given time.
            Returns a float, temperature in [degrees Celsius]
        getPressure(lat,lon,alt,time): request the pressure for the point at the given location at the given time.
            Returns a float, pressure in [millibar]
        getWindDirection(lat,lon,alt,time): request the wind direction for the point at the given location at the given time.
            Returns a float, direction in [degrees] clockwise from north
        getWindSpeed(lat,lon,alt,time): request the wind speed for the point at the given location at the given time.
            Returns a float, speed in [knots]
        getDensity(lat,lon,alt,time): request the density for the point at the given location at the given time.
            Returns a float, density in [kg/m3]
        getViscosity(lat,lon,alt,time): request the viscosity for the point at the given location at the given time.
            Returns a float, viscosity in [Pa s]

        Note: for all the methods above, lat [deg], lon [deg], alt [m], time [datetime.datetime].
    """

    def __init__(self, debugging=False, log_to_file=False):
        global logger

        # Initialize the object

        # COMMON INTERFACE

        # Variables
        self.inflationTemperature = 0.0
        self.launchSiteLat = 0.0
        self.launchSiteLon = 0.0
        self.launchSiteElev = 0.0
        self.dateAndTime = None
        self.UTC_offset = 0
        # Calculated parameters
        self.getTemperature = None
        self.getPressure = None
        self.getWindDirection = None
        self.getWindSpeed = None
        self.getDensity = None
        self.getViscosity = None
        self._UTC_time = None
        # Monte Carlo parameters
        self.getMCWindDirection = []
        self.getMCWindSpeed = []

        self.debugging = debugging
        self.file_logging = log_to_file

        self._weatherLoaded = False

        # SETUP ERROR LOGGING AND DEBUGGING

        logger = logging.getLogger('Environment')

        if debugging:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.WARNING)

        if log_to_file:
            log_handler = logging.FileHandler(filename='error.log')
            log_formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s')
        else:
            log_handler = logging.StreamHandler()
            log_formatter = logging.Formatter('%(levelname)s - %(name)s - %(message)s')

        if debugging:
            log_handler.setLevel(logging.DEBUG)
        else:
            log_handler.setLevel(logging.WARNING)

        log_handler.setFormatter(log_formatter)
        logger.addHandler(log_handler)


class soundingEnvironment(environment):
    """
    The soundingEnvironment class is responsible of generating a sounding-based atmospheric model.

    USAGE
    --------------

    Instantiate a new soundingEnvironment object if you are providing a sounding file to generate an atmospheric model.
    Accepted formats for sounding files are .ftr and .sounding .

    The new object can be instantiated and initialized using the following syntax:

        import weather
        my_sounding_atmosphere = weather.soundingEnvironment()

    Upon instantiation, two optional parameters can be passed: debugging and log_to_file. By default they are both FALSE.
    If log_to_file is set to TRUE, all error and debug logs will be stored in a file called error.log.
    If FALSE, all error and debug logs will be directly displayed on the terminal or command line.
    If debugging is TRUE, all the information available will be logged. Otherwise, only errors will be logged.

    Once the object has been initialized, it needs to be configured with all the environmental data required.
    The following parameters need to be defined (parameters in [brackets] are optional):
        launchSiteLat: latitude of the launch site [deg] (type: float)
        launchSiteLon: longitude of the launch site [deg] (type: float)
        launchSiteElev: elevation of the launch site above Mean Sea Level [m] (type: float)
        dateAndTime: the launch time (type: datetime.datetime)
        UTC_offset: the offset in hours between the current time zone and UTC (for example, Florida in winter has a
            UTC_offset = -5) (type: float)
        inflationTemperature: the ambient temperature during the balloon inflation [degC] (type: float)
        distanceFromSounding: the distance between the launch site and the location where the sounding file was
            recorded [km] (type: float)
        timeFromSounding: the duration between the balloon launch and the time when the sounding file was recorded
            [hours] (type: float)
        [maxAltitude]: the upper altitude limit below which the atmosphere is modelled. Default is 50999m [m] (type: float)

    Once all the parameters have been set, the sounding file can then be loaded and the atmosphere modelled.
    The standard syntax to load the sounding file is the following:

        my_sounding_atmosphere.loadSounding(filePath)

    where filePath is a string containing either the relative or the absolute path to the sounding file.

    If the loadSounding(filePath) method completes successfully, no errors are thrown and the soundingForecast object
    is now ready for used.

    Other methods available:
        perturbWind()   Perform a wind perturbation for the purpose of Monte Carlo simulations.

    Refer to the error.log file for any errors.


    EXAMPLE
    --------------

    A typical usage of the soundingEnvironment object would be as follows:

        from datetime import datetime
        import weather

        my_sounding_atmosphere = weather.soundingEnvironment()
        my_sounding_atmosphere.launchSiteLat = 50.2245
        my_sounding_atmosphere.launchSiteLon = -5.3069
        my_sounding_atmosphere.launchSiteElev = 60
        my_sounding_atmosphere.dateAndTime = datetime.now()
        my_sounding_atmosphere.UTC_offset = 3
        my_sounding_atmosphere.inflationTemperature = 10.5
        my_sounding_atmosphere.distanceFromSounding = 36
        my_sounding_atmosphere.timeFromSounding = 2.8

        my_sounding_atmosphere.loadSounding('/path/to/mySoundingFile.ftr')

    """

    def __init__(self, debugging=False, log_to_file=False):
        """
        Initialize the soundingEnvironment object
        """

        # Run the environment class initialization first
        environment.__init__(self, debugging, log_to_file)

        # Initialize extra sounding-specific variables
        self.distanceFromSounding = 0.0
        self.timeFromSounding = 0.0
        self.maxAltitude = 50999

        self._interpolationPrecision = 200


    def loadSounding(self, soundingFile):
        """
        Load the sounding file and generate the atmospheric model.

        The argument soundingFile should be a string containing either the relative or the absolute path to the
        sounding file.

        This method first reads the file, loads all the data, validates it, and stores it.
        Once this is done, all the data is interpolated (cubic spline interpolation) and the interpolators are then
        stored in the getTemperature, getPressure, getWindDirection, getWindSpeed, getDensity and getViscosity variables.
        This allows a user to request the value of any of the variables by calling the get... methods.

        Note: since sounding files only provide information for one location and at one point in time, the lat, lon and
        time arguments for the get... methods are optional.
        """

        def process_sounding_data(PRESS, HGHT, TEMP, DRCT, SKNT):
            """
            This is an internal method responsible for processing raw data fetched from the sounding file.

            It requires pressure, altitude, temperature, wind speed and direction as inputs.
            All the inputs are cleaned up, prepared for interpolation, interpolated and then stored in the appropriate
            variables.

            This method is independent of the sounding file format
            """

            # Convert to Numpy arrays
            HGHT = numpy.array(HGHT)
            PRESS = numpy.array(PRESS)
            TEMP = numpy.array(TEMP)
            DRCT = numpy.array(DRCT)
            SKNT = numpy.array(SKNT)

            # Remove duplicate altitudes and trim the other data accordingly
            HGHT, indexes = numpy.unique(HGHT, return_index=True)
            PRESS = PRESS[indexes]
            TEMP = TEMP[indexes]
            DRCT = DRCT[indexes]
            SKNT = SKNT[indexes]

            logger.debug('Duplicate altitudes removed.')

            # Remove NaN entries, if any
            nanValidator = numpy.array([HGHT, PRESS, TEMP, DRCT, SKNT])
            nanValidator.transpose()
            nanFree = nanValidator[~numpy.isnan(nanValidator).any(1)]
            nanFree.transpose()
            HGHT = numpy.array(nanFree[0])
            PRESS = numpy.array(nanFree[1])
            TEMP = numpy.array(nanFree[2])
            DRCT = numpy.array(nanFree[3])
            SKNT = numpy.array(nanFree[4])

            logger.debug('NaN entries removed.')

            # _______________________________________________________________________________________________ #
            # Add missing data if launch site elevation or max altitude are out of sounding bounds

            # If the elevation is lower than the lower bound of data, copy the lowest data available and use it for
            # the launch site elevation.
            if self.launchSiteElev < HGHT[0]:
                HGHT = numpy.insert(HGHT, 0, self.launchSiteElev)
                PRESS = numpy.insert(PRESS, 0, PRESS[0])
                TEMP = numpy.insert(TEMP, 0, TEMP[0])
                DRCT = numpy.insert(DRCT, 0, DRCT[0])
                SKNT = numpy.insert(SKNT, 0, SKNT[0])
                logger.debug('Launch site elevation out of bounds. Low altitude data generated.')

            # If the maxAltitude is higher than the upper bound of data, fill the missing altitude levels in with ISA data
            if self.maxAltitude > HGHT[-1]:
                newHeights = numpy.arange(HGHT[-1] + 1, self.maxAltitude + 1, (self.maxAltitude - HGHT[-1] - 1) / 20.)
                HGHTTEMP = numpy.append(HGHT, newHeights[3:])
                for newTempHeight in newHeights[3:]:
                    # Calculate new temperatures and pressures.
                    _, newTemp, _, newPress, _ = tools.ISAatmosphere(altitude=tools.m2feet(newTempHeight))
                    TEMP = numpy.append(TEMP, newTemp)
                    PRESS = numpy.append(PRESS, newPress)

                if self.maxAltitude > 25000 and HGHT[-1] < 25000:
                    HGHTSKNT = numpy.append(HGHT, [25000, self.maxAltitude])
                    SKNT = numpy.append(SKNT, [0, 0])
                else:
                    HGHTSKNT = numpy.append(HGHT, [self.maxAltitude])
                    SKNT = numpy.append(SKNT, [0])

                HGHT = numpy.append(HGHT, self.maxAltitude)
                DRCT = numpy.append(DRCT, DRCT[-1])

                logger.debug('Max altitude out of bounds. High altitude data generated.')
            else:
                HGHTTEMP = HGHT
                HGHTSKNT = HGHT

            # _______________________________________________________________________________________________ #

            # Check whether all fields have data, otherwise fill up vectors with NaN
            if numpy.size(HGHT) == 0 or numpy.size(TEMP) == 0 or numpy.size(PRESS) == 0 or numpy.size(
                    DRCT) == 0 or numpy.size(SKNT) == 0:
                logger.error('There was a problem while processing the sounding.')
                return False

            # Interpolate data
            logger.debug('Beginning interpolation...')

            temperatureInterpolation = UnivariateSpline(HGHTTEMP, TEMP, s=self._interpolationPrecision)
            pressureInterpolation = UnivariateSpline(HGHTTEMP, PRESS, s=self._interpolationPrecision)
            windDirectionInterpolation = UnivariateSpline(HGHT, DRCT, s=self._interpolationPrecision)
            windSpeedInterpolation = UnivariateSpline(HGHTSKNT, SKNT, s=self._interpolationPrecision)

            def getTemperature(*args):
                if len(args) == 1:
                    return temperatureInterpolation(args[0])
                elif len(args) == 4:
                    return temperatureInterpolation(args[2])
                else:
                    return numpy.nan

            def getPressure(*args):
                if len(args) == 1:
                    return pressureInterpolation(args[0])
                elif len(args) == 4:
                    return pressureInterpolation(args[2])
                else:
                    return numpy.nan

            def getWindDirection(*args):
                if len(args) == 1:
                    return windDirectionInterpolation(args[0])
                elif len(args) == 4:
                    return windDirectionInterpolation(args[2])
                else:
                    return numpy.nan

            def getWindSpeed(*args):
                if len(args) == 1:
                    return windSpeedInterpolation(args[0])
                elif len(args) == 4:
                    return windSpeedInterpolation(args[2])
                else:
                    return numpy.nan

            # Store the interpolators in the correct variables
            self.getTemperature = getTemperature
            self.getPressure = getPressure
            self.getWindDirection = getWindDirection
            self.getWindSpeed = getWindSpeed

            logger.debug('Interpolation completed. Preparing derived functions...')

            # Initialize derived functions
            AirMolecMass = 0.02896
            GasConstant = 8.31447
            standardTempRankine = tools.c2kel(15) * (9. / 5)
            Mu0 = 0.01827 # Mu 0 (15 deg) [cP]
            C = 120 # Sutherland's Constant

            def getDensity(*args):
                if len(args) == 1:
                    return self.getPressure(args[0]) * 100 * AirMolecMass / (
                        GasConstant * tools.c2kel(self.getTemperature(args[0])))
                elif len(args) == 4:
                    return self.getPressure(args[0], args[1], args[2], args[3]) * 100 * AirMolecMass / (
                        GasConstant * tools.c2kel(self.getTemperature(args[0], args[1], args[2], args[3])))
                else:
                    return numpy.nan

            def getViscosity(*args):
                if len(args) == 1:
                    tempRankine = tools.c2kel(self.getTemperature(args[0])) * (9. / 5)
                    TTO = (tempRankine / standardTempRankine) ** 1.5 # T/TO [Rankine/Rankine]
                    TR = ((0.555 * standardTempRankine) + C) / ((0.555 * tempRankine) + C)
                    vcP = Mu0 * TTO * TR
                    return vcP / 1000.
                elif len(args) == 4:
                    tempRankine = tools.c2kel(self.getTemperature(args[0], args[1], args[2], args[3])) * (9. / 5)
                    TTO = (tempRankine / standardTempRankine) ** 1.5 # T/TO [Rankine/Rankine]
                    TR = ((0.555 * standardTempRankine) + C) / ((0.555 * tempRankine) + C)
                    vcP = Mu0 * TTO * TR
                    return vcP / 1000.
                else:
                    return numpy.nan

            self.getDensity = getDensity
            self.getViscosity = getViscosity

            logger.debug('All derived functions ready.')

            return True

        if self.UTC_offset == 0:
            self.UTC_offset = tools.getUTCOffset(self.launchSiteLat,self.launchSiteLon,self.dateAndTime)
            logger.debug('Fetched time zone data about the launch site: UTC offset is %f hours' % self.UTC_offset)

        self._UTC_time = self.dateAndTime - timedelta(seconds=self.UTC_offset * 3600)
        logger.debug('Using UTC time %s' % self._UTC_time.strftime('%d/%m/%y %H:%M'))

        if soundingFile[-3:] == 'ftr':
            # READ THE .FTR INPUT FILE AND PASS ITS DATA TO THE process_sounding_data(...) METHOD.

            # Try to open the file
            try:
                datastream = open(soundingFile, 'r')

                logger.debug('Sounding file successfully opened.')

            except IOError:
                logger.error('The sounding file you specified cannot be opened.')
                return
                # Read it and split it in lines
            filelines = datastream.readlines()
            datastream.close()

            # Check that data exists
            if not filelines:
                logger.error('The sounding file you specified is empty.')
                return

            logger.debug('File reading completed.')

            # Initialize variables
            PRESS = []
            HGHT = []
            TEMP = []
            DRCT = []
            SKNT = []

            # Import data into variables
            filelines = filelines[1:] # Remove header
            for line in filelines:
                lineEntries = line.split()
                if len(lineEntries) == 12:
                    PRESS.append(float(lineEntries[2]))
                    HGHT.append(float(lineEntries[1]))
                    TEMP.append(float(lineEntries[3]))
                    DRCT.append(float(lineEntries[7]))
                    SKNT.append(float(lineEntries[6]))

            logger.debug('Data imported.')

            if process_sounding_data(PRESS, HGHT, TEMP, DRCT, SKNT) is True:
                self._weatherLoaded = True
                logger.debug('Weather successfully loaded.')
            else:
                self._weatherLoaded = False
                return


        elif soundingFile[-8:] == 'sounding':
            # READ THE .SOUNDING INPUT FILE AND PASS ITS DATA TO THE process_sounding_data(...) METHOD.

            # Try to open the file
            try:
                datastream = open(soundingFile, 'r')

                logger.debug('Sounding file successfully opened.')

            except IOError:
                logger.error('The sounding file you specified cannot be opened.')
                return
                # Read it and split it in lines
            filelines = datastream.readlines()
            datastream.close()

            # Check that data exists
            if not filelines:
                logger.error('The sounding file you specified is empty.')
                return

            logger.debug('File reading completed.')

            # Initialize variables
            PRESS = []
            HGHT = []
            TEMP = []
            DRCT = []
            SKNT = []

            # Detect first line of useful data
            startLine = 0
            for i in xrange(len(filelines)):
                if filelines[i][0] == '-' and i < 25:
                    startLine = i + 1

            if startLine == 0:
                logger.error('Sounding file format does not follow format standards!')
                self._weatherLoaded = False

            # Import data into variables
            for line in filelines[startLine:]:
                try:
                    data = line.split()
                    if len(data) == 11:
                        PRESS.append(float(data[0]))
                        HGHT.append(float(data[1]))
                        TEMP.append(float(data[2]))
                        DRCT.append(float(data[6]))
                        SKNT.append(float(data[7]))
                except ValueError:
                    # End of useful data
                    break

            logger.debug('Data imported.')

            if process_sounding_data(PRESS, HGHT, TEMP, DRCT, SKNT):
                self._weatherLoaded = True
                logger.debug('Weather successfully loaded.')
            else:
                self._weatherLoaded = False
                return

        else:
            logger.error('Unknown sounding file format. File extension should be .ftr or .sounding')
            self._weatherLoaded = False


    def perturbWind(self, numberOfFlights):
        """
        Perturb the wind profiles for the purpose of Monte Carlo simulations.

        Given the numberOfFlights, this method generates N different perturbed wind profiles, where N = numberOfFlights,
        and stores them in the getMCWindDirection and getMCWindSpeed lists.
        The perturbation is based on the distance and time from sounding and is performed by picking a random experimentally-
        measured perturbation and applying it to each wind profile.

        Wind data should then be requested using the following syntax:

            getMCWindDirection[flightNumber].getTemperature(lat,lon,alt,time)

        where flightNumber should be in the range 0...numberOfFlights-1.

        """

        # Before anything, check that the weather is loaded.
        if not self._weatherLoaded:
            logger.error(
                'Weather not loaded! You need to load a sounding or download weather before perturbing the wind!')
            return

        self.getMCWindDirection = []
        self.getMCWindSpeed = []

        import wind_time_perturbation
        import wind_space_perturbation

        def make_perturbedWind(iDevTime, iDevSpace, randomChance, resultType=None):
            """
            Constructor function responsible of applying perturbations to wind profiles and returning closure functions.
            """

            # Interpolate the deviations
            timeDeviationU = UnivariateSpline(wind_time_perturbation.AltitudesM[:31],
                                              wind_time_perturbation.M_U_Kts[iDevTime], s=5)
            timeDeviationV = UnivariateSpline(wind_time_perturbation.AltitudesM[:31],
                                              wind_time_perturbation.M_V_Kts[iDevTime], s=5)
            spaceDeviationU = UnivariateSpline(wind_space_perturbation.AltitudesM[:31],
                                               wind_space_perturbation.M_U_Kts[iDevSpace], s=5)
            spaceDeviationV = UnivariateSpline(wind_space_perturbation.AltitudesM[:31],
                                               wind_space_perturbation.M_V_Kts[iDevSpace], s=5)

            def perturbedWind(*args):

                if len(args) == 1:
                    altitude = args[0]
                elif len(args) == 4:
                    altitude = args[2]
                else:
                    return numpy.nan

                # Convert non-perturbed wind direction and speed to u- and v-components
                UKts, VKts = tools.dirspeed2uv(self.getWindDirection(altitude), self.getWindSpeed(altitude))

                # Apply the time perturbation
                if randomChance[0] < 0.5:
                    UKts -= timeDeviationU(altitude) * self.timeFromSounding
                else:
                    UKts += timeDeviationU(altitude) * self.timeFromSounding

                if randomChance[1] < 0.5:
                    VKts -= timeDeviationV(altitude) * self.timeFromSounding
                else:
                    VKts += timeDeviationV(altitude) * self.timeFromSounding

                # Apply the space perturbation
                if randomChance[2] < 0.5:
                    UKts -= spaceDeviationU(altitude) * self.distanceFromSounding
                else:
                    UKts += spaceDeviationU(altitude) * self.distanceFromSounding

                if randomChance[3] < 0.5:
                    VKts -= spaceDeviationV(altitude) * self.distanceFromSounding
                else:
                    VKts += spaceDeviationV(altitude) * self.distanceFromSounding

                # Reconvert wind u- and v-components to direction and speed
                newDir, newSpd = tools.uv2dirspeed(UKts, VKts)

                # Deliver the results
                if resultType is None:
                    return newDir, newSpd
                elif resultType == 'direction':
                    return newDir
                elif resultType == 'speed':
                    return newSpd
                else:
                    return

            return perturbedWind

        # Apply perturbations to flights
        for _ in xrange(numberOfFlights):
            # Pick a random deviation out of the available ones.
            devTime = numpy.random.random_integers(0, 1758)
            devSpace = numpy.random.random_integers(0, 896)
            randChance = numpy.random.random(4)

            # Perturb and store
            self.getMCWindDirection.append(make_perturbedWind(devTime, devSpace, randChance, 'direction'))
            self.getMCWindSpeed.append(make_perturbedWind(devTime, devSpace, randChance, 'speed'))


class forecastEnvironment(environment):
    """
    The forecastEnvironment class is responsible of downloading weather forecast data from the GFS and generating a
    forecast-based atmospheric model.

    USAGE
    --------------

    Instantiate a new forecastEnvironment object if you wish to download weather information from the NOAA's Global
    Forecast System (more information can be found on http://nomads.ncep.noaa.gov/ )

    The new object can be instantiated and initialized using the following syntax:

    import weather
    my_forecast_atmosphere = weather.forecastEnvironment(debugging=False)

    Upon instantiation, two optional parameters can be passed: debugging and log_to_file. By default they are both FALSE.
    If log_to_file is set to TRUE, all error and debug logs will be stored in a file called error.log.
    If FALSE, all error and debug logs will be directly displayed on the terminal or command line.
    If debugging is TRUE, all the information available will be logged

    Once the object has been initialized, it needs to be configured with all the environmental data required.
    The following parameters need to be defined (there are no optional parameters):
        launchSiteLat: latitude of the launch site [deg] (type: float)
        launchSiteLon: longitude of the launch site [deg] (type: float)
        launchSiteElev: elevation of the launch site above Mean Sea Level [m] (type: float)
        dateAndTime: the launch time (type: datetime.datetime)
        UTC_offset: the offset in hours between the current time zone and UTC (for example, Florida in winter has a
            UTC_offset = -5) (type: float)
        forceNonHD: if TRUE, the weather forecast download will be forced to be not HD (i.e. 1deg x 1deg). Default FALSE

    Once all the parameters have been set, the forecast data can be downloaded from the GFS and the atmosphere modelled.

    The methods available to perform these operations are the following:
        loadForecast()      create a link to the Global Forecast System and download the required atmospheric data.
             Once the data has been downloaded, altitude data is generated and interpolated. The processed data is then
             prepared to be used with standard environment get... methods. This is enough for the environment to be
             ready for a deterministic flight simulation. Returns TRUE if successful.
        perturbWind()       perform a wind perturbation for the purpose of Monte Carlo simulations. Currently not
             implemented. Note: this method should only be called AFTER loadForecast() has been performed.

    The standard syntax to create the GFS link and download the data is the following:

        my_forecast_atmosphere.loadForecast()

    This downloads all the required atmospheric data and processes it for it to be ready to use.
    Refer to the GFS Module documentation for more information about the GFS link and the features and limitations of
    this system.
    loadForecast() returns TRUE is the download and processing was successful, FALSE otherwise.

    Refer to the error.log file for any errors.


    NOTE
    --------------
    forceNonHD limits the forecast download to non-HD only.
    This is NOT recommended, since the definition of the weather forecast is much better in HD.
    It's important to note that the simulator automatically switches to non-HD forecast if the amount of data requested
    is too high.


    EXAMPLE
    --------------

    A typical usage of the forecastEnvironment object would be as follows:

    from datetime import datetime, timedelta
    import weather

    my_forecast_atmosphere = weather.soundingEnvironment()
    my_forecast_atmosphere.launchSiteLat = 50.2245
    my_forecast_atmosphere.launchSiteLon = -5.3069
    my_forecast_atmosphere.launchSiteElev = 60
    my_forecast_atmosphere.dateAndTime = datetime.now() + timedelta(days=1)
    my_forecast_atmosphere.UTC_offset = 3


    my_forecast_atmosphere.loadForecast()

    """

    def __init__(self, debugging=False, log_to_file=False):
        """
        Initialize the forecastEnvironment object
        """

        # Run the environment class initialization first
        environment.__init__(self, debugging, log_to_file)

        # Initialize extra forecast-specific variables
        self.forceNonHD = False
        self.maxFlightTime = 18000
        self._GFSmodule = None

    def loadForecast(self, progressHandler=None):
        """
        Create a link to the Global Forecast System and download the required atmospheric data.
        Once the data has been downloaded, altitude data is generated and interpolated.
        The processed data is then prepared to be used with standard environment get... methods.

        Returns TRUE if successful.
        """

        # Data validation
        if self._weatherLoaded:
            logger.warning('The weather was already loaded. All data will be overwritten.')

        if self.launchSiteLat == 0.0:
            logger.debug('The launch site latitude is set to 0! Are you sure this is correct?')
        if self.launchSiteLon == 0.0:
            logger.debug('The launch site longitude is set to 0! Are you sure this is correct?')
        if self.dateAndTime is None:
            logger.error('The flight date and time has not been set and is required!')
            return

        if self.UTC_offset == 0:
            self.UTC_offset = tools.getUTCOffset(self.launchSiteLat,self.launchSiteLon,self.dateAndTime)
            logger.debug('Fetched time zone data about the launch site: UTC offset is %f hours' % self.UTC_offset)

        self._UTC_time = self.dateAndTime - timedelta(seconds=self.UTC_offset * 3600)
        logger.debug('Using UTC time %s' % self._UTC_time.strftime('%d/%m/%y %H:%M'))

        # Setup the GFS link
        self.handler = GFS.GFS_Handler(self.launchSiteLat, self.launchSiteLon, self._UTC_time,
                                       HD=self.forceNonHD is False,
                                       forecast_duration=int(self.maxFlightTime / 3600.), debugging=self.debugging,
                                       log_to_file=self.file_logging)
        self._GFSmodule = self.handler

        # Connect to the GFS and download data
        if self._GFSmodule.downloadForecast(progressHandler):
            logger.debug('GFS data successfully downloaded.')
        else:
            logger.error('Error while downloading GFS data.')
            return


        # Setup the standard environment data access functions


        # Linearly interpolate all data downloaded from the GFS
        pressureInterpolation, temperatureInterpolation, windDirectionInterpolation, windSpeedInterpolation = self._GFSmodule.interpolateData(
            'press', 'temp', 'windrct', 'windspd')

        self.getPressure = lambda lat, lon, alt, time: float(pressureInterpolation(lat, lon, alt,
                                                                                   self._GFSmodule.getGFStime(
                                                                                       time - timedelta(
                                                                                           seconds=self.UTC_offset * 3600))))
        self.getTemperature = lambda lat, lon, alt, time: float(temperatureInterpolation(lat, lon, alt,
                                                                                         self._GFSmodule.getGFStime(
                                                                                             time - timedelta(
                                                                                                 seconds=self.UTC_offset * 3600))))
        self.getWindDirection = lambda lat, lon, alt, time: float(windDirectionInterpolation(lat, lon, alt,
                                                                                             self._GFSmodule.getGFStime(
                                                                                                 time - timedelta(
                                                                                                     seconds=self.UTC_offset * 3600))))
        self.getWindSpeed = lambda lat, lon, alt, time: float(windSpeedInterpolation(lat, lon, alt,
                                                                                     self._GFSmodule.getGFStime(
                                                                                         time - timedelta(
                                                                                             seconds=self.UTC_offset * 3600))))

        # Extra definitions for derived quantities (density and viscosity)
        AirMolecMass = 0.02896
        GasConstant = 8.31447
        standardTempRankine = tools.c2kel(15) * (9. / 5)
        Mu0 = 0.01827 # Mu 0 (15 deg) [cP]
        C = 120 # Sutherland's Constant

        self.getDensity = lambda lat, lon, alt, time: self.getPressure(lat, lon, alt, time) * 100 * AirMolecMass / (
            GasConstant * tools.c2kel(self.getTemperature(lat, lon, alt, time)))

        def viscosity(lat, lon, alt, time):
            tempRankine = tools.c2kel(self.getTemperature(lat, lon, alt, time)) * (9. / 5)
            TTO = (tempRankine / standardTempRankine) ** 1.5 # T/TO [Rankine/Rankine]
            TR = ((0.555 * standardTempRankine) + C) / ((0.555 * tempRankine) + C)
            vcP = Mu0 * TTO * TR
            return vcP / 1000.

        self.getViscosity = viscosity

        self._weatherLoaded = True


    def perturbWind(self, numberOfFlights):
        """
        Perform a wind perturbation for the purpose of Monte Carlo simulations.

        Note: this method should only be called AFTER loadForecast() has been performed.
        """

        # For the time being, just return the standard un-perturbed forecast.
        # ### MODIFY HERE TO ADD WIND PERTURBATION TO FORECAST ###

        if not self._weatherLoaded:
            logger.error('Error: the weather has not been loaded yet! Wind cannot be perturbed.')
            return

        self.getMCWindDirection = []
        self.getMCWindSpeed = []

        def perturbedWindDirection(lat, lon, alt, time):
            return self.getWindDirection(lat, lon, alt, time)

        def perturbedWindSpeed(lat, lon, alt, time):
            return self.getWindSpeed(lat, lon, alt, time)

        for _ in numpy.arange(numberOfFlights):
            self.getMCWindDirection.append(perturbedWindDirection)
            self.getMCWindSpeed.append(perturbedWindSpeed)