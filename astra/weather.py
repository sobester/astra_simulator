# coding=utf-8

"""
This module contains classes for modelling the weather and atmosphere.

The environment classes providing the simulator module with all the atmospheric
data throughout a simulation. There are two environment
subclasses defined here:
    forecastEnvironment : GFS-based atmospheric model
    soundingEnvironment : atmospheric model generated from a sounding file

Refer to the individual classes and subclasses for details on how to use them.

University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""
from datetime import timedelta
import logging
from six.moves import range, builtins
import numpy
from scipy.interpolate import UnivariateSpline
from . import wind_time_perturbation
from. import wind_space_perturbation
from . import global_tools as tools
from . import GFS
from types import MethodType
import os


# SETUP ERROR LOGGING AND DEBUGGING
logger = logging.getLogger(__name__)

# Pass through the @profile decorator if line profiler (kernprof) is not in use
try:
    builtins.profile
except AttributeError:
    def profile(func):
        return func


class environment(object):
    """
    Defines a common interface for the Simulator module.

    This is a meta class that should not be instantiated directly, and is
    provided mainly for code reuse and a uniform API for the simulator class,
    regardless of the environment data used.
    
    Parameters
    ----------
    launchSiteLat : float
        latitude of the launch site [deg]
    launchSiteLon : float
        longitude of the launch site [deg]
    launchSiteElev : float
        elevation of the launch site above Mean Sea Level [m]
    dateAndTime : :obj:`datetime.datetime`
        Date and time of launch
    inflationTemperature : float
        the ambient temperature during the balloon inflation [degC]
    UTC_offset : float
        the offset in hours between the current time zone and UTC
        (for example, Florida in winter has a UTC_offset = -5)

    Attributes
    ----------
    launchSiteLat : float
        latitude of the launch site [deg]
    launchSiteLon : float
        longitude of the launch site [deg]
    launchSiteElev : float
        elevation of the launch site above Mean Sea Level [m]
    dateAndTime : :obj:`datetime.datetime`
        Date and time of launch
    UTC_offset : float
        The offset in hours between the current time zone and UTC. NOTE: If
        zero, UTC offset is AUTOMATICALLY retrieved using the launch site GPS
        location

    Notes
    -----
    The primary base class methods that should be overridden are the 'getter'
    functions for Temperature, Pressure, WindDirection, WindSpeed, Density,
    Viscosity

    See Also
    --------
    astra.global_tools.getUTCOffset
    """

    def __init__(self,
                 launchSiteLat,
                 launchSiteLon,
                 launchSiteElev,
                 dateAndTime,
                 inflationTemperature=0.0,
                 UTC_offset=0.0,
                 debugging=False,
                 load_on_init=False):

        # COMMON INTERFACE

        # Variables
        # Set all kwargs as attributes - could move this to Base class
        self.inflationTemperature = inflationTemperature
        self.launchSiteLat = launchSiteLat
        self.launchSiteLon = launchSiteLon
        self.launchSiteElev = launchSiteElev
        self.dateAndTime = dateAndTime
        self.UTC_offset = UTC_offset
        self.debugging = debugging

        self._UTC_time = None

        # Monte Carlo parameters
        self.getMCWindDirection = []
        self.getMCWindSpeed = []

        self._weatherLoaded = False

        if debugging:
            log_lev = logging.DEBUG
        else:
            log_lev = logging.WARNING

        logger.setLevel(log_lev)

    def getTemperature(self, lat, lon, alt, time):
        """Request the temperature for an input location and time.

        Returns
        -------
        temperature : float
            temperature in degrees Celsius
        """
        raise NotImplementedError(
            "getTemperature method must be implemented by class {}".format(
                type(self).__name__))

    def getPressure(self, lat, lon, alt, time):
        """request the pressure for the point at the given location at the
        given time.
        Returns a float, pressure in [millibar]"""
        raise NotImplementedError(
            "getPressure method must be implemented by class {}".format(
                type(self).__name__))

    def getDensity(self, lat, lon, alt, time):
        """request the density for the point at the given location at the
        given time.
        Returns a float, density in [kg/m3]"""
        raise NotImplementedError(
            "getDensity method must be implemented by class {}".format(
            type(self).__name__))

    def getViscosity(self, lat, lon, alt, time):
        """getViscosity(): request the viscosity for the point
        at the given location at the given time.
        Returns a float, viscosity in [Pa s]"""
        raise NotImplementedError(
            "getViscosity method must be implemented by class {}".format(
            type(self).__name__))

    def getWindSpeed(self, *args):
        """request the wind speed for the point at the given location at the
        given time.
        Returns a float, speed in [knots]"""
        raise NotImplementedError(
            "getWindSpeed method must be implemented by class {}".format(
                type(self).__name__))

    def getWindDirection(self, *args):
        """getWindDirection(lat,lon,alt,time): request the wind direction for
        the point at the given location at the given time.
        Returns a float, direction in [degrees] clockwise from north"""
        raise NotImplementedError(
            "getWindDirection method must be implemented by class {}".format(
                type(self).__name__))


class soundingEnvironment(environment):
    """
    Class for generating an atmospheric model from sounding data.

    Instantiate a new soundingEnvironment object if you are providing a
    sounding file to generate an atmospheric model. Accepted formats for
    sounding files are .ftr and .sounding .

    Parameters
    ----------
    launchSiteLat : float
        latitude of the launch site [deg]
    launchSiteLon : float
        longitude of the launch site [deg]
    launchSiteElev : float
        elevation of the launch site above Mean Sea Level [m]
    dateAndTime : :obj:`datetime.datetime`
        The launch time
    UTC_offset : float
        the offset in hours between the current time zone and UTC (for example,
        Florida in winter has a UTC_offset = -5)
    inflationTemperature : float
        the ambient temperature during the balloon inflation [degC]
    distanceFromSounding : float
        the distance between the launch site and the location where the 
        sounding file was recorded [km]
    timeFromSounding : float
        the duration between the balloon launch and the time when the sounding
        file was recorded [hours]
    [maxAltitude] : float (default 50999m)
        the upper altitude limit below which the atmosphere is modelled. [m]
    [debug] : bool, optional (default False)
        If TRUE, all the information available will be logged
    [log_to_file] : bool, optional (default False)
         If true, all error and debug logs will be stored in an error.log file

    :Example:

        >>> from datetime import datetime
        >>> import weather

        >>> my_sounding_atmosphere = weather.soundingEnvironment()
        >>> my_sounding_atmosphere.launchSiteLat = 50.2245
        >>> my_sounding_atmosphere.launchSiteLon = -5.3069
        my_sounding_atmosphere.launchSiteElev = 60
        my_sounding_atmosphere.dateAndTime = datetime.now()
        my_sounding_atmosphere.UTC_offset = 3
        my_sounding_atmosphere.inflationTemperature = 10.5
        my_sounding_atmosphere.distanceFromSounding = 36
        my_sounding_atmosphere.timeFromSounding = 2.8

        my_sounding_atmosphere.loadSounding('/path/to/mySoundingFile.ftr')

    """

    def __init__(self,
                 launchSiteLat,
                 launchSiteLon,
                 launchSiteElev,
                 dateAndTime,
                 soundingFile,
                 timeFromSounding,
                 distanceFromSounding,
                 inflationTemperature=0.0,
                 UTC_offset=0.,
                 debugging=False,
                 load_on_init=False):
        """Initialize the soundingEnvironment object.

        See class documentation.
        """
        # Initialize sounding-specific variables
        self.distanceFromSounding = distanceFromSounding
        self.timeFromSounding = timeFromSounding
        self.maxAltitude = 50999
        self.soundingFile = soundingFile

        self._interpolationPrecision = 200

        # Run the environment class initialization first
        super(soundingEnvironment, self).__init__(
            inflationTemperature=inflationTemperature,
            launchSiteLat=launchSiteLat,
            launchSiteLon=launchSiteLon,
            launchSiteElev=launchSiteElev,
            dateAndTime=dateAndTime,
            UTC_offset=UTC_offset,
            debugging=debugging,
            load_on_init=load_on_init) 

    def load(self, progressHandler=None):
        """Load the sounding file and generate the atmospheric model.

        Data is validated upon loading. Once this is done, cubic spline
        (scipy.interpolate.UnivariateSpline) interpolator functions are then
        stored as the getTemperature, getPressure, getWindDirection,
        getWindSpeed, getDensity and getViscosity methods.

        Parameters
        ----------
        soundingFile : string
            the relative or the absolute path to the sounding file.

        
        Notes
        -----
        * since sounding files only provide information for one location and at 
          one point in time, the lat, lon and time arguments for the get...
          methods are optional.
        """
        soundingFile = self.soundingFile

        # create a null handler if input progressHandler is None:
        if not progressHandler:
            def progressHandler(*args):
                return None

        if self.UTC_offset == 0:
            self.UTC_offset = tools.getUTCOffset(self.launchSiteLat,
                self.launchSiteLon,self.dateAndTime)
            logger.debug('Fetched time zone data about the launch site: UTC offset is %f hours' % self.UTC_offset)

        self._UTC_time = self.dateAndTime - timedelta(seconds=self.UTC_offset * 3600)
        logger.debug('Using UTC time %s' % self._UTC_time.strftime('%d/%m/%y %H:%M'))

        ext = os.path.splitext(soundingFile)[1]
        assert(ext.lower() in ['.ftr', '.sounding']),\
            'Unknown sounding file format: {}. Should be .ftr or .sounding'.format(ext)

        # Try to open the file
        try:
            datastream = open(soundingFile, 'r')
            logger.debug('Sounding file successfully opened.')

        except IOError:
            logger.error('The sounding file you specified cannot be opened.')
            raise
            # Read it and split it in lines
        filelines = datastream.readlines()
        datastream.close()
        # Check that data exists
        assert(filelines), 'The sounding file is empty.'

        # Initialize variables
        PRESS = []
        HGHT = []
        TEMP = []
        DRCT = []
        SKNT = []

        if ext == '.ftr':
            # READ THE .FTR INPUT FILE AND PASS ITS DATA TO THE
            # process_sounding_data(...) METHOD.

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
            progressHandler(0.5, 1)

            if self._process_sounding_data(PRESS, HGHT, TEMP, DRCT, SKNT):
                progressHandler(1.0, 1)
                self._weatherLoaded = True
                logger.debug('Weather successfully loaded.')
            else:
                self._weatherLoaded = False
                return


        elif ext == '.sounding':
            # READ THE .SOUNDING INPUT FILE AND PASS ITS DATA TO THE
            # process_sounding_data(...) METHOD.

            # Check that data exists
            if not filelines:
                logger.error('The sounding file you specified is empty.')
                return

            logger.debug('File reading completed.')

            # Detect first line of useful data
            startLine = 0
            for i in range(len(filelines)):
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

            if self._process_sounding_data(PRESS, HGHT, TEMP, DRCT, SKNT):
                self._weatherLoaded = True
                logger.debug('Weather successfully loaded.')
            else:
                self._weatherLoaded = False
                return

    def _process_sounding_data(self, PRESS, HGHT, TEMP, DRCT, SKNT):
        """
        This is an internal method responsible for processing raw data
        fetched from the sounding file.

        It requires pressure, altitude, temperature, wind speed and
        direction as inputs. All the inputs are cleaned up, prepared for
        interpolation, interpolated and then stored in the appropriate
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

        # _______________________________________________________________ #
        # Add missing data if launch site elevation or max altitude are out
        # of sounding bounds

        # If the elevation is lower than the lower bound of data, copy the
        # lowest data available and use it for the launch site elevation.
        if self.launchSiteElev < HGHT[0]:
            HGHT = numpy.insert(HGHT, 0, self.launchSiteElev)
            PRESS = numpy.insert(PRESS, 0, PRESS[0])
            TEMP = numpy.insert(TEMP, 0, TEMP[0])
            DRCT = numpy.insert(DRCT, 0, DRCT[0])
            SKNT = numpy.insert(SKNT, 0, SKNT[0])
            logger.debug('Launch site elevation out of bounds. Low altitude data generated.')

        # If the maxAltitude is higher than the upper bound of data, fill
        # the missing altitude levels in with ISA data
        if self.maxAltitude > HGHT[-1]:
            newHeights = numpy.arange(HGHT[-1] + 1, self.maxAltitude + 1,
                (self.maxAltitude - HGHT[-1] - 1) / 20.)
            HGHTTEMP = numpy.append(HGHT, newHeights[3:])
            for newTempHeight in newHeights[3:]:
                # Calculate new temperatures and pressures.
                _, newTemp, _, newPress, _ = tools.ISAatmosphere(
                    altitude=tools.m2feet(newTempHeight))
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

        # _______________________________________________________________ #

        # Check whether all fields have data, otherwise fill up vectors with NaN
        if numpy.size(HGHT) == 0 or numpy.size(TEMP) == 0 or\
            numpy.size(PRESS) == 0 or numpy.size(DRCT) == 0 or\
            numpy.size(SKNT) == 0:
            logger.error('There was a problem while processing the sounding.')
            return False

        # Interpolate data
        logger.debug('Beginning interpolation...')

        # TODO: Fix this part, it's unreadable and difficult to debug
        temperatureInterpolation = UnivariateSpline(HGHTTEMP, TEMP,
            s=self._interpolationPrecision)
        pressureInterpolation = UnivariateSpline(HGHTTEMP, PRESS,
            s=self._interpolationPrecision)
        windDirectionInterpolation = UnivariateSpline(HGHT, DRCT,
            s=self._interpolationPrecision)
        windSpeedInterpolation = UnivariateSpline(HGHTSKNT, SKNT,
            s=self._interpolationPrecision)

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
                return self.getPressure(args[0], args[1], args[2], args[3])\
                    * 100 * AirMolecMass / (GasConstant * tools.c2kel(
                        self.getTemperature(args[0],
                                            args[1],
                                            args[2],
                                            args[3])))
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
                tempRankine = tools.c2kel(self.getTemperature(args[0],args[1], args[2], args[3])) * (9. / 5)
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

    def make_perturbedWind(self, iDevTime, iDevSpace, randomChance, resultType=None):
            """
            Constructor function responsible for applying perturbations to wind
            profiles and returning closure functions.
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
                UKts, VKts = tools.dirspeed2uv(self.getWindDirection(altitude),
                    self.getWindSpeed(altitude))

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

    def perturbWind(self, numberOfFlights):
        """
        Perturb the wind profiles for the purpose of Monte Carlo simulations.

        Given the numberOfFlights, this method generates N different perturbed
        wind profiles, where N = numberOfFlights, and stores them in the
        getMCWindDirection and getMCWindSpeed lists. The perturbation is based
        on the distance and time from sounding and is performed by picking a
        random experimentally-measured perturbation and applying it to each
        wind profile.

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

        # Apply perturbations to flights
        for _ in range(numberOfFlights):
            # Pick a random deviation out of the available ones.
            devTime = numpy.random.random_integers(0, 1758)
            devSpace = numpy.random.random_integers(0, 896)
            randChance = numpy.random.random(4)

            # Perturb and store
            self.getMCWindDirection.append(self.make_perturbedWind(devTime,
                devSpace, randChance, 'direction'))
            self.getMCWindSpeed.append(self.make_perturbedWind(devTime, devSpace,
                randChance, 'speed'))

        return None


class forecastEnvironment(environment):
    """
    Class responsible for downloading weather forecast data from the Global
    Forecast System (GFS) and generating a forecast-based atmospheric model.

    Parameters
    ----------
    launchSiteLat : float
        latitude of the launch site [deg]
    launchSiteLon : float
        longitude of the launch site [deg]
    launchSiteElev : float
        elevation of the launch site above Mean Sea Level [m]
    dateAndTime : :obj:`datetime.datetime`
        The launch time
    UTC_offset : float
        the offset in hours between the current time zone and UTC (for example,
        Florida in winter has a UTC_offset = -5)
    inflationTemperature : float
        the ambient temperature during the balloon inflation [degC]
    [forceNonHD] : bool (default False)
        if TRUE, the weather forecast download will be forced to a lower
        resolution (i.e. 1deg x 1deg)
    [forecastDuration] : float (default 4)
        The number of hours from dateAndTime for which to download weather data
    [use_async] : bool (default True)
        Use an asynchronous request for downloads. This should speed up the
        download, but may incur larger memory overhead for large forecastDuration.
    [requestSimultaneous] : bool (default True)
        If True, populate a dictionary of responses from the web download
        requests, then process the data. If False, each response will be
        removed once the data has been processed: This has better memory
        overhead for large ForecastDuration, but is slower, and does not work
        with asynchronous requesting (use_async)
    [debug] : bool, optional (default False)
        If TRUE, all the information available will be logged
    [log_to_file]: bool, optional (default False)
         If true, all error and debug logs will be stored in an error.log file
    [progressHandler] : function, or None (default None)
        Progress for each downloaded parameter (in %) will be passed to this
        function, if provided.
    [load_on_init] : bool, optional (default False)
        If True, the forecast will be downloaded when the environment object is
        created. This is set to False by default, as the :obj:`flight` object
        should preferably be used to load the forecast (input validation and
        preflight checks should be done before expensive data download) 

    See Also
    --------
    astra.GFS

    Notes
    -----
    * The simulator automatically switches to non-HD forecast if the amount of
    data requested is too high.

    :Example:
    
        >>> from datetime import datetime, timedelta
        >>> import weather

        >>> my_forecast_atmosphere = weather.forecastEnvironment()
        >>> my_forecast_atmosphere.launchSiteLat = 50.2245
        >>> my_forecast_atmosphere.launchSiteLon = -5.3069
        >>> my_forecast_atmosphere.launchSiteElev = 60
        >>> my_forecast_atmosphere.dateAndTime = datetime.now() + timedelta(days=1)
        >>> my_forecast_atmosphere.UTC_offset = 3
        >>> my_forecast_atmosphere.loadForecast()

    """

    def __init__(self,
                 launchSiteLat,
                 launchSiteLon,
                 launchSiteElev,
                 dateAndTime,
                 UTC_offset=0,
                 inflationTemperature=0.0,
                 forceNonHD=False,
                 forecastDuration=4,
                 use_async=True,
                 requestSimultaneous=True,
                 debugging=False,
                 progressHandler=None,
                 load_on_init=False):
        """
        Initialize the forecastEnvironment object
        """
        # Initialize extra forecast-specific variables
        self.forceNonHD = forceNonHD
        self.forecastDuration = forecastDuration
        self.use_async = use_async
        self.requestSimultaneous = requestSimultaneous

        self._GFSmodule = None

        # This should be the last thing that is called on init, since the base
        # (environment) class init calls self.load (if load_on_init is True)
        super(forecastEnvironment, self).__init__(
            inflationTemperature=inflationTemperature,
            launchSiteLat=launchSiteLat,
            launchSiteLon=launchSiteLon,
            launchSiteElev=launchSiteElev,
            dateAndTime=dateAndTime,
            UTC_offset=UTC_offset,
            debugging=debugging,
            load_on_init=load_on_init)

    @profile
    def load(self, progressHandler=None):
        """
        Create a link to the Global Forecast System and download the required
        atmospheric data.

        Once the data has been downloaded, altitude data is generated and
        interpolated. The processed data is then prepared to be used with
        standard environment get... methods.

        Parameters
        ----------
        progressHandler : function
            Returns progress 

        Returns
        -------
        status : TRUE if successful.
        """

        # Data validation
        if self._weatherLoaded:
            logger.warning(
                'The weather was already loaded. All data will be overwritten.')

        if self.launchSiteLat == 0.0:
            logger.debug(
                'The launch site latitude is set to 0!')
        if self.launchSiteLon == 0.0:
            logger.debug(
                'The launch site longitude is set to 0!')
        if self.dateAndTime is None:
            raise ValueError(
                'The flight date and time has not been set and is required!')

        if self.UTC_offset == 0:
            self.UTC_offset = tools.getUTCOffset(
                self.launchSiteLat,self.launchSiteLon,self.dateAndTime)
            logger.debug('Fetched time zone data about the launch site: UTC offset is %f hours' % self.UTC_offset)

        self._UTC_time = self.dateAndTime - timedelta(seconds=self.UTC_offset * 3600)
        logger.debug('Using UTC time %s' % self._UTC_time.strftime('%d/%m/%y %H:%M'))

        # log the current parameters
        logger.info('Preparing to download weather data for parameters:')
        logger.debug("    Launch site Latitude: {}".format(self.launchSiteLat))
        logger.debug("    Launch site Longitude: {}".format(self.launchSiteLon))
        logger.debug("    Launch time: {}".format(self._UTC_time))

        # Setup the GFS link
        self._GFSmodule = GFS.GFS_Handler(self.launchSiteLat,
                                          self.launchSiteLon,
                                          self._UTC_time,
                                          use_async=self.use_async,
                                          requestSimultaneous=self.requestSimultaneous,
                                          HD=(not self.forceNonHD),
                                          forecastDuration=
                                            self.forecastDuration,
                                          debugging=self.debugging)

        # Connect to the GFS and download data
        if self._GFSmodule.downloadForecast(progressHandler):
            logger.debug('GFS data successfully downloaded.')
        else:
            logger.error('Error while downloading GFS data.')
            return


        # Setup the standard environment data access functions


        # Linearly interpolate all data downloaded from the GFS
        pressureInterpolation, temperatureInterpolation,\
            windDirectionInterpolation, windSpeedInterpolation = \
                self._GFSmodule.interpolateData('press',
                                                'temp',
                                                'windrct',
                                                'windspd')

        self.getPressure = lambda lat, lon, alt, time: float(
            pressureInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600)))
            )
        self.getTemperature = lambda lat, lon, alt, time: float(
            temperatureInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600)))
            )
        self.getWindDirection = lambda lat, lon, alt, time: float(
            windDirectionInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600))))
        self.getWindSpeed = lambda lat, lon, alt, time: float(
            windSpeedInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600))))

        # Extra definitions for derived quantities (density and viscosity)
        AirMolecMass = 0.02896
        GasConstant = 8.31447
        standardTempRankine = tools.c2kel(15) * (9. / 5)
        Mu0 = 0.01827 # Mu 0 (15 deg) [cP]
        C = 120 # Sutherland's Constant

        self.getDensity = lambda lat, lon, alt, time: \
            self.getPressure(lat, lon, alt, time) * 100 * AirMolecMass / (GasConstant * 
                tools.c2kel(self.getTemperature(lat, lon, alt, time))
            )

        def viscosity(lat, lon, alt, time):
            tempRankine = tools.c2kel(self.getTemperature(lat, lon, alt, time)) * (9. / 5)
            TTO = (tempRankine / standardTempRankine) ** 1.5 # T/TO [Rankine/Rankine]
            TR = ((0.555 * standardTempRankine) + C) / ((0.555 * tempRankine) + C)
            vcP = Mu0 * TTO * TR
            return vcP / 1000.

        self.getViscosity = viscosity

        self._weatherLoaded = True

    def loadFromNOAAFiles(self, fileDict):
        """Creates and stores the interpolator functions using data provided
        in the file names described in input fileDict.

        Parameters
        ----------
        fileDict : :obj:`dict`
            A dictionary of noaa_name: filename pairs, indicating which file
            should be used for each noaa variable. The following keys must be
            defined in this dictionary: 'tmpprs' (Temperature), 
            'hgtprs' (Altitude), 'ugrdprs' (U Winds), 'vgrdprs' (V Winds)
        """

        # Data validation
        if self._weatherLoaded:
            logger.warning(
                'The weather was already loaded. All data will be overwritten.')

        if self.UTC_offset == 0:
            self.UTC_offset = tools.getUTCOffset(
                self.launchSiteLat,self.launchSiteLon,self.dateAndTime)
            logger.debug('Fetched time zone data about the launch site: UTC offset is %f hours' % self.UTC_offset)

        self._UTC_time = self.dateAndTime - timedelta(seconds=self.UTC_offset * 3600)
        logger.debug('Using UTC time %s' % self._UTC_time.strftime('%d/%m/%y %H:%M'))

        # log the current parameters
        logger.info('Preparing to download weather data for parameters:')
        logger.debug("    Launch site Latitude: {}".format(self.launchSiteLat))
        logger.debug("    Launch site Longitude: {}".format(self.launchSiteLon))
        logger.debug("    Launch time: {}".format(self._UTC_time))

        # LOAD THE DATA: curently HD is not allowed, since only the named noaa
        # parameters are loaded from a file (different files are required for
        # high altitude data in the case of an HD simulation).
        self._GFSmodule = GFS.GFS_Handler.fromFiles(fileDict,
            lat=self.launchSiteLat, lon=self.launchSiteLon,
            date_time=self._UTC_time, HD=False,
            forecastDuration=self.forecastDuration,
            debugging=self.debugging)

        # Setup the standard environment data access functions
        # Linearly interpolate all data downloaded from the GFS
        pressureInterpolation, temperatureInterpolation,\
            windDirectionInterpolation, windSpeedInterpolation = \
                self._GFSmodule.interpolateData('press',
                                                'temp',
                                                'windrct',
                                                'windspd')

        self.getPressure = lambda lat, lon, alt, time: float(
            pressureInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600)))
            )
        self.getTemperature = lambda lat, lon, alt, time: float(
            temperatureInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600)))
            )
        self.getWindDirection = lambda lat, lon, alt, time: float(
            windDirectionInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600))))
        self.getWindSpeed = lambda lat, lon, alt, time: float(
            windSpeedInterpolation(lat, lon, alt, self._GFSmodule.getGFStime(
                time - timedelta(seconds=self.UTC_offset * 3600))))

        # Extra definitions for derived quantities (density and viscosity)
        AirMolecMass = 0.02896
        GasConstant = 8.31447
        standardTempRankine = tools.c2kel(15) * (9. / 5)
        Mu0 = 0.01827 # Mu 0 (15 deg) [cP]
        C = 120 # Sutherland's Constant

        self.getDensity = lambda lat, lon, alt, time: \
            self.getPressure(lat, lon, alt, time) * 100 * AirMolecMass / (GasConstant * 
                tools.c2kel(self.getTemperature(lat, lon, alt, time))
            )

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

        Note: this method should only be called AFTER loadForecast() has been
        performed.
        """
        # For the time being, just return the standard un-perturbed forecast.
        # TODO: MODIFY HERE TO ADD WIND PERTURBATION TO FORECAST
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