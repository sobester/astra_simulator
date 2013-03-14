# coding=utf-8

"""
GFS.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

GFS Module
Connect to the NOAA's Global Forecast System and download weather forecast data for any location on the planet for
up to 7 days from the last cycle.
The data is downloaded from the HD service, providing data on a  0.5 deg latitude x 0.5 deg longitude x 47 altitude
levels grid.
The NOAA issues new datasets (cycles) at 00:00, 06:00, 12:00 and 18:00 UTC every day and the 7 days of forecast data
refer to the cycle issuing time.

Warning: because it's not possible to know in advance what the latest cycle available is (there is a delay between the
cycle issuing time and its availability), it's recommended to limit the time span of data downloads from the GFS to 6
days.


USAGE
--------------

Instantiate a new GFS_Handler object to deal with data download from the GFS and further processing.
When a new GFS_Handler object is instantiated, the following parameters are REQUIRED for initialization:
    lat                 latitude of the location where forecast is required. Type: float
    lon                 longitude of the location where forecast is required. Type: float
    date_time           UTC time of the required forecast (min 30 days from today, max 7 days from latest cycle issuing
                        time, see warning above). Type: datetime.datetime
    forecast_duration   the duration in hours of the forecast data required. By default this is 4 hours. Type: float

Additionally, the following three optional parameters can be specified as well:
    HD                  if FALSE, non-HD forecast will be downloaded. Otherwise, HD forecast will be downloaded and, if
                        the data size is too high, the GFS_Handler will automatically switch to non-HD. Type: bool
    debugging, log_to_file
        If log_to_file is set to TRUE, all error and debug logs will be stored in a file called error.log.
        If FALSE, all error and debug logs will be directly displayed on the terminal or command line.
        If debugging is TRUE, all the information available will be logged. Otherwise, only errors will be logged.

The standard syntax to instantiate and initialize a GFS_Handler object is:
    myGFSHandler = GFS_Handler(lat,lon,dateAndTime,forecast_duration)

Once the GFS_Handler object has been initialized, all the following methods can be used:

    downloadForecast()  connects to the Global Forecast System and downloads the closest cycle available to the date_time
                        required. The cycle date and time is stored in the object's cycleDateTime variable.
                        Returns TRUE if the download was successful, FALSE otherwise.

    interpolateData(*variables)  sets up a linear 4d interpolation for each variable given and returns it. The returned
                        interpolator can then be used by calling it with standard 4D coordinates (lat, lon, alt, time).
                        Acceptable values for *variables can be a combination of any of these (separated by commas):
                            't'|'temp'|'temperature'        temperature interpolation
                            'p'|'press'|'pressure'          pressure interpolation
                            'h'|'alt'|'altitude'            altitude interpolation
                            'd'|'windrct'|'wind_direction'  wind direction interpolation
                            's'|'windspd'|'wind_speed'      wind speed interpolation

                        The linear interpolation is based on the interpolate.Linear4DInterpolator class.
                        Warning: the time parameter required by the Linear4DInterpolator objects must be in GFS units!
                        Use the getGFStime(time) method for conversion.
                        Returns a list of Linear4DInterpolator objects, in the same order as *variables.

    getGFStime(time)    converts standard datetime.datetime objects to GFS time units. The time parameter should be
                        a datetime.datetime object.
                        Returns a float corresponding to the converted time, or nan if there was no downloaded data found.

Note: the GFS_Map class is only used for internal data management only and should not be used on its own.


EXAMPLE
--------------

Typical usage of the GFS module:
    # Download forecast data for tomorrow at the University of Southampton and process it

    import datetime

    myGFSlink = GFS_Handler(50.93543,-1.39619,datetime.now()+datetime.timedelta(days=1))
    myGFSlink.downloadForecast()

    getTemp,getPress = myGFSlink.interpolateData('t','p')
    # Request temperature at lat = 51.2, lon = 0.46, alt = 32000m, time = 1 hour from now
    getTemp(51.2,0.46,32000,myGFSlink.getGFStime(datetime.now()+datetime.timedelta(days=1,seconds=3600)))


University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 12/02/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from datetime import datetime, timedelta
from math import floor, ceil
import urllib2
import logging

import numpy
from scipy.interpolate import UnivariateSpline

import global_tools as tools
from interpolate import Linear4DInterpolator


# Error and warning logger
logger = None

earthRadius = 6371009 # m


class GFS_Handler:
    """
    Core of the GFS module.
    Use this to perform any operation related to the Global Forecast System.
    See module's documentation for more information.
    """

    def __init__(self, lat, lon, date_time, HD=True, forecast_duration=4, debugging=False, log_to_file=False):
        global logger

        # Initialize Parameters
        self.launchDateTime = date_time
        self.lat = lat
        self.lon = lon
        self.maxFlightTime = forecast_duration
        self.HD = HD
        self.cycleDateTime = None
        self.firstAvailableTime = None

        # These are the 4D data matrices with all the information needed.
        self.altitudeData = None
        self.temperatureData = None
        self.windDirData = None
        self.windSpeedData = None

        # These are variables storing the mapping criteria of the 4D matrices.
        # They are arrays of dictionaries, one per dimension (i.e. axis) of the matrix.
        # The dictionary keys are the latitude, longitude, pressure and time, while the values are the
        # corresponding matrix data indices.
        self.altitudeMap = None
        self.temperatureMap = None
        self.windsMap = None

        # Grid size setup
        # NOTE: Grid sizes are defined as the difference between the highest and the lowest lat/lon requested, NOT the
        # difference between the highest and the requested point or the lowest and the requested point!
        if forecast_duration == 5:
            # Standard flight
            self.latGridSize = 6
            self.lonGridSize = 12
        else:
            # Non-standard flight (either floating or customized from command line)
            self.latGridSize = 2 * ceil(0.1 * self.maxFlightTime + 0.6) + 3
            self.lonGridSize = 2 * ceil(0.9 * self.maxFlightTime) + 3

        # Automatically switch to non HD if close to the poles, as all the longitudes will be downloaded
        if lat < -80 or lat > 80:
            self.HD = False

        # Force non HD weather if requested data is too high
        if self.latGridSize * self.lonGridSize * ceil(self.maxFlightTime / 3.) > 250:
            self.HD = False

        # HD/SD forecast setup
        if self.HD:
            self.latStep = 0.5
            self.lonStep = 0.5
        else:
            self.latStep = 1.0
            self.lonStep = 1.0
            # Prepare download of high altitude HD data
            self._highAltitudeGFS = GFS_High_Altitude_Handler(lat, lon, date_time, forecast_duration, debugging,
                                                              log_to_file)
            self._highAltitudePressure = None

        # SETUP ERROR LOGGING AND DEBUGGING

        logger = logging.getLogger('GFS')

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


    def downloadForecast(self, progressHandler=None):
        """
        Connect to the Global Forecast System and download the closest cycle available to the date_time
        required. The cycle date and time is stored in the object's cycleDateTime variable.
        Returns TRUE if the download was successful, FALSE otherwise.
        """

        #############################################################################################################
        # INITIALIZE DATA, OR RESET IT IF THEY ALREADY EXISTS
        temperatureMatrix = None
        geopotentialMatrix = None
        uWindsMatrix = None
        vWindsMatrix = None
        temperatureMap = None
        geopotentialMap = None
        uWindsMap = None


        #############################################################################################################
        # EXTRA PARAMETERS
        # The base URL depends on whether the HD service has been requested or not.
        weatherParameters = ['tmpprs', 'hgtprs', 'ugrdprs', 'vgrdprs']
        baseURL = {
            True: 'http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/',
            False: 'http://nomads.ncep.noaa.gov:9090/dods/gfs/'
        }[self.HD]


        #############################################################################################################
        # INITIALIZE TIME HANDLERS AND CALCULATE MOST RECENT CYCLE POSSIBLY AVAILABLE

        simulationDateTime = self.launchDateTime
        currentDateTime = datetime.now()

        if simulationDateTime < currentDateTime:
            # Simulating a flight in the past. In this case, use the simulation date time to calculate the cycles
            # that need to be downloaded.
            currentDateTime = simulationDateTime

        # These are the hours at which new forecast cycles are issued
        dailyCycles = [0, 6, 12, 18]
        # Determine which cycle issuing hour is the closest to the current one
        cycleTime = dailyCycles[numpy.digitize([currentDateTime.hour], dailyCycles)[0] - 1]
        latestCycleDateTime = datetime(currentDateTime.year, currentDateTime.month, currentDateTime.day, cycleTime)


        #############################################################################################################
        # INITIALIZE COMMON REQUEST PARAMETERS. These are lat, lon and alt and do not vary with each cycle.
        # Perform some validation checks to see if they need to be limited (close to the poles) or two requests need to
        # be made (close to the Greenwich meridian).
        requestAllLongitudes = False
        multipleRequests = False
        # This is the grid size converted to GFS units. It basically returns the same number if the gridSize is even,
        # or it returns gridSize-1 if gridSize is odd. This is to get the right grid around the launch point.
        latGridStep = (int(self.latGridSize) / 2) / self.latStep
        lonGridStep = (int(self.lonGridSize) / 2) / self.lonStep

        # ALTITUDE (Download ALL altitude levels available)
        requestAltitude = {
            True: [0, 46],
            False: [0, 25]
        }[self.HD]

        # LATITUDE
        targetLatitude = (round(self.lat) + 90) / self.latStep
        requestLatitude = [targetLatitude - latGridStep, targetLatitude + latGridStep]

        if self.HD:
            # Limit latitude to +/- 90 degrees
            if requestLatitude[0] < 0:
                requestLatitude[0] = 0
            if requestLatitude[1] > 360:
                requestLatitude[1] = 360

            # Request all longitudes if close to the poles
            if requestLatitude[0] < 20 or requestLatitude[1] > 340:
                requestAllLongitudes = True
        else:
        # Limit latitude to +/- 90 degrees
            if requestLatitude[0] < 0:
                requestLatitude[0] = 0
            if requestLatitude[1] > 180:
                requestLatitude[1] = 180

            # Request all longitudes if close to the poles or if simulation is beyond 2 days
            if requestLatitude[0] < 10 or requestLatitude[1] > 170 or self.maxFlightTime > 48:
                requestAllLongitudes = True

        # LONGITUDE
        # Note: There is no need to check if the grid size is higher than the whole world, causing overlapping regions
        # to only request data in the overlapped region, since a grid size approximately equal to half of the world is
        # enough to hit the poles and therefore request worldwide data (longitudes are very close at the poles, hence
        # it's worth keeping data for all of them)

        if requestAllLongitudes:
            requestLongitude = {
                True: [0, 719],
                False: [0, 359]
            }[self.HD]
        else:
            if self.lon >= 0:
                targetLongitude = round(self.lon) / self.lonStep
            else:
                targetLongitude = (360 - abs(round(self.lon))) / self.lonStep

            requestLongitude = [targetLongitude - lonGridStep, targetLongitude + lonGridStep]

            # Check if the values are within the bounds and correct if needed
            if self.HD:
                if requestLongitude[0] < 0:
                    requestLongitude[0] += 720
                if requestLongitude[1] > 719:
                    requestLongitude[1] -= 720
            else:
                if requestLongitude[0] < 0:
                    requestLongitude[0] += 360
                if requestLongitude[1] > 359:
                    requestLongitude[1] -= 360

            # Check if crossing the Greenwich meridian and split the requests
            # If the Greenwich meridian is being crossed, the left bound of the longitude interval will have a
            # higher value than the right one (in HD: Western hemisphere [360:719], Eastern hemisphere:[0:359]), so if
            # the difference between the right and the left one is negative, the Greenwich meridian is being crossed
            if requestLongitude[1] - requestLongitude[0] < 0:
                # SPLIT
                requestLongitude = {
                    True: [[0, requestLongitude[1]], [requestLongitude[0], 719]],
                    False: [[0, requestLongitude[1]], [requestLongitude[0], 359]]
                }[self.HD]
                multipleRequests = True


        #############################################################################################################
        # TRY TO DOWNLOAD DATA WITH THE LATEST CYCLE. IF NOT AVAILABLE, TRY WITH AN EARLIER ONE

        # pastCycle is a variable that indicates how many cycles in the past we're downloading data from.
        # It first tries with 0, indicating the most recent cycle: this is calculated, knowing that cycles are issued
        # every six hours. Sometimes, however, it takes a few hours for data to become available, or a cycle is skipped.
        # This means that data for the latest cycle is not guaranteed to be available.
        # If data is not available, it tries with 1 cycle older, until one is found with data available. If no cycles
        # are found, the method returns FALSE.
        pastCycle = 0
        while True:
            thisCycle = latestCycleDateTime - timedelta(hours=pastCycle * 6)
            self.cycleDateTime = thisCycle

            # This is just a flag to make sure that if no data is found, it stops the whole download and re-starts
            # the loop with an older cycle.
            thisCycleNotAvailable = False

            # Initialize time parameter
            timeFromForecast = simulationDateTime - thisCycle
            hoursFromForecast = timeFromForecast.total_seconds() / 3600.

            # GFS time index for the first dataset to be requested
            requestTime = floor(hoursFromForecast / 3.)

            # This stores the actual time of the first dataset downloaded. It's going to be used to convert real time
            # to GFS "time coordinates" (see getGFStime(time) function)
            self.firstAvailableTime = self.cycleDateTime + timedelta(hours=requestTime * 3)

            # Always download an extra time dataset
            requestTime = [requestTime, requestTime + ceil(self.maxFlightTime / 3.) + 1]

            #########################################################################################################
            # DOWNLOAD DATA FOR ALL PARAMETERS

            for requestVariable in weatherParameters:

                if multipleRequests:
                    numberOfRequests = 2
                else:
                    numberOfRequests = 1

                dataResults = []

                # Check if we need more than 1 request (ie if we are crossing the Greenwich meridian)
                for req_number in xrange(numberOfRequests):

                    if multipleRequests:
                        thisRequestLongitude = requestLongitude[req_number]
                    else:
                        thisRequestLongitude = requestLongitude

                    requestURL = '%sgfs%s%d%02d%02d/gfs%s_%02dz.ascii?%s[%d:%d][%d:%d][%d:%d][%d:%d]' % (
                        baseURL,
                        {True: '_hd', False: ''}[self.HD],
                        thisCycle.year,
                        thisCycle.month,
                        thisCycle.day,
                        {True: '_hd', False: ''}[self.HD],
                        thisCycle.hour,
                        requestVariable,
                        requestTime[0], requestTime[1],
                        requestAltitude[0], requestAltitude[1],
                        requestLatitude[0], requestLatitude[1],
                        thisRequestLongitude[0], thisRequestLongitude[1]
                    )

                    logger.debug('Requesting URL: %s' % requestURL)

                    try:
                        HTTPresponse = urllib2.urlopen(requestURL)
                    except:
                        logger.error('Error while connecting to the GFS server.')
                        logger.error('URL: %s' % requestURL)
                        return False
                    dataResults.append(HTTPresponse.read())

                    if dataResults[-1][0] == "<":
                        # Data from this cycle is not available. Try with the earlier one.
                        logger.debug('Cycle not available.')
                        thisCycleNotAvailable = True
                        break

                if thisCycleNotAvailable:
                    break

                if requestVariable == 'tmpprs':
                    temperatureMatrix, temperatureMap = self._generate_matrix(dataResults)
                    logger.debug('Temperature data downloaded and processed.')
                    if progressHandler is not None:
                        progressHandler(0.25, 1)
                elif requestVariable == 'hgtprs':
                    geopotentialMatrix, geopotentialMap = self._generate_matrix(dataResults)
                    logger.debug('Altitude data downloaded and processed.')
                    if progressHandler is not None:
                        progressHandler(0.5, 1)
                elif requestVariable == 'ugrdprs':
                    uWindsMatrix, uWindsMap = self._generate_matrix(dataResults)
                    logger.debug('U Winds data downloaded and processed.')
                    if progressHandler is not None:
                        progressHandler(0.75, 1)
                elif requestVariable == 'vgrdprs':
                    logger.debug('V Winds data downloaded and processed.')
                    vWindsMatrix, vWindsMap = self._generate_matrix(dataResults)
                    if progressHandler is not None:
                        progressHandler(0.95, 1)




            # Restart the loop if the data wasn't available.
            if thisCycleNotAvailable:
                if pastCycle == 24:
                    logger.error('No GFS cycles available found!')
                    return False
                else:
                    pastCycle += 1
                    logger.debug('Moving to next cycle...')
                    continue

            # If it didn't break in the previous if statement, data was available and has been downloaded.
            # End of the loop.
            break

        #############################################################################################################
        # PROCESS DATA AND PERFORM CONVERSIONS AS REQUIRED

        # Convert temperatures from Kelvin to Celsius
        temperatureMatrix -= 273.15

        # Convert geopotential height to geometric altitude
        altitudeMatrix = geopotentialMatrix * earthRadius / ( earthRadius - geopotentialMatrix )

        # Convert u and v winds to wind direction and wind speed matrices

        # Store the current shape of the 4D matrices
        matrixShape = uWindsMatrix.shape
        # Convert to KNOTS and the turn into direction and speed
        dirspeedWinds = map(tools.uv2dirspeed, (uWindsMatrix * 1.9438445).ravel(), (vWindsMatrix * 1.9438445).ravel())
        # Extract results
        windDirectionMatrix = numpy.array([dirspeed[0] for dirspeed in dirspeedWinds]).reshape(matrixShape)
        windSpeedMatrix = numpy.array([dirspeed[1] for dirspeed in dirspeedWinds]).reshape(matrixShape)



        # Store results
        self.temperatureData = temperatureMatrix
        self.altitudeData = altitudeMatrix
        self.windDirData = windDirectionMatrix
        self.windSpeedData = windSpeedMatrix

        self.temperatureMap = temperatureMap
        self.altitudeMap = geopotentialMap
        self.windsMap = uWindsMap

        # IF NOT HD, DOWNLOAD HIGH ALTITUDE HD DATA
        if not self.HD:
            self._highAltitudeGFS.downloadForecast()
            self._highAltitudePressure = self._highAltitudeGFS.interpolateData('p')

        logger.debug('Forecast successfully downloaded!')

        return True


    def interpolateData(self, *variables):
        """
        Set up a linear 4d interpolation for each variable given and returns it. The returned
        interpolator can then be used by calling it with standard 4D coordinates (lat, lon, alt, time).
        Acceptable values for *variables can be a combination of any of these (separated by commas):
            't'|'temp'|'temperature'        temperature interpolation
            'p'|'press'|'pressure'          pressure interpolation
            'd'|'windrct'|'wind_direction'  wind direction interpolation
            's'|'windspd'|'wind_speed'      wind speed interpolation

        The linear interpolation is based on the interpolate.Linear4DInterpolator class.

        Warning: the time parameter required by the Linear4DInterpolator objects must be in GFS units!
        Use the getGFStime(time) method for conversion.

        Returns a list of Linear4DInterpolator objects, in the same order as *variables.
        """

        results = []

        # Cycle through each variable requested and append the appropriate interpolator to the results list.
        for variable in variables:
            if variable in ('temp', 't', 'temperature'):
                # Interpolate temperature
                if self.HD:
                    results.append(
                        GFS_data_interpolator(self, self.temperatureData, self.temperatureMap.mappingCoordinates))
                else:
                    results.append(
                        GFS_data_interpolator(self, self.temperatureData, self.temperatureMap.mappingCoordinates,
                                              self._highAltitudeGFS.interpolateData('t')))

            elif variable in ('press', 'p', 'pressure'):
                # Interpolate pressure
                results.append(self._pressure_interpolator)

            elif variable in ('windrct', 'd', 'wind_direction'):
                # Interpolate wind direction
                if self.HD:
                    results.append(GFS_data_interpolator(self, self.windDirData, self.windsMap.mappingCoordinates))
                else:
                    results.append(GFS_data_interpolator(self, self.windDirData, self.windsMap.mappingCoordinates,
                                                         self._highAltitudeGFS.interpolateData('d')))

            elif variable in ('windspd', 's', 'wind_speed'):
                # Interpolate wind speed
                if self.HD:
                    results.append(GFS_data_interpolator(self, self.windSpeedData, self.windsMap.mappingCoordinates))
                else:
                    results.append(GFS_data_interpolator(self, self.windSpeedData, self.windsMap.mappingCoordinates,
                                                         self._highAltitudeGFS.interpolateData('s')))

            else:
                logger.error('A wrong interpolation parameter (%s) was passed to the interpolator.' % variable)

        # Return the results.
        if len(results) == 1:
            return results[0]
        else:
            return results


    def getGFStime(self, time):
        """
        Convert standard datetime.datetime objects to GFS time units.
        The time parameter should be a datetime.datetime object.

        Returns a float corresponding to the converted time, or nan if there was no downloaded data found.
        """

        # Check if data is available. If it isn't return NaN.
        if self.firstAvailableTime is None or self.altitudeMap is None:
            return float('nan')
        else:
            try:
                # Try to define the time from the first dataset.
                # If it fails for any reasons, return NaN.
                timeFromFirstDataset = time - self.firstAvailableTime
                return self.altitudeMap.fwdTime[
                           0] + timeFromFirstDataset.days + timeFromFirstDataset.seconds / 3600. / 24.
            except TypeError:
                return float('nan')


    def _generate_matrix(self, dataStreams):
        """
        This is a private method used to generate data matrices and coordinates mapping.
        It's called by the downloadForecast() method when data is downloaded and it supports joining together two
        datasets.

        It should NOT be used alone!
        """

        overallResults = []
        overallMaps = []

        # Run this either once or twice, according to how many datasets have been downloaded (Greenwich meridian
        # crossing.
        for dataStream in dataStreams:

            dataLines = dataStream.split('\n')

            # Count how many latitude, longitude, pressure and time points are available in the datastream.
            # This is used to initialize the results matrix.
            timePoints = int(dataLines[0].split()[1].split('][')[0][1:])
            pressurePoints = int(dataLines[0].split()[1].split('][')[1])
            latitudePoints = int(dataLines[0].split()[1].split('][')[2])
            longitudePoints = int(dataLines[0].split()[1].split('][')[3][:-1])
            totalPoints = timePoints * pressurePoints * latitudePoints * longitudePoints

            # Initialize the matrix
            results = numpy.zeros(totalPoints).reshape((latitudePoints, longitudePoints, pressurePoints, timePoints))

            # Populate the results matrix
            for line in dataLines[1:-12]:
                # Skip empty lines
                if line == '': continue

                # Find the indexes related to this particular latitude line
                #
                # WARNING: THIS IS LIKELY TO CAUSE ISSUES IF THE GFS FORMAT CHANGES!
                #
                # If the GFS data format changes, modify it here!
                # This is VERY format-dependent!
                timeIndex = int(line.split(',')[0].split('][')[0][1:])
                pressureIndex = int(line.split(',')[0].split('][')[1])
                latitudeIndex = int(line.split(',')[0].split('][')[2][:-1])

                # Store values
                results[latitudeIndex, :, pressureIndex, timeIndex] = [float(x) for x in line.split(',')[1:]]

            # Generate the mapping. This is an object containing mapping information between GFS indices for lat,lon,press,
            # time and actual values, and viceversa. "Forward" maps are GFS indices -> real values and are LISTS.
            # "Reverse" maps are real values -> GFS indices and are DICTIONARIES.
            # These are used to be able to find a data value given real coordinates and for interpolation.

            resultsMap = GFS_Map()

            resultsMap.fwdLatitude = [float(lat) for lat in dataLines[-4].split(',')]
            resultsMap.revLatitude = {lat: ind for (lat, ind) in
                                      zip(resultsMap.fwdLatitude, xrange(len(resultsMap.fwdLatitude)))}

            resultsMap.fwdLongitude = [float(lon) - 360 if float(lon) > 180 else float(lon) for lon in
                                       dataLines[-2].split(',')]
            resultsMap.revLongitude = {lon: ind for (lon, ind) in
                                       zip(resultsMap.fwdLongitude, xrange(len(resultsMap.fwdLongitude)))}

            resultsMap.fwdPressure = [float(press) for press in dataLines[-6].split(',')]
            resultsMap.revPressure = {press: ind for (press, ind) in
                                      zip(resultsMap.fwdPressure, xrange(len(resultsMap.fwdPressure)))}

            resultsMap.fwdTime = [float(time) for time in dataLines[-8].split(',')]
            resultsMap.revTime = {time: ind for (time, ind) in zip(resultsMap.fwdTime, xrange(len(resultsMap.fwdTime)))}

            overallResults.append(results)
            overallMaps.append(resultsMap)

        # Join results together, if needed
        if len(overallResults) > 1:
            # Identify which one is to the right
            isFirstOneEast = overallMaps[0].fwdLongitude[0] == 0
            if isFirstOneEast:
                results = numpy.concatenate((overallResults[1], overallResults[0]), axis=1)
                resultsMap = overallMaps[1].rjoin(overallMaps[0])
            else:
                results = numpy.concatenate((overallResults[0], overallResults[1]), axis=1)
                resultsMap = overallMaps[0].rjoin(overallMaps[1])
        else:
            results = overallResults[0]
            resultsMap = overallMaps[0]

        # Store all the forward maps to a single list called mappingCoordinates.
        # This is used by interpolateData(*variables) to correctly initialize the interpolators.
        resultsMap.mapCoordinates()

        return results, resultsMap


    def _pressure_interpolator(self, lat, lon, alt, time):
        """
        Method used to extract pressure from (lat,lon,alt,time) coordinates.
        This is essentially a 3D interpolation of the column of altitudes (ie not just a single value) at the point
        defined just by latitude, longitude and time.
        Once the column of altitudes is extracted, a UnivariateSpline (1D) interpolation is run between that column of
        altitudes and the corresponding pressures (which are always constant).
        Then, root at which the altitude corresponds to the target (the one passed as a parameter) is found.

        This corresponds to the pressure at that point, which is then returned.

        If the requested point is outside latitude, longitude or time bounds, the nearest value available for the out-
        of-bounds axis is used.
        If the requested point is outside altitude bounds, ISA conditions are used above the upper limit and the nearest
        value is used below the lower limit.

        Note: This is not using the standard interpolate.Linear4DInterpolator because it's NOT a 4D interpolation.
        """

        # Clip out-of-bounds coordinates and limit them
        lat = numpy.clip(lat, self.altitudeMap.fwdLatitude[0], self.altitudeMap.fwdLatitude[-1])
        lon = numpy.clip(lon, self.altitudeMap.fwdLongitude[0], self.altitudeMap.fwdLongitude[-1])
        time = numpy.clip(time, self.altitudeMap.fwdTime[0], self.altitudeMap.fwdTime[-1])

        # Find closest indices and subtract 1 if the upper limit is being reached, to avoid a KeyError
        i = numpy.digitize([lat], self.altitudeMap.fwdLatitude)[0]
        if i == len(self.altitudeMap.fwdLatitude):
            i -= 1
        idxLat = [i - 1, i]
        i = numpy.digitize([time], self.altitudeMap.fwdTime)[0]
        if i == len(self.altitudeMap.fwdTime):
            i -= 1
        idxTime = [i - 1, i]

        # Reverse mapping for longitude, to account for crossing the 180th meridian
        # Note: this method is used for two reasons:
        #   1. It's faster than numpy.digitize
        #   2. numpy.digitize fails at longitudes between -180 and 179.5 (for HD), since there isn't an entry for -180.
        lonGrid = [floor(lon / self.lonStep) * self.lonStep, ceil(lon / self.lonStep) * self.lonStep]
        if lonGrid[0] == -180:
            lonGrid[0] = 180

        if lonGrid[0] == lonGrid[1]:
            try:
                idxLon = [self.altitudeMap.revLongitude[lonGrid[0]],
                          self.altitudeMap.revLongitude[lonGrid[1] + self.lonStep]]
            except KeyError:
                idxLon = [self.altitudeMap.revLongitude[lonGrid[0] - self.lonStep],
                          self.altitudeMap.revLongitude[lonGrid[1]]]
        else:
            idxLon = [self.altitudeMap.revLongitude[lonGrid[0]], self.altitudeMap.revLongitude[lonGrid[1]]]

        try:
            fracLat = 1 - abs((lat - self.altitudeMap.fwdLatitude[idxLat[0]]) / (
                self.altitudeMap.fwdLatitude[idxLat[1]] - self.altitudeMap.fwdLatitude[idxLat[0]]))
            fracLon = 1 - abs((lon - self.altitudeMap.fwdLongitude[idxLon[0]]) / (
                self.altitudeMap.fwdLongitude[idxLon[1]] - self.altitudeMap.fwdLongitude[idxLon[0]]))
            fracTime = 1 - abs((time - self.altitudeMap.fwdTime[idxTime[0]]) / (
                self.altitudeMap.fwdTime[idxTime[1]] - self.altitudeMap.fwdTime[idxTime[0]]))
        except IndexError:
            # Requested point outside bounds

            logger.error('An error occurred while clipping coordinate points.')
            raise

        alt00 = fracLat * self.altitudeData[idxLat[0], idxLon[0], :, idxTime[0]] + (1 - fracLat) * self.altitudeData[
                                                                                                   idxLat[1], idxLon[0],
                                                                                                   :, idxTime[0]]
        alt01 = fracLat * self.altitudeData[idxLat[0], idxLon[0], :, idxTime[1]] + (1 - fracLat) * self.altitudeData[
                                                                                                   idxLat[1], idxLon[0],
                                                                                                   :, idxTime[1]]
        alt10 = fracLat * self.altitudeData[idxLat[0], idxLon[1], :, idxTime[0]] + (1 - fracLat) * self.altitudeData[
                                                                                                   idxLat[1], idxLon[1],
                                                                                                   :, idxTime[0]]
        alt11 = fracLat * self.altitudeData[idxLat[0], idxLon[1], :, idxTime[1]] + (1 - fracLat) * self.altitudeData[
                                                                                                   idxLat[1], idxLon[1],
                                                                                                   :, idxTime[1]]
        alt_time0 = fracLon * alt00 + (1 - fracLon) * alt10
        alt_time1 = fracLon * alt01 + (1 - fracLon) * alt11
        altitude = fracTime * alt_time0 + (1 - fracTime) * alt_time1 - alt

        if altitude.min() > 0:
            # NEAREST NEIGHBOR: if requested point is below minimum altitude, return lowest point available
            return self.altitudeMap.fwdPressure[0]
        elif altitude.max() < 0:
            # NEAREST NEIGHBOR: if requested point is above maximum altitude, return highest point available
            if self.HD:
                return self.altitudeMap.fwdPressure[-1]
            else:
                return self._highAltitudePressure(lat, lon, alt, time)
        else:
            # LINEAR INTERPOLATION
            # (see method documentation for details)
            f = UnivariateSpline(self.altitudeMap.fwdPressure[::-1], altitude[::-1], s=0)
            return f.roots()[0]


class GFS_High_Altitude_Handler(GFS_Handler):
    def __init__(self, lat, lon, date_time, forecast_duration=4, debugging=False, log_to_file=False):
        global logger

        # Initialize Parameters
        self.launchDateTime = date_time
        self.lat = lat
        self.lon = lon
        self.maxFlightTime = forecast_duration
        self.cycleDateTime = None
        self.firstAvailableTime = None

        # These are the 4D data matrices with all the information needed.
        self.altitudeData = None
        self.temperatureData = None
        self.windDirData = None
        self.windSpeedData = None

        # These are variables storing the mapping criteria of the 4D matrices.
        # They are arrays of dictionaries, one per dimension (i.e. axis) of the matrix.
        # The dictionary keys are the latitude, longitude, pressure and time, while the values are the
        # corresponding matrix data indices.
        self.altitudeMap = None
        self.temperatureMap = None
        self.windsMap = None

        # Grid size setup
        # NOTE: Grid sizes are defined as the difference between the highest and the lowest lat/lon requested, NOT the
        # difference between the highest and the requested point or the lowest and the requested point!
        if forecast_duration == 4:
            # Standard flight
            self.latGridSize = 6
            self.lonGridSize = 12
        else:
            # Non-standard flight (either floating or customized from command line)
            self.latGridSize = 2 * ceil(0.1 * self.maxFlightTime + 0.6) + 3
            self.lonGridSize = 2 * ceil(0.9 * self.maxFlightTime) + 3

        self.HD = True
        self.latStep = 0.5
        self.lonStep = 0.5

        # SETUP ERROR LOGGING AND DEBUGGING

        logger = logging.getLogger('GFS_High_Alt')

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


    def downloadForecast(self, progressHandler=None):
        """
        Connect to the Global Forecast System and download the closest cycle available to the date_time
        required. The cycle date and time is stored in the object's cycleDateTime variable.
        Returns TRUE if the download was successful, FALSE otherwise.
        """

        #############################################################################################################
        # INITIALIZE DATA, OR RESET IT IF THEY ALREADY EXISTS
        temperatureMatrix = None
        geopotentialMatrix = None
        uWindsMatrix = None
        vWindsMatrix = None
        temperatureMap = None
        geopotentialMap = None
        uWindsMap = None


        #############################################################################################################
        # EXTRA PARAMETERS
        # The base URL depends on whether the HD service has been requested or not.
        weatherParameters = ['tmpprs', 'hgtprs', 'ugrdprs', 'vgrdprs']
        baseURL = 'http://nomads.ncep.noaa.gov:9090/dods/gfs_hd/'


        #############################################################################################################
        # INITIALIZE TIME HANDLERS AND CALCULATE MOST RECENT CYCLE POSSIBLY AVAILABLE

        simulationDateTime = self.launchDateTime
        currentDateTime = datetime.now()

        if simulationDateTime < currentDateTime:
            # Simulating a flight in the past. In this case, use the simulation date time to calculate the cycles
            # that need to be downloaded.
            currentDateTime = simulationDateTime

        # These are the hours at which new forecast cycles are issued
        dailyCycles = [0, 6, 12, 18]
        # Determine which cycle issuing hour is the closest to the current one
        cycleTime = dailyCycles[numpy.digitize([currentDateTime.hour], dailyCycles)[0] - 1]
        latestCycleDateTime = datetime(currentDateTime.year, currentDateTime.month, currentDateTime.day, cycleTime)


        #############################################################################################################
        # INITIALIZE COMMON REQUEST PARAMETERS. These are lat, lon and alt and do not vary with each cycle.
        # Perform some validation checks to see if they need to be limited (close to the poles) or two requests need to
        # be made (close to the Greenwich meridian).
        requestAllLongitudes = False
        multipleRequests = False
        # This is the grid size converted to GFS units. It basically returns the same number if the gridSize is even,
        # or it returns gridSize-1 if gridSize is odd. This is to get the right grid around the launch point.
        latGridStep = (int(self.latGridSize) / 2) / self.latStep
        lonGridStep = (int(self.lonGridSize) / 2) / self.lonStep

        # ALTITUDE (Download ALL altitude levels available)
        requestAltitude = [41, 46]

        # LATITUDE
        targetLatitude = (round(self.lat) + 90) / self.latStep
        requestLatitude = [targetLatitude - latGridStep, targetLatitude + latGridStep]

        if requestLatitude[0] < 0:
            requestLatitude[0] = 0
        if requestLatitude[1] > 360:
            requestLatitude[1] = 360

        # Request all longitudes if close to the poles
        if requestLatitude[0] < 6 or requestLatitude[1] > 354:
            requestAllLongitudes = True

        # LONGITUDE
        # Note: There is no need to check if the grid size is higher than the whole world, causing overlapping regions
        # to only request data in the overlapped region, since a grid size approximately equal to half of the world is
        # enough to hit the poles and therefore request worldwide data (longitudes are very close at the poles, hence
        # it's worth keeping data for all of them)

        if requestAllLongitudes:
            requestLongitude = [0, 719]
        else:
            if self.lon >= 0:
                targetLongitude = round(self.lon) / self.lonStep
            else:
                targetLongitude = (360 - abs(round(self.lon))) / self.lonStep

            requestLongitude = [targetLongitude - lonGridStep, targetLongitude + lonGridStep]

            # Check if the values are within the bounds and correct if needed
            if requestLongitude[0] < 0:
                requestLongitude[0] += 720
            if requestLongitude[1] > 719:
                requestLongitude[1] -= 720

            # Check if crossing the Greenwich meridian and split the requests
            # If the Greenwich meridian is being crossed, the left bound of the longitude interval will have a
            # higher value than the right one (in HD: Western hemisphere [360:719], Eastern hemisphere:[0:359]), so if
            # the difference between the right and the left one is negative, the Greenwich meridian is being crossed
            if requestLongitude[1] - requestLongitude[0] < 0:
                # SPLIT
                requestLongitude = [[0, requestLongitude[1]], [requestLongitude[0], 719]]
                multipleRequests = True


        #############################################################################################################
        # TRY TO DOWNLOAD DATA WITH THE LATEST CYCLE. IF NOT AVAILABLE, TRY WITH AN EARLIER ONE

        # pastCycle is a variable that indicates how many cycles in the past we're downloading data from.
        # It first tries with 0, indicating the most recent cycle: this is calculated, knowing that cycles are issued
        # every six hours. Sometimes, however, it takes a few hours for data to become available, or a cycle is skipped.
        # This means that data for the latest cycle is not guaranteed to be available.
        # If data is not available, it tries with 1 cycle older, until one is found with data available. If no cycles
        # are found, the method returns FALSE.
        pastCycle = 0
        while True:
            thisCycle = latestCycleDateTime - timedelta(hours=pastCycle * 6)
            self.cycleDateTime = thisCycle

            # This is just a flag to make sure that if no data is found, it stops the whole download and re-starts
            # the loop with an older cycle.
            thisCycleNotAvailable = False

            # Initialize time parameter
            timeFromForecast = simulationDateTime - thisCycle
            hoursFromForecast = timeFromForecast.total_seconds() / 3600.

            # GFS time index for the first dataset to be requested
            requestTime = floor(hoursFromForecast / 3.)

            # This stores the actual time of the first dataset downloaded. It's going to be used to convert real time
            # to GFS "time coordinates" (see getGFStime(time) function)
            self.firstAvailableTime = self.cycleDateTime + timedelta(hours=requestTime * 3)

            # Always download an extra time dataset
            requestTime = [requestTime, requestTime + ceil(self.maxFlightTime / 3.) + 1]

            #########################################################################################################
            # DOWNLOAD DATA FOR ALL PARAMETERS

            for requestVariable in weatherParameters:

                if multipleRequests:
                    numberOfRequests = 2
                else:
                    numberOfRequests = 1

                dataResults = []

                # Check if we need more than 1 request (ie if we are crossing the Greenwich meridian)
                for req_number in xrange(numberOfRequests):

                    if multipleRequests:
                        thisRequestLongitude = requestLongitude[req_number]
                    else:
                        thisRequestLongitude = requestLongitude

                    requestURL = '%sgfs%s%d%02d%02d/gfs%s_%02dz.ascii?%s[%d:%d][%d:%d][%d:%d][%d:%d]' % (
                        baseURL,
                        '_hd',
                        thisCycle.year,
                        thisCycle.month,
                        thisCycle.day,
                        '_hd',
                        thisCycle.hour,
                        requestVariable,
                        requestTime[0], requestTime[1],
                        requestAltitude[0], requestAltitude[1],
                        requestLatitude[0], requestLatitude[1],
                        thisRequestLongitude[0], thisRequestLongitude[1]
                    )

                    logger.debug('Requesting URL: %s' % requestURL)

                    try:
                        HTTPresponse = urllib2.urlopen(requestURL)
                    except:
                        logger.error('Error while connecting to the GFS server.')
                        logger.error('URL: %s' % requestURL)
                        return False
                    dataResults.append(HTTPresponse.read())

                    if dataResults[-1][0] == "<":
                        # Data from this cycle is not available. Try with the earlier one.
                        logger.debug('Cycle not available.')
                        thisCycleNotAvailable = True
                        break

                if thisCycleNotAvailable:
                    break

                if requestVariable == 'tmpprs':
                    temperatureMatrix, temperatureMap = self._generate_matrix(dataResults)
                    logger.debug('Temperature data downloaded and processed.')
                elif requestVariable == 'hgtprs':
                    geopotentialMatrix, geopotentialMap = self._generate_matrix(dataResults)
                    logger.debug('Altitude data downloaded and processed.')
                elif requestVariable == 'ugrdprs':
                    uWindsMatrix, uWindsMap = self._generate_matrix(dataResults)
                    logger.debug('U Winds data downloaded and processed.')
                elif requestVariable == 'vgrdprs':
                    logger.debug('V Winds data downloaded and processed.')
                    vWindsMatrix, vWindsMap = self._generate_matrix(dataResults)




            # Restart the loop if the data wasn't available.
            if thisCycleNotAvailable:
                if pastCycle == 24:
                    logger.error('No GFS cycles available found!')
                    return False
                else:
                    pastCycle += 1
                    logger.debug('Moving to next cycle...')
                    continue

            # If it didn't break in the previous if statement, data was available and has been downloaded.
            # End of the loop.
            break

        #############################################################################################################
        # PROCESS DATA AND PERFORM CONVERSIONS AS REQUIRED

        # Convert temperatures from Kelvin to Celsius
        temperatureMatrix -= 273.15

        # Convert geopotential height to geometric altitude
        altitudeMatrix = geopotentialMatrix * earthRadius / ( earthRadius - geopotentialMatrix )

        # Convert u and v winds to wind direction and wind speed matrices

        # Store the current shape of the 4D matrices
        matrixShape = uWindsMatrix.shape
        # Convert to KNOTS and the turn into direction and speed
        dirspeedWinds = map(tools.uv2dirspeed, (uWindsMatrix * 1.9438445).ravel(), (vWindsMatrix * 1.9438445).ravel())
        # Extract results
        windDirectionMatrix = numpy.array([dirspeed[0] for dirspeed in dirspeedWinds]).reshape(matrixShape)
        windSpeedMatrix = numpy.array([dirspeed[1] for dirspeed in dirspeedWinds]).reshape(matrixShape)



        # Store results
        self.temperatureData = temperatureMatrix
        self.altitudeData = altitudeMatrix
        self.windDirData = windDirectionMatrix
        self.windSpeedData = windSpeedMatrix

        self.temperatureMap = temperatureMap
        self.altitudeMap = geopotentialMap
        self.windsMap = uWindsMap

        logger.debug('High altitude HD forecast successfully downloaded!')

        return True


class GFS_Map:
    """
    Private class used to store 4D mapping data for a specific parameter.

    Methods available:
        rjoin(map)          Joins together this map and the one passed in the parameter, adding it to the right. This
                            method returns the joint GFS_Map
        ljoin(map)          Joins together this map and the one passed in the parameter, adding it to the left. This
                            method returns the joint GFS_Map.
        mapCoordinates()    Prepares the mappingCoordinates variable by putting together all forward and reverse maps.
                            Used by the Linear4DInterpolator to access interpolation data maps.
    """

    def __init__(self):
        # Forward lists map GFS indices to real world coordinates.
        self.fwdLatitude = []
        self.fwdLongitude = []
        self.fwdPressure = []
        self.fwdTime = []

        # Reverse dictionaries map real world coordinates to GFS indices.
        # (keys: real world coordinates, values: GFS indices)
        self.revLatitude = {}
        self.revLongitude = {}
        self.revPressure = {}
        self.revTime = {}

        # Interpolation coordinates is a list of 4D real world coordinates for each value in the data matrix
        self.mappingCoordinates = []

    def rjoin(self, data_map):
        # This method joins together this map and the one passed in the parameter, concatenating it to the RIGHT.
        # Returns the joint GFS_Map

        lonLen = len(self.fwdLongitude)

        if self.fwdLatitude == data_map.fwdLatitude and self.fwdPressure == data_map.fwdPressure and self.fwdTime == data_map.fwdTime:
            # Join
            self.fwdLongitude = list(numpy.concatenate((self.fwdLongitude, data_map.fwdLongitude)))

        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        if self.revLatitude == data_map.revLatitude and self.revPressure == data_map.revPressure and self.revTime == data_map.revTime:
            # Join
            self.revLongitude.update({key: value + lonLen for (key, value) in data_map.revLongitude.iteritems()})

        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        return self

    def ljoin(self, data_map):
        # This method joins together this map and the one passed in the parameter, concatenating it to the LEFT.
        # Returns the joint GFS_Map

        lonLen = len(data_map.fwdLongitude)

        if self.fwdLatitude == data_map.fwdLatitude and self.fwdPressure == data_map.fwdPressure and self.fwdTime == data_map.fwdTime:
            self.fwdLongitude = list(numpy.concatenate((data_map.fwdLongitude, self.fwdLongitude)))
        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        if self.revLatitude == data_map.revLatitude and self.revPressure == data_map.revPressure and self.revTime == data_map.revTime:
            self.revLongitude = {key: value + lonLen for (key, value) in self.revLongitude.iteritems()}.update(
                data_map.revLongitude)
        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        return self

    def mapCoordinates(self):
        # Prepare the mappingCoordinates variable by putting together all forward and reverse maps.
        # Used by the Linear4DInterpolator to access interpolation data maps.

        self.mappingCoordinates = [
            self.fwdLatitude,
            self.fwdLongitude,
            self.fwdPressure,
            self.fwdTime,
            self.revLatitude,
            self.revLongitude,
            self.revPressure,
            self.revTime
        ]


class GFS_data_interpolator:
    """
    Private class used by GFS_Handler to interpolate data.

    This class acts as an interface to be able to use the correct coordinates (lat,lon,alt,time) when requesting
    interpolated data.

    Upon initialization, the following parameters need to be passed:
        GFS_Handler     The GFS handler object. This is used to fetch pressure data upon request
        data            The 4D matrix containing all the data that needs to be interpolated
        map             The map of the data. This should be the mappingCoordinates variable generated by the
                        GFS_Map.mapCoordinates() method.

    Once initialized, the object can be called with the requested coordinates.
    At this point, the object converts the altitude to pressure using the GFS_Handler's _pressure_interpolator(...)
    method and then uses that to interpolate the required data and return it.

    """


    def __init__(self, GFS_Handler, data, dmap, high_alt_interpolator=None, min_pressure=None):
        # Store data
        self.GFS = GFS_Handler
        self._interpolator = Linear4DInterpolator(data, dmap)
        self._high_alt_interpolator = high_alt_interpolator
        self._min_press = min_pressure

    def __call__(self, lat, lon, alt, time):
        # Get the pressure with a pressure interpolator and use it for the generic one.

        if type(time) is not float and type(time) is not int:
            logger.error('The time passed is not in GFS coordinates! Use the getGFStime() method.')
            return

        pressure = self.GFS._pressure_interpolator(lat, lon, alt, time)

        if self._min_press is not None and pressure < self._min_press:
            return self._high_alt_interpolator(lat, lon, alt, time)
        else:
            return self._interpolator(
                lat,
                lon,
                pressure,
                time
            )