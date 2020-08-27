# coding=utf-8

"""
Module containing the classes for interacting with the NOAA's Global Forecast
System. 

GFS_Handler is the primary class used outside of this module. See class
documentation for use

University of Southampton
"""
from datetime import datetime, timedelta
from math import floor, ceil
from six.moves import range, builtins
from six.moves.urllib.request import urlopen
import logging
import grequests
import itertools
import numpy
from scipy.interpolate import UnivariateSpline
import urllib.parse

from . import global_tools as tools
from .interpolate import Linear4DInterpolator

# Error and warning logger
logger = logging.getLogger(__name__)

earthRadius = 6371009  # m


# Pass through the @profile decorator if line profiler (kernprof) is not in use
try:
    builtins.profile
except AttributeError:
    def profile(func):
        return func


def get_urldict_async(urls_dict, hooks_dict=None):
    """An asynchronous helper function for making multiple url requests per
    key in url_dict

    Uses grequests to produce the requests concurrently.

    Parameters
    ----------
    url_dict : dict of lists
        (key: [url1, url2...]) pairs for logically related urls.

    hooks_dict : dict
        The post request hook for response processing (see requests lib for
        examples)

    Returns
    -------
    results_dict: dict of list
        (key: [reponse_text1, reponse_text2, ...]) for all keys in url_dict
    """
    # first get all urls as a list
    urls_list = list(itertools.chain.from_iterable(urls_dict.values()))
    reqs = (grequests.get(u, hooks=hooks_dict) for u in urls_list)
    responses = grequests.map(reqs)
    results = {urllib.parse.unquote(r.url): r.text for r in responses}
    if any(result[0] == "<" for result in results.values()):
        logger.debug("GFS cycle not found.")
        return
    else:
        # get the url expected at each location for each key in url_dict, then
        # insert the retrieved data at this location
        results_dict = {key: [results[url] for url in url_list] for key, url_list in urls_dict.items()}
        return results_dict



class GFS_Handler(object):
    """
    Use this to perform any operation related to the Global Forecast System.

    Downloads weather forecast data for any location on the planet for up to 10
    days from the closest cycle. The data is downloaded from either the HD
    service, providing data on a 0.25 deg latitude x 0.25 deg longitude x 26
    altitude levels grid, or from the SD service, with a coarser grid. Note
    that even when the HD service is used, the high altitude data is always
    requested in SD since the HD coverage is only up to 10mB pressure alt.

    Parameters
    ----------
    lat : float
        latitude of the location where forecast is required.
    lon : float
        longitude of the location where forecast is required.
    date_time : :obj:`datetime.datetime.`
        UTC time of the required forecast (min 30 days from today, max 10 days
        from latest cycle issuing time, see warning above).
    [HD] : bool (default True)
        if FALSE, non-HD forecast will be downloaded. Otherwise, HD forecast
        will be downloaded and, if the data size is too high, the GFS_Handler
        will automatically switch to non-HD.
    [forecastDuration] : scalar (default 4)
        the duration in hours of the forecast data required.
    [use_async] : bool (default True)
        Use an asynchronous request for downloads. This should speed up the
        download, but may incur larger memory overhead for large forecastDuration.
    [requestSimultaneous] : bool (default True)
        If True, populate a dictionary of responses from the web download
        requests, then process the data. If False, each response will be
        removed once the data has been processed: This has better memory
        overhead for large ForecastDuration, but is slower, and does not work
        with asynchronous requesting (use_async)
    [debugging] : bool (default False)
        If debugging is TRUE, all the information available will be logged.
        Otherwise, only errors will be logged.
    [log_to_file] : bool (default False)
        If log_to_file is set to TRUE, all error and debug logs will be stored
        in a file called error.log. If FALSE, all error and debug logs will be
        directly displayed on the terminal or command line.
    [progressHandler] : function (default None)
        Progress for each downloaded parameter (in %) will be passed to this
        function, if provided.

    Notes
    -----
    * The NOAA issues new datasets (cycles) at 00:00, 06:00, 12:00 and 18:00 UTC
    every day and the 10 days of forecast data. Refer to the cycle issuing
    time.

    * As it's not possible to know in advance what the latest cycle available
    is (there is a delay between the cycle issuing time and its
    availability), it's recommended to limit the time span of data downloads
    from the GFS to 8 days.

    :Example:
    >>> # Download forecast data for tomorrow at the University of Southampton

    >>> import datetime

    >>> myGFSlink = GFS_Handler(50.93543,-1.39619,datetime.now()+datetime.timedelta(days=1))
    >>> myGFSlink.downloadForecast()

    >>> getTemp,getPress = myGFSlink.interpolateData('t','p')
    >>> # Request temperature at 51.2 lat, 0.46 lon, 32000m alt, 1 hr from now
    >>> getTemp(51.2,0.46,32000,myGFSlink.getGFStime(
    >>>    datetime.now()+datetime.timedelta(days=1,seconds=3600)))
    """
    weatherParameters = {'tmpprs': 'Temperature',
                          'hgtprs': 'Altitude',
                          'ugrdprs': 'U Winds',
                          'vgrdprs': 'V Winds'}

    def __init__(self, lat, lon, date_time, HD=True, forecastDuration=4,
        use_async=True, requestSimultaneous=True, debugging=False,
        progressHandler=None):
        # Initialize Parameters
        self.launchDateTime = date_time
        self.lat = lat
        self.lon = lon
        self.forecastDuration = forecastDuration
        self.HD = HD
        self.use_async = use_async
        self.requestSimultaneous = requestSimultaneous
        self.cycleDateTime = None
        self.firstAvailableTime = None

        # These are the 4D data matrices with all the information needed.
        self.altitudeData = None
        self.temperatureData = None
        self.windDirData = None
        self.windSpeedData = None

        # These are variables storing the mapping criteria of the 4D matrices.
        # They are arrays of dictionaries, one per dimension (i.e. axis) of the
        # matrix. The dictionary keys are the latitude, longitude, pressure and
        # time, while the values are the corresponding matrix data indices.
        self.altitudeMap = None
        self.temperatureMap = None
        self.windsMap = None

        requestAllLongitudes = False
        multipleRequests = False

        # Grid size setup
        # NOTE: Grid sizes are defined as the difference between the highest
        # and the lowest lat/lon requested, NOT the difference between the
        # highest and the requested point or the lowest and the requested point
        if forecastDuration == 4:
            # Standard flight
            self.latGridSize = 6
            self.lonGridSize = 12
        else:
            # Non-standard flight (floating or customized from command line)
            self.latGridSize = 2 * ceil(0.1 * self.forecastDuration + 0.6) + 3
            self.lonGridSize = 2 * ceil(0.9 * self.forecastDuration) + 3

        # Automatically switch to non HD if close to the poles, as all the
        # longitudes will be downloaded
        if lat < -80 or lat > 80:
            self.HD = False

        # Force non HD weather if requested data is too high
        if self.latGridSize * self.lonGridSize * ceil(self.forecastDuration / 3.) > 250:
            self.HD = False

        # HD/SD forecast setup
        if self.HD:
            self.latStep = 0.25
            self.lonStep = 0.25
            # TODO: Using a derived class and storing it as an attribute for a
            #   special case is a bad design pattern in general: This could be
            #   avoided by using two gfs handlers in the weather class, or
            #   downloading data for specific altitudes and lat/lonStep if
            #   self.HD is True
            # Prepare download of high altitude SD data
            self._highAltitudeGFS = GFS_High_Altitude_Handler(lat,
                lon, date_time, forecastDuration, debugging)
            self._highAltitudePressure = None
        else:
            self.latStep = 0.5
            self.lonStep = 0.5

        # This is the grid size converted to GFS units. It basically returns
        # the same number if the gridSize is even, or it returns gridSize-1 if
        # gridSize is odd. This is to get the right grid around the launch
        # point.
        latGridStep = (int(self.latGridSize) / 2) / self.latStep
        lonGridStep = (int(self.lonGridSize) / 2) / self.lonStep

        # ALTITUDE (Download ALL altitude levels available)
        self.requestAltitude = {
            True: [0, 25],
            False: [0, 46]
        }[self.HD];

        # LATITUDE
        targetLatitude = (round(self.lat) + 90) / self.latStep
        self.requestLatitude = [targetLatitude - latGridStep, targetLatitude + latGridStep]

        if self.HD:
            # Limit latitude to +/- 90 degrees
            if self.requestLatitude[0] < 0:
                self.requestLatitude[0] = 0
            if self.requestLatitude[1] > 720:
                self.requestLatitude[1] = 720

            # Request all longitudes if close to the poles
            if self.requestLatitude[0] < 40 or self.requestLatitude[1] > 680:
                requestAllLongitudes = True
        else:
        # Limit latitude to +/- 90 degrees
            if self.requestLatitude[0] < 0:
                self.requestLatitude[0] = 0
            if self.requestLatitude[1] > 360:
                self.requestLatitude[1] = 360

            # Request all longitudes if close to the poles or if simulation is beyond 2 days
            if self.requestLatitude[0] < 20 or self.requestLatitude[1] > 340 or self.forecastDuration > 48:
                requestAllLongitudes = True

        # LONGITUDE
        # Note: There is no need to check if the grid size is higher than the
        # whole world, causing overlapping regions to only request data in the
        # overlapped region, since a grid size approximately equal to half of
        # the world is enough to hit the poles and therefore request worldwide
        # data (longitudes are very close at the poles, hence it's worth
        # keeping data for all of them)

        if requestAllLongitudes:
            self.requestLongitudes = {
                True: [[0, 1439]],
                False: [[0, 719]]
            }[self.HD]
        else:
            if self.lon >= 0:
                targetLongitude = round(self.lon) / self.lonStep
            else:
                targetLongitude = (360 - abs(round(self.lon))) / self.lonStep

            self.requestLongitudes = [[targetLongitude - lonGridStep, targetLongitude + lonGridStep]]

            # Check if the values are within the bounds and correct if needed
            if self.HD:
                if self.requestLongitudes[0][0] < 0:
                    self.requestLongitudes[0][0] += 1440
                if self.requestLongitudes[0][1] > 1439:
                    self.requestLongitudes[0][1] -= 1440
            else:
                if self.requestLongitudes[0][0] < 0:
                    self.requestLongitudes[0][0] += 720
                if self.requestLongitudes[0][1] > 719:
                    self.requestLongitudes[0][1] -= 720

            # Check if crossing the Greenwich meridian and split the requests
            # If the Greenwich meridian is being crossed, the left bound of the
            # longitude interval will have a higher value than the right one
            # (in HD: Western hemisphere [360:719], Eastern hemisphere 
            # [0:359]), so if the difference between the right and the left one
            # is negative, the Greenwich meridian is being crossed
            if self.requestLongitudes[0][1] - self.requestLongitudes[0][0] < 0:
                # SPLIT
                self.requestLongitudes = {
                    True: [[0, self.requestLongitudes[0][1]], [self.requestLongitudes[0][0], 1439]],
                    False: [[0, self.requestLongitudes[0][1]], [self.requestLongitudes[0][0], 719]]
                }[self.HD]

        # The base URL depends on whether the HD service has been requested.
        self.baseURL = {
            True: 'https://nomads.ncep.noaa.gov/dods/gfs_0p25/',
            False: 'https://nomads.ncep.noaa.gov/dods/gfs_0p50/'
        }[self.HD]

        if debugging:
            log_lev = logging.DEBUG
        else:
            log_lev = logging.WARNING

        logger.setLevel(log_lev)


    def _get_NOAA_REST_url(self, requestVar, requestLongitude, cycle, requestTime):
        """
        Parameters
        ----------
        requestVar : string
            noaa identifier of the variable name:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        requestLongitude : list of int, length 2
            The [start, end] window of longitude for which to get data
            (GFS units)
        cycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast

        returns
        -------
        requestURL : string
            The noaa API request url
        """
        requestURL = '%sgfs%d%02d%02d/gfs_%s_%02dz.ascii?%s[%d:%d][%d:%d][%d:%d][%d:%d]' % (
                self.baseURL,
                cycle.year,
                cycle.month,
                cycle.day,
                {True: '0p25', False: '0p50'}[self.HD],
                cycle.hour,
                requestVar,
                requestTime[0], requestTime[1],
                self.requestAltitude[0], self.requestAltitude[1],
                self.requestLatitude[0], self.requestLatitude[1],
                requestLongitude[0], requestLongitude[1]
            )
        return requestURL

    def _NOAA_request(self, requestVar, cycle, requestTime):
        """
        Parameters
        ----------
        requestVar : string
            noaa identifier of the variable name:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        cycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast
        
        Returns
        -------
        dataResults : list
            list of responses for each request longitude
        """
        dataResults = []

        # Check if we need more than 1 request (ie if we are crossing
        # the Greenwich meridian)
        for requestLongitude in self.requestLongitudes:

            requestURL = self._get_NOAA_REST_url(requestVar, requestLongitude, cycle, requestTime)

            logger.debug('Requesting URL: %s' % requestURL)

            try:
                HTTPresponse = urlopen(requestURL)
                response = HTTPresponse.read().decode('utf-8')
            except:
                logger.exception(
                    'Error while connecting to the GFS server.')
                return
            if response[0] == "<":
                logger.debug("GFS cycle not found.")
                return
            else:
                dataResults.append(response)
        return dataResults

    @profile
    def processNOAARequest(self, requestVar, cycle, requestTime):
        """Downloads data from a NOAA request url, before generating the
        matrix and interpolator

        Allows the response to go out of scope after creating the matrices,
        hence it may decrease the memory overhead compared with downloading
        all items simultaneously

        Parameters
        ----------
        requestVar : string
            noaa identifier of the variable name:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        cycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast

        Returns
        -------
        data_matrix : numpy array
            The downloaded 4D data (empty if the response failed)
        data_map : dict
            the (empty if the response failed)

        """
        response = self._NOAA_request(requestVar, cycle, requestTime)
        if response:
            requestReadableName = self.weatherParameters[requestVar]
            logger.debug('{} data downloaded'.format(
                requestReadableName))

            data_matrix, data_map = self._generate_matrix(response)
        else:
            data_matrix, data_map = {}, {}
        return data_matrix, data_map


    def _NOAA_request_all(self, cycle, requestTime, progressHandler):
        """Requests temperature, altitude and U-V wind direction data from the
        NOAA GFS system for the ranges specified in the class, and for the
        input cycle time and requestTime.

        Parameters
        ----------
        requestLongitudes : list
            list of the longitudes for which to make the request (usually one,
            but two are required around the greenwich median)
        cycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast
        progressHandler : function
            Function that handles the download progress. This is usually
            the member function astra.flight.updateProgress.

        Returns
        -------
        results : dict
            ('noaa_name' : response) pairs, where response is the returned
            data string. Keys:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        """
        results = {}
        progressHandler(0, 1)
        for ivar, requestVar in enumerate(self.weatherParameters.keys()):
            dataResults = self._NOAA_request(requestVar,
                                                 cycle,
                                                 requestTime
                                                 )
            if dataResults:
                progressHandler(1. / len(self.weatherParameters) * (ivar + 1), 1)
                requestReadableName = self.weatherParameters[requestVar]
                logger.debug('{} data downloaded'.format(
                    requestReadableName))
            else:
                return
            results[requestVar] = dataResults
        return results

    def _NOAA_request_all_async(self, cycle, requestTime, progressHandler):
        """Collects all urls for noaa data requests for parameters in
        self.weatherParameters, then submits a combined request asynchronously.

        Offers a concurrent, asynchronous alternative to the standard approach
        for getting urls one at a time, as is done with self._NOAA_request.
        This method uses asyncio. See notes for benchmarking.

        Parameters
        ----------
        requestLongitudes : list
            list of the longitudes for which to make the request (usually one,
            but two are required around the greenwich median)
        cycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast
        progressHandler : function
            Function that handles the download progress. This is usually
            the member function astra.flight.updateProgress.

        Returns
        -------
        result : list of strings
            The (utf-8) decoded response string from the request


        Notes
        -----
        * This method has been benchmarked, and should be approximately 80% than
        the standard method in self._NOAA_request.
        * asyncio was also trialled, for which it was found that the
        performance difference was about the same, on average. asyncio also
        required code that was less readable, and that also hid log messages.
        """
        urls = {}
        progressHandler(0, 1)

        logger.debug('Requesting weather urls asynchronously: status will be sent to requests logger')

        for var in self.weatherParameters.keys():
            urls[var] = [self._get_NOAA_REST_url(var, reqLon, cycle, requestTime)
                for reqLon in self.requestLongitudes]

        progress_increment = 1./sum([len(url_list) for url_list in urls.values()])
        # Have to keep a mutable progress dict to be updated by the hook
        # function:
        progress = [0]
        def increment_progress(reponse, *args, **kwargs):
            # Hook function that updates a progress file via progressHandler.
            # args do nothing, but will be passed by the request hook anyway,
            # so need to do a type of hack that will ignore those inputs
            progress[0] += progress_increment
            logger.debug('Updating Download progress: {}%% complete'.format(100*progress[0]))
            progressHandler(progress[0], 1)

        # Note: currently omiting the hook function, as it doesn't seem to
        # write to the file after the first update
        results = get_urldict_async(urls, hooks_dict={'response': increment_progress})

        return results

    def getNOAAMatricesMapsCycle(self, thisCycle, requestTime, progressHandler):
        """For an input cycle, this function will load all weather variables.

        Multiple attempts will be made for each variable, if the response fails.
        This is, therefore, not the best function to use if it is not known
        beforehand whether the data cycle exists: See Notes.

        Parameters
        ----------
        requestLongitudes : list
            list of the longitudes for which to make the request (usually one,
            but two are required around the greenwich median)
        thisCycle : :obj:`datetime.datetime`
            The cycle datetime for which to obtain the forecast
        requestTime : :obj:`datetime.datetime`
            The launch datetime for which to obtain the forecast
        progressHandler : function
            Function that handles the download progress. This is usually
            the member function astra.flight.updateProgress.

        Returns
        -------
        data_matrices : dict
            ('noaa_name' : data_matrix) pairs, where response is the returned
            data string. Keys:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        data_maps : dict
            keys as in data_matrices, but values instead contain the data_maps
            reutrned by GFS_Handler.processNOAARequest

        Notes
        -----
        GFS_Handler.getNOAAData is the primary download function, and will
        make multiple attempts to find a cycle before referring the download
        of all other variables to this function.

        See Also
        --------
        GFS_Handler.getNOAAData, GFS_Handler.downloadForecast
        """
        data_matrices = {}
        data_maps = {}
        # This is just a flag to make sure that if no data is found, it
        # stops the whole download and re-starts the loop with an older
        # cycle.
        progressHandler(0, 1)

        if self.requestSimultaneous:
            ###############################################################
            # DOWNLOAD DATA FOR ALL PARAMETERS
            if self.use_async:
                results = self._NOAA_request_all_async(thisCycle,
                                                       requestTime,
                                                       progressHandler)

            else:
                results = self._NOAA_request_all(thisCycle,
                                                 requestTime,
                                                 progressHandler)
            if results:
                for ivar, (requestVar, dataResults) in enumerate(results.items()):
                    # Convert the data to matrix and map
                    data_matrices[requestVar], data_maps[requestVar] =\
                        self._generate_matrix(dataResults)

        else:
            for ivar, requestVar in enumerate(self.weatherParameters.keys()):
                # Convert the data to matrix and map and store progress
                for i in range(5):
                    if requestVar not in data_matrices:
                        data_matrix, data_map = \
                            self.processNOAARequest(requestVar, thisCycle,
                                                    requestTime)
                        if data_map:
                            progressHandler(1. / len(self.weatherParameters) * (ivar+1), 1)
                            data_matrices[requestVar], data_maps[requestVar] =\
                                data_matrix, data_map
                        else:
                            raise RuntimeError("'{}' data failed to download for this cycle. Cannot proceed.".format(requestVar))
                            
        # If it got here, the GFS was found: break the outer loop
        return data_matrices, data_maps

    def getNOAAData(self, simulationDateTime, latestCycleDateTime, progressHandler):
        """Makes multiple attempts to find an available data set (cycle) from
        the noaa web service, before downloading all results with
        getNOAAMatricesMapsCycle.

        Parameters
        ----------
        simulationDateTime : :obj:`datetime.datetime`
        latestCycleDateTime : :obj:`datetime.datetime
            The (assumed) earliest available weather forecast date and time.
            Requests will begin searching for valid datasets in reversing
            6 hour intervals from this date time.
        progressHandler : function 
            See simulator.updateProgress

        Returns
        -------
        data_matrices : dict
            ('noaa_name' : data_matrix) pairs, where response is the returned
            data string. Keys:
            'tmpprs': Temperature,
            'hgtprs': 'Altitude',
            'ugrdprs': 'U Winds',
            'vgrdprs': 'V Winds'
        data_maps : dict
            keys as in data_matrices, but values instead contain the data_maps
        """
        #######################################################################
        # TRY TO DOWNLOAD DATA WITH THE LATEST CYCLE. IF NOT AVAILABLE, TRY
        # WITH AN EARLIER ONE

        # pastCycle is a variable that indicates how many cycles in the past
        # we're downloading data from. It first tries with 0, indicating the
        # most recent cycle: this is calculated, knowing that cycles are issued
        # every six hours. Sometimes, however, it takes a few hours for data to
        # become available, or a cycle is skipped. This means that data for the
        # latest cycle is not guaranteed to be available. If data is not
        # available, it tries with 1 cycle older, until one is found with data
        # available. If no cycles are found, the method raises a runtime error.

        for pastCycle in range(25):
            logger.debug('Attempting to download cycle data.')
            thisCycle = latestCycleDateTime - timedelta(hours=pastCycle * 6)
            self.cycleDateTime = thisCycle

            # Initialize time parameter
            timeFromForecast = simulationDateTime - thisCycle
            hoursFromForecast = timeFromForecast.total_seconds() / 3600.

            # GFS time index for the first dataset to be requested
            # (1 GFS index = three hours)
            requestTime = floor(hoursFromForecast / 3.)

            # This stores the actual time of the first dataset downloaded. It's
            # going to be used to convert real time to GFS "time coordinates"
            # (see getGFStime(time) function)
            self.firstAvailableTime = self.cycleDateTime + timedelta(hours=requestTime * 3)

            # Always download an extra time dataset
            # PChambers note: probably +1 because of the index slicing system
            # used on the GFS servers, i.e., times 0:5 will get times
            # 0, 1, 2, 3 and 4
            requestTime = [requestTime, requestTime + ceil(
                self.forecastDuration / 3.) + 1]
            thisCycle = latestCycleDateTime - timedelta(hours=pastCycle * 6)
            self.cycleDateTime = thisCycle

            # Probe the system to see if data is available for this cycle:
            dataResults = self._NOAA_request('tmpprs', thisCycle,
                [requestTime[0], requestTime[0] + 1])

            if (dataResults):
                break
            else:
                logger.debug("Moving to next cycle")

        # Main download
        data_matrices, data_maps = self.getNOAAMatricesMapsCycle(
            thisCycle, requestTime, progressHandler)

        if not (data_matrices and data_maps):
            raise RuntimeError('No available GFS cycles found!')
        return data_matrices, data_maps

    @profile
    def downloadForecast(self, progressHandler=None):
        """Connect to the Global Forecast System and download the closest cycle
        available to the date_time required, for the lon/lat window described
        in this GFS_Handler object.

        The cycle date and time is stored in the object's cycleDateTime
        variable, while the data and interpolation maps are stored in
        [self.]temperatureData, altitudeData, windDirData, windSpeedData,
        temperatureMap, etc.

        Parameters
        ----------
        [progressHandler] : function (default None)
            Function that handles the download progress. This is usually
            the member function astra.flight.updateProgress.

        [use_async] : bool (default True)
            If true, this function will build a series of asynchronous requests

        Returns
        -------
        status : bool
            TRUE if the download was successful, FALSE otherwise.
        """
        # create a null handler if input progressHandler is None:
        if not progressHandler:
            def progressHandler(*args):
                return None

        #######################################################################
        # INITIALIZE TIME HANDLERS AND CALCULATE MOST RECENT CYCLE AVAILABLE

        simulationDateTime = self.launchDateTime
        currentDateTime = datetime.now()

        if simulationDateTime < currentDateTime:
            # Simulating a flight in the past. In this case, use the simulation
            # date time to calculate the cycles
            # that need to be downloaded.
            currentDateTime = simulationDateTime

        # These are the hours at which new forecast cycles are issued
        dailyCycles = [0, 6, 12, 18]
        # Determine which cycle issuing hour is the closest to the current one
        # TODO: Fix this unreadable line: Is it trying to get the earliest
        # before? What about 2300: does it go later on earlier?
        cycleTime = dailyCycles[numpy.digitize([currentDateTime.hour], dailyCycles)[0] - 1]
        latestCycleDateTime = datetime(currentDateTime.year,
                                       currentDateTime.month,
                                       currentDateTime.day,
                                       cycleTime)

        data_matrices, data_maps = self.getNOAAData(simulationDateTime,
            latestCycleDateTime, progressHandler)

#######################################################################
        # PROCESS DATA AND PERFORM CONVERSIONS AS REQUIRED

        # Convert temperatures from Kelvin to Celsius
        data_matrices['tmpprs'] -= 273.15

        # Convert geopotential height to geometric altitude
        altitudeMatrix = (data_matrices['hgtprs'] * earthRadius /
                          (earthRadius - data_matrices['hgtprs']))

        # Convert u and v winds to wind direction and wind speed matrices

        # Convert to KNOTS and the turn into direction and speed
        self.windDirData, self.windSpeedData = tools.uv2dirspeed(
            data_matrices['ugrdprs'],
            data_matrices['vgrdprs'])

        # Store results
        self.temperatureData = data_matrices['tmpprs']
        self.altitudeData = data_matrices['hgtprs']
        self.windSpeedData *= 1.9438445    # Convert to knots

        self.temperatureMap = data_maps['tmpprs']
        self.altitudeMap = data_maps['hgtprs']
        self.windsMap = data_maps['ugrdprs']

        # DOWNLOAD HIGH ALTITUDE 0.5 x 0.5 DATA
        if (self.HD):
            logger.debug('Preparing to download high altitude forecast...')
            self._highAltitudeGFS.downloadForecast()
            self._highAltitudePressure = self._highAltitudeGFS.interpolateData('p')

        logger.debug('Forecast successfully downloaded!')

        return True

    @classmethod
    def fromFiles(cls, fileDict, lat, lon, date_time, HD=False, **kwargs):
        """
        Parameters
        ----------
        fileDict : :obj:`dict`
            A dictionary of noaa_name: filename pairs, indicating which file
            should be used for each noaa variable. The following keys must be
            defined in this dictionary: 'tmpprs' (Temperature), 
            'hgtprs' (Altitude), 'ugrdprs' (U Winds), 'vgrdprs' (V Winds)

        lat : float
            the latitude that the files refer to
        lon : float
            the longitude that the files refer to
        date_time : :obj:`datetime.datetime`
            The datetime that the files refer to
        [HD] : bool (default False)
            Passed to the class cls.

        **kwargs : dict
            The dictionary of keyword: value parameters to be passed to the
            new GFS_Handler. It is the users responsibility to ensure that the
            latitude, longitude, and all other parameters passed to this method
            match those intended for the input data file; there is no obvious
            way to obtain launch site from the noaa files.

        Notes
        -----
        This method should be primarily used for testing and benchmarking

        See Also
        --------
        astra.weather.forecastEnvironment.loadFromNOAAFiles
        """
        # Check that all weather parameters appear in the input dict
        assert(all(k in fileDict for k in GFS_Handler.weatherParameters.keys())),\
            'All NOAA parameter names must appear in input fileDict. See (self).weatherParameters'
        module = cls(lat=lat, lon=lon, date_time=date_time, HD=HD, **kwargs)
        data_matrices = {}
        data_maps = {}

        module.firstAvailableTime = date_time

        for name, fname in fileDict.items():
            requestReadableName = module.weatherParameters[name]
            with open(fname, 'rb') as f:
                dataResults = f.read().decode('utf-8')
            # Convert the data to matrix and map and store progress. Note that
            # this has to be passed as a list to generate_matrix, since it
            # ordinarily allows for multiple longitude data streams if near
            # the greenwich meridian.
            data_matrices[name], data_maps[name] =\
                module._generate_matrix([dataResults])
            logger.debug('{} data processed.'.format(
                requestReadableName))

        #######################################################################
        # PROCESS DATA AND PERFORM CONVERSIONS AS REQUIRED

                # Convert temperatures from Kelvin to Celsius
        data_matrices['tmpprs'] -= 273.15

        # Convert geopotential height to geometric altitude
        altitudeMatrix = (data_matrices['hgtprs'] * earthRadius /
                          (earthRadius - data_matrices['hgtprs']))

        # Convert u and v winds to wind direction and wind speed matrices

        # Convert to KNOTS and the turn into direction and speed
        module.windDirData, module.windSpeedData = tools.uv2dirspeed(
            data_matrices['ugrdprs'],
            data_matrices['vgrdprs'])

        # Store results
        module.temperatureData = data_matrices['tmpprs']
        module.altitudeData = data_matrices['hgtprs']
        module.windSpeedData *= 1.9438445    # Convert to knots

        module.temperatureMap = data_maps['tmpprs']
        module.altitudeMap = data_maps['hgtprs']
        module.windsMap = data_maps['ugrdprs']

        return module

    @profile
    def interpolateData(self, *variables):
        """
        Set up a linear 4d interpolation for each variable given and returns
        it. The returned interpolator can then be used by calling it with
        standard 4D coordinates (lat, lon, alt, time).

        Parameters
        ----------
        variables : tuple of strings
            Acceptable values for *variables can be a combination of any of
                't'|'temp'|'temperature'        temperature interpolation
                'p'|'press'|'pressure'          pressure interpolation
                'd'|'windrct'|'wind_direction'  wind direction interpolation
                's'|'windspd'|'wind_speed'      wind speed interpolation

        Returns
        -------
        list of :obj:`Linear4DInterpolator` objects, in the same order as
        *variables.

        Notes
        -----
        * The linear interpolation is based on the
        interpolate.Linear4DInterpolator class.

        * Warning: the time parameter required by the Linear4DInterpolator
        objects must be in GFS units! Use the getGFStime(time) method for
        conversion.
        """

        results = []

        # Cycle through each variable requested and append the appropriate
        # interpolator to the results list.
        for variable in variables:
            if variable in ('temp', 't', 'temperature'):
                # Interpolate temperature
                if not self.HD:
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
                if not self.HD:
                    results.append(GFS_data_interpolator(self, self.windDirData, self.windsMap.mappingCoordinates))
                else:
                    results.append(GFS_data_interpolator(self, self.windDirData, self.windsMap.mappingCoordinates,
                                                         self._highAltitudeGFS.interpolateData('d')))

            elif variable in ('windspd', 's', 'wind_speed'):
                # Interpolate wind speed
                if not self.HD:
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

        Returns
        -------
        GFStime (scalar) or NaN:
            float corresponding to the converted time, or nan if there
            was no downloaded data found.
        """

        # Check if data is available. If it isn't return NaN.
        if self.firstAvailableTime is None or self.altitudeMap is None:
            return float('nan')
        else:
            try:
                # Try to define the time from the first dataset.
                # If it fails for any reasons, return NaN.
                timeFromFirstDataset = time - self.firstAvailableTime
                return self.altitudeMap.fwdTime[0] + timeFromFirstDataset.days\
                    + timeFromFirstDataset.seconds / 3600. / 24.
            except TypeError:
                return float('nan')


    def _generate_matrix(self, dataStreams):
        """
        Generates data matrices and coordinates mapping.

        Called by the downloadForecast() method when data is downloaded, and
        supports joining together two datasets. This method should NOT be used
        alone.
        """

        overallResults = []
        overallMaps = []

        # Run this either once or twice, according to how many datasets have
        # been downloaded (Greenwich meridian crossing.
        for dataStream in dataStreams:
            dataLines = dataStream.split('\n')

            # Count how many latitude, longitude, pressure and time points are
            # available in the datastream. This is used to initialize the
            # results matrix.
            timePoints = int(dataLines[0].split()[1].split('][')[0][1:])
            pressurePoints = int(dataLines[0].split()[1].split('][')[1])
            latitudePoints = int(dataLines[0].split()[1].split('][')[2])
            longitudePoints = int(dataLines[0].split()[1].split('][')[3][:-1])
            totalPoints = timePoints * pressurePoints * latitudePoints * longitudePoints

            # Initialize the matrix
            results = numpy.zeros(totalPoints).reshape((latitudePoints,
                longitudePoints, pressurePoints, timePoints))

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
                results[latitudeIndex, :, pressureIndex, timeIndex] = [
                    float(x) if float(x) < 1e8 else 0 for x in line.split(',')[1:]]

            # Generate the mapping. This is an object containing mapping
            # information between GFS indices for lat,lon,press, time and
            # actual values, and viceversa. "Forward" maps are GFS indices ->
            # real values and are LISTS. "Reverse" maps are real values -> GFS
            # indices and are DICTIONARIES. These are used to be able to find a
            # data value given real coordinates and for interpolation.

            resultsMap = GFS_Map()

            resultsMap.fwdLatitude = [float(lat) for lat in dataLines[-4].split(',')]
            resultsMap.revLatitude = {lat: ind for (lat, ind) in
                                      zip(resultsMap.fwdLatitude,
                                          range(len(resultsMap.fwdLatitude)))}

            resultsMap.fwdLongitude = [float(lon) - 360 if float(lon) > 180 else float(lon) for lon in
                                       dataLines[-2].split(',')]
            resultsMap.revLongitude = {lon: ind for (lon, ind) in
                                       zip(resultsMap.fwdLongitude, 
                                           range(len(resultsMap.fwdLongitude)))}

            resultsMap.fwdPressure = [float(press) for press in dataLines[-6].split(',')]
            resultsMap.revPressure = {press: ind for (press, ind) in
                                      zip(resultsMap.fwdPressure, range(
                                          len(resultsMap.fwdPressure)))}

            resultsMap.fwdTime = [float(time) for time in dataLines[-8].split(',')]
            resultsMap.revTime = {time: ind for (time, ind) in zip(
                resultsMap.fwdTime, range(len(resultsMap.fwdTime)))}

            overallResults.append(results)
            overallMaps.append(resultsMap)

        # Join results together, if needed
        if len(overallResults) > 1:
            # Identify which one is to the right
            isFirstOneEast = overallMaps[0].fwdLongitude[0] == 0
            if isFirstOneEast:
                results = numpy.concatenate((overallResults[1],
                    overallResults[0]), axis=1)
                resultsMap = overallMaps[1].rjoin(overallMaps[0])
            else:
                results = numpy.concatenate((overallResults[0],
                    overallResults[1]), axis=1)
                resultsMap = overallMaps[0].rjoin(overallMaps[1])
        else:
            results = overallResults[0]
            resultsMap = overallMaps[0]

        # Store all the forward maps to a single list called mappingCoordinates.
        # This is used by interpolateData(*variables) to correctly initialize
        # the interpolators.
        resultsMap.mapCoordinates()

        return results, resultsMap


    def _pressure_interpolator(self, lat, lon, alt, time):
        """
        Extracts pressure from (lat,lon,alt,time) coordinates.

        This is essentially a 3D interpolation of the column of altitudes
        (ie not just a single value) at the point defined just by latitude,
        longitude and time. Once the column of altitudes is extracted, a
        UnivariateSpline (1D) interpolation is run between that column of
        altitudes and the corresponding pressures (which are always constant).
        Then, root at which the altitude corresponds to the inputs is found.

        If the requested point is outside latitude, longitude or time bounds,
        the nearest value available for the out-of-bounds axis is used.
        If the requested point is outside altitude bounds, ISA conditions are
        used above the upper limit and the nearest value is used below the
        lower limit.

        Notes
        -----
        * This is not using the standard interpolate.Linear4DInterpolator
        because it's NOT a 4D interpolation.
        """

        # Clip out-of-bounds coordinates and limit them
        lat = numpy.clip(lat, self.altitudeMap.fwdLatitude[0],
            self.altitudeMap.fwdLatitude[-1])
        lon = numpy.clip(lon, self.altitudeMap.fwdLongitude[0],
            self.altitudeMap.fwdLongitude[-1])
        time = numpy.clip(time, self.altitudeMap.fwdTime[0],
            self.altitudeMap.fwdTime[-1])

        # Find closest indices and subtract 1 if the upper limit is being
        # reached, to avoid a KeyError
        i = numpy.digitize([lat], self.altitudeMap.fwdLatitude)[0]
        if i == len(self.altitudeMap.fwdLatitude):
            i -= 1
        idxLat = [i - 1, i]
        i = numpy.digitize([time], self.altitudeMap.fwdTime)[0]
        if i == len(self.altitudeMap.fwdTime):
            i -= 1
        idxTime = [i - 1, i]

        # Reverse mapping for longitude, to account for crossing the 180th
        # meridian
        # Note: this method is used for two reasons:
        #   1. It's faster than numpy.digitize
        #   2. numpy.digitize fails at longitudes between -180 and 179.5
        #      (for HD), since there isn't an entry for -180.
        lonGrid = [floor(lon / self.lonStep) * self.lonStep,
                   ceil(lon / self.lonStep) * self.lonStep]
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
            idxLon = [self.altitudeMap.revLongitude[lonGrid[0]],
                      self.altitudeMap.revLongitude[lonGrid[1]]]

        try:
            fracLat = 1 - (abs((lat - self.altitudeMap.fwdLatitude[idxLat[0]]) /
                (self.altitudeMap.fwdLatitude[idxLat[1]] -
                    self.altitudeMap.fwdLatitude[idxLat[0]]))
            )
            fracLon = 1 - (abs((lon - self.altitudeMap.fwdLongitude[idxLon[0]]) /
                (self.altitudeMap.fwdLongitude[idxLon[1]] - 
                    self.altitudeMap.fwdLongitude[idxLon[0]]))
            )
            fracTime = 1 - (abs((time - self.altitudeMap.fwdTime[idxTime[0]]) /
                (self.altitudeMap.fwdTime[idxTime[1]] -
                    self.altitudeMap.fwdTime[idxTime[0]]))
            )
        except IndexError:
            # Requested point outside bounds

            logger.error('An error occurred while clipping coordinate points.')
            raise

        alt00 = fracLat * self.altitudeData[idxLat[0], idxLon[0], :, idxTime[0]]\
            + (1 - fracLat) * self.altitudeData[idxLat[1], idxLon[0], :, idxTime[0]]
        alt01 = fracLat * self.altitudeData[idxLat[0], idxLon[0], :, idxTime[1]]\
            + (1 - fracLat) * self.altitudeData[idxLat[1], idxLon[0], :, idxTime[1]]
        alt10 = fracLat * self.altitudeData[idxLat[0], idxLon[1], :, idxTime[0]]\
            + (1 - fracLat) * self.altitudeData[idxLat[1], idxLon[1], :, idxTime[0]]
        alt11 = fracLat * self.altitudeData[idxLat[0], idxLon[1], :, idxTime[1]]\
            + (1 - fracLat) * self.altitudeData[idxLat[1], idxLon[1], :, idxTime[1]]
        alt_time0 = fracLon * alt00 + (1 - fracLon) * alt10
        alt_time1 = fracLon * alt01 + (1 - fracLon) * alt11
        altitude = fracTime * alt_time0 + (1 - fracTime) * alt_time1 - alt

        if altitude.min() > 0:
            # NEAREST NEIGHBOR: if requested point is below minimum altitude,
            # return lowest point available
            return self.altitudeMap.fwdPressure[0]
        elif altitude.max() < 0:
            # NEAREST NEIGHBOR: if requested point is above maximum altitude,
            # return highest point available
            if not self.HD:
                return self.altitudeMap.fwdPressure[-1]
            else:
                return self._highAltitudePressure(lat, lon, alt, time)
        else:
            # LINEAR INTERPOLATION
            # (see method documentation for details)
            f = UnivariateSpline(self.altitudeMap.fwdPressure[::-1],
                                 altitude[::-1], s=0)
            return f.roots()[0]


# TODO: Migrate this class back to GFS_Handler - there is no need for an
# additional class here, we just need to switch on standard definition settings
# for the high altitude segment of the GFS_Handler downloadForecast
class GFS_High_Altitude_Handler(GFS_Handler):
    def __init__(self, lat, lon, date_time, forecastDuration=4,
        debugging=False, log_to_file=False):

        super(GFS_High_Altitude_Handler, self).__init__(lat=lat,
            lon=lon,
            date_time=date_time,
            forecastDuration=forecastDuration,
            debugging=debugging,
            HD=False)

        # Set the altitude to the higher levels, since the lower altitude
        # handler covers most of the data
        self.requestAltitude = [41, 46]


class GFS_Map(object):
    """
    Private class used to store 4D mapping data for a specific parameter.

    Note: the GFS_Map class is only used for internal data management only and
    should not be used on its own.

    Attributes
    ----------
    fwdLatitude : list
        Forward list for mapping GFS indices to real world coordinates 
    fwdLongitude : list
        Forward list for mapping GFS indices to real world coordinates 
    fwdPressure : list
        Forward list for mapping GFS indices to real world coordinates 
    fwdTime : list
        Forward list for mapping GFS indices to real world coordinates 
    revLatitude : dict
        Forward list for mapping GFS indices to real world coordinates 
    revLongitude : dict
        Forward list for mapping GFS indices to real world coordinates 
    revPressure : dict
        Forward list for mapping GFS indices to real world coordinates 
    revTime : dict
        Forward list for mapping GFS indices to real world coordinates 
    mappingCoordinates : list
        4D real world Interpolation coordinates for each value in the data
        matrix
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

        # Interpolation coordinates is a list of 4D real world coordinates for
        # each value in the data matrix
        self.mappingCoordinates = []

    def rjoin(self, data_map):
        """This method joins and returns this map and the one passed in
        data_map, concatenating it to the RIGHT."""
        lonLen = len(self.fwdLongitude)

        if self.fwdLatitude == data_map.fwdLatitude and self.fwdPressure == data_map.fwdPressure and self.fwdTime == data_map.fwdTime:
            # Join
            self.fwdLongitude = list(numpy.concatenate((self.fwdLongitude, data_map.fwdLongitude)))

        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        if self.revLatitude == data_map.revLatitude and self.revPressure == data_map.revPressure and self.revTime == data_map.revTime:
            # Join
            self.revLongitude.update({key: value + lonLen for (key, value) in data_map.revLongitude.items()})

        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        return self

    def ljoin(self, data_map):
        """This method joins and returns this map and the one passed in
        data_map, concatenating it to the LEFT."""
        lonLen = len(data_map.fwdLongitude)

        if self.fwdLatitude == data_map.fwdLatitude and self.fwdPressure == data_map.fwdPressure and self.fwdTime == data_map.fwdTime:
            self.fwdLongitude = list(numpy.concatenate((data_map.fwdLongitude, self.fwdLongitude)))
        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        if self.revLatitude == data_map.revLatitude and self.revPressure == data_map.revPressure and self.revTime == data_map.revTime:
            self.revLongitude = {key: value + lonLen for (key, value) in self.revLongitude.items()}.update(
                data_map.revLongitude)
        else:
            logger.error("Map joining failed. Latitudes, Pressures and/or Times don't match between the two maps!")

        return self

    def mapCoordinates(self):
        """Prepare the mappingCoordinates variable by putting together all
        forward and reverse maps. 

        Used by the Linear4DInterpolator to access interpolation data maps.
        """
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


class GFS_data_interpolator(object):
    """
    Private class used by GFS_Handler to interpolate data.

    This class acts as an interface to be able to use the correct coordinates
    (lat,lon,alt,time) when requesting interpolated data.

    Once initialized, the object can be called with the requested coordinates.
    At this point, the object converts the altitude to pressure using the
    GFS_Handler's _pressure_interpolator(...)
    method and then uses that to interpolate the required data and return it.

    Parameters
    ----------
    GFS_Handler : :obj:`GFS_Handler
        This is used to fetch pressure data upon request
    data : numpy array (4D)
        matrix containing all the data that needs to be interpolated
    map : list
        The map of the data. This should be the mappingCoordinates variable
        generated by the GFS_Map.mapCoordinates() method.
    """

    @profile
    def __init__(self, GFS_Handler, data, dmap, high_alt_interpolator=None,
        min_pressure=None):
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
