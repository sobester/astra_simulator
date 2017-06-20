# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-06-19 16:45:06
from .simulator import flight, flightProfile
from .weather import forecastEnvironment
from datetime import datetime, timedelta
import logging
import operator
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from bisect import bisect_right
from copy import deepcopy
from operator import eq
import functools
from . import global_tools as tools
from .flight_tools import nozzleLiftFixedAscent
import matplotlib.dates as mdates


logger = logging.getLogger(__name__)


class HallOfFame(object):
    """Class copied and edited from the deap python library.

    The hall of fame contains the best individual that ever lived in the
    population during the evolution. It is lexicographically sorted at all
    time so that the first element of the hall of fame is the individual that
    has the best first distanceFromTarget value ever seen, according to the weights
    provided to the distanceFromTarget at creation time.
    
    The insertion is made so that old individuals have priority on new
    individuals. A single copy of each individual is kept at all time, the
    equivalence between two individuals is made by the operator passed to the
    *similar* argument.
    :param maxsize: The maximum number of individual to keep in the hall of
                    fame.
    :param similar: An equivalence operator between two individuals, optional.
                    It defaults to operator :func:`operator.eq`.
    
    The class :class:`HallOfFame` provides an interface similar to a list
    (without being one completely). It is possible to retrieve its length, to
    iterate on it forward and backward and to get an item or a slice from it.
    """
    def __init__(self, maxsize=10, similar=eq):
        self.maxsize = maxsize
        self.keys = list()
        self.items = list()
        self.similar = similar
    
    def update(self, profile):
        """Update the hall of fame with the *population* by replacing the
        worst individuals in it by the best individuals present in
        *population* (if they are better). The size of the hall of fame is
        kept constant.
        
        :param population: A list of individual with a distanceFromTarget attribute to
                           update the hall of fame with.
        """
        if len(self) == 0 and self.maxsize !=0:
            # Working on an empty hall of fame is problematic for the
            # "for else"
            self.insert(profile)
        
        if profile.objective < self[-1].objective or len(self) < self.maxsize:
            for hofer in self:
                # Loop through the hall of fame to check for any
                # similar individual
                if self.similar(profile, hofer):
                    break
            else:
                # The profile is unique and strictly better than
                # the worst
                if len(self) >= self.maxsize:
                    self.remove(-1)
                self.insert(profile)
    
    def insert(self, item):
        """Insert a new individual in the hall of fame using the
        :func:`~bisect.bisect_right` function. The inserted individual is
        inserted on the right side of an equal individual. Inserting a new 
        individual in the hall of fame also preserve the hall of fame's order.
        This method **does not** check for the size of the hall of fame, in a
        way that inserting a new individual in a full hall of fame will not
        remove the worst individual to maintain a constant size.
        
        :param item: The individual with a distanceFromTarget attribute to insert in the
                     hall of fame.
        """
        item = deepcopy(item)
        i = bisect_right(self.keys, item.objective)
        self.items.insert(len(self) - i, item)
        self.keys.insert(i, item.objective)
    
    def remove(self, index):
        """Remove the specified *index* from the hall of fame.
        
        :param index: An integer giving which item to remove.
        """
        del self.keys[len(self) - (index % len(self) + 1)]
        del self.items[index]
    
    def clear(self):
        """Clear the hall of fame."""
        del self.items[:]
        del self.keys[:]

    def __len__(self):
        return len(self.items)

    def __getitem__(self, i):
        return self.items[i]

    def __iter__(self):
        return iter(self.items)

    def __reversed__(self):
        return reversed(self.items)
    
    def __str__(self):
        return str(self.items)


class targetFlight(flight):
    """Class for estimating the required launch parameters (date/time, nozzle
    lift, floating vs standard, cutdown altitude, etc.) for minimal distance
    from a target landing site.

    Parameters
    ----------
    start_dateTime : :obj:`datetime.datetime`
        The date and time of the start of the window for which to test flights
    targetLat, targetLon : scalar
        The target landing latitude/longitude
    launchSites : list of lat, lon, elev triplets
        lat, lon and elev triplets for each launch site 
    nozzleLift : scalar
        The nozzle lift (in kg). This class currently uses a static nozzle
        lift, and searches time only. Future versions will also search the
        nozzle lift parameter.
    windowDuration : int
        The duration for which to download weather data (in hours). Flights
        will be sampled within this range.
    N_Pareto : int
        The number of 

    See Also
    --------
    :obj:`astra.simulator.flight`
    """

    def __init__(self,
                 start_dateTime,
                 targetLat,
                 targetLon,
                 targetElev,
                 launchSites,
                 balloonGasType,
                 balloonModel,
                 nozzleLift,
                 payloadTrainWeight,
                 inflationTemperature,
                 N_Pareto=10,
                 windowDuration=24,
                 requestSimultaneous=False,
                 HD=False,
                 maxFlightTime=18000,
                 parachuteModel=None,
                 trainEquivSphereDiam=0.1,
                 floatingFlight=False,
                 floatingAltitude=None,
                 ventingStart=1000,
                 excessPressureCoeff=1.,
                 outputFile='',
                 debugging=False,
                 log_to_file=False,
                 progress_to_file=False,
                 launchSiteForecasts=[]):
        self.targetLat = targetLat
        self.targetLon = targetLon
        self.targetElev = targetElev
        self.start_dateTime = start_dateTime
        self.windowDuration = windowDuration
        self.launchSites = launchSites

        self.end_dateTime = self.start_dateTime + timedelta(hours=windowDuration)
        self.bestProfile = None
        self.launchSiteForecasts = []

        super(targetFlight, self).__init__(environment=None,
                                balloonGasType=balloonGasType,
                                balloonModel=balloonModel,
                                nozzleLift=nozzleLift,
                                payloadTrainWeight=payloadTrainWeight, 
                                maxFlightTime=maxFlightTime,
                                parachuteModel=parachuteModel,
                                numberOfSimRuns=1,      # This uses mean params (not MC)
                                trainEquivSphereDiam=trainEquivSphereDiam,
                                floatingFlight=False,   # For now, only allow bursting flights
                                ventingStart=ventingStart,
                                excessPressureCoeff=excessPressureCoeff,
                                debugging=False,
                                log_to_file=False,
                                progress_to_file=False)

        # Use a maximum lateral upper atmospheric speed to determine if this
        # flight is feasible, given an 'as the crow flies' flight entirely in
        # the jet stream, assuming max speed of 200 km/hr
        self.maxLateralSpeed = 200    # m/s 
        self.cutoffDistance = self.maxLateralSpeed * self.maxFlightTime / 3600.

        # This section could be made faster if launch sites are clustered into
        # single forecasts (and therefore, less download requests)
        if launchSiteForecasts:
            self.launchSiteForecasts = launchSiteForecasts
        else:
            self.launchSiteForecasts = []
            for site in launchSites:
                if tools.haversine(site[0], site[1], targetLat, targetLon) < self.cutoffDistance:
                    self.launchSiteForecasts.append(forecastEnvironment(site[0], site[1], site[2], start_dateTime,
                                                inflationTemperature=inflationTemperature,
                                                forceNonHD=(not HD),
                                                forecastDuration=windowDuration,
                                                debugging=debugging,
                                                progressHandler=None,
                                                load_on_init=False,
                                                requestSimultaneous=requestSimultaneous,
                                                )
                    )
                else:
                    logger.warning("Launch site lat {}, lon {} is too far from the landing site.".format(
                        site[0], site[1]))
                    # Need to log progress to file for the web interface here

        # Module level loggers seem to cause issues for derived classes spread
        # accross modules: trying to fix this here
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

        # Set some max a min bounds on ascent rate
        self.minAscentRate = 1.5
        self.maxAscentRate = 6

        # Define a list of input vectors Xs (to the target distance function),
        # and a list of objective scores Ys for plotting later:
        self.Xs = []


    # @functools.wraps(flight.fly)
    # def fly(self, flightNumber, launchDateTime, runPreflight=True,
    #     profileClass=targetFlightProfile):
        
    #     See flight class fly method for docs.

    #     Default profile class is changed in this method, to allow self.results
    #     to be populated with profiles that contain information about their
    #     distance from the target, and the input vector X passed to
    #     self.targetDistance.
        
    #     super(targetFlight, self).fly(flightNumber, launchDateTime,
    #         runPreflight=True, profileClass=profileClass)

    def targetDistance(self, X, appendResult=False):
        """
        Parameters
        ----------
        X : iterable
            The input vector. Should contain elements,
                x[0] :  time, in hours, after the start of the optimisation
                window (self.startTime)

        appendResult : 
        
        returns
        -------
        distance : scalar
            The euclidean norm of [target_lon - lon, target_lat - lat] (degrees)
        """
        t = X[0]
        nozzleLift = X[1]

        logger.debug("Running flight for datetime {}, nozzleLift={}kg".format(
                        self.start_dateTime + timedelta(hours=t), nozzleLift))

        # nozzle lift is an access controlled variable, which may raise an error if
        # lower than the payload lift: return a penalty if this is the case
        try:
            self.nozzleLift = nozzleLift
        except ValueError:
            # Pay a penalty larger than the maximum possible distance (Earth circumference)
            return 5e6

        launchDateTime = self.start_dateTime + timedelta(hours=t)

        resultProfile, solution = self.fly(0, launchDateTime)
        
        landing_lat = resultProfile.latitudeProfile[-1]
        landing_lon = resultProfile.longitudeProfile[-1]
        dist = tools.haversine(landing_lat, landing_lon, self.targetLat, self.targetLon)

        # Storing attributes on the fly is bad from a software engineering
        # standpoint, but it works:
        resultProfile.X = X
        resultProfile.distanceFromTarget = dist

        # Store normalised objective (negative) as this is a minimization. Used
        # in the hall of fame production and sorting.
        resultProfile.objective = - dist

        # Use the hall of fame to store this profile
        self.results.update(resultProfile)

        return dist

    def _callbackStoreResult(self, xk, convergence):
        self.Xs.append(xk)

    def bruteForce(self):
        """ Sample the parameter space and form a map of the landing site
        distance landscape

        Currently only considers time
        """
        # Number of x-y points to use
        Nx = int(self.windowDuration / 4) + 1    # run 4 hourly calculations
        Ny = 21

        # Keep all paths for now, to visualize the landscape of nozzle lift
        self.results = HallOfFame(maxsize=Nx * Ny)

        dts = np.linspace(0, self.windowDuration, Nx)

        distances = np.zeros([Nx, Ny])

        dateTimeVec = [self.start_dateTime + timedelta(hours=dt) for dt in dts]


        for launchSiteForecast in self.launchSiteForecasts:
            self.environment = launchSiteForecast
            # Estimated maximum bound of nozzle lift for target ascent rates
            self.nozzleLiftLowerBound = nozzleLiftFixedAscent(self.minAscentRate,
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.)
            self.nozzleLiftUpperBound = nozzleLiftFixedAscent(self.maxAscentRate,
                    self._balloonWeight, self.payloadTrainWeight,
                    self.environment.inflationTemperature,
                    self.environment.getPressure(self.launchSiteLat,
                                                      self.launchSiteLon,
                                                      self.launchSiteElev,
                                                      self.start_dateTime),
                    self._gasMolecularMass, self.excessPressureCoeff,
                    CD=(0.225 + 0.425)/2.)


            nozzleLiftVec = np.linspace(self.nozzleLiftLowerBound, self.nozzleLiftUpperBound, Ny)
            X, Y = np.meshgrid(dts, nozzleLiftVec)

            logger.info("Date range: [{}, {}], Nx={} points".format(
                dateTimeVec[0], dateTimeVec[-1], Nx))
            logger.info("Nozzle Lift range: [{}, {}] (kg), Ny={} points".format(
                self.nozzleLiftLowerBound, self.nozzleLiftUpperBound, Ny))

            logger.debug("Running brute force calculation")
            for i, t in enumerate(dts):
                for j, nozzleLift in enumerate(nozzleLiftVec):
                    # brute force approach
                    distance = self.targetDistance([t, nozzleLift])
                    # distance_lift_vec[k] = distance
                    distances[i, j] = distance

        bestProfile = self.results[0]
        self.bestProfile = bestProfile
        return bestProfile, dateTimeVec, nozzleLiftVec, distances

    def optimizeTargetLandingSite(self, method='Nelder-Mead', maxsize=10, **kwargs):
        """
        Parameters
        ----------
        method : string
            'brute-force' : Sample flight every hour within the time window
                *Note: this probably isn't suitable given more than two params
            'Nelder-Mead' : Use scipy's Nelder mead pattern search. This has
                proven to be sensitive to initial guesses for this application.
            'DE' : Use scipy differential evolution.

        [maxsize] : int (default 10)
            The number of (lexographically sorted by distance from target site)
            profiles to keep in self.results

        **kwargs - extra args to pass to scipy
        """
        # run the simulation every hour over the time window. Note that we need
        # to also get weather to cover a flight of duration <maxFlightTime>,
        # starting at the end of the window.
        # scipy.optimize.fmin_ 

        self.results = HallOfFame(maxsize=10)
        self.Xs = []

        # distance_map = {}
        # # For now, assume the first is the only interesting launch site
        # Estimated maximum bound of nozzle lift for target ascent rates
        self.environment = self.launchSiteForecasts[0]

        self.nozzleLiftLowerBound = nozzleLiftFixedAscent(self.minAscentRate,
            self._balloonWeight, self.payloadTrainWeight,
            self.environment.inflationTemperature,
            self.environment.getPressure(self.launchSiteLat,
                                              self.launchSiteLon,
                                              self.launchSiteElev,
                                              self.start_dateTime),
            self._gasMolecularMass, self.excessPressureCoeff,
            CD=(0.225 + 0.425)/2.)
        self.nozzleLiftUpperBound = nozzleLiftFixedAscent(self.maxAscentRate,
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.)

        if method in ('Nelder-Mead', 'L-BFGS-B'):
            try:
                x0 = kwargs.pop('x0')
            except KeyError:
                logger.info('No Initial guess x0 provided. Using half of the windowDuration.')
                x0 = 0.5 * self.windowDuration
            res = opt.minimize(self.targetDistance, x0=x0,
                                          method=method,
                                          callback=(lambda x: self._callbackStoreResult(x, convergence=None)),
                                          **kwargs)
            bestProfile = self.results[0]
            self.bestProfile = bestProfile
        elif method == 'DE':
            res = opt.differential_evolution(self.targetDistance, 
                bounds=[(0, self.windowDuration - self.maxFlightTime/3600.),
                        (self.nozzleLiftLowerBound, self.nozzleLiftUpperBound)],
                                             callback=self._callbackStoreResult,
                                             **kwargs)
            bestProfile = self.results[0]
            self.bestProfile = bestProfile

        else:
            raise ValueError('No known optimization method for {}'.format(method))
        return res

    # def plotBruteForceComparison(self):
    #     """
    #     """
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)

    #     plt.axis('equal')
    #     dateTimes = []
    #     distances = []


    #     raise NotImplementedError()


    def plotPaths3D(self, fig=None, ax=None):
        """Plots the resulting trajectories contained in self.results, along
        with target 

        This method should be used after calls to self.targetDistance, or
        any method in optimizeTargetLandingSite, since these will populate
        self.results
        """
        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111, projection='3d')

        plt.axis('equal')

        for profile in self.results:
            lat_arr, lon_arr, alt_arr = profile.latitudeProfile, profile.longitudeProfile, profile.altitudeProfile
            ax.plot(lat_arr, lon_arr, alt_arr, 'b-')

        best_lat, best_lon, best_alt = self.bestProfile.latitudeProfile, self.bestProfile.longitudeProfile,self.bestProfile.altitudeProfile
        ax.plot(best_lat, best_lon, best_alt, 'k-', label='Best')
        ax.plot([self.targetLat], [self.targetLon], [self.targetElev], 'gx', label='target')
        ax.legend(loc='lower left')
        ax.set_xlabel('Lat (deg)')
        ax.set_ylabel('Lon (deg)')
        ax.set_zlabel('Alt (km)')
        fig.show()
        return fig, ax

    def plotLandingSites(self, fig=None, ax=None, landingMarker='bx'):
        """
        """
        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111, aspect='equal')
            ax.set_xlabel('Lat (deg)')
            ax.set_ylabel('Lon (deg)')


        lats = [prof.latitudeProfile[-1] for prof in self.results]
        lons = [prof.longitudeProfile[-1] for prof in self.results]

        ax.plot(lats, lons, landingMarker, label='Trajectory simulations')
            
        best_lat, best_lon, best_alt = self.bestProfile.latitudeProfile, self.bestProfile.longitudeProfile,self.bestProfile.altitudeProfile
        m = ax.plot(self.targetLat, self.targetLon, 'gx', label='Target')
        ax.plot(self.launchSites[0][0], self.launchSites[0][1], 'go', label='Launch Site')
        ax.plot(best_lat[-1], best_lon[-1], 'kx', label='Best')
            
        legend = ax.legend(loc='lower left')

        fig.show()
        return fig, ax

    def plotLandingSiteDistanceContours(self, dateTimeVector, nozzleLiftVector, distances,
        fig=None, ax=None, bestMarker='b*', appendLabel=''):
        """Plots the ground distance of the landing sites contained in
        self.results from the target landing site.

        Requires a call to any of the optimisation of brute force calculations
        to populate self.results.
        """
        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111)
            ax.set_xlabel('Date and Time')
            ax.set_ylabel('Nozzle Lift')
            xtick_locator = mdates.AutoDateLocator()
            xtick_formatter = mdates.AutoDateFormatter(xtick_locator)
            ax.xaxis.set_major_locator(xtick_locator)
            ax.xaxis.set_major_formatter(xtick_formatter)
            # ax.xaxis.set_minor_locator(hours)
            date_start = self.start_dateTime
            date_end = (self.start_dateTime + timedelta(hours=self.windowDuration))
            ax.set_xlim(date_start, date_end)
            fig.autofmt_xdate()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0,
                             box.width, box.height * 0.9])
            
        CS = ax.contour(dateTimeVector, nozzleLiftVector, distances, label='Landing sites'+appendLabel)
        ax.plot(self.bestProfile.launchDateTime, self.bestProfile.X[1], bestMarker, markersize=8, label='min' + appendLabel)

        fig.show()
        legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2)
        ax.clabel(CS, inline=1, fontsize=12)
        return fig, ax

    def plotLandingSiteDistances(self, fig=None, ax=None, marker='bx', bestMarker='b*', appendLabel=''):
        """Plots the ground distance of the landing sites contained in
        self.results from the target landing site.

        Requires a call to any of the optimisation of brute force calculations
        to populate self.results.
        """
        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111)
            ax.set_xlabel('Date and Time')
            ax.set_ylabel('Nozzle Lift')
            xtick_locator = mdates.AutoDateLocator()
            xtick_formatter = mdates.AutoDateFormatter(xtick_locator)
            ax.xaxis.set_major_locator(xtick_locator)
            ax.xaxis.set_major_formatter(xtick_formatter)
            # ax.xaxis.set_minor_locator(hours)
            date_start = self.start_dateTime
            date_end = (self.start_dateTime + timedelta(hours=self.windowDuration))
            ax.set_xlim(date_start, date_end)
            fig.autofmt_xdate()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0,
                             box.width, box.height * 0.9])
 
        dts, nozzleLifts = zip(*self.Xs)
        dateTimes = [self.start_dateTime + timedelta(hours=dt) for dt in dts]
        ax.plot(dateTimes, nozzleLifts, marker, label='Landing sites'+appendLabel)
        ax.plot(self.bestProfile.launchDateTime, self.bestProfile.X[1], bestMarker, markersize=8, label='min' + appendLabel)

        fig.show()
        legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2)
        return fig, ax
