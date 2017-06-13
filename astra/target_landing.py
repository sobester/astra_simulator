# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   Paul Chambers
# @Last Modified time: 2017-06-13 14:27:12
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

        # This section could be made faster if launch sites are clustered into
        # single forecasts (and therefore, less download requests)
        if launchSiteForecasts:
            self.launchSiteForecasts = launchSiteForecasts
        else:
            self.launchSiteForecasts = [forecastEnvironment(site[0], site[1], site[2], start_dateTime,
                                        inflationTemperature=inflationTemperature,
                                        forceNonHD=(not HD),
                                        forecastDuration=windowDuration,
                                        debugging=False,
                                        progressHandler=None,
                                        load_on_init=False
                                        ) for site in launchSites]

    # self.simflight = flight(balloonGasType,
    #                         balloonModel,
    #                         nozzleLift,
    #                         payloadTrainWeight, 
    #                         HD=HD,
    #                         maxFlightTime=maxFlightTime,
    #                         parachuteModel=parachuteModel,
    #                         numberOfSimRuns=1,      # This uses mean params (not MC)
    #                         trainEquivSphereDiam=trainEquivSphereDiam,
    #                         floatingFlight=False,   # For now, only allow bursting flights
    #                         ventingStart=ventingStart,
    #                         excessPressureCoeff=excessPressureCoeff,
    #                         debugging=False,
    #                         log_to_file=False,
    #                         progress_to_file=False)
                            
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
        t = int(X[0])
        launchDateTime = self.start_dateTime + timedelta(hours=t)
        resultProfile, solution = self.fly(0, launchDateTime)
        resultProfile.X = X

        landing_lat = resultProfile.latitudeProfile[-1]
        landing_lon = resultProfile.longitudeProfile[-1]
        dist = np.linalg.norm(np.array([self.targetLat, self.targetLon]) - np.array([landing_lat, landing_lon]))
        resultProfile.distanceFromTarget = dist
        # Store normalised objective (negative) as this is a minimization. Used
        # in the hall of fame production and sorting.
        resultProfile.objective = - dist

        # Use the hall of fame to store this profile
        self.results.update(resultProfile)

        return dist

    def bruteForce(self):
        """ Sample the parameter space and form a map of the landing site
        distance landscape

        Currently only considers time
        """
        dts = np.arange(self.windowDuration)

        # ensure the hall of fame remembers all individuals profiles
        self.results = HallOfFame(maxsize=len(dts))

        dateTimeVec = self.start_dateTime + dts * timedelta(hours=1)

        for i, launchSiteForecast in enumerate(self.launchSiteForecasts):
            self.environment = launchSiteForecast

            for j, t in enumerate(dts):
                # distance_lift_vec = np.zeros(np.length(nozzelLift_Vector))
                # brute force approach
                dtime = dateTimeVec[j]
                self.targetDistance([t])
                # distance_lift_vec[k] = distance

        bestProfile = self.results[0]
        self.bestProfile = bestProfile
        return bestProfile

    def optimizeTargetLandingSite(self, method='Nelder-Mead', **kwargs):
        """
        Parameters
        ----------
        method : string
            'brute-force' : Sample flight every hour within the time window
                *Note: this probably isn't suitable given more than two params
            'Nelder-Mead' : Use scipy's Nelder mead pattern search. This has
                proven to be sensitive to initial guesses for this application.
            'DE' : Use scipy differential evolution.

        **kwargs - extra args to pass to scipy
        """
        # run the simulation every hour over the time window. Note that we need
        # to also get weather to cover a flight of duration <maxFlightTime>,
        # starting at the end of the window.
        # scipy.optimize.fmin_ 

        # # Estimated maximum bound of nozzle lift is that required for an ascent
        # # rate of 6 m/s

        # nozzleLiftLowerBound = 
        # nozzleLiftUpperBound = 
        # nozzelLift_Vector = []
        self.results = []

        # distance_map = {}

        if method == 'Nelder-Mead':
            try:
                x0 = kwargs.pop('x0')
            except KeyError:
                logger.info('No Initial guess x0 provided. Using half of the windowDuration.')
                x0 = 0.5 * self.windowDuration
            res = scipy.optimize.minimize(self.targetDistance, x0=x0,
                                          method='Nelder-Mead',
                                          **kwargs)
            best_x = self.start_dateTime + timedelta(hours=res.x[0])

        elif method == 'DE':
            raise NotImplementedError('{} method not yet implemented'.format(method))

        else:
            raise ValueError('No known optimization method for {}'.format(method))

    # def plotBruteForceComparison(self):
    #     """
    #     """
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)

    #     plt.axis('equal')
    #     dateTimes = []
    #     distances = []


    #     raise NotImplementedError()


    def plotPaths3D(self):
        """Plots the resulting trajectories contained in self.results, along
        with target 

        This method should be used after calls to self.targetDistance, or
        any method in optimizeTargetLandingSite, since these will populate
        self.results
        """
        fig = plt.figure()
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
        return fig, ax

    def plotLandingSites(self):
        """
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)

        plt.axis('equal')

        for profile in self.results:
            lat_arr, lon_arr, alt_arr = profile.latitudeProfile, profile.longitudeProfile, profile.altitudeProfile
            ax.plot(lat_arr[-1], lon_arr[-1], 'bx')
            
        best_lat, best_lon, best_alt = self.bestProfile.latitudeProfile, self.bestProfile.longitudeProfile,self.bestProfile.altitudeProfile
        ax.plot(self.targetLat, self.targetLon, 'gx', label='Target')
        ax.plot(self.launchSites[0][0], self.launchSites[0][1], 'rx', label='Launch Site')
        ax.plot(best_lat[-1], best_lon[-1], 'kx', label='Best')
            
        ax.legend(loc='lower left')
        ax.set_xlabel('Lat (deg)')
        ax.set_ylabel('Lon (deg)')
        return fig, ax