# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-05-26 15:57:13
from .simulator import flight
from .weather import forecastEnvironment
from datetime import datetime, timedelta
import logging
import operator
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



logger = logging.getLogger(__name__)

class flightProfile(object):
    def __init__(self, X, result_path, distance):
        self.X = X
        self.result_path = result_path
        self.distance = distance


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
                 progress_to_file=False):

        self.targetLat = targetLat
        self.targetLon = targetLon
        self.targetElev = targetElev
        self.start_dateTime = start_dateTime
        self.windowDuration = windowDuration
        self.launchSites = launchSites

        # This section could be made faster if launch sites are clustered into
        # single forecasts (and therefore, less download requests)
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

    def targetDistance(self, X, best=[], appendResult=True):
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
        launchDateTime = self.start_dateTime + timedelta(t)
        result, solution = self.fly(0, launchDateTime)
        landing_lat = result['result'][2][-1]
        landing_lon = result['result'][3][-1]
        if appendResult:
            self.results.append(result)
        best = []

        return np.linalg.norm(np.array([self.targetLat, self.targetLon]) - np.array([landing_lat, landing_lon]))

    def bruteForce(self):
        """ Sample the parameter space and form a map of the landing site
        distance landscape

        Currently only considers time
        """
        self.results = []

        distance_map = {}
        for i, launchSiteForecast in enumerate(self.launchSiteForecasts):
            self.environment = launchSiteForecast

            for j, t in enumerate(datetimeVector):
                # distance_lift_vec = np.zeros(np.length(nozzelLift_Vector))
                # brute force approach
                distance = self.targetDistance(t)
                # distance_lift_vec[k] = distance

                distance_map[t] = distance

        best_datetime = min(distance_map.items(), key=operator.itemgetter(1))[0]

        # This seemed like the best way to extract the path from the
        # key obtained from the distance map. Note that the best_result
        # dict should only contain one entry, unless simulations are
        # duplicated
        best_result = [result_dict for result_dict in self.results if result_dict['launchDateTime'] == best_datetime]
        self.best_result = best_result[0]
        return best_result[0], distance_map


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
         
        datetimeVector = [self.start_dateTime + timedelta(hours=t)
                            for t in range(self.windowDuration)]


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

        for result_dict in self.results:
            result = result_dict['result']
            lat_arr, lon_arr, alt_arr = result[2:5]
            ax.plot(lat_arr, lon_arr, alt_arr, 'b-')

        best_lat, best_lon, best_alt = self.best_result['result'][2:5]
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

        for result_dict in self.results:
            result = result_dict['result']
            lat_arr, lon_arr, alt_arr = result[2:5]
            ax.plot(lat_arr[-1], lon_arr[-1], 'bx')
            
        best_lat, best_lon, best_alt = self.best_result['result'][2:5]
        ax.plot(self.targetLat, self.targetLon, 'gx', label='Target')
        ax.plot(self.launchSites[0][0], self.launchSites[0][1], 'rx', label='Launch Site')
        ax.plot(best_lat[-1], best_lon[-1], 'kx', label='Best')
            
        ax.legend(loc='lower left')
        ax.set_xlabel('Lat (deg)')
        ax.set_ylabel('Lon (deg)')
        return fig, ax