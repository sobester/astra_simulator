# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   Paul Chambers
# @Last Modified time: 2017-05-11 00:19:08
from .simulator import flight
from .weather import forecastEnvironment
from datetime import datetime, timedelta
import logging
import operator
import numpy as np


logger = logging.getLogger(__name__)


class targetFlight(flight):
    """

    Parameters
    ----------
    launchSites : list of lat, lon, elev triplets
        lat, lon and elev triplets for each launch site 

    windowDuration : int
        number of hours for which to run flights
    """

    def __init__(self,
                 start_dateTime,
                 targetLat,
                 targetLon,
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

    def targetDistance(self, launchDateTime, appendResult=True):
        result = self.fly(0, launchDateTime)
        landing_lat = result['result'][2][-1]
        landing_lon = result['result'][3][-1]
        self.results.append(result)
        return np.linalg.norm(np.array([self.targetLat, self.targetLon]) - np.array([landing_lat, landing_lat]))

    def optimizeTargetLandingSite(self):
        """
        """
        # run the simulation every hour over the time window
        datetimeVector = [self.start_dateTime + timedelta(hours=t)
                            for t in range(self.windowDuration)]
        self.results = {}

        distance_map = {}
        for i, launchSiteForecast in enumerate(self.launchSiteForecasts):
            self.environment = launchSiteForecast

            try:
                self.preflight(self.environment.dateAndTime)
            except:
                logger.exception(
                    "Error during preflight validations and calculations:")
                raise
            for j, t in enumerate(datetimeVector):
                # brute force approach
                distance = self.targetDistance(t)
                distance_map[t] = distance

        best_datetime = min(distance_map.items(), key=operator.itemgetter(1))[0]

        best_result = [result_dict for result_dict in self.results if result_dict['launchDateTime'] == best_datetime]
        return best_result, distance_map
