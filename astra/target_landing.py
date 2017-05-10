# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-05-10 16:35:33
from .simulator import flight
from .weather import forecastEnvironment


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
    self.start_dateTime
    self.windowDuration = windowDuration
    self.launchSites = launchSites

    # This section could be made faster if launch sites are clustered into
    # single forecasts (and therefore, less download requests)
    self.launchSiteForecasts = [forecastEnvironment(lat, lon, elev, dateAndTime,
                                    inflationTemperature=inflationTemperature,
                                    forceNonHD=(not HD),
                                    forecastDuration=windowDuration,
                                    debugging=False,
                                    progressHandler=None,
                                    load_on_init=False
                                    ) for (lat, lon, elev) in launchSites]

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
    #                         
    super(targetFlight, self).__init__(balloonGasType,
                            balloonModel,
                            nozzleLift,
                            payloadTrainWeight, 
                            HD=HD,
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


    def optimizeTargetLandingSite(self):
        """
        """
        # run the simulation every hour over the time window
        timeVector = numpy.arange(0, self.maxFlightTime + self.samplingTime,
            self.samplingTime)
        datetimeVector = [self.start_dateTime + datetime.timedelta(hours=t)
                            for t in timeVector]

        for i, landingSite in self.landingSites:
            launchSiteForecast = launchSiteForecasts[i]
            self.simflight.environment = launchSiteForecast
            # brute force approach
               
