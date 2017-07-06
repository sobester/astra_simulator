# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-05-08 11:36:23
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-07-06 18:12:32
from .simulator import flight, flightProfile
from .weather import forecastEnvironment
from .available_balloons_parachutes import balloons
from .target_landing import targetProfile

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
from .flight_tools import nozzleLiftFixedAscent, liftingGasMass
import matplotlib.dates as mdates
import inspect
from deap import tools as dptools
from deap import creator
from deap import base
from deap import algorithms
import matplotlib.ticker as ticker
import matplotlib.tri as mtri
import itertools
import random
import functools


# Helper functions for deap
#############################################################################
def interpIndividual(bounds, individual):
    interp = []
    for i in range(len(individual)):
        bmin, bmax = bounds[i]
        interp.append(bmin + (bmax-bmin) * individual[i])
    return interp

def checkBounds(bmin, bmax):
    def decorator(func):
        def wrappper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    if child[i] > bmax:
                        child[i] = bmax
                    elif child[i] < bmin:
                        child[i] = bmin
            return offspring
        return wrappper
    return decorator
#############################################################################


# Create the class that will be used for assessing multi-objective fitness. By
# default, we want to apply equal minimizing weights to all objectives, but
# this will be overwritten by the targetFlight.optimizeTargetLanding function
creator.create("flightFitness", base.Fitness, weights=(-1, -1, -1))

logger = logging.getLogger(__name__)


class targetProfile(flightProfile):
    def __init__(self, launchDateTime,
                 nozzleLift,
                 flightNumber,
                 timeVector,
                 latitudeProfile,
                 longitudeProfile,
                 altitudeProfile,
                 highestAltIndex,
                 highestAltitude,
                 hasBurst,
                 balloonModel,
                 fitness,
                 X
                 ):
        super(targetProfile, self).__init__(launchDateTime=launchDateTime,
            nozzleLift=nozzleLift,
            flightNumber=flightNumber,
            timeVector=timeVector,
            latitudeProfile=latitudeProfile,
            longitudeProfile=longitudeProfile,
            altitudeProfile=altitudeProfile,
            highestAltIndex=highestAltIndex,
            highestAltitude=highestAltitude,
            hasBurst=hasBurst,
            balloonModel=balloonModel)
        self.fitness = fitness
        self.X = X


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
    weights : tuple of signed int
        The weights for each of the objectives assessed by the targetDistance
        function. 

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
                 weights=(-1, -1, -1),
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

        self.targetLat = targetLat
        self.targetLon = targetLon
        self.targetElev = targetElev
        self.start_dateTime = start_dateTime
        self.windowDuration = windowDuration
        self.launchSites = launchSites

        self.end_dateTime = self.start_dateTime + timedelta(hours=windowDuration)
        self.bestProfile = None
        self.launchSiteForecasts = []
        self.balloonsSelected = [balloonModel]

        # Use a maximum lateral upper atmospheric speed to determine if this
        # flight is feasible, given an 'as the crow flies' flight entirely in
        # the jet stream, assuming max speed of 200 km/hr
        self.maxLateralSpeed = 200    # m/s 
        self.cutoffDistance = self.maxLateralSpeed * self.maxFlightTime / 3600.

        # Set some allowable modes of flight (these are used by the 
        # targetDistance function)
        self.flightModes = ['standard', 'cutdown', 'floating']


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

        # Set some max a min bounds on ascent rate
        self.minAscentRate = 1.5
        self.maxAscentRate = 6

        # Define a list of input vectors Xs (to the target distance function),
        # and a list of objective scores Ys for plotting later:
        self.Xs = []

        self.flightFitness = deepcopy(creator.flightFitness)

        # This will update the weights used by creator.flightFitness, and hence
        # the order of the individuals added to the ParetoFront by the
        # targetDistance function.
        self.weights = weights

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, newWeights):
        """Takes a tuple of weights, corresponding to the multiobjective
        values of fitness evaluated for each flight profile.

        Parameters
        ----------
        newWeights : tuple
            Weightings used to multiply values in the objective function
            newWeights[0] : normalised distance from target
        """
        self.flightFitness.weights = newWeights
        self._weights = newWeights

    @property
    def balloonsSelected(self):
        return self._balloonsSelected

    @balloonsSelected.setter
    def balloonsSelected(self, newBalloonModels):
        """Takes an input list of selected balloon models, and stores it as an
        attribute. Also updates the self.largestBalloon, which is used to
        normalise the cost (gas mass) objective.
        """
        self._balloonsSelected = {k: balloons[k] for k in newBalloonModels}
        self.largestBalloon = max(self._balloonsSelected,
            key=lambda k: abs(self._balloonsSelected.get(k)[1]))

    def targetDistanceFactory(self, constantsDict):
        """Constructs a N dimensional objective function, where N is the number
        of arguments in self.targetDistance minus the number of arguments in
        constantsDict

        Effectively, removes the parameters in constantsDict from the
        targetDistance function, and creates a new objective where the input
        is a vector of the remaining args.

        Parameters
        ----------
        constantsDict : :obj:`dict`
            The mapping of key: value pairs that will be made constant in the
            targetDistance master objective function
        """
        args = inspect.signature(self.targetDistance).parameters

        remaining_args = [arg for arg in args if arg not in constantsDict]

        logger.debug('Arguments in the objective function are: {}'.format(remaining_args))

        # Store the list of ordered variable names that represent data in the 
        # input X vector
        self.variables = remaining_args

        def objective(X):
            X_kwargs = dict(zip(remaining_args, X))
            kwargs = {**constantsDict, **X_kwargs}

            # Check that Floating flight is within the bounds of zero and 3,
            # then round to give an integer (Note: this doesn't work with Nelder
            # Mead
            print(kwargs['flightMode'])
            if not isinstance(kwargs['flightMode'], str):
                assert(kwargs['flightMode'] >= 0 and kwargs['flightMode'] <= len(self.flightModes)-1),\
                    "flightMode value {} out of bounds (0, {})".format(
                        kwargs['floatingFlight'], len(self.flightModes))
                # Convert integer to a mode (1: standard, 2: cutdown, 3: floating)
                kwargs['flightMode'] = self.flightModes[int(round(kwargs['flightMode']))]

            return self.targetDistance(**kwargs)

        return objective

    def createObjectiveAndBounds(self, flightModes=['standard'],
        flexibleBalloon=False, deviceActivationAltitudeBounds=(),
        balloonModels=[], returnWeightedSum=True):
        """

        Returns
        -------
        objective : function

        boundsDict : list
            The bounds of each parameter in the returned objective function,
            in the order they appear in the arguments.
        """
        # Build the dictionary of constants that the objective function will
        # ignore:
        constantsDict = {}

        self.flightModes = flightModes

        if len(flightModes) == 1:
            mode = flightModes[0]
            constantsDict['flightMode'] = mode
        if 'floating' not in flightModes:
            constantsDict['floatDuration'] = np.inf
        if 'cutdown' not in flightModes and 'floating' not in flightModes:
            constantsDict['deviceActivationAltitude'] = np.inf

        if not flexibleBalloon:
            constantsDict['balloonNominalBurstDia'] = balloons[self.balloonModel][1]

        constantsDict['returnWeightedSum'] = returnWeightedSum

        objective = self.targetDistanceFactory(constantsDict)

        nozzleLiftLowerBound = nozzleLiftFixedAscent(self.minAscentRate,
            self._balloonWeight, self.payloadTrainWeight,
            self.environment.inflationTemperature,
            self.environment.getPressure(self.launchSiteLat,
                                              self.launchSiteLon,
                                              self.launchSiteElev,
                                              self.start_dateTime),
            self._gasMolecularMass, self.excessPressureCoeff,
            CD=(0.225 + 0.425)/2.)
        nozzleLiftUpperBound = nozzleLiftFixedAscent(self.maxAscentRate,
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.)

        # Set up the dictionary of bounds on all variables: note that some may
        # not be passed to scipy, as this depends on the constants used in the
        # objective function
        boundsDict = {}
        boundsDict['t'] = (0, self.windowDuration)
        boundsDict['targetAscentRate'] = (1.5, 6.0)
        boundsDict['floatingFlight'] = (0, 1)
        # Limit float duration to a continuous venting (0) to the entire flight
        # (The penalty imposed in targetDistance function should avoid this area)
        boundsDict['floatDuration'] = (0, self.maxFlightTime)

        # device altitude parameter bounds
        boundsDict['deviceActivationAltitude'] = deviceActivationAltitudeBounds

        # Variable balloon parameter bounds:
        if balloonModels:
            self.balloonsSelected = balloonModels
        else:
            self.balloonsSelected = [self.balloonModel]
        balloonRads = [bal[1] for bal in self.balloonsSelected.values()]
        boundsDict['balloonNominalBurstDia'] = (min(balloonRads), max(balloonRads))

        # Need to subtract 1 from the selected flightModes, since this will be
        # converted to an index (starting from 0):
        boundsDict['flightMode'] = (0, len(self.flightModes) - 1)

        # Get the maximum gas mass expected, for normalising the cost
        # (gas mass) objective. ISA pressure is used, since the difference
        # in max pressure vs isa pressure over this time window should not
        # have a great effect on the maximum estimated gas mass
        self.maxGasMass = liftingGasMass(nozzleLiftUpperBound,
            balloons[self.largestBalloon][0], self.environment.inflationTemperature,
            1013.25, self._gasMolecularMass, self.excessPressureCoeff)[0]

        # self.variables was updated by self.targetDistanceFactory, and includes
        # the name of all variables that will be used for optimisation, in order
        bounds = [boundsDict[var] for var in self.variables]

        return objective, bounds

    def targetDistance(self, t, targetAscentRate, flightMode, deviceActivationAltitude,
        floatDuration, balloonNominalBurstDia, returnWeightedSum=True):
        """
        Parameters
        ----------
        X : iterable
            The input vector. Should contain elements,
                x[0] :  time, in hours, after the start of the optimisation
                window (self.startTime)
        
        returns
        -------
        distance : scalar
            The euclidean norm of [target_lon - lon, target_lat - lat] (degrees)
        """
        # t = X[0]
        # nozzleLift = X[1]

        # Find the balloon model with nearest burst diameter to that of the
        # input
        if balloonNominalBurstDia:
            self.balloonModel = min(self.balloonsSelected,
                key=lambda k: abs(self.balloonsSelected.get(k)[1]-balloonNominalBurstDia))

        # Convert ascent rate to nozzle lift for the nearest balloon model
        # giving the input burst diameter
        nozzleLift = nozzleLiftFixedAscent(targetAscentRate,
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.)

        log_msg = "Running flight for datetime {}, nozzleLift={}kg, balloon {}".format(
                        self.start_dateTime + timedelta(hours=t), nozzleLift, self.balloonModel)

        if flightMode == 'floating':
            log_msg += ", Floating Altitude {}m, Duration {} seconds".format(deviceActivationAltitude, floatDuration)
            self.floatingFlight = True
            self.floatingAltitude = deviceActivationAltitude
            self.floatDuration = floatDuration

        elif flightMode == 'cutdown':
            log_msg += ", cutdown Altitude {}m".format(deviceActivationAltitude)
            self.cutdown = True
            self.cutdownAltitude = deviceActivationAltitude

        logger.debug(log_msg)

        Nobjs = 3
        fitnessArr = np.zeros(Nobjs)

        # nozzle lift is an access controlled variable, which may raise an error if
        # lower than the payload lift: return a penalty if this is the case
        try:
            self.nozzleLift = nozzleLift
        except ValueError:
            # Pay a penalty larger than the maximum possible values
            fitnessArr += [5e6] * Nobjs
            if returnWeightedSum:
                return sum(fitnessArr)
            else:
                return fitnessArr

        launchDateTime = self.start_dateTime + timedelta(hours=t)

        resultProfile, solution = self.fly(0, launchDateTime)

        # Penalty for flights not bursting within the time limit
        if not resultProfile.hasBurst:
            fitnessArr += [5e6] * Nobjs
        
        # Distance objective (normalised) - CURRENTLY REMOVING NORMALISATION
        landing_lat = resultProfile.latitudeProfile[-1]
        landing_lon = resultProfile.longitudeProfile[-1]
        distNorm = (tools.haversine(landing_lat, landing_lon, self.targetLat, self.targetLon)) # /
            #self.cutoffDistance)

        # Cost related objective (currently evaluates the gas mass required to
        # achieve the nozzle lift used for this profile, normalised by the max
        # )
        if self.balloonsSelected:
            balloonsSelected = self.balloonsSelected
        else:
            # Just use the full dict from available_balloons_parachutes.py
            balloonsSelected = balloons

        gasMassNorm = self._gasMassAtInflation # / self.maxGasMass

        # Time related objective (could be useful for minimizing cold soak time)
        timeNorm = resultProfile.flightDurationSecs  #/ self.maxFlightTime

        fitnessArr += [distNorm, gasMassNorm, timeNorm]
        fitness = self.flightFitness(fitnessArr)
        self.fitnesses.append(fitness)
        X = [t, targetAscentRate, flightMode, deviceActivationAltitude,
            floatDuration, balloonNominalBurstDia]
        resultProfile = targetProfile.fromProfile(resultProfile, fitness=fitness, X=X)

        # Use the ParetoFront update method to store this profile lexographically
        self.results.update([resultProfile])

        if returnWeightedSum:
            return - sum(fitness.wvalues)
        else:
            return fitness.values

    def _callbackStoreResult(self, xk, convergence):
        xk[1] = nozzleLiftFixedAscent(xk[1],
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.)
        self.Xs.append(xk)

    def bruteForce(self, Nx, Ny, balloonModel, flightMode='standard', deviceActivationAltitude=None, floatDuration=None):
        """ Sample the parameter space and form a map of the landing site
        distance landscape
        Currently only considers time
        """
        # Set up the arrays where results will be stored. Also storing the 
        # multiple objectives 
        self.results = dptools.ParetoFront()
        scores = np.zeros([Nx, Ny])
        self.fitnesses = []

        dts = np.linspace(0, self.windowDuration, Nx)
        dateTimeVec = [self.start_dateTime + timedelta(hours=dt) for dt in dts]

        self.balloonModel = balloonModel
        self.balloonsSelected = [balloonModel]

        for launchSiteForecast in self.launchSiteForecasts:
            self.environment = launchSiteForecast
            # Estimated maximum bound of nozzle lift for target ascent rates

            ascentRateVec = np.linspace(self.minAscentRate, self.maxAscentRate, Ny)
            X, Y = np.meshgrid(dts, ascentRateVec)

            nozzleLiftVec = np.array([nozzleLiftFixedAscent(ascRate,
                    self._balloonWeight, self.payloadTrainWeight,
                    self.environment.inflationTemperature,
                    self.environment.getPressure(self.launchSiteLat,
                                                      self.launchSiteLon,
                                                      self.launchSiteElev,
                                                      self.start_dateTime),
                    self._gasMolecularMass, self.excessPressureCoeff,
                    CD=(0.225 + 0.425)/2.) for ascRate in ascentRateVec])

            # Get the maximum gas mass expected, for normalising the cost
            # (gas mass) objective. ISA pressure is used, since the difference
            # in max pressure vs isa pressure over this time window should not
            # have a great effect on the maximum estimated gas mass
            self.maxGasMass = liftingGasMass(nozzleLiftVec[-1],
                balloons[self.largestBalloon][0], self.environment.inflationTemperature,
                1013.25,
                self._gasMolecularMass,
                self.excessPressureCoeff)[0]

            logger.info("Date range: [{}, {}], Nx={} points".format(
                dateTimeVec[0], dateTimeVec[-1], Nx))
            logger.info("Nozzle Lift range: [{}, {}] (kg), Ny={} points".format(
                nozzleLiftVec[0], nozzleLiftVec[-1], Ny))

            logger.debug("Running brute force calculation")
            for i, t in enumerate(dts):
                for j, ascRate in enumerate(ascentRateVec):
                    # brute force approach
                    objective = self.targetDistance(t=t, targetAscentRate=ascRate,
                        flightMode=flightMode, deviceActivationAltitude=deviceActivationAltitude,
                        floatDuration=floatDuration,
                        balloonNominalBurstDia=balloons[self.balloonModel][1])
                    # distance_lift_vec[k] = distance
                    scores[i, j] = objective

        bestProfile = max(self.results, key=lambda prof: sum(prof.fitness.wvalues))
        self.bestProfile = bestProfile
        return bestProfile, dateTimeVec, nozzleLiftVec, scores

    def bruteForceSlice(self, Nx=21, Ny=21, balloonModel=None, flightMode='standard',
        sliceParam='', Nslices=None, sliceBounds=(), sliceParam_subset=[]):
        """Runs a brute force (discrete grid) of the targetDistance objective
        function, using Nx

        Parameters
        ----------
        Nx, Ny : int
            Number of sampling points in Time and nozzle Lift respectively,
            between the bounds of [0, self.windowDuration], and nozzle lift
            calculated to achieve ascent rates between the bounds of
            self.minAscentRate and self.maxAscentRate
        sliceParam : string
            variable name for which to slice the distance vs nozzle lift and
            time landscape. Expected names are any of
                'floatingAltitude', 'cutdownAltitude'
        Nslices : int
            Number of slices to use between sliceBounds to use for sliceParam
        sliceBounds : tuple
            the lower limit, upper limit between which to sample sliceParam

        Notes
        -----
        Whichever parameter is chosen for slicing will be overwritten. E.g.:
        if balloonModel='TA100', but the 

        :Example:
            >>> 
        """

        # Keep all paths for now, to visualize the landscape of nozzle lift
        self.results = dptools.ParetoFront()
        self.fitnesses = []

        dts = np.linspace(0, self.windowDuration, Nx)
        distances = np.zeros(Nx)
        dateTimeVec = [self.start_dateTime + timedelta(hours=dt) for dt in dts]

        if balloonModel:
            self.balloonModel = balloonModel
            self.balloonsSelected = [balloonModel]

        assert(Nslices), "Argument 'Nslices' is required to slice the landscape through {}".format(sliceParam)
        assert(sliceBounds), "Argument 'sliceBounds' is required to slice the landscape through {}".format(sliceParam)

        # Start with a 'vanilla' objective function, where only t and targetAscentRate
        # are used, and edit the single parameter which defines the slice:
        objectivebounds_kwargs = {'flightModes': [flightMode],
                                  'flexibleBalloon': False,
                                  'deviceActivationAltitudeBounds': (),
                                  'returnWeightedSum': True,
                                  'balloonModels': [balloonModel]
                                  }

        if sliceParam == 'balloonModel':
            objectivebounds_kwargs['balloonModels'] = sliceParam_subset
            objectivebounds_kwargs['flexibleBalloon'] = True

        elif sliceParam == 'flightMode':
            assert(sliceParam_subset),\
                'A deviceActivationAltitudeBounds tuple is required to slice through deviceActivationAltitude'
            objectivebounds_kwargs['flightModes'] = sliceParam_subset

        elif sliceParam == 'deviceActivationAltitude':
            assert(sliceBounds),\
                'A sliceBounds tuple is required to slice through deviceActivationAltitude'
            objectivebounds_kwargs['deviceActivationAltitude'] = sliceBounds


        objectiveFunction, bounds = self.createObjectiveAndBounds(flightModes=flightModes,
            flexibleBalloon=flexibleBalloon,
            deviceActivationAltitudeBounds=deviceActivationAltitudeBounds,
            balloonModels=balloonModels,
            returnWeightedSum=returnWeightedSum)

        assert(len(bounds) == 3),\
            '{} parameters appear in the objective function: not specified which to slice through.'.format(len(bounds))

        if sliceParam in ('balloonModel', 'flightMode'):
            # These are discrete parameters: ensure that only the correct number
            # of slices is used
            logger.warning('A discrete slicing parameter was selected: ignoring input Nslices')
            Nslices = len(sliceParam_subset)

        # The bounds in the last position belong to the slicing parameter
        sliceVec = np.linspace(min(bounds[-1]), max(bounds[-1]), Nslices)

        # Build the objective function for these slices:
        # if sliceParam == 'floatingAltitude':
        #     sliceVec = np.linspace(min(sliceBounds), max(sliceBounds), Nslices)
        #     constantsDict = {'floatingFlight': True, 'cutdown': False,
        #         'cutdownAltitude':np.inf, 'floatDuration': np.inf, 
        #         'balloonNominalBurstDia': balloons[self.balloonModel][1]}

        # elif sliceParam == 'cutdownAltitude':
        #     cutdownFlight = True
        #     sliceVec = np.linspace(min(sliceBounds), max(sliceBounds), Nslices)
        #     constantsDict = {'cutdown': True, 'floatingFlight': False,
        #         'floatingAltitude':0, 'floatDuration': np.inf,
        #         'balloonNominalBurstDia': balloons[self.balloonModel][1]}

        # else:
        #     logger.debug("{} is an unknown slicing parameter. See targetFlight.targetDistance parameter names".format(sliceParam))


        self.environment = launchSiteForecast




        ascentRateVec = np.linspace(self.minAscentRate, self.maxAscentRate, Ny)

        nozzleLiftVec = np.array([nozzleLiftFixedAscent(ascRate,
                self._balloonWeight, self.payloadTrainWeight,
                self.environment.inflationTemperature,
                self.environment.getPressure(self.launchSiteLat,
                                                  self.launchSiteLon,
                                                  self.launchSiteElev,
                                                  self.start_dateTime),
                self._gasMolecularMass, self.excessPressureCoeff,
                CD=(0.225 + 0.425)/2.) for ascRate in ascentRateVec])


        # Get the maximum gas mass expected, for normalising the cost
        # (gas mass) objective. ISA pressure is used, since the difference
        # in max pressure vs isa pressure over this time window should not
        # have a great effect on the maximum estimated gas mass
        self.maxGasMass = liftingGasMass(nozzleLiftVec[-1],
            balloons[self.largestBalloon][0], self.environment.inflationTemperature,
            1013.25,
            self._gasMolecularMass,
            self.excessPressureCoeff)[0]

        X, Y, Z = np.meshgrid(dts, ascentRateVec, sliceVec, indexing="ij")

        sliceUnits = {'floatingAltitude': 'm', 'cutdownAltitude': 'm'}

        logger.info("Date range: [{}, {}], Nx={} points".format(
            dateTimeVec[0], dateTimeVec[-1], Nx))
        logger.info("Nozzle Lift range: [{}, {}] (kg), Ny={} points".format(
            self.nozzleLiftLowerBound, self.nozzleLiftUpperBound, Ny))
        logger.info("{} range: [{}, {}] {}, Nz={} points]".format(
            sliceParam, min(sliceBounds), max(sliceBounds),
            sliceUnits[sliceParam], Nslices))
        logger.debug("Running brute force calculation")
        
        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        # Main calculation here
        objectiveVals = [objectiveFunction(point) for point in points]

        distances = np.reshape(objectiveVals, (np.size(dts), np.size(nozzleLiftVec), np.size(sliceVec)))

        bestProfile = max(self.results, key=lambda prof: sum(prof.fitness.wvalues))
        self.bestProfile = bestProfile
        return bestProfile, dateTimeVec, nozzleLiftVec, sliceVec, distances

    def optimizeTargetLandingSite(self, flightModes=['standard'],
        flexibleBalloon=False, deviceActivationAltitudeBounds=(), balloonModels=[],
        method='Nelder-Mead', weights=(), seed=None, **kwargs):
        """
        Parameters
        ----------
        flightMode : list of string
            defines the flight modes to use in the optimization:
            'standard' = standard ascent, burst, descent
            'floating' = ascent with venting, float at deviceActivationAltitude,
                burst after floatDuration seconds, then descend
            'cutdown' : ascend up to deviceActivationAltitude, force burst, 
                then descend



        params : :obj:`dict`
            The variables which will be used in the optimization

        method : string
            'brute-force' : Sample flight every hour within the time window
                *Note: this probably isn't suitable given more than two params
            'Nelder-Mead' : Use scipy's Nelder mead pattern search. This has
                proven to be sensitive to initial guesses for this application.
            'DE' : Use scipy differential evolution.
            'GA' : Use a mu plus lambda genetic algorithm, from DEAP

        [maxsize] : int (default 10)
            The number of (lexographically sorted by distance from target site)
            profiles to keep in self.results

        [weights] : tuple of int
            if provided, the value of self.weights will be overridden, changing
            the weightings used in the objective function and Pareto Front
            ordering. See self.weights.

        **kwargs - extra args to pass to scipy
        
        See Also
        --------
        targetFlight.targetDistance
        
        """
        # run the simulation every hour over the time window. Note that we need
        # to also get weather to cover a flight of duration <maxFlightTime>,
        # starting at the end of the window.
        # scipy.optimize.fmin_ 

        self.results = dptools.ParetoFront()
        # Store all objective scores, for Pareto plotting
        self.fitnesses = []

        self.Xs = []

        if weights:
            self.weights = weights
 
         # Need to return a tuple of fitness for ga, or a weighted sum for scipy:
        returnWeightedSum = (method.lower() != 'ga')

        # # For now, assume the first is the only interesting launch site
        # Estimated maximum bound of nozzle lift for target ascent rates
        self.environment = self.launchSiteForecasts[0]

        objective, bounds = self.createObjectiveAndBounds(flightModes=flightModes,
            flexibleBalloon=flexibleBalloon,
            deviceActivationAltitudeBounds=deviceActivationAltitudeBounds,
            balloonModels=balloonModels,
            returnWeightedSum=returnWeightedSum)

        if method.lower() in ('nelder-mead', 'l-bfgs-b'):
            try:
                x0 = kwargs.pop('x0')
            except KeyError:
                logger.exception('An initial guess x0 is required for method {}.'.format(method))

            try:
                res = opt.minimize(objective, x0=x0, method=method,
                                   callback=(lambda x: self._callbackStoreResult(x, convergence=None)),
                                   bounds=bounds, args=(), **kwargs)
            except TypeError:
                # Likely that this method does not support bounded optimisation,
                # so try without that argument
                res = opt.minimize(objective, x0=x0, method=method,
                   callback=(lambda x: self._callbackStoreResult(x, convergence=None)),
                   args=(), **kwargs)

            bestProfile = self.results[0]
            self.bestProfile = bestProfile

        elif method.lower() == 'de':
            res = opt.differential_evolution(objective, bounds=bounds,
                callback=self._callbackStoreResult, seed=seed, **kwargs)

            bestProfile = self.results[0]
            self.bestProfile = bestProfile

        elif method.lower() == 'ga':
            def evaluateIndividualTarget(individual):
                X = interpIndividual(bounds, individual)
                return objective(X)
            creator.create("FitnessMin", base.Fitness, weights=weights)
            creator.create("Individual", list, fitness=simulator.FitnessMin)

            toolbox = base.Toolbox()

            # Attribute generator
            toolbox.register("attr_float", random.uniform, 0., 1.)

            # Structure initializers
            toolbox.register("individual", dptools.initRepeat, creator.Individual, toolbox.attr_float, len(self.variables))
            toolbox.register("population", dptools.initRepeat, list, toolbox.individual)

            # define the population to be a list of individuals
            toolbox.register("population", dptools.initRepeat, list, toolbox.individual)

            pop = toolbox.population(n=150)


            toolbox.register("evaluate", evaluateIndividualTarget)
            toolbox.register("mate", dptools.cxBlend, alpha=1.5)
            toolbox.register("mutate", dptools.mutGaussian, mu=0, sigma=3, indpb=0.3)
            toolbox.register("select", dptools.selNSGA2)

            toolbox.decorate("mate", checkBounds(0, 1))
            toolbox.decorate("mutate", checkBounds(0, 1))
            
            if seed:
                random.seed(seed)

            MU, LAMBDA = 50, 100
            pop = toolbox.population(n=MU)
            stats = tools.Statistics(lambda ind: ind.fitness.values)
            stats.register("avg", np.mean, axis=0)
            stats.register("std", np.std, axis=0)
            stats.register("min", np.min, axis=0)
            stats.register("max", np.max, axis=0)
            
            # Limit to 1500 evals (about 25 minutes)
            maxevals = 1000
            cxpb = 0.5
            mutpb = 0.2
            
            # Max number of evaluations is (on average) the probability of replacement
            # times by the population size, for each generation: inverse
            ngen = int(round(maxevals / (MU * (cxpb + mutpb))))
                              
            algorithms.eaMuPlusLambda(pop, toolbox, mu=MU, lambda_=LAMBDA, 
                                      cxpb=cxpb, mutpb=mutpb, ngen=ngen, 
                                      stats=stats)

        else:
            raise ValueError('No known optimization method for {}'.format(method))
        return res

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
        ax.set_zlabel('Alt (m)')
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

    def plotObjectiveContours(self, dateTimeVector, nozzleLiftVector, scores,
        fig=None, ax=None, bestMarker='b*', appendLabel='', bestProfile=None,  **kwargs):
        """Plots the ground distance of the landing sites contained in
        self.results from the target landing site.

        Requires a call to any of the optimisation of brute force calculations
        to populate self.results.

        Parameters
        ----------
        kwargs :
            Additional named args will be passed to ax.contour for plot
            settings
        """
        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111)
            ax.set_xlabel('Date and Time')
            ax.set_ylabel('Nozzle Lift (kg)')
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
            
        CS = ax.contour(dateTimeVector, nozzleLiftVector, scores, label='Landing sites'+appendLabel, **kwargs)

        if bestProfile:
            ax.plot(bestProfile.launchDateTime, bestProfile.nozzleLift, bestMarker, markersize=8, label='min' + appendLabel)
            legend = ax.legend(loc='upper right', ncol=2)

        fig.show()
        ax.clabel(CS, inline=1, fontsize=12, fmt='%.2f')
        return fig, ax

    def plotObjectiveContours3D(self, dateTimeVector, nozzleLiftVector, scores,
        fig=None, ax=None, bestMarker='b*', appendLabel='', bestProfile=None, **kwargs):
        """Plots the ground distance of the landing sites contained in
        self.results from the target landing site.

        Requires a call to any of the optimisation of brute force calculations
        to populate self.results.

        Parameters
        ----------
        kwargs :
            Additional named args will be passed to ax.contour for plot
            settings
        """
        # Find the rows that are outside the expected range, and scaled ranges
        # of the plot
        scores_filtered = scores[scores < 5]
        if np.size(scores_filtered) < np.size(scores):
            logger.warning("Some scores are poorly scaled. Plot values will be cut off at +5")

        zmin, zmax = np.min(scores_filtered), np.max(scores_filtered)

        if not fig:
            fig = plt.figure()
        if not ax:
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlabel('Date and Time')
            ax.set_ylabel('Nozzle Lift (kg)')
            ax.set_zlabel('Objective Score')

            ax.set_zlim([zmin, zmax])

            def format_date(x, pos=None):
                return datetime.fromtimestamp(x).strftime('%m-%d %H') #use FuncFormatter to format dates

            # Use 6 x ticks, and rotate them
            xticks = [(self.start_dateTime + timedelta(hours=h)).timestamp() for h in np.linspace(0, self.windowDuration, 6)]
            ax.w_xaxis.set_major_locator(ticker.FixedLocator(xticks))
            ax.w_xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
             # re-create what autofmt_xdate but with w_xaxis
            for tl in ax.w_xaxis.get_ticklabels():
                tl.set_ha('right')
                tl.set_rotation(30)

            # ax.set_xlim(date_start, date_end)
            # fig.autofmt_xdate()
        
        # Need to convert to unix time for 3d plotting, since matplotlib tries
        # to autoscale_xyz, which doesn't work on non float types
        tstamps = [dtime.timestamp() for dtime in dateTimeVector]
        T, L = np.meshgrid(tstamps, nozzleLiftVector)

        CS = ax.plot_surface(T, L, scores, label='Landing sites'+appendLabel,
            vmin=zmin, vmax=zmax, **kwargs)

        if bestProfile:
            ax.plot([bestProfile.launchDateTime.timestamp()], [self.bestProfile.nozzleLift],
                [-sum(bestProfile.fitness.wvalues)], bestMarker, markersize=8,
                label='min' + appendLabel)

        # Manually set up x ticks, since 
        # ax.set_xticks()

        # Add a color bar which maps values to colors.
        fig.colorbar(CS, shrink=0.5, aspect=5)

        fig.show()
        # legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2)
        return fig, ax

    def plotObjectiveLocations(self, fig=None, ax=None, marker='bx',
        bestMarker='b*', appendLabel='', bestProfile=None, **kwargs):
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
            ax.set_ylabel('Nozzle Lift (kg)')
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

        if bestProfile:
            ax.plot(bestProfile.launchDateTime, bestProfile.nozzleLift, bestMarker,
                markersize=8, label='min' + appendLabel)


        fig.show()
        legend = ax.legend(loc='upper right', ncol=2)
        return fig, ax

    def plotParetoFront(self):
        """Plots the points stored in the simulator.ParetoFront (containing the non-dominated
        Pareto efficient individuals).
        
        Notes
        -----
        * Only works for 2 or 3 objectives. 
        * Assumes that all fitnesses have the same number of values and weightings
        """
        fitnesses = self.fitnesses

        nonzero_indices = np.nonzero(fitnesses[0].weights)[0]
        assert(len(nonzero_indices) >= 2),\
            "Input fitness array has less than 2 weighted objectives: no useful Pareto solutions."
            
        # Dictionary of weight_index : label pairs. Extract the axes labels for the non zero weighted vars
        axlabelsAvailable = {0: r'$\bar{\Delta}$',
                           1: r'$\bar{m}_{gas}$',
                           2: r'$\bar{t} \cdot$'}
        axlabels = [axlabelsAvailable[idx] for idx in nonzero_indices]
        
        if len(nonzero_indices) == 2:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            fs = list(zip(*[fitness.values for fitness in fitnesses]))
            f1 = np.array(fs[nonzero_indices[0]])
            f2 = np.array(fs[nonzero_indices[1]])
            ax.plot(f1, f2, 'kx')
            pareto_vars = list(zip(*[profile.fitness.values for profile in self.results]))
            pareto_array = np.array([pareto_vars[nonzero_indices[0]], pareto_vars[nonzero_indices[1]]]).T
            pareto_array = pareto_array[np.argsort(pareto_array[:, 1])]
            ax.plot(pareto_array[:, 0], pareto_array[:, 1], 'ko--', mfc='none',
                    markersize=16, label='Pareto front')

            # Find best (assumes that the weights are tuned to how import each parameter is)
            iBest = np.argmin(pareto_array[:,0] + pareto_array[:, 1])
            ax.plot(pareto_array[iBest, 0], pareto_array[iBest, 1], 'k*', markersize=12, label='Best')

            # Formatting
            ax.set_xlabel(axlabels[0])
            ax.set_ylabel(axlabels[1])
            ax.legend(loc='upper right')
            
        elif len(nonzero_indices) == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            f1, f2, f3 = zip(*[fitness.values for fitness in fitnesses])
            f1 = np.array(f1)
            f2 = np.array(f2)
            f3 = np.array(f3)
            ax.scatter(f1, f2, f3, label='Brute Force Points')
            pareto_x, pareto_y, pareto_z = zip(*[profile.fitness.values for profile in self.results])
            pareto_array = np.array([pareto_x, pareto_y, pareto_z]).T
        
            ax.scatter(pareto_array[:, 0], pareto_array[:, 1], pareto_array[:, 2])#, facecolors='k',
    #                    label='Pareto points')
        
            # Plot the Pareto surface triangulation
            triang = mtri.Triangulation(pareto_array[:,0],pareto_array[:,1])
            ax.plot_trisurf(triang,pareto_array[:,2],color='red', alpha=0.5)

            # Find best (assumes that the weights are tuned to how important each parameter is)
            iBest = np.argmin(pareto_array[:,0] + pareto_array[:, 1])
            ax.plot([pareto_array[iBest, 0]], [pareto_array[iBest, 1]], [pareto_array[iBest, 2]], 'k*', markersize=12, label='Best')

            # Formatting
            ax.set_xlabel(axlabels[0])
            ax.set_ylabel(axlabels[1])
            ax.set_zlabel(axlabels[2])
            ax.legend(loc='upper right')        
        
        else:
            raise ValueError("Can't support plotting {}D Pareto front. Adjust weightings.".format(len(nonzero_indices)))
        
            return None