# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-04-25 12:40:25
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-04-25 16:28:06
import pytest
from astra.weather import soundingEnvironment
from astra import flight
import os
import tempfile
from datetime import datetime, timedelta
import astra

def test_soundingEnvironment():
    """This roughly runs the sounding example and performs some primitive
    assertions"""
    launch_datetime = datetime.now() + timedelta(days=1)
    soundingFile = os.path.join(os.path.dirname(astra.__file__),
        '../examples/sp.sounding')

    simEnvironment = soundingEnvironment(launchSiteLat=50.2245,         # deg
                                         launchSiteLon=-5.3069,         # deg
                                         launchSiteElev=60,             # m
                                         distanceFromSounding=0,        # km
                                         timeFromSounding=3,            # hours
                                         inflationTemperature=10.5,     # degC
                                         dateAndTime=datetime.now(),
                                         soundingFile=soundingFile,
                                         debugging=True)

    # The forecast download is automatically called by the flight object below.
    # However, if you'd like to download it in advance, uncomment the following
    # line. The flight object will automatically detect it and won't download
    # the forecast twice.
    # simEnvironment.loadForecast()

    # Launch setup
    output_dir = os.path.join(tempfile.gettempdir(), 'astra_output')

    simFlight = flight(environment=simEnvironment,
                       balloonGasType='Helium',
                       balloonModel='TA800',
                       nozzleLift=1,                                # kg
                       payloadTrainWeight=0.433,                    # kg
                       parachuteModel='SPH36',
                       numberOfSimRuns=2,
                       trainEquivSphereDiam=0.1,                    # m
                       excessPressureCoeff=1,
                       outputFile=output_dir,
                       debugging=True,
                       log_to_file=True)

    simFlight.run()

    # Check that the output was created
    test_fname = os.path.join(output_dir, 'out.kml')
    assert(os.path.isfile(test_fname))

