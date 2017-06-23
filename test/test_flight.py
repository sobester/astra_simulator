# -*- coding: utf-8 -*-
# @Author: p-chambers
# @Date:   2017-04-21 17:21:23
# @Last Modified by:   p-chambers
# @Last Modified time: 2017-06-22 15:01:12
import pytest
from datetime import datetime, timedelta
from astra.weather import forecastEnvironment
from astra import flight
import astra
from astra.GFS import GFS_Handler
import tempfile
import os
import numpy as np
import filecmp


@pytest.fixture(params=[
  # These parameters are designed to break the code - hopefully no one makes
  # balloons with these names, or uses a heavier than air lifting gas (SF6)
    ('balloonModel', 'MYMADEUPBALLOON'),
    ('baloonGasType', 'SF6'),
    ('parachuteModel', 'RandomChuteModelX'),
    ('nozzleLift', 0.01),    # This is much smaller than the Payload weight
    ('simEnvironment', None),

])
def example_inputs(request):
    """Initializes a standard set of inputs to the flight class and injects
    a requested incorrect value input to check that errors are raised
    appropriately"""
    launch_datetime = datetime.now() + timedelta(days=1)
    # simEnvironment = forecastEnvironment(launchSiteLat=29.2108,      # deg
    #                                      launchSiteLon=-81.0228,     # deg
    #                                      launchSiteElev=4,           # m
    #                                      dateAndTime=launch_datetime,
    #                                      forceNonHD=True,
    #                                      debugging=True)
    output_dir = os.path.join(tempfile.gettempdir(), 'astra_output')

    inputs = {'environment': None,
              'balloonGasType': 'Helium',
              'balloonModel': 'TA800',
              'nozzleLift': 1,                            # kg
              'payloadTrainWeight': 0.433,                    # kg
              'parachuteModel': 'SPH36',
              'trainEquivSphereDiam': 0.1,
              'outputFile': output_dir}
    inputs[request.param[0]] = request.param[1] 
    return inputs


def test_forecast_example_inputs():
    """Note: this test is extremely minimal, and only checks if a solution is
    similar to a previous result.

    This test may also fail if parallelisation is used (which may be in the
    case in MKL versions of numpy), since the order of random number accessing
    from different processes/threads could be different.
    """
    np.random.seed(42)     # This should not be changed
    launch_datetime = datetime(2017, 4, 24,hour=12,minute=15)
    simEnvironment = forecastEnvironment(launchSiteLat=29.2108,      # deg
                                         launchSiteLon=-81.0228,     # deg
                                         launchSiteElev=4,           # m
                                         dateAndTime=launch_datetime,
                                         forceNonHD=True,
                                         debugging=True)

    # Set up the example input data files (from 24/04/2017, Daytona Beach)
    fileDict = {}
    for paramName in GFS_Handler.weatherParameters.keys():
        fileDict[paramName] = os.path.join(os.path.dirname(astra.__file__),
            '../test/example_data',
            'gfs_0p50_06z.ascii?{}[12:15][0:46][231:245][545:571]'.format(paramName))


    simEnvironment.loadFromNOAAFiles(fileDict)

    output_dir = tempfile.mktemp(suffix='')

    inputs = {'environment': simEnvironment,
              'balloonGasType': 'Helium',
              'balloonModel': 'TA800',
              'nozzleLift': 1,                            # kg
              'payloadTrainWeight': 0.433,                    # kg
              'parachuteModel': 'SPH36',
              'trainEquivSphereDiam': 0.1,
              'outputFile': output_dir}
    simFlight = flight(**inputs)
    simFlight.run()
    # TODO: Add more checks here to compare the path obtained with a reference
    # solution. This will require some saved forecast files
    test_fname = os.path.join(output_dir, 'out.kml')
    assert(os.path.isfile(test_fname))
    # check that the kml is a close match to the reference solution
    reference_fname = os.path.join(os.path.dirname(__file__), 'example_data/expected_output/out.kml')
    assert(filecmp.cmp(reference_fname, test_fname))


def test_invalid_inputs(example_inputs):
    inputs = example_inputs
    with pytest.raises(Exception) as e_info:
        simFlight = flight(**inputs)
        simFlight.preflight()
