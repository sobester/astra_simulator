astra API reference
===================

The main classes for user interaction here should be the `astra.simulator.flight`_ class, for running Monte Carlo simulations of flights for a given set of inputs, and the `astra.target_landing.targetFlight <#astra.target_landing.targetFlight>`_ class, for optimizing the flight parameters given a target landing site.

Flight Profiles
---------------

The `flightProfile` is used to store data about the input parameters and flight path obtained. Additionally, `targetProfile` contains data about the multiple objectives for optimizing the flight (distance from a landing site, cost, and flight duration) via the `astra.target_landing.targetFlight <#astra.target_landing.targetFlight>`_.

.. autoclass:: astra.simulator.flightProfile
	:members:
	:undoc-members:
	:show-inheritance:

.. autoclass:: astra.target_landing.targetProfile
	:members:
	:undoc-members:
	:show-inheritance:

Flight 
------

.. autoclass:: astra.simulator.flight
	:members:
	:undoc-members:
	:show-inheritance:

Target Landing Optimizer
------------------------

.. autoclass:: astra.target_landing.targetFlight
	:members:
	:undoc-members:
	:show-inheritance:

Environments
------------

.. autoclass:: astra.weather.environment
	:members:
	:undoc-members:
	:show-inheritance:

.. autoclass:: astra.weather.forecastEnvironment
	:members:
	:undoc-members:
	:show-inheritance:

.. autoclass:: astra.weather.soundingEnvironment
	:members:
	:undoc-members:
	:show-inheritance:


GFS_Handler
-----------
.. autoclass:: astra.GFS.GFS_Handler
	:members:
	:undoc-members:
	:show-inheritance:

Linear 4D interpolator
----------------------
.. autoclass:: astra.interpolate.Linear4DInterpolator
	:members:
	:undoc-members:
	:show-inheritance:


astra.flight_tools module
-------------------------

.. automodule:: astra.flight_tools
	:members:
	:undoc-members:


astra.global_tools module
-------------------------
.. automodule:: astra.global_tools
	:members:
	:undoc-members:

Available Balloons and Parachutes
---------------------------------
.. automodule:: astra.available_balloons_parachutes
	:members:
	:undoc-members:
