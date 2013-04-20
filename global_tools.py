# coding=utf-8

"""
global_tools.py
ASTRA High Altitude Balloon Flight Planner

DESCRIPTION
--------------

Global Tools
A collection of small global tools.
These tools generally perform basic unit conversions, provide simple data and analyze lists.


USAGE
--------------

The functions available are:

    feet2m(length)
        Convert a length in feet to meters.
        Returns a float with the length [m]

        Parameters:
            length: length in feet to be converted

    m2feet(length)
        Convert a length in meters to feet.
        Returns a float with the length [ft]

        Parameters:
            length: length in meters to be converted

    kel2c(temperature)
        Convert a temperature in Kelvin to degrees Celsius.
        Returns a float with the temperature [degC]

        Parameters:
            temperature: temperature in Kelvin to be converted

    c2kel(temperature)
        Convert a temperature in degrees Celsius to Kelvin.
        Returns a float with the temperature [K]

        Parameters:
            temperature: temperature in Celsius to be converted

    pa2mbar(pressure)
        Convert a pressure in Pascal to millibar
        Returns a float with the pressure [mbar]

        Parameters:
            pressure: pressure in Pa to be converted

    dirspeed2uv(windDirection,windSpeed,resultType=None)
        Convert wind direction and speed to u and v components [same units as windSpeed].
        If resultType is 'u', this function returns only the u component.
        If resultType is 'v', this function returns only the v component.
        If resultType is None, both u and v components are returned

        Parameters:
            windDirection: direction of the wind in degrees clockwise from the north. Direction indicates where the
                           wind is blowing FROM.
            windSpeed: speed of the wind [any units]
            resultType: 'u'|'v'|None, see above.

    uv2dirspeed(u,v)
        Convert u and v wind components to wind direction and speed.
        Returns a list: [windDirection, windSpeed], where wind direction is in degrees clockwise from the north,
        wind speed has same units as u and v components.

        Parameters:
           u: wind's u-component [any units]
           v: wind's v-component [any units]

    m2deg(dLat,dLon,latitude)
        Converts meters to degrees of latitude and longitude.
        Returns a list: [deltaLat, deltaLon], both in degrees. Sign convention: deltaLat positive to the north,
        deltaLon positive to the west.

        Parameters:
           dLat: distance in the direction of meridians [m]. Sign convention: positive to the north
           dLon: distance in the direction of parallels [m]. Sign convention: positive to the east
           latitude: current latitude [deg], used to calculate the distance between meridians.

    deg2m(degLat,degLon,latitude)
        Convert degrees of latitude and longitude to meters.
        Returns a list: [dLat,dLon], both in meters, representing the distance in the direction of the meridians and of
        the parallels respectively.

        Parameters:
            degLat: distance in degrees of latitude
            degLon: distance in degrees of longitude
            latitude: current latitude [deg], used to calculate the distance between meridians.

    prettySeconds(seconds)
        Convert seconds in hours, minutes and seconds.
        Returns a list: [hours, minutes, seconds]

        Parameters:
            seconds: number of seconds to be converted

    getUTCOffset(latitude, longitude, timestamp)
        Use the Google Maps API Time Zone service to obtain UTC offset information about
        the given location.
        Returns a string in the format '+HHMM'

        Parameters:
            latitude: latitude of the location
            longitude: longitude of the location
            date_time: datetime object

    find_nearest_index(array,value)
        Find the index of the entry in the array which is nearest to value.
        Returns an integer with the index.

        Parameters:
            array: a 1D numpy array of values
            value: the target value

    ISAatmosphere(altitude=None,temperature=None,density=None,pressure=None,speedOfSound=None)
        Return ISA atmospheric conditions for the input parameters given.
        Returns a list: [altitude,temperature,density,pressure,speedOfSound], with altitude [m], temperature [degC],
        density [kg/m3], pressure [mbar], speedOfSound [m/s].

        Parameters:
            (Note: Either the altitude or the temperature MUST be given!)

            altitude: [m]
            temperature: [degC]
            density: [kg/m3]
            pressure: [mbar]
            speedOfSound: [m/s]

        Note: if only the temperature is given as an argument and this is higher than 489.8 degrees Celsius, two
        possible ISA solutions exist. In this case, the function will return a list as follows:
        [[altitude0, altitude1],temperature,[density0,density1],[pressure0,pressure1],[speedOfSound0, speedOfSound1]].


University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from math import sqrt, exp, sin, cos, radians, atan, atan2, tan, pi
import numpy


def feet2m(lengthFeet):
    """
    feet2m(length)
        Convert a length in feet to meters.
        Returns a float with the length [m]

        Parameters:
            length: length in feet to be converted
    """
    return lengthFeet * 0.3048


def m2feet(lengthM):
    """
    m2feet(length)
        Convert a length in meters to feet.
        Returns a float with the length [ft]

        Parameters:
            length: length in meters to be converted
    """
    return lengthM / 0.3048


def kel2c(tempK):
    """
    kel2c(temperature)
        Convert a temperature in Kelvin to degrees Celsius.
        Returns a float with the temperature [degC]

        Parameters:
            temperature: temperature in Kelvin to be converted
    """
    return tempK - 273.15


def c2kel(tempC):
    """
    c2kel(temperature)
        Convert a temperature in degrees Celsius to Kelvin.
        Returns a float with the temperature [K]

        Parameters:
            temperature: temperature in Celsius to be converted
    """
    return tempC + 273.15


def pa2mbar(pressPa):
    """
    pa2mbar(pressure)
        Convert a pressure in Pascal to millibar
        Returns a float with the pressure [mbar]

        Parameters:
            pressure: pressure in Pa to be converted
    """
    return pressPa * 0.01


def dirspeed2uv(windDirection, windSpeed, resultType=None):
    """
    dirspeed2uv(windDirection,windSpeed,resultType=None)
        Convert wind direction and speed to u and v components [same units as windSpeed].
        If resultType is 'u', this function returns only the u component.
        If resultType is 'v', this function returns only the v component.
        If resultType is None, both u and v components are returned

        Parameters:
            windDirection: direction of the wind in degrees clockwise from the north. Direction indicates where the
                           wind is blowing FROM.
            windSpeed: speed of the wind [any units]
            resultType: 'u'|'v'|None, see above.
    """
    v = -windSpeed * cos(radians(windDirection))
    u = -windSpeed * sin(radians(windDirection))

    if resultType == 'u':
        return u
    elif resultType == 'v':
        return v
    elif resultType == 'uv' or resultType is None:
        return u, v
    else:
        print 'Unexpected argument format for resultType. Please refer to documentation.'


def uv2dirspeed(u, v):
    """
    uv2dirspeed(u,v)
        Convert u and v wind components to wind direction and speed.
        Returns a list: [windDirection, windSpeed], where wind direction is in degrees clockwise from the north,
        wind speed has same units as u and v components.

        Parameters:
           u: wind's u-component [any units]
           v: wind's v-component [any units]
    """

    windSpeed = sqrt(u ** 2 + v ** 2)
    windDirDeg = (180 / pi) * atan2(-u, -v)
    if windDirDeg < 0:
        windDirDeg += 360

    return windDirDeg, windSpeed


def m2deg(dLat, dLon, latitude):
    """
    m2deg(dLat,dLon,latitude)
        Converts meters to degrees of latitude and longitude.
        Returns a list: [deltaLat, deltaLon], both in degrees. Sign convention: deltaLat positive to the north,
        deltaLon positive to the west.

        Parameters:
           dLat: distance in the direction of meridians [m]. Sign convention: positive to the north
           dLon: distance in the direction of parallels [m]. Sign convention: positive to the east
           latitude: current latitude [deg], used to calculate the distance between meridians.
    """

    # Equatorial radius
    R = 6378137 #[m]

    # One degree of latitude and longitude
    beta = atan(0.99664719 * tan(radians(abs(latitude))))
    oneDegLon = (pi / 180) * R * cos(beta)
    oneDegLat = 111100

    return dLat / oneDegLat, dLon / oneDegLon


def deg2m(degLat, degLon, latitude):
    """
    deg2m(degLat,degLon,latitude)
        Convert degrees of latitude and longitude to meters.
        Returns a list: [dLat,dLon], both in meters, representing the distance in the direction of the meridians and of
        the parallels respectively.

        Parameters:
            degLat: distance in degrees of latitude
            degLon: distance in degrees of longitude
            latitude: current latitude [deg], used to calculate the distance between meridians.
    """

    # Equatorial radius
    R = 6378137 #[m]

    # One degree of latitude and longitude
    beta = atan(0.99664719 * tan(radians(abs(latitude))))

    # One deg of Lon in meters
    oneDegLon = (pi / 180) * R * cos(beta)

    # One deg of Lat in meters
    oneDegLat = 111100

    return degLat * oneDegLat, degLon * oneDegLon


def prettySeconds(seconds):
    """
    prettySeconds(seconds)
        Convert seconds in hours, minutes and seconds.
        Returns a list: [hours, minutes, seconds]

        Parameters:
            seconds: number of seconds to be converted
    """

    is_negative = seconds < 0

    seconds = abs(seconds)
    hours = seconds // 3600
    minutes = (seconds // 60) % 60
    secs = seconds - minutes * 60 - hours * 3600

    if is_negative:
        if hours == 0:
            if minutes == 0:
                secs *= -1
            minutes *= -1
        hours *= -1

    return hours, minutes, secs


def find_nearest_index(array, value):
    """
    find_nearest_index(array,value)
        Find the index of the entry in the array which is nearest to value.
        Returns an integer with the index.

        Parameters:
            array: a 1D numpy array of values
            value: the target value
    """

    return numpy.abs(array - value).argmin()


def getUTCOffset(latitude, longitude, date_time):
    """
    getUTCOffset(latitude, longitude, timestamp)
        Use the Google Maps API Time Zone service to obtain UTC offset information about
        the given location.
        Returns a float with the UTC offset in hours

        Parameters:
            latitude: latitude of the location
            longitude: longitude of the location
            date_time: datetime object
    """

    import urllib2, time, json

    timestamp = time.mktime(date_time.timetuple())

    requestURL = "https://maps.googleapis.com/maps/api/timezone/json?location=%f,%f&timestamp=%d&sensor=true" % (
        latitude,
        longitude,
        timestamp
    )


    try:
        HTTPresponse = urllib2.urlopen(requestURL)
    except:
        return 0

    data = json.JSONDecoder().decode(HTTPresponse.read())
    if str(data['status']) == 'OK':
        total_seconds = data['dstOffset'] + data['rawOffset']
        return total_seconds / 3600.
    else:
        return 0


def ISAatmosphere(altitude=None, temperature=None, density=None, pressure=None, speedOfSound=None):
    """
    ISAatmosphere(altitude=None,temperature=None,density=None,pressure=None,speedOfSound=None)
        Return ISA atmospheric conditions for the input parameters given.
        Returns a list: [altitude,temperature,density,pressure,speedOfSound], with altitude [ft], temperature [degC],
        density [kg/m3], pressure [mbar], speedOfSound [m/s].

        Parameters:
            (Note: Either the altitude or the temperature MUST be given!)

            altitude: [ft]
            temperature: [degC]
            density: [kg/m3]
            pressure: [mbar]
            speedOfSound: [m/s]

        Note: if only the temperature is given as an argument and this is higher than 489.8 degrees Celsius, two
        possible ISA solutions exist. In this case, the function will return a list as follows:
        [[altitude0, altitude1],temperature,[density0,density1],[pressure0,pressure1],[speedOfSound0, speedOfSound1]].
    """

    # Constants

    # Level limits
    Level1 = 11000 # m
    Level2 = 20000 # m
    Level3 = 32000 # m
    Level4 = 47000 # m
    Level5 = 51000 # m
    #Level6 = 71000 # m

    # Level 1
    A1 = 288.15
    B1 = -6.5e-3
    C1 = 8.9619638
    D1 = -0.20216125e-3
    E1 = 5.2558797
    I1 = 1.048840
    J1 = -23.659414e-6
    L1 = 4.2558797

    # Level 2
    A2 = 216.65
    B2 = 0
    F2 = 128244.5
    G2 = -0.15768852e-3
    M2 = 2.0621400
    N2 = -0.15768852e-3

    # Level 3
    A3 = 196.65
    B3 = 1e-3
    C3 = 0.70551848
    D3 = 3.5876861e-6
    E3 = -34.163218
    I3 = 0.9726309
    J3 = 4.94600e-6
    L3 = -35.163218

    # Level 4
    A4 = 139.05
    B4 = 2.8e-3
    C4 = 0.34926867
    D4 = 7.0330980e-6
    E4 = -12.201149
    I4 = 0.84392929
    J4 = 16.993902e-6
    L4 = -13.201149

    # Level 5
    A5 = 270.65
    B5 = 0
    F5 = 41828.42
    G5 = -0.12622656e-3
    M5 = 0.53839563
    N5 = -0.12622656e-3

    R = 287.05287


    # CASE 1: All parameters known. Problem over-defined.
    # Arguments returned unchanged.
    if altitude is not None and temperature is not None and density is not None and pressure is not None and speedOfSound is not None:
        print 'Overconstrained ISA problem. Arguments returned unchanged.'
        return altitude, temperature, density, pressure, speedOfSound

    # CASE 2: All parameters unknown. Problem under-defined.
    # Arguments returned unchanged.
    if altitude is None and temperature is None and density is None and pressure is None and speedOfSound is None:
        print 'Underconstrained ISA problem. Arguments returned unchanged.'
        return altitude, temperature, density, pressure, speedOfSound

    # CASE 3: Only altitude known.
    if altitude is not None and temperature is None and density is None and pressure is None and speedOfSound is None:

        if altitude < 0:
            print 'Altitude out of bounds, set to 0.'
            altitude = 0
        elif feet2m(altitude) > Level5:
            print 'Altitude out of bounds, set to 32km.'
            altitude = m2feet(Level5)

        altitudeM = feet2m(altitude)
        # Determine the appropriate layer of the atmosphere
        if altitudeM < Level1:
            # Troposphere, temperature decreases linearly
            temperature = kel2c(A1 + B1 * altitudeM)
            pressure = pa2mbar((C1 + D1 * altitudeM) ** E1)
            density = (I1 + J1 * altitudeM) ** L1
        elif altitudeM < Level2:
            # Lower stratosphere, temperature is constant
            temperature = kel2c(A2 + B2 * altitudeM)
            pressure = pa2mbar(F2 * exp(G2 * altitudeM))
            density = M2 * exp(N2 * altitudeM)
        elif altitudeM < Level3:
            # Upper stratosphere, temperature is increasing
            temperature = kel2c(A3 + B3 * altitudeM)
            pressure = pa2mbar((C3 + D3 * altitudeM) ** E3)
            density = (I3 + J3 * altitudeM) ** L3
        elif altitudeM < Level4:
            # Between 32 and 47 km
            temperature = kel2c(A4 + B4 * altitudeM)
            pressure = pa2mbar((C4 + D4 * altitudeM) ** E4)
            density = (I4 + J4 * altitudeM) ** L4
        else:
            # Between 47 and 51 km
            temperature = kel2c(A5 + B5 * altitudeM)
            pressure = pa2mbar(F5 * exp(G5 * altitudeM))
            density = M5 * exp(N5 * altitudeM)
        speedOfSound = sqrt(1.4 * R * c2kel(temperature))

        return altitude, temperature, density, pressure, speedOfSound

    # CASE 4: Only temperature known.
    if altitude is None and temperature is not None and density is None and pressure is None and speedOfSound is None:
        # There is either one solution...
        thresholdTemperature = kel2c(A3 + B3 * Level3)
        if temperature > thresholdTemperature:
            altitude = m2feet((c2kel(temperature) - A1) / B1)
            speedOfSound = sqrt(1.4 * R * c2kel(temperature))
            pressure = pa2mbar((C1 + D1 * feet2m(altitude)) ** E1)
            density = (I1 + J1 * feet2m(altitude)) ** L1

            return altitude, temperature, density, pressure, speedOfSound
        # ... or two
        elif c2kel(temperature) >= A2:
            # Initialize arrays
            altitude = [0.0, 0.0]
            pressure = [0.0, 0.0]
            density = [0.0, 0.0]

            # Lower altitude solution
            altitude[0] = m2feet((c2kel(temperature) - A1) / B1)
            pressure[0] = pa2mbar((C1 + D1 * feet2m(altitude[0])) ** E1)
            density[0] = (I1 + J1 * feet2m(altitude[0])) ** L1
            # Higher altitude solution
            altitude[1] = m2feet((c2kel(temperature) - A3) / B3)
            pressure[1] = pa2mbar((C3 + D3 * feet2m(altitude[1])) ** E3)
            density[1] = (I3 + J3 * feet2m(altitude[1])) ** L3
            speedOfSound = sqrt(1.4 * R * c2kel(temperature))

            return altitude, temperature, density, pressure, speedOfSound
        else:
            print 'Temperature out of bounds.'
            return altitude, temperature, density, pressure, speedOfSound