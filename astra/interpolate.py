# coding=utf-8

"""
This module defines the class Linear4DInterpolator used to interpolate data
from the Global Forecast System.

The Linear4DInterpolator is a simple 4D linear interpolator, designed to work
on rectangular grids only (they don't need to be necessarily uniform).
Linear4DInterpolator essentially isolates the 4D hyper-tetrahedron inside
which the requested point lies and performs a quadrilinear interpolation using
the 16 vertices of that hyper-tetrahedron only.

This method has proven much faster than scipy's
interpolate.LinearNDInterpolator, which is based on Qhull's triangulation
method.


University of Southampton
Niccolo' Zapponi, nz1g10@soton.ac.uk, 22/04/2013
"""

__author__ = "Niccolo' Zapponi, University of Southampton, nz1g10@soton.ac.uk"

from math import floor, ceil
import numpy
from six.moves import builtins

# Pass through the @profile decorator if line profiler (kernprof) is not in use
try:
    builtins.profile
except AttributeError:
    def profile(func):
        return func


class Linear4DInterpolator(object):
    """
    The Linear4DInterpolator is a simple 4D linear interpolator, designed to
    work on rectangular grids only (they don't need to be necessarily uniform).

    Linear4DInterpolator essentially isolates the 4D hyper-tetrahedron inside
    which the requested point lies and performs a quadrilinear interpolation
    using the 16 vertices of that hyper-tetrahedron only. This method has
    proven much faster than scipy's interpolate.LinearNDInterpolator, which is
    based on Qhull's triangulation method, very slow on rectangular grids.

    Parameters
    ----------
    data : numpy array (4D)
        the 4D matrix containing data to be interpolated
    data_map : list
        contains four lists and four dictionaries, formed of the real world
        values corresponding to each data matrix axis and their reverse mapping
        dictionaries. For example, for a 3x2x4x1 matrix, coordinates should be
        a list of four lists and four dictionaries, the first one being a list
        of the coordinates corresponding to the 3 'rows' of the first axis, the
        second one a list of the coordinates corresponding to the 2 'rows' of
        the second axis, etc. Then, the first dictionary should contain a map
        between real world values and 'row' number for the first axis (real
        world value is the dictionary's key and the row number is the
        corresponding value), and so on.
        Note: to make sure this is formatted correctly, use the GFS.GFS_Map to
        prepare the mapping and then use its mapCoordinates() method to
        generate the variable to be used here as the map input.

    Notes
    -----
    * Although this implementation of the interpolator has been specifically 
    designed for use with the Global Forecast System (and therefore expects
    latitude,longitude,pressure,time coordinates), it can actually be used with
    any 4D rectangular grid of data, provided that the appropriate data map is
    passed (see below).

    * For more information about quadrilinear interpolation, see
    http://en.wikipedia.org/wiki/Trilinear_interpolation. The same concept can
    be extended to 4 dimensions.

    :Example:
        >>> myInterpolator = Linear4DInterpolator(my4DMatrix,
            [axis0map,axis1map,axis2map,axis3map,axis0dict,axis1dict,axis2dict,
            axis3dict])

        >>> # Request data (requires latitude, longitude, pressure and time):
        >>> myInterpolator(lat, lon, press, time)
    """

    @profile
    def __init__(self, data, data_map):
        # Store data
        self.data = data
        self.dmap = data_map

        # Extract data bounds from map
        self.min = numpy.array([
            numpy.array(self.dmap[0]).min(),
            numpy.array(self.dmap[1]).min(),
            numpy.array(self.dmap[2]).min(),
            numpy.array(self.dmap[3]).min()
        ])
        self.max = numpy.array([
            numpy.array(self.dmap[0]).max(),
            numpy.array(self.dmap[1]).max(),
            numpy.array(self.dmap[2]).max(),
            numpy.array(self.dmap[3]).max()
        ])

        # Extend minimum longitude bound if the worldwide data has been passed
        if self.max[1] == 180 and self.min[1] == -179.5:
            self.min[1] = -180

        # Calculate longitude step to account for HD and SD data
        self.lonStep = self.dmap[1][1] - self.dmap[1][0]

    def __call__(self, lat, lon, press, time):
        # Check if within bounds and switch to nearest neighbour if not
        lat = numpy.clip(lat, self.min[0], self.max[0])
        lon = numpy.clip(lon, self.min[1], self.max[1])
        press = numpy.clip(press, self.min[2], self.max[2])
        time = numpy.clip(time, self.min[3], self.max[3])

        # Find closest indices and subtract 1 if the upper limit is being reached, to avoid a KeyError
        i = numpy.digitize([lat], self.dmap[0])[0]
        if i == len(self.dmap[0]):
            i -= 1
        idx0 = [i - 1, i]
        i = numpy.digitize([press], self.dmap[2])[0]
        if i == len(self.dmap[1]):
            i -= 1
        idx2 = [i - 1, i]
        i = numpy.digitize([time], self.dmap[3])[0]
        if i == len(self.dmap[3]):
            i -= 1
        idx3 = [i - 1, i]

        # Reverse mapping for longitude, to account for crossing the 180th meridian
        lonGrid = [floor(lon / self.lonStep) * self.lonStep, ceil(lon / self.lonStep) * self.lonStep]
        if lonGrid[0] == -180:
            lonGrid[0] = 180

        if lonGrid[0] == lonGrid[1]:
            try:
                idx1 = [self.dmap[5][lonGrid[0]], self.dmap[5][lonGrid[1] + self.lonStep]]
            except KeyError:
                idx1 = [self.dmap[5][lonGrid[0] - self.lonStep], self.dmap[5][lonGrid[1]]]
        else:
            idx1 = [self.dmap[5][lonGrid[0]], self.dmap[5][lonGrid[1]]]

        # Calculate the normalized distance of the requested point from the lower vertex, along all 4 axes
        frac0 = 1 - abs((lat - self.dmap[0][idx0[0]]) / (self.dmap[0][idx0[1]] - self.dmap[0][idx0[0]]))
        frac1 = 1 - abs((lon - self.dmap[1][idx1[0]]) / (self.dmap[1][idx1[1]] - self.dmap[1][idx1[0]]))
        frac2 = 1 - abs((press - self.dmap[2][idx2[0]]) / (self.dmap[2][idx2[1]] - self.dmap[2][idx2[0]]))
        frac3 = 1 - abs((time - self.dmap[3][idx3[0]]) / (self.dmap[3][idx3[1]] - self.dmap[3][idx3[0]]))

        # Interpolate (one dimension at a time)
        # 1st dimension
        tx000 = frac0 * self.data[idx0[0], idx1[0], idx2[0], idx3[0]] + (1 - frac0) * self.data[
            idx0[1], idx1[0], idx2[0], idx3[0]]
        tx001 = frac0 * self.data[idx0[0], idx1[0], idx2[0], idx3[1]] + (1 - frac0) * self.data[
            idx0[1], idx1[0], idx2[0], idx3[1]]
        tx010 = frac0 * self.data[idx0[0], idx1[0], idx2[1], idx3[0]] + (1 - frac0) * self.data[
            idx0[1], idx1[0], idx2[1], idx3[0]]
        tx011 = frac0 * self.data[idx0[0], idx1[0], idx2[1], idx3[1]] + (1 - frac0) * self.data[
            idx0[1], idx1[0], idx2[1], idx3[1]]
        tx100 = frac0 * self.data[idx0[0], idx1[1], idx2[0], idx3[0]] + (1 - frac0) * self.data[
            idx0[1], idx1[1], idx2[0], idx3[0]]
        tx101 = frac0 * self.data[idx0[0], idx1[1], idx2[0], idx3[1]] + (1 - frac0) * self.data[
            idx0[1], idx1[1], idx2[0], idx3[1]]
        tx110 = frac0 * self.data[idx0[0], idx1[1], idx2[1], idx3[0]] + (1 - frac0) * self.data[
            idx0[1], idx1[1], idx2[1], idx3[0]]
        tx111 = frac0 * self.data[idx0[0], idx1[1], idx2[1], idx3[1]] + (1 - frac0) * self.data[
            idx0[1], idx1[1], idx2[1], idx3[1]]

        # 2nd dimension
        txy00 = frac1 * tx000 + (1 - frac1) * tx100
        txy01 = frac1 * tx001 + (1 - frac1) * tx101
        txy10 = frac1 * tx010 + (1 - frac1) * tx110
        txy11 = frac1 * tx011 + (1 - frac1) * tx111

        # 3rd dimension
        txyz0 = frac2 * txy00 + (1 - frac2) * txy10
        txyz1 = frac2 * txy01 + (1 - frac2) * txy11

        # 4th dimension
        result = frac3 * txyz0 + (1 - frac3) * txyz1

        return result