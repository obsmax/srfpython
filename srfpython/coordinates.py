import numpy as np


def haversine(loni, lati, lonj, latj, R=6371.):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    source https://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points
    lone, late : coordinates of the first point, in degrees
    lons, lats : coordinates of the second point, in degrees
    :return distance in km

    consistent with Ll2DGAB
    """
    # convert decimal degrees to radians
    q = np.pi / 180.

    # haversine formula
    dlon = (loni - lonj)
    dlat = (lati - latj)
    a = np.sin(q * dlat / 2.) ** 2. + np.cos(q * latj) * np.cos(q * lati) * np.sin(q * dlon / 2.) ** 2.
    c = 2. * np.arcsin(np.sqrt(a))
    km = R * c
    return km