#!/usr/bin/python
__author__="morganlnance"



def calc_mean_degrees( data ):
    """
    Given a list of <data> in degrees from -180 to 180, return the mean
    Since I can't use NumPy on Jazz
    :param data: list( data points )
    :return: float( mean )
    """
    '''
    # imports
    Found from https://rosettacode.org/wiki/Averages/Mean_angle
    from cmath import rect, phase
    from math import radians, degrees
    return degrees( phase( sum( rect( 1, radians( deg ) ) for deg in data ) / len( data ) ) )
    '''
    return sum( data ) / len( data )



def calc_stddev_degrees( data ):
    """
    Given a list of <data> in degrees from -180 to 180, return the sample standard deviation
    Since I can't use NumPy on Jazz
    :param data: list( data points )
    :return: float( standard deviation )
    """
    # imports
    from math import sqrt

    # calculate the mean
    N = len( data )
    mean = calc_mean_degrees( data )

    # calculate the sample standard deviation
    stddev = sqrt( sum( [ ( x - mean ) ** 2 for x in data ] ) / N - 1 )
    
    '''
    # log is ln here
    from math import log
    Found from https://rdrr.io/cran/circular/man/sd.circular.html
    stddev = -2 * log( mean )
    '''

    return stddev
