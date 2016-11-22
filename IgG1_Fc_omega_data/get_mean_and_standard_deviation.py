#!/usr/bin/python
__author__="morganlnance"


from util import calc_mean_degrees, calc_stddev_degrees
from omega1_and_omega2_data_N_linked_GlcNAc import omega1_data, \
    omega2_data


print "\nMean omega1:", calc_mean_degrees( omega1_data )
print "Standard deviation omega1:", calc_stddev_degrees( omega1_data )

print "\nMean omega2:", calc_mean_degrees( omega2_data )
print "Standard deviation omega2:", calc_stddev_degrees( omega2_data ), "\n"
