#!/usr/bin/python
__author__="morganlnance"


import matplotlib
try:
    # this might only be a Mac thing
    matplotlib.use( "TKAgg" )
except:
    pass
import matplotlib.pyplot as plt
from omega1_and_omega2_data_N_linked_GlcNAc import omega1_data, \
    omega2_data


plt.scatter( omega1_data, omega2_data )
plt.xlim( [ -180, 180 ] )
plt.xlabel( "omega1" )
plt.ylim( [ -180, 180 ] )
plt.ylabel( "omega2" )
plt.grid( True )
plt.plot( [ 0 for x in range( 360 + 1 ) ], range( -180, 180 + 1 ), color="black" )
plt.plot( range( -180, 180 + 1 ), [ 0 for x in range( 360 + 1 ) ], color="black" )
plt.title( "Compiled IgG Fc statistics for ASN-GlcNAc omegas" )

plt.show()
