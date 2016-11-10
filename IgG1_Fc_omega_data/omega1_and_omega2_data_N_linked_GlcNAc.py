#!/usr/bin/python
__author__="morganlnance"


# data for structures collected via PyMol on native crystal structures
# omega1 data for N-linked GlcNAc
omega1_data = [ 
    #-154.566, -162.504,   # 3ay4 - IgG1 Fc G2 to FcgRIIIa
    -163.850, 176.923,    # 3ave - IgG1 Fc no paper, but it comes out as G0F2
     -164.387, -146.038,  # 5d4q - IgG1 Fc G2F1 ( check paper )
     -146.449, -146.340,  # 5d6d - IgG1 Fc G2F2 to FcgRIIIa ( check paper )
     -166.996, -171.113,  # 1h3x - IgG1 Fc G0F2
     -150.421, -155.169,  # 4w4n - IgG1 Fc G0F1 ( in crystal at least, check paper )
     -155.977, -141.183,  # 4wi2 - IgG1 Fc G0F2 ( in crstyal at least, paper to be published )
     -177.420, 164.920    # 3do3 - IgG1 Fc G2F2 ( in crystal at least, no paper )
     ]
# for standardization, make everything negative
omega1_data = [ deg - 360 if deg > 0 else deg for deg in omega1_data ]

# omega2 data for N-linked GlcNAc
omega2_data = [ 
    #59.198, 59.055,  # 3ay4
    53.823, 47.082,  # 3ave
    48.590, 63.976,  # 5d4q
    64.005, 63.988,  # 5d6d
    45.997, 59.196,  # 1h3x
    48.583, 56.529,  # 4w4n
    54.787, 46.017,  # 4wi2
    50.729, 66.251   # 3do3
    ]
# for standardization, make everything negative
#omega2_data = [ deg - 360 if deg > 0 else deg for deg in omega2_data ]
