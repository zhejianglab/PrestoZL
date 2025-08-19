
from __future__ import print_function
from builtins import zip
from builtins import range
import os

def myexecute(cmd):
    print("'%s'"%cmd)
    os.system(cmd)

# By default, do not output subbands
outsubs = False

nsub = 1024

basename = 'Dec+0640_00_05_arcdrift-M04_0642'
rawfiles = '/presto_data/test/Dec/Dec+0640_00_05/Dec+0640_00_05_arcdrift-M04_0642.fits'

# dDM steps from DDplan.py
dDMs        = [0.1, 0.2, 0.3, 0.5]
# dsubDM steps
dsubDMs     = [75.8, 151.6, 245.39999999999998, 409.0]
# downsample factors
downsamps   = [1, 2, 4, 8]
# number of calls per set of subbands
subcalls    = [9, 4, 4, 2]
# The low DM for each set of DMs
startDMs    = [59.4, 741.6, 1348.0, 2329.6]
# DMs/call
dmspercalls = [758, 758, 818, 818]


# Loop over the DDplan plans
for dDM, dsubDM, dmspercall, downsamp, subcall, startDM in zip(dDMs, dsubDMs, dmspercalls, downsamps, subcalls, startDMs):
    # Loop over the number of calls
    for ii in range(subcall):
        subDM = startDM + (ii+0.5)*dsubDM
        loDM = startDM + ii*dsubDM
        if outsubs:
            # Get our downsampling right
            subdownsamp = downsamp // 2
            datdownsamp = 2
            if downsamp < 2: subdownsamp = datdownsamp = 1
            # First create the subbands
            myexecute("prepsubband -sub -subdm %.2f -nsub %d -downsamp %d -o %s %s" %
                      (subDM, nsub, subdownsamp, basename, rawfiles))
            # And now create the time series
            subnames = basename+"_DM%.2f.sub[0-9]*"%subDM
            myexecute("prepsubband -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (loDM, dDM, dmspercall, datdownsamp, basename, subnames))
        else:
            myexecute("prepsubband -nsub %d -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s" %
                      (nsub, loDM, dDM, dmspercall, downsamp, basename, rawfiles))
