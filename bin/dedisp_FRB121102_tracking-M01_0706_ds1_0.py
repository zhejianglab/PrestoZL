
from __future__ import print_function
from builtins import zip
from builtins import range
import os

def myexecute(cmd):
    print("'%s'"%cmd)
    os.system(cmd)

# By default, do not output subbands
outsubs = False

nsub = 3280

basename = 'FRB121102_tracking-M01_0706_ds1_0'
rawfiles = '/presto_data/test/FRB/FRB121102_tracking-M01_0706_ds1_0.fits'

# dDM steps from DDplan.py
dDMs        = [0.1, 0.2, 0.5]
# dsubDM steps
dsubDMs     = [262.2, 524.4, 1212.0]
# downsample factors
downsamps   = [1, 2, 4]
# number of calls per set of subbands
subcalls    = [3, 1, 1]
# The low DM for each set of DMs
startDMs    = [59.4, 846.0, 1370.4]
# DMs/call
dmspercalls = [2622, 2622, 2424]


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
