#! /usr/bin/env python
from numpy import *
import sys

class dmxrange:
    def __init__(self, lofreqs, hifreqs):
        self.los = lofreqs
        self.his = hifreqs
        self.min = min(lofreqs + hifreqs)
        self.max = max(lofreqs + hifreqs)
    def sum_print(self):
        print "%8.2f-%8.2f (%6.2f days): " % (self.min, self.max, self.max-self.min),
        print "lofs:", self.los,
        print "hifs:", self.his

divide_freq = 1000.0 # MHz
offset = 0.5 # days (for making DMX ranges)
max_diff = 15.0 # days (if no freq match, throw out these dates)

fname = "oneday.resids"
if (len(sys.argv)>1): 
    fname = sys.argv[1]
MJDs, freqs = loadtxt(fname, usecols = (0,1), unpack=True)

loMJDs = MJDs[freqs < divide_freq]
hiMJDs = MJDs[freqs > divide_freq]
loMJDs.sort()
hiMJDs.sort()
print "There are %d dates with freqs > %.0f MHz" % (len(hiMJDs), divide_freq)
print "There are %d dates with freqs < %.0f MHz\n" % (len(loMJDs), divide_freq)

DMXs = []

good_his = set([])
bad_los = []
# Walk through all of the low freq obs
for ii, loMJD in enumerate(loMJDs):
    # find all the high freq obs within max_diff days
    # of the low freq obs
    hi_close = hiMJDs[fabs(hiMJDs - loMJD) < max_diff]
    # and where they are closer to this loMJD compared to the
    # other nearby ones
    if (ii > 0):
        diffs = fabs(hi_close - loMJD)
        lodiffs = fabs(hi_close - loMJDs[ii-1])
        hi_close = hi_close[diffs < lodiffs]
    if (ii < len(loMJDs)-1):
        diffs = fabs(hi_close - loMJD)
        hidiffs = fabs(hi_close - loMJDs[ii+1])
        hi_close = hi_close[diffs < hidiffs]
    if len(hi_close):  # add a DMXrange
        DMXs.append(dmxrange([loMJD], list(hi_close)))
        good_his = good_his.union(set(hi_close))
    else:
        bad_los.append(loMJD)

bad_los = set(bad_los)
saved_los = []
#print bad_los
# Now walk through the DMXs and see if we can't fit a bad_lo freq in
for bad_lo in bad_los:
    absmindiff = 2 * max_diff
    ind = 0
    for ii, DMX in enumerate(DMXs):
        if (fabs(bad_lo - DMX.min) < max_diff and
            fabs(bad_lo - DMX.max) < max_diff):
            mindiff = min(fabs(bad_lo - DMX.min), fabs(bad_lo - DMX.max))
            if mindiff < absmindiff:
                absmindiff = mindiff
                ind = ii
    if (absmindiff < max_diff):
        # print DMXs[ind].min, DMXs[ind].max, bad_lo
        DMXs[ind].los.append(bad_lo)
        # update the min and max vals
        DMXs[ind].min = min(DMXs[ind].los + DMXs[ind].his)
        DMXs[ind].max = max(DMXs[ind].los + DMXs[ind].his)
        saved_los.append(bad_lo)

# These are the low-freq obs we can't save
bad_los -= set(saved_los)
bad_los = sorted(list(bad_los))

# These are the high-freq obs we can't save
bad_his = set(hiMJDs) - good_his
bad_his = sorted(list(bad_his))

print "\n These are the 'good' ranges for DMX and days are low/high freq:"
for DMX in DMXs:
    DMX.sum_print()

print "\nRemove high-frequency data from these days:"
for hibad in bad_his:
    print "%8.2f" % hibad
print "\nRemove low-frequency data from these days:"
for lobad in bad_los:
    print "%8.2f" % lobad

print "\n Enter the following in your parfile"
print "-------------------------------------"
print "DMX         %.2f" % max_diff
oldmax = 0.0
for ii, DMX in enumerate(DMXs):
    print "DMX_%04d      0.0       %d" % (ii+1, 1 if ii else 0)
    print "DMXR1_%04d      %10.4f" % (ii+1, DMX.min-offset)
    print "DMXR2_%04d      %10.4f" % (ii+1, DMX.max+offset)
    if DMX.min < oldmax:
        print "Ack!  This shouldn't be happening!"
    oldmax = DMX.max
