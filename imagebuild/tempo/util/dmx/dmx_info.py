#!/usr/bin/env python

import sys
import argparse, struct

try:
  import tempo_utils
except:
  print "tmepo_utils is required to run this script"
  print "make sure it is in your python path"
  print "or install it from https://github.com/demorest/tempo_utils"
  sys.exit()

parser = argparse.ArgumentParser(description="List infomration about TOAs in each DMX range")
parser.add_argument("parfile",  help="name of .par file")
parser.add_argument("timfile",  help="name of .tim file")
parser.add_argument("-u","--unsorted", action="store_true",
                     help="print list DMX number order (don't sort by MJD)")
 
args = parser.parse_args()
parfile  = args.parfile
timfile  = args.timfile
unsorted = args.unsorted

# read in the DMX ranges

fpar = open(parfile,"r")
dmxr = {}
for s in fpar:
  if s.startswith("DMXR"):
    ss = s.split()
    parname = ss[0]
    if s.startswith("DMXR"):
      irange = int(ss[0].split("_")[-1])
      if not irange in dmxr:
        dmxr[irange] = [0.,0.]
      parval  = float(ss[1].replace("D","E"))
      if s.startswith("DMXR1"):
        dmxr[irange][0] = parval
      elif s.startswith("DMXR2"):
        dmxr[irange][1] = parval
fpar.close()


# read in the MJDs and INFO flags

toas = tempo_utils.read_toa_file(timfile)

# loop through the MJDs and INFO flags,
# develop a list of INFO strings,
# and accumulate the number of TOAs in
# each DMX range segretaged by INFO string.

infolist = []
orphans = []
norange = []

dmxn = {}
for k in dmxr:
  dmxn[k] = []

for t in toas:
  if not t.is_toa():
    continue
  t.parse_line()
  mjd = t.mjd
  info = t.flags['f']

  if not info in infolist:
    infolist.append(info)
    for k in dmxn:
      dmxn[k].append(0)
    norange.append(0)
  iinfo = infolist.index(info) 

  idx = -1
  for k in dmxr:
    if mjd>=dmxr[k][0] and mjd<=dmxr[k][1]:
      dmxn[k][iinfo] += 1
      break
  else:
    norange[iinfo] += 1
    orphans.append(mjd)  

print "# information about DMX ranges"
print "#"
print "# parfile:  ",parfile
print "# timfile:  ",timfile 
print "#"
print "# Shorthand used in table below:"
for i,info in zip(range(len(infolist)),infolist):
  print "# I%2.2d ==  %s" % (i,info)
print "#"

if len(orphans)>0:
  print "# TOAs not in DMX ranges:"
  print ("# %4s %8s  %8s  %6s  "+("   I%2.2d ")*len(infolist)) % \
      tuple( [ "" , "","","" ] + range(len(infolist)) )
  print ("  %3s  %8s  %8s  %6s  "+(" %5d ")*len(norange)) % tuple( ["","","",""] + norange)
  print "# A full list of these TOAs is below the table"

print "#"


print ("# %4s %10s  %10s  %10s  "+("   I%2.2d ")*len(infolist)) % \
    tuple( [ "DMXn" , "  start","   end"," length" ] + range(len(infolist)) )

klist = dmxr.keys()
klist.sort()

for k in klist:
  print ("  %3d  %10.4f  %10.4f  %10.4f  "+(" %5d ")*len(dmxn[k])) % \
            tuple( [k,dmxr[k][0],dmxr[k][1],dmxr[k][1]-dmxr[k][0]] + dmxn[k] )

if len(orphans)>0:
  print "#"
  print "# Dates of TOAs not in any DMX range:"
  for mjd in orphans:
    print "#    %f" %(mjd)


