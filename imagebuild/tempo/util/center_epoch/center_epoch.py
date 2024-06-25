#!/usr/bin/env python

import sys
import argparse, struct
import subprocess
import os

parser = argparse.ArgumentParser(description="Center period, position, and orbital epoch")
parser.add_argument("infile",   help="input .par file")
parser.add_argument("timfile",  help="input .tim file")
parser.add_argument("outfile",  help="output .par file")
parser.add_argument("-t","--tpoflag",  help="tempo flag(s)", default="")

args = parser.parse_args()
infile   = args.infile
timfile  = args.timfile
outfile  = args.outfile
tpoflag  = args.tpoflag

tmpfile  =    "_tmp.par_"
npulsefile = "_npulse.tmp_"

# make the pulse number file

if tpoflag=="":
  subprocess.call(["tempo","-f",infile,timfile,"-no",npulsefile])
else:
  subprocess.call(["tempo",tpoflag,"-f",infile,timfile,"-no",npulsefile])

# read the old .par file and make a new .par file with epochs moved
#  and a meager attempt to shift f0 and (where present) t0 or tasc near
#  its new value

fin  = open(infile,"r")
ftmp = open(tmpfile,"w")

uset0 = False
usetasc = False
usepb = False
usefb = False

for s in fin:
  ftmp.write(s)
  ss = s.split()
  if len(ss)==0:
    continue
  if len(ss)>1:
    par, val = ss[0:2]
  if par=="START":
    start = float(val)
  elif ss[0]=="FINISH":
    finish = float(val)
  elif ss[0]=="F0":
    f0 = float(val.replace("D","E"))
  elif ss[0]=="F1":
    f1 = float(val.replace("D","E"))
  elif ss[0]=="PEPOCH":
    pepoch = float(val.replace("D","E"))
  elif ss[0]=="PB":
    usepb = True
    pb = float(val.replace("D","E"))
  elif ss[0]=="FB0":
    usefb = True
    fb = float(val.replace("D","E"))
  elif ss[0]=="T0":
    t0 = float(val.replace("D","E"))
    uset0 = True
  elif ss[0]=="TASC":
    tasc = float(val.replace("D","E"))
    usetasc = True
  elif ss[0]=="PSR":
    psr = val


newepoch = int(0.5*(start+finish))
newf0 = f0 + f1*(newepoch-pepoch)*86400.
ftmp.write("POSEPOCH %20.10f\n" % newepoch)
ftmp.write("PEPOCH   %20.10f\n" % newepoch)
ftmp.write("F0       %20.15e 1\n" % newf0)

if usefb:
  pb = 1./(86400.*fb)  # fb is in s^-1
  usepb = True

if uset0:
  if usepb:
    newt0 = t0 + pb*int((newepoch-t0)/pb+0.5)
  else:
    print "Error: T0 is present but neither PB nor FB0 is set"
    sys.exit(1)
  ftmp.write("T0     %20.15e 1\n" % newt0)

if usetasc:
  if usepb:
    newtasc = tasc + pb*int((newepoch-tasc)/pb+0.5)
  else:
    print "Error: T0 is present but neither PB nor FB0 is set"
    sys.exit(1)
  ftmp.write("TASC   %20.15e 1\n" % newtasc)

ftmp.close()


# run tempo to get better new parameter values
# iterate (run twice) before moving to the desired outfile location

tpooutfile = psr+".par"
if tpoflag == "":
  subprocess.call(["tempo","-f",tmpfile,timfile,"-ni",npulsefile])
  subprocess.call(["tempo","-f",tpooutfile,timfile])
else:
  subprocess.call(["tempo",tpoflag,"-f",tmpfile,timfile,"-ni",npulsefile])
  subprocess.call(["tempo",tpoflag,"-f",tpooutfile,timfile])

os.rename(tpooutfile,outfile)


