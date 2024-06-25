#!/usr/bin/env python

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Adjust DM and DMX_xxxx so median DMX value is zero")
parser.add_argument("infile",   help="input .par file")
parser.add_argument("outfile",  help="output .par file")

args = parser.parse_args()
infile   = args.infile
outfile  = args.outfile

with open(infile,"r") as fin:
  p=fin.readlines()

dm = False
dmxlist = []

for s in p:
  if s.startswith("DMX_"):
    dmxlist.append(float(s.split()[1].replace("D","E")))

dmxmed = np.round(np.median(dmxlist),decimals=6)

with open(outfile,"w") as fout:
  for s in p:
    if s.startswith("DM "):
      ss = s.split()
      dm = float(s.split()[1].replace("D","E"))+dmxmed
      dms = ("%23.6f" % dm).replace("e","D")
      # s = "DM "+dms+s[26:]
      s = "DM"
      s = s + ("%24.6f" % dm)
      if len(ss)>2:
        s = s + ("%3s" % ss[2])
      if len(ss)>3:
        s = s + ("%20s" % ss[3])
      if len(ss)>4:
        s = s + " ".join(ss[4:])
      s = s + "\n"
    if s.startswith("DMX_"):
      ss = s.split()
      dmx = float(ss[1].replace("D","E"))-dmxmed
      s = "%-8s" % (ss[0])
      s = s + ("%18.8e" % dmx).replace("e","D")
      if len(ss)>2:
        s = s + ("%3s" % ss[2])
      if len(ss)>3:
        s = s + ("%20s" % ss[3])
      if len(ss)>4:
        s = s + " ".join(ss[4:])
      s = s + "\n"
    fout.write(s)


