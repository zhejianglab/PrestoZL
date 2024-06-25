import numpy as np
import string as st
import sys

# Initialize option flags
flags = {'-h':False, '-l':False, '-n':False, '-p':False, '-s':False, '-t':False}

def printhelp():
   print ""
   print "  Usage:"
   print "    Command line:  python pubpar.py [options] <parfile>"
   print "         IPython:  run pubpar.py [options] <parfile>"
   print ""
   print ""
   print "  This script takes a parfile output by a Tempo fit and displays"
   print "  the fitted parameters with their uncertainties in parentheses"
   print "  in a format appropriate for publication."
   print ""
   print ""
   print "  Options:"
   print "    -h  Show this useful help text."
   print ""
   print "    -l  Output exponential notation in LaTeX X$\\times 10^{Y}$"
   print "        format rather than simple XeY format."
   print ""
   print "    -n <precision>"
   print "        Force number of decimals in quoted error to <precision>,"
   print "        rather than the default, which is 1 digit, or 2 digits"
   print "        when the first digit is a 1."
   print ""
   print "    -p  Convert rotation frequencies and frequency derivatives to"
   print "        periods and period derivatives. F0 must be included."
   print "          F0 --> PER0, F1 --> PER1, etc."
   print ""
   print "    -s  Force scientific notation."
   print ""
   print "    -t  Leave Tempo errors as reported in parfile rather than"
   print "        doubling them."
   print ""

# Function to get exponent for a number's exponential form
def getexp(number):
   return np.floor(np.log10(np.abs(number)))

nargs = len(sys.argv)
# Skip first argument (the name of the program) and last argument (presumably
# the name of the parfile) in the following loop, which sets the options
precision = 1
j = 1
while j < (nargs-1):
   flags[sys.argv[j]] = True
   if sys.argv[j] == '-n':
      try:
         precision = np.int(sys.argv[j+1])
         if precision <= 0: sys.exit(1)
      except:
         print "Option '-n' must be followed by a positive integer. Use '-h' option for help."
         sys.exit(0)
      j += 1
   j += 1

# Show help text and exit if necessary (allow '-h' to be last or only argument)
if flags['-h'] or nargs == 1 or sys.argv[nargs-1] == '-h':
   printhelp()
   sys.exit(0)

# The parfile must be the last argument
parfile = sys.argv[nargs-1]

# Read in par file and pull out lines that appear to have been fit
# (in other words, those with four columns)
try:
   par = np.loadtxt(parfile, dtype='str', delimiter='$')
except:
   print "Could not open parfile. Use '-h' option for help."
   sys.exit(0)

# Assume 4 entries on a line mean it was potentially fit--if third entry
# was '1', it was presumably fit. If it's got colons, deal with it.
params = []
for row in par:
   line = st.split(row)
   if len(line) == 4:
      if line[2] == '1':
         colon_from_end = st.find(line[1][::-1], ':')
         if colon_from_end >= 0:
            start_i = len(line[1]) - colon_from_end
            if line[1][start_i] == '0':
               lookahead = 1
               while line[1][start_i + lookahead] == '0':
                  lookahead += 1
               if line[1][start_i + lookahead] == '.': start_i += lookahead - 1
               else: start_i += lookahead
            precolon = line[1][:start_i]
            line[1] = line[1][start_i:]
         else: precolon = ""
         params.append([line[0], np.float(st.replace(line[1], 'D', 'e')), \
            np.float(st.replace(line[3], 'D', 'e')), precolon])

# Convert frequency and frequency derivatives to period and period
# derivatives if flags['-p'] = True
already_gave_F0_warning = False
f0 = 0.0
deriv = -1
if flags['-p']:
   for i in range(0, len(params)):
      success = True
      freq = True
      name = params[i][0]
      ndigits = len(name) - 1
      if name[0] == 'F' and (ndigits == 1 or ndigits == 2):
         if ndigits == 1: deriv = np.int(name[1])
         else: deriv = np.int(name[1:3])
      else: freq = False
      if deriv == 0: f0 = params[i][1]
      try: per = -params[i][1] / (f0*f0)
      except: success = False
      if (freq and success):
         params[i][0] = 'PER' + np.str(deriv)
         if deriv == 0: params[i][1] = np.abs(per)
         else: params[i][1] = per
         params[i][2] /= (f0*f0)
      elif (freq and not success and not already_gave_F0_warning):
         print "F0 must come before F(>0) in order to convert frequencies to periods."
         already_gave_F0_warning = True

# Define dictionary of units
units = {\
'F0':'s^-1',\
'F1':'s^-2',\
'F2':'s^-3',\
'F3':'s^-4',\
'F4':'s^-5',\
'F5':'s^-6',\
'F6':'s^-7',\
'F7':'s^-8',\
'F8':'s^-9',\
'F9':'s^-10',\
'F10':'s^-11',\
'F11':'s^-12',\
'F12':'s^-13',\
'PER0':'s',\
'PER1':'',\
'PER2':'s^-1',\
'PER3':'s^-2',\
'PER4':'s^-3',\
'PER5':'s^-4',\
'PER6':'s^-5',\
'PER7':'s^-6',\
'PER8':'s^-7',\
'PER9':'s^-8',\
'PER10':'s^-9',\
'PER11':'s^-10',\
'PER12':'s^-11',\
'RA':'hrs, min, sec',\
'RAJ':'hrs, min, sec',\
'DEC':'deg, min, sec',\
'DECJ':'deg, min, sec',\
'PMRA':'mas yr^-1',\
'PMDEC':'mas yr^-1',\
'PMRV':'mas yr^-1',\
'BETA':'degrees',\
'LAMBDA':'degrees',\
'PMBETA':'mas yr^-1',\
'PMLAMBDA':'mas yr^-1',\
'PX':'mas',\
'DM':'pc cm^-3',\
'A1':'lightsec',\
'E':'',\
'EDOT':'10^-12 s^-1',\
'T0':'MJD',\
'TASC':'MJD',\
'PB':'days',\
'PBDOT':'10^-12',\
'XPBDOT':'10^-12',\
'GAMMA':'s',\
'PPNGAMMA':'',\
'SINI':'',\
'MTOT':'M_sun',\
'M2':'M_sun',\
'OM':'degrees',\
'OMDOT':'deg yr^-1',\
'OM2DOT':'rad s^-2',\
'XOMDOT':'deg yr^-1',\
'XDOT':'10^-12',\
'X2DOT':'s^-1',\
'DTHETA':'10^-6',\
'AFAC':'sin(eta)/sin(lambda)'\
}
# NOTE: units for DMn, DMX_???? and DMX1_???? are handled in an if-statement
# later on

# Format parameters with uncertainties in parentheses
output = []
for i in range(0, len(params)):
   prec = precision

   name = params[i][0]
   dat = params[i][1]
 
   # double Tempo 1-sigma errors unless flags['-t']=True
   if flags['-t']: sig = params[i][2]
   else: sig = 2*params[i][2]
   
   # if there were colons in the number, they'll be slapped back in at
   # the end
   precolon = params[i][3]

   exp_dat = getexp(dat)
   exp_sig = getexp(sig)

   # check for first sigma digit 1...
   sig_first_1 = np.int(sig / (10**exp_sig))
   sig_first_2 = np.int(sig / (10**(exp_sig-1)))
   sig_first_3 = np.int(sig / (10**(exp_sig-2)))
   if ((sig_first_1 == 1 and sig_first_3 - 100 < 95) or sig_first_2 >= 95) and not flags['-n']:
      prec = 2

   exp_shift = exp_dat

   dat /= 10**exp_shift
   sig /= 10**exp_shift

   exp_dat -= exp_shift
   exp_sig -= exp_shift

   # position of least significant digit
   exp_LSD = exp_sig - prec + 1

   dat_str = round(dat, np.int(np.abs(exp_LSD)))
   sig_str = np.int(round(sig / (10**exp_LSD)))

   exp_dat_str = getexp(dat_str)
   exp_sig_str = getexp(sig_str)

   fix_dat = (exp_dat_str > exp_dat)
   fix_sig = (exp_sig_str > exp_sig - exp_LSD)

   if fix_dat:
      exp_shift += 1
      dat_str /= 10
      exp_sig -= 1
      exp_sig_str -= 1
      exp_LSD -= 1

   if fix_sig:
      sig_str /= 10
      exp_sig += 1
      exp_LSD += 1

   if exp_sig >= 0:
      sig_str = sig_str * 10**exp_LSD

   nchars = np.int(np.abs(exp_LSD)) + 2

   # deal with uncertainty in first digit
   if exp_sig == exp_dat and getexp(dat_str) >= prec - 1:
      dat_str = np.int(dat_str)
      sig_str = np.int(sig_str)
      nchars = 1

   # add one to character count for negative sign
   if dat_str < 0: nchars += 1

   # remove scientific notation if sensible (and desired)
   if ((np.abs(exp_LSD) >= exp_shift and exp_shift > 0) or exp_shift == -1.0) and not flags['-s']:
      dat_str *= 10**exp_shift
      exp_sig += exp_shift
      exp_LSD += exp_shift
      # if we're making this a decimal < 1, there's a zero on the left, which
      # adds 1 to the number of characters
      if exp_shift < 0: nchars += 1
      exp_shift = 0

   # turn data and sigma into strings (conversion to string rounds, which can
   # be helpful, but not if we want really precise numbers, so if nchars > 13,
   # we jump through some hoops to get things to come out ok... maybe there
   # is a cleaner solution!)
   if nchars > 13:
      arbitrary_cutoff = 10
      s = repr(dat_str)
      part1 = s[0:arbitrary_cutoff]
      part2_flt = np.float(s[arbitrary_cutoff:])
      part2_exp = getexp(part2_flt)
      part2_flt /= 10**part2_exp
      part2 = np.str(np.int(round(part2_flt, nchars-1-arbitrary_cutoff) \
         * (10**(nchars-1-arbitrary_cutoff))))
      dat_str = part1 + part2
   else: dat_str = np.str(dat_str)
   sig_str = np.str(sig_str)

   # add 0's to end of data string if necessary
   dec_pos = st.find(dat_str, '.')
   if dec_pos >= 0:
      while len(dat_str) - dec_pos <= np.abs(exp_LSD):
         dat_str += '0'

   # create strings of info to be displayed
   if exp_shift >= 0:
      exp_str_sign = "+"
   else:
      exp_str_sign = "-"
   #
   if exp_shift == 0.0 and not flags['-s']:
      exp_str = ""
   else:
      if flags['-l']:
         exp_str_num = np.str(np.int(exp_shift))
         exp_str = "$\\times 10^{" + exp_str_num + "}$"
      else:
         exp_str_num = np.str(np.int(np.abs(exp_shift)))
         if np.abs(exp_shift) < 10:
            exp_str = "e" + exp_str_sign + "0" + exp_str_num
         else:
            exp_str = "e" + exp_str_sign + exp_str_num
   #
   str_a = name
   str_b = precolon + dat_str + "(" + sig_str + ")" + exp_str
   try: str_c = units[name]
   except:
      # check for some DM-related units before giving up
      if name[0:4] == "DMX_": str_c = "pc cm^-3"
      elif name[0:5] == "DMX1_": str_c = "pc cm^-3 yr^-1"
      elif name[0:2] == "DM":
         try:
            dm_num = np.int(name[2:])
            str_c = "pc cm^-3 yr^-" + np.str(dm_num)
         except: str_c = "(Units unknown)"
      else: str_c = "(Units unknown)"
   #
   output.append([str_a, str_b, str_c])

# Find minimum column widths necessary
wcol1 = wcol2 = wcol3 = 0
for i in range(0, len(output)):
   if len(output[i][0]) > wcol1: wcol1 = len(output[i][0])
   if len(output[i][1]) > wcol2: wcol2 = len(output[i][1])
   if len(output[i][2]) > wcol3: wcol3 = len(output[i][2])

# Print results!
print ""
print "File '" + parfile + "'."
if flags['-t']: print "Fit parameters reported with Tempo error in parfile."
else: print "Fit parameters reported with Tempo error doubled."
print ""
for row in output:
   print row[0].rjust(wcol1), row[1].rjust(wcol2), row[2].ljust(wcol3)
print ""

