#! /usr/bin/env python
import argparse, sys, struct, math, numpy

class tempo_cov:

    def __init__(self,filename='matrix.tmp'):
        self.read_matrix(filename)

    def reset(self,npar):
        self.npar = npar
        self.params = [''] * npar
        self.errs = [0.0] * npar
        self.cov = [[]] * npar

    def add_row(self,raw):
        stuff = struct.unpack("=ii5sdd" + "d"*self.npar, raw)
        idx = stuff[1] - 1
        self.params[idx] = stuff[2]
        self.errs[idx] = stuff[4]
        self.cov[idx] = list(stuff[5:])

    def get_cov(self,par1,par2):
        idx1 = self.params.index(par1)
        idx2 = self.params.index(par2)
        return self.cov[idx1][idx2]

    def read_matrix(self,filename='matrix.tmp'):
        """Read a tempo matrix.tmp file."""
        f = open(filename,'r')
        run = True
        first = True
        nbhdr = 4 + 4 + 5 + 8 + 8
        while run:
            try:
                l1 = struct.unpack("=i", f.read(4))[0]
                ncov = (l1 - nbhdr)/8
                if first:
                    self.reset(ncov)
                    first = False
                self.add_row(f.read(l1))
                l2 = struct.unpack("=i", f.read(4))[0]
                # TODO could put more error handling here
            except:
                run = False

class dmxparse:
    def __init__(self,parfile,xmx=False):
        self.val = {}
        self.err = {}
        self.epoch = {}
        self.r1 = {}
        self.r2 = {}
        self.f1 = {}
        self.f2 = {}
        if xmx:
          d = 'X'
        else:
          d = 'D'
        for l in open(parfile).readlines():
            if not l.startswith(d+'MX'): continue
            fields = l.split()
            key = fields[0]
            val = fields[1].replace('D','e')
            if len(fields)==4: 
                flag = int(fields[2])
                err = fields[3].replace('D','e')
            pfx=None
            idx=None
            (pfx,idx,newkey) = (None,None,None)
            if '_' in key: 
                (pfx, idx) = key.split('_')
                newkey = d+"X%03d" % int(idx)
            if l.startswith(d+'MX_') and flag==1:
                self.val[newkey] = float(val)
                self.err[newkey] = float(err)
            if l.startswith(d+'MXEP_'): self.epoch[newkey] = float(val)
            if l.startswith(d+'MXR1_'): self.r1[newkey] = float(val)
            if l.startswith(d+'MXR2_'): self.r2[newkey] = float(val)
            if l.startswith(d+'MXF1_'): self.f1[newkey] = float(val)
            if l.startswith(d+'MXF2_'): self.f2[newkey] = float(val)
        # the following are intended for xmx, in which epochs and frequency
        #   ranges aren't defined (as of 2018-Jun-18, when I am adding this)
        for key in self.r1.keys():
            if not key in self.epoch:
                self.epoch[key] = 0.5*(self.r1[key]+self.r2[key])
            if not key in self.f1:
                self.f1[key] = 0.
            if not key in self.f2:
                self.f1[key] = 0.
                self.f2[key] = 0.

    def fix_errs(self,tmp_cov):
        self.verr = {}
        n = len(self.err)
        cc = numpy.zeros((n,n))
        k = sorted(self.err.keys())
        for i in range(n):
            for j in range(n):
                cc[i,j] = tmp_cov.get_cov(k[i],k[j]) \
                        * self.err[k[i]] * self.err[k[j]]

        # Find error in mean DM
        self.mean = numpy.mean(self.val.values())
        self.mean_err = numpy.sqrt(cc.sum())/float(n)

        # Do the correction for varying DM
        m = numpy.identity(n) - numpy.ones((n,n))/float(n)
        cc = numpy.dot(numpy.dot(m,cc),m)
        for i in range(n):
            self.verr[k[i]] = math.sqrt(cc[i,i])

def get_base_dm(parfile):
    dm = 0.
    for l in open(parfile).readlines():
        if not l.startswith('DM '): continue
        dm = float(l.split()[1])
    return dm
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Parse DMX (or XMX) values from tempo .par file.\nSubtract mean DM value (default, may be over-ridden by flags.")
    parser.add_argument('parfile',    metavar='par_file', help="Parameter file (output from tempo)")
    parser.add_argument('matrixfile', metavar='cov_matrix_file',nargs="?",default="matrix.tmp",help="Covariance file (output from tempo)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-d','--dmfull',      action="store_true", help="Report full DM value at each epoch (no offset)")
    group.add_argument('-o','--originaldmx', action="store_true", help="Report original DMX values (no mean subtracttion)")
    group.add_argument('-x','--xmx',         action="store_true", help="Analyze XMX coefficient values instead of DMX values")

    args = parser.parse_args()

    dmfull = args.dmfull
    originaldmx = args.originaldmx
    parfile = args.parfile
    matrixfile = args.matrixfile
    xmx = args.xmx

    d = dmxparse(parfile,xmx)
    c = tempo_cov(matrixfile)

    if xmx:
      dx = 'XX'
      dmx = 'XMX'
      dm = 'XMX'
    else:
      dx = 'DX'
      dmx = 'DMX'
      dm = 'DM'
    # Check that DMX's match up
    ndmx1 = [dx in p for p in c.params].count(True)
    ndmx2 = len(d.err)
    if ndmx1 != ndmx2:
        print >>sys.stderr, "Error: "+dmx+" entries do not match up"
        print >>sys.stderr, "  parfile N%s = %d, matrix N%s = %d" % (dmx, ndmx1, dmx, ndmx2)
        sys.exit(1)

    d.fix_errs(c)
    print "# Mean %s value = %+.6e" % (dmx,d.mean)
    print "# Uncertainty in average %s = %.5e" % (dm,d.mean_err)
    print "# Columns: %sEP %s_value %s_var_err %sR1 %sR2 %sF1 %sF2 %s_bin" % (dmx,dmx,dmx,dmx,dmx,dmx,dmx,dmx)
    doffset = 0.
    if dmfull:
        if xmx:
            doffset = 0.
        else: 
            doffset = get_base_dm(parfile)
    elif not originaldmx:
        doffset = -d.mean
    for k in sorted(d.val.keys()):
        # print d.epoch[k], d.val[k], d.err[k], d.verr[k], k
        print "%.4f %+.7e %.3e %.4f %.4f %7.2f %7.2f %s" % (d.epoch[k], 
                d.val[k]+doffset, d.verr[k], 
                d.r1[k], d.r2[k], d.f1[k], d.f2[k], k)
