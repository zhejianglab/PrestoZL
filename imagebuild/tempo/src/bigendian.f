c     Endian check now done by configure
      logical function bigendian()
      include 'config.h'
      bigendian = ISBIGENDIAN
      return
      end
