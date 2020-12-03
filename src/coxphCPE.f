      subroutine cpesub(n, p, xmat, xbeta, bw, cpe, scpe, vderiv,
     1     ursum, ussq)
      integer n, p
      double precision xmat(n,p), xbeta(n), bw, cpe, scpe, vderiv(p),
     1     ursum(n), ussq

      integer i, j, k
      double precision xbij, pij, pji, dij, dji, eij, eji, uij, uji,
     1     xijk

      double precision fpnorm, fdnorm
      external fpnorm, fdnorm

      do 100 i = 1, n-1
         do 90 j = i+1, n
            xbij = xbeta(i) - xbeta(j)
            pij = fpnorm(-xbij/bw)
            pji = 1.0 - pij
            dij = fdnorm(-xbij/bw)
            dji = dij
            eij = 1.0+exp(xbij)
            eji = 1.0+exp(-xbij)
            uij = pij/eij
            uji = pji/eji
            if (xbij .le. 0.0) then
               cpe = cpe + 1.0/eij
            else
               cpe = cpe + 1.0/eji
            endif
            scpe = scpe + uij + uji
            ussq = ussq + 2.0*(uij + uji)**2
            ursum(i) = ursum(i) + uij + uji
            ursum(j) = ursum(j) + uij + uji
            do 10 k = 1, p
               xijk = xmat(i,k) - xmat(j,k)
               if (xijk .ne. 0.0) then
                  vderiv(k) = vderiv(k) + (xijk/bw)*(dij/eij - dji/eji)
     1                 + xijk*(pij*(eij-1)/eij**2 - pji*(eji-1)/eji**2)
               endif
 10         continue
 90      continue
 100  continue

      return
      end
