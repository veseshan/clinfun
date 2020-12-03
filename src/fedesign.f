c     Computes Fisher's exact 2-sided rejection region for margin n1,n2
      subroutine ferej(nmx,n1,n2,alpha,fcl,lgamma)
c------------------------------------------------------------------------------
      implicit none
c------------------------------------------------------------------------------
      integer nmx,n1,n2,fcl(nmx,2)
c     fcl stores the lower (lcp) and upper (ucp) critical values
      double precision alpha,lgamma(nmx)
c     lgamma stores log(n!) for 0,1,...,n1+n2 sent from S-plus
      integer n,k,lcp,ucp,ii
      double precision hgpc,hgp,hgps,alpha2
c     hgp  - hypergeometric probabilities of lower and upper points
c     hgps - sum of hypergeometric probabilities of rejection region
c     hgpc - common factor of hypergeometric probabilities for n1,n2,k

      alpha2 = alpha/2.0d0
      n = n1 + n2

c     store conditional critical regions
      do 100 k = 0,n
         lcp = max(0,k-n2)
         ucp = min(n1,k)
         hgpc = lgamma(n1+1) + lgamma(n2+1) + lgamma(k+1)
     1        + lgamma(n-k+1) - lgamma(n+1)
         ii = lcp
         hgps = 0.0d0
 20      hgp = exp(hgpc - lgamma(ii+1) - lgamma(n1-ii+1)
     1        - lgamma(k-ii+1) - lgamma(n2-k+ii+1))
         hgps = hgps + hgp
         if (hgps.gt.alpha2) go to 25
         ii = ii+1
         go to 20
 25      fcl(k+1,1) = ii
         ii = ucp
         hgps = 0.0d0
 30      hgp = exp(hgpc - lgamma(ii+1) - lgamma(n1-ii+1)
     1        - lgamma(k-ii+1) - lgamma(n2-k+ii+1))
         hgps = hgps + hgp
         if (hgps.gt.alpha2) go to 35
         ii = ii-1
         go to 30
 35      fcl(k+1,2) = ii
 100  continue

      return
      end

c     Computes unconditional power for Fisher's exact rejection region
      subroutine fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
      integer nmx,n1,n2,fcl(nmx,2)
      double precision p1,p2,lgamma(nmx),upow
      integer n,i,j,k,lcp,ucp
      double precision bprob1,bprob2,lp1,l1p1,lp2,l1p2

      n = n1 + n2
      upow = 0.0d0
      lp1 = log(p1)
      l1p1 = log(1.0d0-p1)
      lp2 = log(p2)
      l1p2 = log(1.0d0-p2)
      do 20 k = 0,n
         lcp = max(0,k-n2)
         ucp = min(n1,k)
         do 10 i = lcp,ucp
            j = k - i
            if ((i.lt.fcl(k+1,1)).or.(i.gt.fcl(k+1,2))) then
               bprob1 = exp(lgamma(n1+1) - lgamma(i+1) - lgamma(n1-i+1)
     1              + dfloat(i)*lp1 + dfloat(n1-i)*l1p1)
               bprob2 = exp(lgamma(n2+1) - lgamma(j+1) - lgamma(n2-j+1)
     1              + dfloat(j)*lp2 + dfloat(n2-j)*l1p2)
               upow = upow + bprob1*bprob2
            endif
 10      continue
 20   continue

      return
      end

      subroutine femdor(nmx,n1,n2,p1,alpha,power,fcl,lgamma,omdor)
      integer nmx,n1,n2,fcl(nmx,2)
      double precision p1,alpha,power,lgamma(nmx),omdor(3,2)
c     omdor contains the Schlesselman and FE odd ratios and power
      integer n,i
      double precision d,smdor,p2,upow,p2l,p2r

      n = n1 + n2
      call ferej(nmx,n1,n2,alpha,fcl,lgamma)
c     Exact power for the Schlesselman mdor
      smdor = omdor(1,1)
      d = (smdor - 1.0d0)*p1*(1.0d0-p1)/(1.0d0 + (smdor - 1.0d0)*p1)
      p2 = p1 + d
      call fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
      omdor(1,2) = upow
c     Exact power for the CPS mdor
      smdor = omdor(2,1)
      d = (smdor - 1.0d0)*p1*(1.0d0-p1)/(1.0d0 + (smdor - 1.0d0)*p1)
      p2 = p1 + d
      call fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
      omdor(2,2) = upow


      p2l = p1
      p2r = 1.0d0
      do 20 i = 1,20
         p2 = (p2l + p2r)/2.0d0
         call fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
         if (upow.le.power) then
            p2l = p2
         else
            p2r = p2
         endif
 20   continue
      p2 = (p2l + p2r)/2.0d0
      call fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
      omdor(3,1) = (p2/(1.0d0-p2))/(p1/(1.0d0-p1))
      omdor(3,2) = upow

      return
      end

      subroutine fessiz(nmx,p1,p2,r,alpha,power,npm,fcl,lgamma,ossiz)
      integer nmx,npm,fcl(nmx,2)
      double precision p1,p2,r,alpha,power,lgamma(nmx),ossiz(2,3)
c     ossiz contains the Fleiss and FE sample size and power
      integer n1,n2,i,n1f,n2f,fen1,fen2
      double precision upow

      n1f = int(ossiz(1,1)) + 1
      n2f = int(ossiz(1,2)) + 1
      call ferej(nmx,n1f,n2f,alpha,fcl,lgamma)
      call fepow(nmx,n1f,n2f,p1,p2,fcl,lgamma,upow)
      ossiz(1,3) = upow

      fen1 = int(ossiz(1,1)-dfloat(npm)) + 1
      fen2 = int(ossiz(1,2)-r*dfloat(npm)) + 1
      do 20 i = -npm,npm
         n1 = int(ossiz(1,1)+dfloat(i)) + 1
         n2 = int(ossiz(1,2)+r*dfloat(i)) + 1
         call ferej(nmx,n1,n2,alpha,fcl,lgamma)
         call fepow(nmx,n1,n2,p1,p2,fcl,lgamma,upow)
         if (upow.lt.power) then
            fen1 = int(ossiz(1,1)+dfloat(i+1)) + 1
            fen2 = int(ossiz(1,2)+r*dfloat(i+1)) + 1
         endif
 20   continue
      call ferej(nmx,fen1,fen2,alpha,fcl,lgamma)
      call fepow(nmx,fen1,fen2,p1,p2,fcl,lgamma,upow)
      ossiz(2,1) = dfloat(fen1)
      ossiz(2,2) = dfloat(fen2)
      ossiz(2,3) = upow

      return
      end
