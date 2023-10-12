c     nv, nn, nd are number of markers, normal and diseased subjects, n = nn+nd
c     calculates the empirical ROC curve and 
c     Delong, Delong, Clarke-Pearson variance for area under ROC
      subroutine roccurve(n, nn, nd, marker, status, nu, tpr, fpr)
      integer n, nn, nd, nu, status(n)
      double precision marker(n), tpr(nu), fpr(nu)

      integer i, l
      double precision rn, rd, normgt, disgt

      rn = dble(nn)
      rd = dble(nd)

c     now the contribution of each unique observation
      l = nu
      tpr(l) = 1
      fpr(l) = 1
c     normgt is #normal > x  &  disgt is #diseased > x
      normgt = rn
      disgt = rd
      do 10 i = 1, n-1
         if (marker(i) .eq. marker(i+1)) then
            if (status(i) .eq. 0) then
               normgt = normgt - 1
            else
               disgt = disgt - 1
            endif
         else
            if (status(i) .eq. 0) then
               normgt = normgt - 1
            else
               disgt = disgt - 1
            endif
            l = l-1
            tpr(l) = disgt/rd
            fpr(l) = normgt/rn
         endif
 10   continue
      tpr(1) = 0
      fpr(1) = 0

      return
      end
