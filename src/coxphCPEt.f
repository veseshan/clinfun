c     CPE for discrete risk score; finite number of categories and lots of ties
c     there is no smoothing involved; variance in Mo & Heller 2016 paper
c     we will start with data sorted (in decreasing order) by xbeta

      subroutine cpesubt(n, p, xmat, xbeta, nPI, cpe, vderiv, ursum,
     1     ussq)
      integer n, p
      double precision xmat(n,p), xbeta(n), nPI, cpe, vderiv(p),
     1     ursum(n), ussq

      integer i, j, k
      double precision xbji, exbji, eij, uij, xjik

c     number of pairs with beta'Xi > beta'Xj
      nPI = 0
c     calculate CPW with ties out
      cpe = 0
      do 20 i = 1, n-1
         do 10 j = i+1, n
            xbji = xbeta(j) - xbeta(i)
c     if pair is not tied
            if (xbji .ne. 0.0) then
c     increment pair count
               nPI = nPI + 1
c     eij is I(beta'Xij > 0)/{1+exp(beta'Xji)}
               eij = 1.0+exp(xbji)
c     pair's contribution to CPE
               cpe = cpe + 1.0/eij
            endif
 10      continue
 20   continue
c     estimated CPE
      cpe = cpe/nPI

c     variance calculation
      do 100 i = 1, n-1
         do 90 j = i+1, n
            xbji = xbeta(j) - xbeta(i)
c     if pair is not tied
            if (xbji .ne. 0.0) then
c     eij is I(beta'Xij > 0)/{1+exp(beta'Xji)}
               exbji = exp(xbji)
               eij = 1.0/(1.0+exbji)
c     uij is centered by estimated cpe
               uij = eij - cpe
c     eij can now be used for exp(beta'Xji)/{1+exp(beta'Xji)}^2
               eij = exbji*eij*eij
c     the variance has both uij and uji; need to use it in both order
               ursum(i) = ursum(i) + uij
               ursum(j) = ursum(j) + uij
               ussq = ussq + 2*uij**2
c     derivative term bit
               do 80 k = 1, p
                  xjik = xmat(j,k) - xmat(i,k)
                  vderiv(k) = vderiv(k) - xjik*eij
 80            continue
            endif
 90      continue
 100  continue

      return
      end
