c     subroutine to calculate the matrices needed for nested AUC inference
c     x is the data matrix p covariates (full model has p+1 covariates)
c     xb (xb0) is the X*beta from the full (reduced) model
c     data are ordered such that first n0 rows are y=0 and rest y=1
c     dmat, wmat (matrices) and valt (scalar) objects that are returned
      subroutine daucmats(n, p, n0, x, xb, xb0, delta, dmat, wmat, valt)
      integer n, p, n0
      double precision x(n, p), xb(n), xb0(n), delta, dmat(p, p),
     1     wmat(p, p), valt

      double precision, allocatable :: uidot(:,:), udotj(:,:)
      double precision, allocatable :: uij(:), xij(:)
      double precision, allocatable :: w1(:,:), w2(:,:), wsq(:,:)
      double precision, allocatable :: eidot(:), edotj(:)

      integer n1, i, i1, j, k, l
      double precision xbij, dxbij, xbdxbij, dcontr, nw1, nw2, n0n1
      double precision xb0ij, eij, esq, sig1, sig2, bwratio

      double precision fpnorm, fdnorm
      external fpnorm, fdnorm

      n1 = n - n0
      n0n1 = dfloat(n0)*dfloat(n1)

      allocate(uidot(n1, p), udotj(n0,p))
      allocate(uij(p), xij(p), w1(p,p), w2(p,p), wsq(p,p))

c     initialize variables
      do 10 i = 1, n1
         do 9 k = 1, p
            uidot(i, k) = 0.0
 9       continue
 10   continue
      do 20 j = 1, n0
         do 19 k = 1, p
            udotj(j, k) = 0.0
 19      continue
 20   continue

      do 30 k = 1, p
         do 29 l = 1, p
            wsq(k,l) = 0.0
 29      continue
 30   continue

c     start computing the matrices
      do 40 i = n0+1, n
         i1 = i-n0
         do 39 j = 1, n0
c     X*beta of ith subject - X*beta of jth subject (scaled by bandwidth)
            xbij = xb(i) - xb(j)
c     standard normal density (phi) of xbij
            dxbij = fdnorm(xbij)
c     product xbij*phi(xbij) needed for second derivative matrix
            xbdxbij = xbij*dxbij
            do 38 k = 1, p
c     X_ik - X_jk 
               xij(k) = x(i,k) - x(j,k)
c     first derivative Uij is phi(xbij)*xij
               uij(k) = dxbij*xij(k)
c     because of symmetry only the lower triangle is computed
               do 37 l = 1, k
c     second derivative is xbij*phi(xbij)*xij*t(xij)
                  dcontr = xbdxbij * xij(k) * xij(l)
                  dmat(k, l) = dmat(k, l) - dcontr
                  wsq(k, l) = wsq(k, l) - uij(k)*uij(l)
 37            continue
c     need sum_i{Uij} and sum_j{Uij} for W matrix
               uidot(i1,k) = uidot(i1,k) + uij(k)
               udotj(j,k) = udotj(j,k) + uij(k)
 38         continue
c            write(25,*) i, j, uij
 39      continue
 40   continue

c     calculate w1 and w2
c     W1 = sum_i {(sum_j Uij) * t(sum_j Uij)} - sum_ij {Uij * t(Uij)}
c     W2 = sum_j {(sum_i Uij) * t(sum_i Uij)} - sum_ij {Uij * t(Uij)}

c     sum_ij {Uij * t(Uij)} part of both W1 and W2
      do 50 k = 1, p
         do 49 l = 1, p
            w1(k,l) = wsq(k,l)
            w2(k,l) = wsq(k,l)
 49      continue
 50   continue

c     sum_i {(sum_j Uij) * t(sum_j Uij)} part of W1
      do 60 i = 1, n1
         do 59 k = 1, p
            do 58 l = 1, k
               w1(k, l) = w1(k, l) + uidot(i, k)*uidot(i, l)
 58         continue
 59      continue
 60   continue

c     sum_j {(sum_i Uij) * t(sum_i Uij)} part of W2
      do 70 j = 1, n0
         do 69 k = 1, p
            do 68 l = 1, k
               w2(k, l) = w2(k, l) + udotj(j, k)*udotj(j, l)
 68         continue
 69      continue
 70   continue

c     symmetrize dmat, w1 and w2
      if (p .gt. 1) then
         do 80 k = 2, p
            do 79 l = 1, k-1
               dmat(l, k) = dmat(k, l)
               w1(l, k) = w1(k, l)
               w2(l, k) = w2(k, l)
 79         continue
 80      continue
      endif

c     scale dmat by n0*n1
c     scale w1 by n1*n0*(n0-1) and w2 by n1*(n1-1)*n0 and calculate W
      nw1 = dfloat(n)/(dfloat(n1)**2 * dfloat(n0) * dfloat(n0-1))
      nw2 = dfloat(n)/(dfloat(n0)**2 * dfloat(n1) * dfloat(n1-1))
      do 90 k = 1, p
         do 89 l = 1, p
            wmat(k,l) = nw1*w1(k,l) + nw2*w2(k,l)
            dmat(k,l) = dmat(k,l)/n0n1
 89      continue
 90   continue

      deallocate(uidot, udotj, uij, xij, w1, w2)

c     calculating the variance under the alternative
      allocate(eidot(n1), edotj(n0))
      do 100 i = 1, n1
         eidot(i) = 0.0
 100  continue
      do 101 j = 1, n0
         edotj(j) = 0.0
 101  continue

c     xbeta is scaled by bandwidth for density need scaling for cdf
c     bwratio is the bandwidth ratio n^1/3 for cdf and n^1/5 for pdf
      bwratio = dfloat(n)**(-2.0d0/15.0d0)
      do 102 i = 1,n
         xb(i) = xb(i)*bwratio
 102  continue
c     calculation of variance for the CI under non-null deltaauc
      do 120 i = n0+1, n
         i1 = i-n0
         do 110 j = 1, n0
            xbij = xb(i) - xb(j)
            xb0ij = xb0(i) - xb0(j)
c     eij = Phi(xbij) - Phi(xb0ij); but variance uses eij - delta
            eij = fpnorm(xbij) - fpnorm(xb0ij) - delta
c     sum_i eij
            eidot(i1) = eidot(i1) + eij
c     sum_j eij
            edotj(j) = edotj(j) + eij
c     sum_ij eij**2
            esq = esq - eij**2
 110     continue
 120  continue

c     sigma1^2 = sum_i {(sum_j eij)**2} - sum_ij {eij**2}
      sig1 = esq
      do 130 i = 1, n1
         sig1 = sig1 + eidot(i)**2
 130  continue

c     sigma1^2 = sum_j {(sum_i eij)**2} - sum_ij {eij**2}
      sig2 = esq
      do 140 j = 1, n0
         sig2 = sig2 + edotj(j)**2
 140   continue

c     variance is scaled appropriately
      valt = nw1*sig1 + nw2*sig2
      deallocate(eidot, edotj)

      return
      end

c     calculates the nonparametric area under the ROC curves
c     n, nn: #total #normal subjects; marker sorted by disease status
      subroutine rocauc(n, nn, marker, area)
      integer n, nn
      double precision marker(n), area

      integer i, j, j0, k, nties0, nties1, nties, nd
      double precision rnd, disgt

c     storage for x, loc used to store ith marker and order
      integer, allocatable :: loc(:)
      double precision, allocatable :: x(:)
      allocate(x(n), loc(n))

      nd = n - nn
      rnd = dfloat(nn)*dfloat(nd)

      do 10 i = 1, n
         x(i) = marker(i)
         loc(i) = i
 10   continue

c     sort the markers
      call qsort4(x, loc, 1, n)

      area = 0
c     now the contribution of each observation
      nties = 0
      nties1 = 0
c     disgt is #diseased with marker value > x 
      disgt = dfloat(nd)
      do 20 i = 1, n-1
         if (x(i) .eq. x(i+1)) then
c     if next marker value is same as current go to the next i
c     and count the total and diseased subjects at tied marker value 
            nties = nties + 1
            if (loc(i) .gt. nn) nties1 = nties1 + 1
         else
c     if not update area and reset ties counters to zero
            nties = nties + 1
            if (loc(i) .gt. nn) nties1 = nties1 + 1
            nties0 = nties - nties1
            disgt = disgt - dfloat(nties1)
c     increment count by #dis > x and half #dis == x for each normal at x
            area = area + dfloat(nties0)*(disgt + 0.5*dfloat(nties1))
            nties = 0
            nties1 = 0
         endif
 20   continue
c     update number of obsns with value of marker(n) and update area
      nties = nties + 1
      if (loc(n) .gt. nn) nties1 = nties1 + 1
      nties0 = nties - nties1
      disgt = disgt - dfloat(nties1)
c     increment count by #dis > x and half #dis == x for each normal at x
      area = area + dfloat(nties0)*(disgt + 0.5*dfloat(nties1))
c     calculate the area
      area = area/rnd

      deallocate(x, loc)
      return
      end

c     calculates the smoothed area under the ROC curves
c     n, nn: #total #normal subjects; marker sorted by disease status
c     marker has been scaled by the bandwidth
      subroutine smrocauc(n, nn, marker, area)
      integer n, nn
      double precision marker(n), area

      integer i, j
      double precision rnd

      double precision fpnorm
      external fpnorm

      nd = n - nn
      rnd = dfloat(nn)*dfloat(nd)

      area = 0.0d0
      do 20 i = 1, nn
         do 10 j = nn+1, n
            area = area + fpnorm(marker(j) - marker(i))
 10      continue
 20   continue
c     calculate the area
      area = area/rnd

      return
      end
