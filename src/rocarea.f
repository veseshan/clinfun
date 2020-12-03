c     calculates the Delong, Delong, Clarke-Pearson variance for area under ROC
c     nv, nn, nd are number of markers, normal and diseased subjects, n = nn+nd
      subroutine rocarea(n, nv, nn, nd, markers, area, jkarea)
      integer n, nv, nn, nd
      double precision markers(n,nv), area(nv), jkarea(n,nv)

      integer i, j, j0, k, nties0, nties1, nties
      double precision rn, rd, rnd, normlt, disgt

c     storage for x, loc used to store ith marker and order
      integer, allocatable :: loc(:)
      double precision, allocatable :: x(:)
      allocate(x(n), loc(n))

      rn = dfloat(nn-1)*dfloat(nd)
      rd = dfloat(nn)*dfloat(nd-1)
      rnd = dfloat(nn)*dfloat(nd)

      do 100 k = 1,nv
         do 10 i = 1, n
            x(i) = markers(i,k)
            loc(i) = i
 10      continue

c     sort the markers
         call qsort4(x, loc, 1, n)

         area(k) = 0
c     now the contribution of each observation
         nties = 0
         nties0 = 0
c     normlt is #normal < x  &  disgt is #diseased > x 
         normlt = 0
         disgt = dfloat(nd)
         do 30 i = 1, n-1
            if (x(i) .eq. x(i+1)) then
               nties = nties + 1
               if (loc(i) .le. nn) nties0 = nties0 + 1
            else
               nties = nties + 1
               if (loc(i) .le. nn) nties0 = nties0 + 1
               nties1 = nties - nties0
               disgt = disgt - dfloat(nties1)
               do 20 j = (i-nties+1),i
                  j0 = loc(j)
                  if (j0 .le. nn) then
                     jkarea(j0,k) = disgt + 0.5*dfloat(nties1)
                     area(k) = area(k) + jkarea(j0,k)
                  else
                     jkarea(j0,k) = normlt + 0.5*dfloat(nties0)
                  endif
 20            continue
               normlt = normlt + dfloat(nties0)
               nties = 0
               nties0 = 0
            endif
 30      continue
         nties = nties + 1
         if (loc(n) .le. nn) nties0 = nties0 + 1
         nties1 = nties - nties0
         disgt = disgt - dfloat(nties1)
         do 40 j = (n-nties+1),n
            j0 = loc(j)
            if (j0 .le. nn) then
               jkarea(j0,k) = disgt + 0.5*dfloat(nties1)
               area(k) = area(k) + jkarea(j0,k)
            else
               jkarea(j0,k) = normlt + 0.5*dfloat(nties0)
            endif
 40      continue

c     calculate the area if observation i (normal) is left out
         do 50 i = 1, nn
            jkarea(i,k) = (area(k) - jkarea(i,k))/rn
 50      continue
c     calculate the area if observation j (diseased) is left out
         do 60 j = nn+1, n
            jkarea(j,k) = (area(k) - jkarea(j,k))/rd
 60      continue
c     calculate the area
         area(k) = area(k)/rnd
 100  continue

      deallocate(x, loc)
      return
      end
