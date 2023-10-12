c     written by Venkatraman E Seshan 4/20/2011
c     x, y are ordered using order(x,y)
      subroutine ktau(n, x, y, tau)
      integer n
      double precision x(n), y(n), tau

      integer i, n0, nties
      double precision dn, dnx, dny, ties

c     storage for idx used to store indices for blocks of x
      integer, allocatable :: idx(:)
      allocate(idx(n))

      dn = dble(n)*dble(n-1)/2

      dnx = 0.0d0
      nties = 1
      n0 = 0
c     find the ties in x, calculate correction for ties
c     n0 is the number of unique x values
      do 10 i = 1, n-1
         if (x(i) .eq. x(i+1)) then
            nties = nties + 1
         else
            n0 = n0+1
c     here idx[i] = sum(x == ux[i]) where ux is sort(unique(x))
            idx(n0) = nties
            ties = dble(nties)
            dnx = dnx + ties*(ties-1.0d0)/2.0d0
            nties = 1
         endif
 10   continue
      n0 = n0+1
      idx(n0) = nties
      if (x(n-1) .eq. x(n)) then 
         ties = dble(nties)
         dnx = dnx + ties*(ties-1.0d0)/2.0d0
      endif
c     now idx[i] = sum(x <= ux[i]) where ux is sort(unique(x))
      do 15 i = 2,n0
         idx(i) = idx(i-1) + idx(i)
 15   continue

c     call the subroutine to compute the #concordant - #discordant
      call countall(n, y, n0, idx, tau)

c     the call above returns y in sorted order, now adjust for ties in y
      dny = 0
      ties = 1.0d0
      do 20 i = 1, n-1
         if (y(i) .eq. y(i+1)) then
            ties = ties + 1.0d0
         else
            dny = dny + ties*(ties-1.0d0)/2.0d0
            ties = 1.0d0
         endif
 20   continue
      if (y(n-1) .eq. y(n)) then 
         dny = dny + ties*(ties-1.0d0)/2.0d0
      endif

c     calculate Kendall's tau-b
      tau = tau/sqrt((dn-dnx)*(dn-dny))

      deallocate(idx)

      return
      end

c     loop y through the unique values of x, coalescing pair at a time
      subroutine countall(n, y, n0, idx, tau)
      integer n, n0, idx(n0)
      double precision y(n), tau

      integer i, n1, m, m1, m0
      double precision btau

      tau = 0
c     this loop merges consecutive pairs of blocks of x 
c     starting with individual (unique) x values
c     eventually leading to a single block in log2(n0) steps
      do 20 while (n0 .gt. 1)
         n1 = n0/2
         m = idx(2)
         m1 = idx(1)
         m0 = 1
         call blockcount(m, y(m0), m1, btau)
         tau = tau + btau
         idx(1) = idx(2)
         do 10 i = 2,n1
            m0 = idx(2*(i-1)) + 1
            m = idx(2*i) - m0 + 1
            m1 = idx(2*i-1) - m0 + 1
            call blockcount(m, y(m0), m1, btau)
            tau = tau + btau
            idx(i) = idx(2*i)
 10      continue
         if (2*n1 .lt. n0) then
            idx(n1+1) = idx(n0)
            n0 = n1 + 1
         else
            n0 = n1
         endif
 20   continue

      return
      end

c     blocks A = 1,...,m1 and B=m1+1,...,m; x[i] < x[j] for i in A & j in B
c     y's sorted within blocks y[1] <= ... <= y[m1] and y[m1+1] <= ... <= y[m]
c     this subroutine computes {I(y[i] < y[j]) - I(y[i] > y[j])}
c     summed over i = 1,...,m1 and j = m1+1,...,m. also returns sorted y
      subroutine blockcount(m, y, m1, btau)
      integer m, m1
      double precision y(m), btau

      integer i, j, l, mp1
      double precision numlt, numeq, numgt, cury

c     temporary space for y so that sorted y can be returned
      double precision, allocatable :: localy(:)
c     need the m+1th value for the while loop
      mp1 = m+1
      allocate(localy(mp1))


c     copy y into localy
      do 10 i = 1,m
         localy(i) = y(i)
 10   continue
c     need localy[m+1] larger than y[m1] (make it larger than all y's)
      localy(mp1) = max(y(m),y(m1))+1

c     initialize number lt, eq, and gt as well as btau
c     at the start all y[m1+1] to y[m] are assumed larger
      numlt = 0
      numgt = dble(m-m1)
      numeq = 0

c     initialize block count
      btau = 0
c     start with current y smaller than y[1]
      cury = localy(1) - 1
      j = m1+1
      l = 0
c     loops through y[1] ... y[m1]
      do 50 i = 1, m1
         if (cury .lt. localy(i)) then
            cury = localy(i)
            numlt = numlt + numeq
            numeq = 0
c     find y[m1+1 ... m] less than current y[i]
            do 20 while ((localy(j) .lt. cury) .and. (j .le. m))
c     add the y's to the sorted list
               l = l+1
               y(l) = localy(j)
c     adjuest the number greater and less accordingly
               numlt = numlt + 1
               numgt = numgt - 1
               j = j+1
 20         continue
c     now find y[m1+1 ... m] equal to current y[i]
            do 30 while ((localy(j) .eq. cury) .and. (j .le. m))
c     add the y's to the sorted list
               l = l+1
               y(l) = localy(j)
               numeq = numeq + 1
               j = j+1
 30         continue
c     adjust number greater than current y[i]
            numgt = numgt - numeq
            l = l+1
            y(l) = localy(i)
            btau = btau + numgt - numlt
         else
c     contribution to tau when y[i] = y[i-1] 
            l = l+1
            y(l) = localy(i)
            btau = btau + numgt - numlt
         endif
 50   continue
 
      deallocate(localy)

      return
      end
