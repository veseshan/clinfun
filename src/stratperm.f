c     code for stratified permutation index vector
c     where ii is the orginal index vector is 1,...,n
c     ns1 is the number of strata plus 1
c     ns2 is 0 followed by the cumulatice totals of strata size
c     uu is a vector of n uniform(0,1) random numbers

      subroutine strperm1(n, ii, ns1, ns2, uu)
      integer n, ii(n), ns1, ns2(ns1)
      double precision uu(n)

      integer i, j, k, l, tmp

      l = 0
      do 100 i = 1, ns1-1
         j = ns2(i+1) - ns2(i)
         do 50 while(j .gt. 1) 
            l = l+1
            k = int(uu(l)*dble(j))
            tmp = ii(l)
            ii(l) = ii(l+k)
            ii(l+k) = tmp
            j = j-1
 50      continue
         l = l+1
 100  continue

      return
      end
