c     Probability generating function for the Jonckheere-Terpstra statistic
c     Based on the algorithm by Harding (JRSS-C, pp1-6, 1988).
c     The algorithm is reproduced in R mailing list on Oct 6, 2003
c
c     the ranksum for J-T statistic goes from 0 to M (= mxsum-1)
c
      subroutine djonck(mxsum, jrsum, ng, cgsize)
      integer mxsum, ng, cgsize(ng)
      double precision jrsum(mxsum)

      integer mm, m, n, g, t, u, s

      mm = mxsum-1
      do 100 g = 1, ng-1
         m = cgsize(g) - cgsize(g+1)
         n = cgsize(g+1)
         do 20 t = n+1, m+n
            do 10 u = mm, t, -1
               jrsum(u+1) = jrsum(u+1) - jrsum(u+1-t)
 10         continue
 20      continue

         do 40 s = 1, m
            do 30 u = s, mm
               jrsum(u+1) = jrsum(u+1) + jrsum(u+1-s)
 30         continue
 40      continue
 100  continue

      return
      end
