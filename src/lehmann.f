c     Mann-Whitney, Kruskal-Wallis, Jonckheere-Terpstra test
c     (monte carlo) power under Lehmann alternative

c     ng     = number of groups
c     gsize  = group sizes
c     mcgsiz = temporary group size vector for monte carlo step
c     oratio = odds ratio
c     gsor   = inner product of group size and odds ratio
c     rsum   = rank-sum vector
c     nrep   = number of replications (for monte carlo)
c     kw     = indicator of Kruskal-Wallis or Jonckheere-Terpstra
c     tstat  = test statistic vector

      subroutine lehman(ng, gsize, mcgsiz, oratio, gsor, rsum, kw, nrep,
     1     tstat)
      integer ng, gsize(ng), nrep
      logical kw
      double precision mcgsiz(ng), oratio(ng), gsor, rsum(ng), 
     1     tstat(nrep)

      double precision jtstat
      external jtstat

      integer nn, nr, i

c     calculate the total sample size
      nn = 0
      do 10 i = 1, ng
         nn = nn + gsize(i)
 10   continue

      call rndstart()
      if (kw) then
         do 50 nr = 1, nrep
c     initialize rank-sum and mcgsiz (= gsize*oratio) vector
            do 20 i = 1, ng
               rsum(i) = 0.0
               mcgsiz(i) = dfloat(gsize(i))*oratio(i)
 20         continue
c     generate a new rank-sum vector
            call kwrsum(nn, ng, mcgsiz, oratio, rsum, gsor)
c     compute Kruskal-Wallis test statistic
            tstat(nr) = 0.0
            do 30 i = 1, ng
               tstat(nr) = tstat(nr) + rsum(i)**2/dfloat(gsize(i))
 30         continue
 50      continue
      else
         do 100 nr = 1, nrep
c     initialize rank-sum and mcgsiz (= gsize*oratio) vector
            do 60 i = 1, ng
               rsum(i) = dfloat(gsize(i))
               mcgsiz(i) = dfloat(gsize(i))*oratio(i)
 60         continue
c     compute Jonckheere-Terpstra statistic for a new rank-sum vector
            tstat(nr) = jtstat(nn, ng, mcgsiz, oratio, rsum, gsor)
 100     continue
      endif
      call rndend()

      return
      end

      subroutine kwrsum(nn, ng, mcgsiz, oratio, rsum, gsor)
      integer nn, ng
      double precision mcgsiz(ng), oratio(ng), rsum(ng), gsor

      double precision dunif
      external dunif

      double precision cc, thrsh, mcgsor
      integer i, j

c     ith obsn is in group j with prob nj*or(j)/sum(nj*or(j))
c     sum(nj*or(j)) at the beginning 
      mcgsor = gsor
      do 20 i = 1, nn
c     scale the uniform random number to avoid division
         cc = mcgsor*dunif()
         j = 0
         thrsh = 0
c     check if the random number is within range for group j
         do 10 while(cc .gt. thrsh)
c     go the next group
            j = j + 1
c     increment threshold by nj*or(j)
            thrsh = thrsh + mcgsiz(j)
 10      enddo
c     decrease nj*or(j) for the jth group
         mcgsiz(j) = mcgsiz(j) - oratio(j)
c     decrease sum(nj*or(j))
         mcgsor = mcgsor - oratio(j)
c     add to rank sum for jth group
         rsum(j) = rsum(j) + dfloat(i)
 20   continue

      return
      end

      double precision function jtstat(nn, ng, mcgsiz, oratio, rsum, 
     1     gsor)
      integer nn, ng
      double precision mcgsiz(ng), oratio(ng), rsum(ng), gsor

      double precision dunif
      external dunif

      double precision cc, thrsh, mcgsor
      integer i, j

c     initialize J-T statistic as n+(n-1)+...+1
      jtstat = dfloat(nn*(nn+1))/2
c     ith obsn is in group j with prob nj*or(j)/sum(nj*or(j))
c     sum(nj*or(j)) at the beginning 
      mcgsor = gsor
      do 20 i = 1, nn
c     scale the uniform random number to avoid division
         cc = mcgsor*dunif()
         j = 0
         thrsh = 0
c     check if the random number is within range for group j
         do 10 while(cc .gt. thrsh)
c     go the next group
            j = j + 1
c     decrease J-T statistic by number of subjects in group j
            jtstat = jtstat - rsum(j)
c     increment threshold by nj*or(j)
            thrsh = thrsh + mcgsiz(j)
 10      enddo
c         call dblepr("jtstat",6,jtstat,1)
c         call intpr(" j",2,j,1)
c     decrease nj*or(j) for the jth group
         mcgsiz(j) = mcgsiz(j) - oratio(j)
c     decrease sum(nj*or(j))
         mcgsor = mcgsor - oratio(j)
c     reduce jth group size
         rsum(j) = rsum(j) - 1.0
 20   continue

      return
      end
