c     This subroutine is adapted from the code for the JRSS-B paper 
c     it does the conditional logrank in a 3 group setting
c     need to change it to a general k-group setting
c     create k indicator variables and the k vector of U-statistics
c
c     n = number of observations,  ng = number of groups
c     p = number of covariates (maximum of 3) -- moved to R wrapper?
c     osts, ogrp, ocov are the status, group and covariate variables
c     a0 and a1 (per grp) vectors correspond to the same in the paper 
c     vectors Vii, Vij, Vji, Vjj are the group specific U-statistics
c     Vidot, Vdotj, Vijij are the row, column and product sums
c     igrp, jgrp are the group indicator for observation i & j
c     lrmn, lrvar are the test statistic vector and variance
c     bb is the bandwidth for the kernel

      subroutine uclrst(n, ng, p, osts, ogrp, ocov, a0, a1, xi, xj,
     1     Vii, Vij, Vji, Vjj, Vidot, Vdotj, Vijij, igrp, jgrp, lrmn,
     2     lrvar, bb)
      integer n, ng, p, ogrp(n)
      double precision osts(n), ocov(n,p), a0(n), a1(ng,n), xi(p),
     1     xj(p), Vii(ng), Vij(ng), Vji(ng), Vjj(ng), Vidot(ng,n),
     2     Vdotj(ng,n), Vijij(ng,ng), igrp(ng), jgrp(ng), lrmn(ng),
     3     lrvar(ng,ng), bb(p)


      integer i,j,ic,g,gr,gc
      double precision rn,kernel,kwtii,kwtij
      external kernel

      rn = dble(n)

c     calculate the a0 and a1 vector for all time points
c     the vij calculations depend only on these
      do 25 i = 1,n
         a0(i) = 0.0
         do 21 g = 1,ng
            a1(g,i) = 0.0
 21      continue
         do 23 j = i,n
            do 22 ic = 1,p
               xi(ic) = ocov(i,ic)
               xj(ic) = ocov(j,ic)
 22         continue
            kwtij = kernel(p,xj,xi,bb)
c     a1 changes only for the group for jth observation
            a0(i) = a0(i) + kwtij
            a1(ogrp(j), i) = a1(ogrp(j), i) + kwtij
 23      continue
         a0(i) = a0(i)/rn
         do 24 g = 1,ng
            a1(g,i) = a1(g,i)/rn
 24      continue
 25   continue

c     now calculate Vij to compute the statistic and its variance
      do 40 i = 1,n
         do 31 ic = 1,p
            xi(ic) = ocov(i,ic)
 31      continue
c     set the indicator of the ith observation
         igrp(ogrp(i)) = 1.0
         do 39 j = i,n
            do 32 ic = 1,p
               xj(ic) = ocov(j,ic)
 32         continue
            kwtii = kernel(p,xi,xi,bb)
            kwtij = kernel(p,xj,xi,bb)
c     set the indicator of the jth observation
            jgrp(ogrp(j)) = 1.0
            do 33 g = 1,ng
c     compute the Vij's
               Vii(g) = osts(i)*(igrp(g) - a1(g,i)/a0(i))*
     1              (1.0-kwtii/a0(i))/rn
               if (j .gt. i) then
                  Vij(g) = osts(i)*(igrp(g) - a1(g,i)/a0(i) -
     1                 kwtij*(jgrp(g) - a1(g,i)/a0(i))/a0(i))/rn
                  Vji(g) = osts(j)*(jgrp(g) - a1(g,j)/a0(j))/rn
                  Vjj(g) = osts(j)*(jgrp(g) - a1(g,j)/a0(j))*
     1                 (1.0-kwtii/a0(j))/rn
               else
                  Vij(g) = Vii(g)
                  Vji(g) = Vii(g)
                  Vjj(g) = Vii(g)
               endif
 33         continue
c     unset the indicator of the jth observation
            jgrp(ogrp(j)) = 0.0
c     compute the row and column sums of Vij's
            do 34 g = 1,ng
               lrmn(g) = lrmn(g) + Vij(g)
               Vidot(g,i) = Vidot(g,i) + Vij(g)
               Vdotj(g,j) = Vdotj(g,j) + Vij(g)
               if(j .gt. i) then
                  lrmn(g) = lrmn(g) + Vji(g)
                  Vidot(g,j) = Vidot(g,j) + Vji(g)
                  Vdotj(g,i) = Vdotj(g,i) + Vji(g)
               endif
 34         continue
c     compute the Vij product contributions to the variance
            if (i.lt.j) then
c     these are the i neq j terms
               do 36 gr = 1,ng
                  do 35 gc = gr,ng
                     if (gr .eq. gc) then
                        Vijij(gr,gr) = Vijij(gr,gr) + 
     1                       (Vij(gr) + Vji(gr))*(Vij(gr) +
     2                       Vji(gr) + 2*Vii(gr) + 2*Vjj(gr))
                     else
                        Vijij(gr,gc) = Vijij(gr,gc) + 
     1                       (Vij(gr) + Vji(gr))*(Vij(gc) + Vji(gc)) + 
     2                       (Vij(gr) + Vji(gr))*(Vii(gc) + Vjj(gc)) + 
     3                       (Vij(gc) + Vji(gc))*(Vii(gr) + Vjj(gr))
     4                       
                     endif
 35               continue
 36            continue
            else
c     these are the ii terms
               do 38 gr = 1,ng
                  do 37 gc = gr,ng
                     if (gr .eq. gc) then
                        Vijij(gr,gr) = Vijij(gr,gr) + 3*Vii(gr)**2
                     else
                        Vijij(gr,gc) = Vijij(gr,gc) + 
     1                                         3*Vii(gr)*Vii(gc)
                     endif
 37               continue
 38            continue
            endif
 39      continue
c     unset the indicator of the ith observation
         igrp(ogrp(i)) = 0.0
 40   continue

c     finish the variance computation
      do 50 gr = 1,ng
         do 49 gc = gr,ng
            do 48 i = 1,n
               lrvar(gr,gc) = lrvar(gr,gc) + (Vidot(gr,i)+Vdotj(gr,i))*
     1              (Vidot(gc,i)+Vdotj(gc,i))
 48         continue
            lrvar(gr,gc) = lrvar(gr,gc) - Vijij(gr,gc)
            if (gc .gt. gr) lrvar(gc,gr) = lrvar(gr,gc)
 49      continue
 50   continue

      do 60 i = 1,n
         if (osts(i) .eq. 1) then
            igrp(ogrp(i)) = igrp(ogrp(i)) + 1.0
            do 55 g = 1,ng
               jgrp(g) = jgrp(g) + a1(g,i)/a0(i)
 55         continue
         endif
 60   continue

      return
      end

      double precision function kernel(p,xx,xi,bb)
      integer i,p
      double precision xx(p),xi(p),bb(p),uu
      kernel = 1.0
      do 10 i = 1,p
         uu = ((xx(i)-xi(i))/bb(i))**2
         kernel=kernel*(1.0/bb(i))*exp(-uu/2)
 10   continue
      return
      end
