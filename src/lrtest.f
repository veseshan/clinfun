c     G-rho family of tests for comparing survival distributions
c     the expected value is #deaths * #risk(i)/ #risk
c     under permutation #risk(i) changes; #deaths and #risk don't
c     wghts is the survprob(t-)**rho doesn't change under permutation
c
c     n=#subjects; nt=#time-points; tfrq=frequency of times;
c     tdth=#deaths at each time-point (by strata);
c     ng=#groups; ns=#strata; sfrq=#unique times per stratum; 
c     grisk=#subjects at risk in each group
c
      subroutine lrtest(n, nt, ng, ns, tfrq, tdth, sfrq, grisk, wghts,
     1     sts, grp, odeath, edeath, lrvar)
      integer n, nt, ng, ns, tfrq(nt), sfrq(ns), grp(n)
      double precision tdth(nt), grisk(ng), wghts(nt), sts(n),
     1     odeath(ng), edeath(ng), lrvar(ng, ng)


      integer itim, istr, i, j, k, jr, jc
      double precision nrisk, efactr, vfactr

      k = n+1
      i = nt + 1
      do 60 istr = ns,1,-1
         nrisk = 0.0
         do 10 j = 1,ng
            grisk(j) = 0.0
 10      continue
         do 50 itim = sfrq(istr), 1, -1
            i = i-1
            do 20 j = 1,tfrq(i)
               k = k-1
               nrisk = nrisk + 1.0
               grisk(grp(k)) = grisk(grp(k)) + 1.0
               odeath(grp(k)) = odeath(grp(k)) + wghts(i)*sts(k)
 20         continue
            if (tdth(i) .gt. 0) then
               efactr = wghts(i)*tdth(i)/nrisk
               if (nrisk .gt. 1.0) then
                  vfactr = (wghts(i)**2)*tdth(i)*(nrisk-tdth(i))/
     1                 ((nrisk**2)*(nrisk - 1.0))
               else
                  vfactr = 0.0
               endif
               do 40 jr = 1,ng
                  edeath(jr) = edeath(jr) + efactr*grisk(jr)
                  lrvar(jr,jr) = lrvar(jr,jr) + vfactr*grisk(jr)*
     1                 (nrisk-grisk(jr))
                  do 30 jc = 1,jr-1
                     lrvar(jr,jc) = lrvar(jr,jc) - vfactr*grisk(jr)*
     1                    grisk(jc)
 30               continue
 40            continue
            endif
 50      continue
 60   continue

      do 80 jr = 1, ng-1
         do 70 jc = jr+1, ng
            lrvar(jr,jc) = lrvar(jc,jr)
 70      continue
 80   continue

      return
      end
