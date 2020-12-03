c
c     compile with Splus SHLIB -o filename.so ph2simon.f ....
c
      subroutine f2bdry(m,nmax,ep1,ep2,p0,p1,cp0,cp1,bdry,peten,
     1     nmax1,bprob0,bprob1)
c------------------------------------------------------------------------------
      implicit none
c------------------------------------------------------------------------------
      integer m,nmax,bdry(nmax,4),nmax1
      double precision ep1,ep2,p0(m),p1(m),cp0(m),cp1(m),peten(nmax,2),
     1     bprob0(nmax1),bprob1(nmax1)

      double precision pet,ess,dn1,dn2,essn
      integer i,n1,n2,n,r1,r,ind1,ind2,ind21,rr

      do 100 n = 2,nmax
         essn = dfloat(n)
         do 90 n1 = 1,n-1
            n2 = n-n1
            dn1 = dfloat(n1)
            dn2 = dfloat(n2)
            ind1 = n1*(n1+3)/2
            ind2 = n2*(n2+3)/2
            do 10 i=1,n+1
               bprob0(i) = 0.0d0
               bprob1(i) = 0.0d0
 10         continue
            pet = 1.0d0
            do 50 r1 = n1,0,-1
               ind21 = ind2
               pet = pet - p0(ind1)
               do 40 r=n2+r1,r1,-1
                  rr = r + 1
                  bprob0(rr) = bprob0(rr) + p0(ind1)*cp0(ind21)
                  bprob1(rr) = bprob1(rr) + p1(ind1)*cp1(ind21)
                  ind21 = ind21 - 1
                  if ((bprob0(rr).lt.ep1).and.(1-bprob1(rr).lt.ep2))
     1                 then
                     ess = dn1 + (1.0d0-pet)*dn2
                     if (ess.lt.essn) then
                        essn = ess
                        peten(n,1) = ess
                        peten(n,2) = pet
                        bdry(n,1)= r1-1
                        bdry(n,2)= n1
                        bdry(n,3)= r-1
                        bdry(n,4)= n
                     endif
                  endif
 40            continue
               rr = r1 + 1
               do 41 r = 1,r1,1
                  bprob0(r) = bprob0(r1+1)
                  bprob1(r) = bprob1(r1+1)
 41            continue
               ind1 = ind1 -1
 50         continue
 90      continue
 100  continue

      return
      end
