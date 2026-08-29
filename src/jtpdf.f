c     proof of independence in Mark van de Wiel dissertation
c     http://alexandria.tue.nl/extra2/200012964.pdf

c     instead of calling dwilcox here values are passed on
c     ***** problem rose in R 4.6.0 because of need for wilcox_free *****
c     *****  this version gets the density in R and passes it here  *****

c     replace old pdf0 and pdf1 vectors with mxsum by ng matrix
c     last column (ng) is for pdf1 convolution calculation
c     columns 1 <= g < ng are Mann-Whitney pdf (old pdf0) for g vs {g+1,...,ng}

      subroutine jtpdf(mxsum, pdf, ng, cgsize, pdfs)
      integer mxsum, ng, cgsize(ng)
      double precision pdf(mxsum),pdfs(mxsum,ng)

      integer i, j, g, mn0, mn1

c     this is the group ng-1 against group ng comparison
c     if only two groups it reduces to Wilcoxon-Mann-Whitney
      g = ng - 1
      m = cgsize(g) - cgsize(g+1)
      n = cgsize(g+1)
      mn1 = m*n
      do 10 i = 0, mn1
         pdf(i+1) = pdfs(i+1,g)
 10   continue

      do 60 g = ng-2,1,-1
c     current pdf1 is in pdf; set it and reset pdf to 0
         do 20 i = 1, mn1+1
            pdfs(i,ng) = pdf(i)
            pdf(i) = 0
 20      continue
c     obtain pdf for group g against {g+1,...,ng}
         m = cgsize(g) - cgsize(g+1)
         n = cgsize(g+1)
         mn0 = m*n
         do 50 i = 0, mn0
            do 40 j = 0, mn1
               pdf(i+j+1) = pdf(i+j+1) + pdfs(i+1,g)*pdfs(j+1,ng)
 40         continue
 50      continue
c     increment the range for next iteration of pdf1 (if needed)
         mn1 = mn0 + mn1
 60   continue

      return
      end
