c
c Program to generate a grid for Klebanoff's flatplate 
c
      parameter (idim = 65, jdim = 65)
c
      implicit real*8(a-h,o-z)
c
      dimension x (idim,jdim),
     1          y (idim,jdim),
     2          s1(jdim)
c
c ytop  +----------------------------+
c       |                            |
c       |                            |
c       |                            |
c       |                            |
c       |                            |
c ybot  +----------------------------+
c       |                            |
c       xinit                        xfinal
c
c Flatplate data:
c
      xinit  = - 0.5d0
      xfinal =   1.5d0
c
      ybot   =   0.0d0
      ytop   =   0.18d0
c
c Set up equal spaced x distribution
c
      ds     = (xfinal - xinit) / float (idim - 1)
c
c Y distribution:
c
      strtch = 1.00025
      call stretch (jdim, strtch, s1)
c
c Calculate i = 1 distribution
c
      do j = 1, jdim
         x(1,j) = xinit
         y(1,j) = ybot + s1(j) * (ytop - ybot)
      end do
c
c Copy i = 1 distribution throughout domain
c
      do i = 2, idim
         xsta = xinit + float (i-1) * ds
         do j = 1, jdim
            x(i,j) = xsta
            y(i,j) = y(1,j)
         end do
      end do
c
c Output grid
c
      open (7,file='flatplate.grd',form='unformatted',status='new')
      call wrp3d (idim, jdim, x, y)
      close (7)
c
c Finished
c
      stop
      end




      subroutine stretch (jdim, strtch, s)
c
c Subroutine to set up a non-dimensionalized stretching function
c
      implicit real*8(a-h,o-z)
c
      dimension s(jdim)
c
      rh(b,psi) = 1.d0 - b + 2.d0 * b /
     1                    ( 1.d0 + ((b+1.d0)/(b-1.d0))**psi )
c
      s(1) = 0.d0
      do 100 j = 2, jdim
         ps   = float (jdim-j) / float(jdim-1)
         s(j) = rh(strtch, ps)
  100 continue
c
c Finished
c
      return
      end



      subroutine wrp3d (idim, jdim, x, y)
c
c Subroutine to write a binary plot3d grid file
c
      implicit real*8(a-h,o-z)
c
      dimension x(idim,jdim),
     1          y(idim,jdim)
c
      write (7) idim, jdim
      write (7) ((x(i,j),i=1,idim),j=1,jdim),
     1          ((y(i,j),i=1,idim),j=1,jdim)
c
      return
      end
