c
c program to read a grid in Charley Swanson's format and 
c output it for ISAAC (plot3d binary)
c
      parameter (idim = 225, jdim = 33)
c
c Double precision
c
      implicit real*8 (a-h,o-z)
c
      dimension x(idim,jdim), y(idim,jdim)
c
c read in Charley's grid
c
      ingrd = 8
      open (ingrd,file='gr224n20',status='old')
      read (ingrd,*) 
      read (ingrd,*) 
      read (ingrd,*) 
      read (ingrd,*) xi, xj, cord
      read (ingrd,*) ni, nj, nte1, nte2, nbot
      if (ni .ne. idim .or. nj .ne. jdim) then
         write (*,'(a,24i5)')' Input grid size differs from specified ',
     1     idim, jdim, ni, nj
         stop
      endif
      read (ingrd,1000) ((x(i,j),i=1,ni),j=1,nj)
      read (ingrd,1000) ((y(i,j),i=1,ni),j=1,nj)
c
c output grid
c
      igdout = 9
      inc = 1
      niout = idim
      njout = jdim
      open (igdout,file='n12_225_33.grd',status='new',
     1      form='unformatted')
      write (igdout) niout, njout
      write (igdout) ((x(i,j),i=1,idim,inc),j=1,jdim,inc),
     1               ((y(i,j),i=1,idim,inc),j=1,jdim,inc)
      close (igdout)
c
      stop
 1000 format (10e13.7)
      end
