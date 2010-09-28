c
c Routine to read in a 2D single block binary PLOT3D grid and
c write it out formatted
c
      parameter (idim = 257)
      parameter (jdim =  97)
c
      dimension r(idim, jdim, 2)
c
c 1.  Read binary grid
c
      open (7,file='rae9_257.grd',form='unformatted',status='old')
      read (7) idimin, jdimin
      if (idimin .ne. idim .or. jdimin .ne. jdim) then
         write (6,'(a)')     ' Input dimensions do not match '
         write (6,'(a,2i5)') ' Dimension IDIM, JDIM = ', idim,   jdim
         write (6,'(a,2i5)') ' Gridfile  IDIM, JDIM = ', idimin, jdimin
         stop
      endif
c
      read (7) (((r(i,j,l),i=1,idim),j=1,jdim),l=1,2)
c
c 2.  Write formatted grid file
c
      open (8,file='rae_257x97.fmt',form='formatted',status='new')
      write (8,1000) idim, jdim
 1000 format (2i10)
      write (8,1005) (((r(i,j,l),i=1,idim),j=1,jdim),l=1,2)
 1005 format (4e23.15)
c
c 3.  Finished
c
      stop
      end
