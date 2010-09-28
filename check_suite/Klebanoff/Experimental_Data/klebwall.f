c
c Routine to convert the u/u_edge vs. y data from Klebanoff to 
c wall coordinates (u+ vs. y+).
c Some values from the report:
c     Re_L  = 4.2 x 10^6
c     L     = 14.2 ft.
c     U_inf = 50 ft./sec. (u_edge)
c     delta = 3 in.
c     y+    = 27.6 at y/delta = 0.01 (p.1148)
c
      Re    = 4.2e06
      uedge = 50.e0
      rnu   = uedge * 14.2e0 / Re
      delta = 3.e0 / 12.e0
c
c     y+ = y * u_tau / nu
c
      utau  = 27.6 * rnu / (0.01 * delta)
c
c Read blank lines
c
      do iblank = 1, 8
         read (5,*)
      end do
c
c Read u/u_edge vs. y(in.)
c
      write (6,1000) 
 1000 format (' ','Klebanoffs flat plate')
      write (6,1010) 
 1010 format (' ','       y+        ','       u+        ',
     1            '     y/delta     ')
      do i = 1, 33
         read (5,*) j, y, u
         y      = y / 12.e0
         u      = u * uedge
         yplus  = y * utau / rnu
         uplus  = u / utau
         ydelta = y / delta
         write (6,'(3e17.7)') yplus, uplus, ydelta
      end do
c
c Finished
c
      stop
      end
