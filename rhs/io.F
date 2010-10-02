C
C
C Portions Copyright (C) 2001 Joseph H. Morrison
C
C This file is part of ISAAC.
C
C This program is distributed under the terms of the ISAAC Public Source
C License. This program is distributed WITHOUT ANY WARRANTY; without
C even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C PURPOSE. 
C
C You should have received a copy of the ISAAC Public Source License
C with this program. If you did not, you may get a copy of the license
C at http://isaac-cfd.sourceforge.net.
C
C     $Revision: 4.4 $
C     $Author: jhmorr $
C     $Date: 2001/10/29 03:25:31 $
C     $State: Exp $
C     $Log: io.F,v $
C     Revision 4.4  2001/10/29 03:25:31  jhmorr
C     Updated license information
C
C     Revision 4.3  2001/06/22 03:47:30  jhmorr
C     Add output of function name file. Modify contact information
C
C     Revision 4.2  2001/06/08 04:56:16  jhmorr
C     Added notice
C
C     Revision 4.1.1.1  1998/05/30 15:53:53  jhmorr
C     Added input for homogeneous shear boundary condition.
C
c Revision 4.1  1998/04/14  20:34:58  jhmorr
c Alpha 4.1.
c
C



      SUBROUTINE CHAPND (STRING1, STRING2)
C
C Routine to append character STRING2 to STRING1.
C
C N.B. This routine assumes that STRING1 is right padded with spaces. 
C      Unpredictable results will be returned if STRING1 is not right padded.
C
C STRING1 : Input base character string (with at least one terminating space)
C STRING2 : Input character string to be appended to STRING1
C
      CHARACTER*(*) STRING1,
     1              STRING2
C
      LEN1 = ICHLNG (STRING1)
C
      STRING1 = STRING1(1:LEN1) // STRING2
C
C Finished with CHAPND
C
      RETURN
      END



      INTEGER FUNCTION ICHLNG (STRING)
C
C Function to return the length of a string defined as the location where
C the first white space is located.
C
C STRING : Input string
C
      CHARACTER*80 STRING
C
      LOC = 1
  100 CONTINUE
      IF (STRING(LOC:LOC) .EQ. ' ') THEN
         GO TO 200
      ENDIF
      LOC = LOC + 1
      IF (LOC .GT. 80) THEN
         GO TO 200
      ENDIF
      GO TO 100
C
C Found space or end of string
C
  200 CONTINUE
      ICHLNG = LOC - 1
C
      RETURN
      END





      SUBROUTINE MASFLO (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                 Q, S, VOL, FMASS, DQP, DQM, QP, QM, F, FM, STEMP,
     2                 IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX, VISCOS)
C
C     Routine to calculate the mass flow through an I = const
C     plane to be used for convergence tests.
C     This routine is based on the FLUXI routine.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                      (used in dimensioning temporary space used for 
C                       calculating q(+,-) and f.
C     Q              : Primitive variables at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     FMASS          : Mass flow through the Ith plane
C     DQP,DQM        : Storage for delta q = q(i) - q(i-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     F,FM           : Fluxes at cell interfaces -
C                      if Flux Difference Splitting then F = Flux
C                      if Flux Vector Splitting then F = Fplus, FM = Fminus
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for I direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C     RKAP           : Kappa to use in MUSCL interpolation
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     2          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION FMASS (IDIM),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (MXSECT,NQ), 
     4          QM    (MXSECT,NQ),
     5          F     (IPSTRT-1:IPEND,2:JDIM,NF), 
     6          FM    (IPSTRT-1:IPEND,2:JDIM,NF), 
     7          STEMP (IPSTRT-1:IPEND,2:JDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C     Loop over K calculating Q(+,-), F(Q(+,-)) and updating the Residual.
C
      DO 300 K = 2, KDIM
C
C     Inviscid Residuals
C
         CALL MUSCLI (IDIM, JDIM, KDIM, K, IPSTRT, IPEND,
     1                Q, DQP, DQM, QP, QM,
     2                IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C     Calculate the fluxes from q(plus) and q(minus).  Pass only
C     the appropriate metrics from the current direction.
C
         DO 120 I = IPSTRT - 1, IPEND
            DO 110 J = 2, JDIM
               DO 100 L = 1, 3
                  STEMP(I,J,L) = S(J,K,I,L,1)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         NPTS = (IPEND - IPSTRT + 2) * (JDIM - 1)
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
	CALL FDS(IFLUX,ITURB,NPTS,QP,QM,STEMP,DQP,DQM,F)
C
            DO 160 I = IPSTRT - 1, IPEND
               DO 140 J = 2, JDIM
                  FMASS(I) = FMASS(I) + F(I,J,1) * S(J,K,I,4,1) 
  140          CONTINUE
  160       CONTINUE
         ENDIF
  300 CONTINUE
C
C     Output the results
C
      DO 500 I = IPSTRT - 1, IPEND
         WRITE (IOUT,1000) I, FMASS(I)
  500 CONTINUE
 1000 FORMAT (' ',' Mass flow through cross section ',I4,' is ',
     1        1PE14.7)
C
C     Finished calculating mass flow through I cross sections
C
      RETURN
      END






      SUBROUTINE MONITR (IDIM, JDIM, KDIM, Q, NPRMON, IPRMON)
C
C     Routine to monitor output during the iteration cycle
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     Q              : Primitive variables at cell centers
C     NPRMON         : Number of monitor print statements for current block
C     IPRMON         : Monitor print control vector for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IPRMON(10,NPRMON)
C
      DIMENSION ISTRT(3), IEND(3), INC(3)
C
      CHARACTER*10 TITLE(3)
      DATA TITLE / 'J, K =    ',
     1             'K, I =    ',
     2             'I, J =    ' /
C
C     Loop over all print statements
C
      do 130 i = 0, idim+2
         do 120 k = 0, kdim+2
            do 110 j = 0, jdim+2
               if (q(j,k,i,4) .ne. 0.e0) then
                  write (iout,'(a,3i5,e15.7)') ' i,j,k,w ',
     1                                           i,j,k,q(j,k,i,4)
               endif
  110       continue
  120    continue
  130 continue
C
C     Output a few blank lines to make output look better
C
      WRITE (IOUT,1015)
C
      RETURN
 1000 FORMAT (//' ',20X,A10,2I5/
     1        ' ','   I   J   K',
     2            6X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/PINF',8X,'T/TINF')
 1010 FORMAT (' ',3I4,7E13.5)
 1015 FORMAT (//' ')
      END






      SUBROUTINE OUTPUT (IDIM, JDIM, KDIM, Q, RMU, RMUT, YPLUS, ITURB, 
     1                   RC, IBLOCK, NPRINT, IPRINT)
C
C     Routine to write output
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     YPLUS          : Y+ at cell centers
C     ITURB          : Turbulence model
C     RC             : Calculated position of grid at cell centers
C     IBLOCK         : Block number to identify output
C     NPRINT         : Number of print statements for current block
C     IPRINT         : Print control vector for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          YPLUS (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IPRINT(10,NPRINT)
C
      DIMENSION ISTRT(3), IEND(3), INC(3)
C
      CHARACTER*10 TITLE(3)
      DATA TITLE / 'J, K =    ',
     1             'K, I =    ',
     2             'I, J =    ' /
C
C     Loop over all print statements
C
      DO 800 IPR = 1, NPRINT
         IDIR     = IPRINT(1,IPR)
         ISTRT(1) = IPRINT(2,IPR)
         ISTRT(2) = IPRINT(3,IPR)
         ISTRT(3) = IPRINT(4,IPR)
         IEND(1)  = IPRINT(5,IPR)
         IEND(2)  = IPRINT(6,IPR)
         IEND(3)  = IPRINT(7,IPR)
         INC(1)   = IPRINT(8,IPR)
         INC(2)   = IPRINT(9,IPR)
         INC(3)   = IPRINT(10,IPR)
         DO 530 I3 = ISTRT(ICY(IDIR,3)), IEND(ICY(IDIR,3)), 
     1               INC(ICY(IDIR,3))
            DO 520 I2 = ISTRT(ICY(IDIR,2)), IEND(ICY(IDIR,2)), 
     1                  INC(ICY(IDIR,2))
C
C     Output Header
C
               IF (ITURB .EQ. ITLMNR) THEN
                  WRITE (IOUT,1003) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITBLMX) THEN
                  WRITE (IOUT,1005) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITKE) THEN
                  WRITE (IOUT,1010) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITKW) THEN
                  WRITE (IOUT,1015) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITRS) THEN
                  WRITE (IOUT,1025) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITLES) THEN
                  WRITE (IOUT,1030) IBLOCK, TITLE(IDIR), I2, I3
               ELSE
                  WRITE (IOUT,1000) IBLOCK, TITLE(IDIR), I2, I3
               ENDIF
C
C     Inner Loop to Print out
C
               DO 510 I1 = ISTRT(ICY(IDIR,1)), IEND(ICY(IDIR,1)), 
     1                     INC(ICY(IDIR,1))
                  II = I1*IKD(ICY(IDIR,1),1) + I2*IKD(ICY(IDIR,2),1) +
     1                 I3*IKD(ICY(IDIR,3),1)
                  JJ = I1*IKD(ICY(IDIR,1),2) + I2*IKD(ICY(IDIR,2),2) +
     1                 I3*IKD(ICY(IDIR,3),2)
                  KK = I1*IKD(ICY(IDIR,1),3) + I2*IKD(ICY(IDIR,2),3) +
     1                 I3*IKD(ICY(IDIR,3),3)
                  PPINF = GAMMA * Q(JJ,KK,II,5)
                  TTINF = GAMMA * Q(JJ,KK,II,5) / Q(JJ,KK,II,1)
                  IF (ITURB .EQ. ITLMNR) THEN
                     WRITE (IOUT,1103) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II)
                  ELSE IF (ITURB .EQ. ITBLMX) THEN
                     WRITE (IOUT,1105) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II)
                  ELSE IF (ITURB. EQ. ITKE) THEN
                     RHO   = Q(JJ,KK,II,1)
                     P     = Q(JJ,KK,II,5)
                     TKE   = Q(JJ,KK,II,6)
                     EPSLN = Q(JJ,KK,II,7)
                     ASQ   = GAMMA * P / RHO
                     RMT   = SQRT (2.E0 * TKE / ASQ)
                     RET   = RHO*TKE*TKE/(EPSLN*RMU(JJ,KK,II)+RSMALL) *
     1                       RE / FSMACH
                     WRITE (IOUT,1110) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  Q(JJ,KK,II,7), 
     5                  YPLUS(JJ,KK,II), RMT, RET
                  ELSE IF (ITURB .EQ. ITKW) THEN
                     RHO   = Q(JJ,KK,II,1)
                     P     = Q(JJ,KK,II,5)
                     TKE   = Q(JJ,KK,II,6)
                     EPSLN = TKE * Q(JJ,KK,II,7) 
                     ASQ   = GAMMA * P / RHO
                     RMT   = SQRT (2.E0 * TKE / ASQ)
                     RET   = RHO*TKE*TKE/(EPSLN*RMU(JJ,KK,II)+RSMALL) *
     1                       RE / FSMACH
                     WRITE (IOUT,1110) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  Q(JJ,KK,II,7), 
     5                  YPLUS(JJ,KK,II), RMT, RET
                  ELSE IF (ITURB. EQ. ITRS) THEN
                     WRITE (IOUT,1125) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  
     5                  Q(JJ,KK,II,7),  Q(JJ,KK,II,8),  Q(JJ,KK,II,9),  
     6                  Q(JJ,KK,II,10), Q(JJ,KK,II,11), Q(JJ,KK,II,12), 
     7                  YPLUS(JJ,KK,II)
                  ELSE IF (ITURB .EQ. ITLES) THEN
                     WRITE (IOUT,1130) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II)
                  ELSE
                     WRITE (IOUT,1100) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF
                  ENDIF
  510          CONTINUE
  520       CONTINUE
  530    CONTINUE
  800 CONTINUE
C
      RETURN
C
 1000 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',9X,'X',16X,'Y',16X,'Z',
     2            15X,'RHO',15X,'U',16X,'V',
     3            16X,'W',13X,'P/Pinf',12X,'T/Tinf')
 1003 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU')
 1005 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT')
 1010 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'k',10X,'EPSILON',8X,'Yplus',
     5            7X,'Mach_t',7X,'Re_t')
 1015 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'k',11X,'OMEGA',8X,'Yplus',
     5            7X,'Mach_t',7X,'Re_t')
 1025 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'TAUXX',8X,'TAUYY',8X,'TAUZZ',8X,
     5            'TAUXY',8X,'TAUXZ',8X,'TAUYZ',7X,'EPSILON',8X,'Yplus')
 1030 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT')
 1100 FORMAT (' ',3I4,10(1PE17.9))
 1103 FORMAT (' ',3I4, 3(1PE13.5), 8(1PE13.5))
 1105 FORMAT (' ',3I4, 3(1PE13.5), 9(1PE13.5))
 1110 FORMAT (' ',3I4, 3(1PE13.5),14(1PE13.5))
 1125 FORMAT (' ',3I4, 3(1PE13.5),17(1PE13.5))
 1130 FORMAT (' ',3I4, 3(1PE13.5), 9(1PE13.5))
C
      END







 
      SUBROUTINE CUT3D (IDIM, JDIM, KDIM, Q, ICUTS)
C
C     This routine averages the value of the interior and boundary 
C     condition flow variable at the cell center for the CUT boundary
C     condition type to the face to coincide with the PLOT3D output.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION ICUTS (21)
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3), 
     2          IEND  (3),
     3          ISGN  (3),
     4          LEN   (3) 
C
C     Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C     Store boundary condition data in mnemonic names
C
      IDIR     = ICUTS( 2)
      INOUT    = ICUTS( 3)
      ISTRT(1) = ICUTS( 4)
      IEND (1) = ICUTS( 5)
      ISTRT(2) = ICUTS( 6)
      IEND (2) = ICUTS( 7)
      ISTRT(3) = ICUTS( 8)
      IEND (3) = ICUTS( 9)
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INOUT = 1) or the outer boundary (INOUT > 1)
C          INR    = interior point of block
C          IBCPT  = boundary point of block
C
      IF (INOUT .EQ. 1) THEN
         INR    = 2
         IBCPT  = 1
      ELSE
         INR    = IJKDIM(IDIR)
         IBCPT  = IJKDIM(IDIR) + 1
      ENDIF
C
C     Find the length of the two outer loops and loop over these using
C     offsets and delta function to set the two cut data points to the
C     two interior data points of block 2.
C
      DO 100 L = 1, 3
C
         LEN(L) = ABS ( IEND(L) - ISTRT(L) )
C
C     Increment/Decrement 
C
         IF ( IEND(L) .GT. ISTRT(L) ) THEN
            ISGN(L) =   1
         ELSE
            ISGN(L) = - 1
         ENDIF
  100 CONTINUE
C
C     II1 first  boundary point of block 1
C     II2 second boundary point of block 1
C     IN1 first  interior point of block 2
C     IN2 second interior point of block 2
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
      DO 800 I3 = 0, LEN(IC3)
         DO 700 I2 = 0, LEN(IC2)
            IBC = IBCPT                       * IKD (IC1, 1) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 1) 
            JBC = IBCPT                       * IKD (IC1, 2) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 2) 
            KBC = IBCPT                       * IKD (IC1, 3) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 3) 
C
            IN  = INR                         * IKD (IC1, 1) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 1) 
            JN  = INR                         * IKD (IC1, 2) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 2) 
            KN  = INR                         * IKD (IC1, 3) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 3) 
C
CDIR$ SHORTLOOP
            DO 600 L = 1, NQ
               Q(JBC,KBC,IBC,L) = 0.5E0 * (Q(JN,KN,IN,L) +
     1                                     Q(JBC,KBC,IBC,L))
 600        CONTINUE
 700     CONTINUE
 800  CONTINUE
C
C     Finished with Cut boundary condition fix to PLOT3D
C
      RETURN
      END






      SUBROUTINE PLOT3D (IDIM, JDIM, KDIM, Q, PROPS, RC, BINARY)
C
C Routine to write PLOT3D output. The mean flow variables are output in
C the PLOT3D solution (Q) file in 'conserved' form (neglecting TKE in the
C equation of state) and the turbulence variables and properties (if NQ > 5) 
C are output in the function file. The turbulence quantities are output as 
C they are calculated and are not converted to 'conserved' form.
C
C This routine is organized as follows:
C     1. Convert mean flow from primitive to 'conserved' variables 
C        NOTE: do not include TKE in conversion from pressure to total energy
C              to match the equation of state that is assumed in FAST/PLOT3D
C     2. Write PLOT3D output files (FORMATTED or UNFORMATTED)
C        a. Three-dimensional output
C           i.  Binary    { }
C           ii. Formatted { }
C        b. Two  -dimensional output
C           i.  Binary    { }
C           ii. Formatted { }
C               { a. Grid file     (G)
C                 b. Solution file (Q)
C                 c. If (NQ > 5) write function file (F) }
C     3. Convert mean flow from 'conserved' to primitive variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C RC             : Calculated position of grid at cell centers
C BINARY         : Logical for PLOT3D output file type
C                  = TRUE  -> binary (unformatted) file
C                    FALSE -> formatted file
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      LOGICAL BINARY
C 
C 1. Convert mean flow from primitive to 'conserved' variables 
C    NOTE: do not include TKE in conversion from pressure to total energy
C          to match the equation of state that is assumed in FAST/PLOT3D
C 
      DO 130 I = 1, IDIM+1 
         DO 120 K = 1, KDIM+1 
            DO 110 J = 1, JDIM+1 
               Q(J,K,I,5) = Q(J,K,I,5) / GAMM1 + 0.5E0*Q(J,K,I,1)*
     1                           ( Q(J,K,I,2)*Q(J,K,I,2)
     2                           + Q(J,K,I,3)*Q(J,K,I,3)
     3                           + Q(J,K,I,4)*Q(J,K,I,4) )
               Q(J,K,I,2) = Q(J,K,I,2) * Q(J,K,I,1)
               Q(J,K,I,3) = Q(J,K,I,3) * Q(J,K,I,1)
               Q(J,K,I,4) = Q(J,K,I,4) * Q(J,K,I,1)
  110       CONTINUE 
  120    CONTINUE 
  130 CONTINUE 
C 
C 2.   Write PLOT3D output files (FORMATTED or UNFORMATTED)
C
C 2-a. Output for three-dimensional case
C
      IF ( THREED ) THEN
         IF ( BINARY ) THEN
C
C 2-a-i.  Binary output
C
            WRITE (IPLT3G) ((((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                     K=1,KDIM+1),L=1,3)
            WRITE (IPLT3Q) ((((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                     K=1,KDIM+1),L=1,5)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F) 
     1                     ((((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                        K=1,KDIM+1),L=6,NQ),
     2                     ((((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     3                                        K=1,KDIM+1),L=1,NP)
            ENDIF
         ELSE
C
C 2-a-ii. Formatted output
C
            WRITE (IPLT3G,1000) ((((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                          K=1,KDIM+1),L=1,3)
            WRITE (IPLT3Q,1000) ((((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                          K=1,KDIM+1),L=1,5)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F,1000) 
     1                       ((((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                          K=1,KDIM+1),L=6,NQ),
     3                       ((((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                          K=1,KDIM+1),L=1,NP)
            ENDIF
         ENDIF
C
C 2-b. Output for two-dimensional case
C      Output as a three-dimensional case with KDIM=1 for FAST
C
      ELSE
         K = 2
         IF ( BINARY ) THEN
C
C 2-b-i.  Binary output
C
            WRITE (IPLT3G) (((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                    L=1,2)
            WRITE (IPLT3Q) (((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                    L=1,3),
     2                      ((Q (J,K,I,5),I=1,IDIM+1),J=1,JDIM+1)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F) 
     1                     (((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                       L=6,NQ),
     3                     (((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                       L=1,NP)
            ENDIF
         ELSE
C
C 2-b-ii. Formatted output
C
            WRITE (IPLT3G,1000) (((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                         L=1,2)
            WRITE (IPLT3Q,1000) (((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                         L=1,3),
     2                           ((Q (J,K,I,5),I=1,IDIM+1),J=1,JDIM+1)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F,1000) 
     1                       (((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                         L=6,NQ),
     3                       (((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                         L=1,NP)
            ENDIF
         ENDIF
      ENDIF
C 
C 3.  Convert mean flow from 'conserved' to primitive variables
C 
      DO 330 I = 1, IDIM+1 
         DO 320 K = 1, KDIM+1 
            DO 310 J = 1, JDIM+1 
               Q(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,1)
               Q(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,1)
               Q(J,K,I,4) = Q(J,K,I,4) / Q(J,K,I,1)
               Q(J,K,I,5) = GAMM1*( Q(J,K,I,5) - 0.5E0*Q(J,K,I,1)*
     1                               ( Q(J,K,I,2)*Q(J,K,I,2)
     2                               + Q(J,K,I,3)*Q(J,K,I,3)
     3                               + Q(J,K,I,4)*Q(J,K,I,4) ) )
  310       CONTINUE 
  320    CONTINUE 
  330 CONTINUE 
C
C Finished with PLOT3D output
C
      RETURN
 1000 FORMAT (5E16.8)
      END





      SUBROUTINE PLOT3DFN (ITURB)
C
C This routine writes the PLOT3D function name file.
C
C ITURB : Turbulence model
C
C This routine is organized as follows:
C     1. Write function names for turbulence quantities
C        a. K-epsilon
C        b. K-omega
C        c. Reynolds stress
C     2. Write function names for properties
C
      include '../header/common.h'
C
C Local storage
C
      INTEGER      L
C
      CHARACTER*25 FUNNAME (5)
C
C 1.   Write function names for turbulence quantities
C 1-a. K-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         WRITE (IPLT3FN,1000) 'K                        ',
     1                        'Epsilon                  '
C
C 1-b. K-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         WRITE (IPLT3FN,1000) 'K                        ',
     1                        'Omega                    '
C
C 1-c. Reynolds stress
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         WRITE (IPLT3FN,1000) 'Tau_xx                   ',
     1                        'Tau_yy                   ',
     2                        'Tau_zz                   ',
     3                        'Tau_xy                   ',
     4                        'Tau_xz                   ',
     5                        'Tau_yz                   ',
     6                        'Epsilon                  '
C
      ENDIF
C
C 2.   Write function names for properties
C
      FUNNAME(1) = 'mu'
      FUNNAME(2) = 'mu_t'
      FUNNAME(3) = 'y+'
      FUNNAME(4) = 'y_normal'
      FUNNAME(5) = 'C_mu*'
C
      DO L = 1, NP
         WRITE (IPLT3FN,1000) FUNNAME(L)
      END DO
C
C Finished PLOT3DFN writing function name file
C
      RETURN
 1000 FORMAT (A25)
      END






      SUBROUTINE RCENTR (IDIM, JDIM, KDIM, R, RC, IBCDAT, NBCS)
C
C Routine to average the positions of the vertices to get a cell centered 
C value of the grid for output purposes.
C
C This routine is organized as follows:
C      1. Calculate cell centers for (1:dim+1) in all directions
C         a. Corners
C         b. Edges
C         c. Faces
C         d. Interior
C      2. For wall functions, loop over all boundary conditions and
C         correct the cell center calculated above to account for
C         wall function spacing
C
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector of grid vertices (x,y,z)
C RC             : Calculated position of grid at cell centers
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend 
C                         10     = profile identifier
C NBCS           : Number of boundary conditions for current block
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3),
     1          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
C R  is the grid input at vertices
C RC is the grid at cell centers to be sent to PLOT3D
C
C 1. Calculate cell-centered grid
C
C 1-a. Corners
C
      DO 10 IDIR = 1, 3
         RC(1,1,1,IDIR)                = R(1,1,1,IDIR)
         RC(JDIM+1,1,1,IDIR)           = R(JDIM,1,1,IDIR)
         RC(1,1,IDIM+1,IDIR)           = R(1,1,IDIM,IDIR)
         RC(JDIM+1,1,IDIM+1,IDIR)      = R(JDIM,1,IDIM,IDIR)
         RC(1,KDIM+1,1,IDIR)           = R(1,KDIM,1,IDIR)
         RC(JDIM+1,KDIM+1,1,IDIR)      = R(JDIM,KDIM,1,IDIR)
         RC(1,KDIM+1,IDIM+1,IDIR)      = R(1,KDIM,IDIM,IDIR)
         RC(JDIM+1,KDIM+1,IDIM+1,IDIR) = R(JDIM,KDIM,IDIM,IDIR)
   10 CONTINUE
C
C 1-b. Edges
C
      DO 60 IDIR = 1, 3
         DO 30 I = 2, IDIM
            RC(1,1,I,IDIR)           = 0.5E0 * ( R(1,1,I-1,IDIR)
     1                                         + R(1,1,I,IDIR) )
            RC(JDIM+1,1,I,IDIR)      = 0.5E0 * ( R(JDIM,1,I-1,IDIR)
     1                                         + R(JDIM,1,I,IDIR) )
            RC(1,KDIM+1,I,IDIR)      = 0.5E0 * ( R(1,KDIM,I-1,IDIR)
     1                                         + R(1,KDIM,I,IDIR) )
            RC(JDIM+1,KDIM+1,I,IDIR) = 0.5E0 * ( R(JDIM,KDIM,I-1,IDIR)
     1                                         + R(JDIM,KDIM,I,IDIR) )
   30    CONTINUE
         DO 40 J = 2, JDIM
            RC(J,1,1,IDIR)           = 0.5E0 * ( R(J-1,1,1,IDIR)
     1                                         + R(J,1,1,IDIR) )
            RC(J,1,IDIM+1,IDIR)      = 0.5E0 * ( R(J-1,1,IDIM,IDIR)
     1                                         + R(J,1,IDIM,IDIR) )
            RC(J,KDIM+1,1,IDIR)      = 0.5E0 * ( R(J-1,KDIM,1,IDIR)
     1                                         + R(J,KDIM,1,IDIR) )
            RC(J,KDIM+1,IDIM+1,IDIR) = 0.5E0 * ( R(J-1,KDIM,IDIM,IDIR)
     1                                         + R(J,KDIM,IDIM,IDIR) )
   40    CONTINUE
         DO 50 K = 2, KDIM
            RC(1,K,1,IDIR)           = 0.5E0 * ( R(1,K-1,1,IDIR)
     1                                         + R(1,K,1,IDIR) )
            RC(1,K,IDIM+1,IDIR)      = 0.5E0 * ( R(1,K-1,IDIM,IDIR)
     1                                         + R(1,K,IDIM,IDIR) )
            RC(JDIM+1,K,1,IDIR)      = 0.5E0 * ( R(JDIM,K-1,1,IDIR)
     1                                         + R(JDIM,K,1,IDIR) )
            RC(JDIM+1,K,IDIM+1,IDIR) = 0.5E0 * ( R(JDIM,K-1,IDIM,IDIR)
     1                                         + R(JDIM,K,IDIM,IDIR) )
   50    CONTINUE
   60 CONTINUE
C
C 1-c. J = 1, J = JDIM+1 Planes
C
      DO 120 IDIR = 1, 3
         DO 110 I = 2, IDIM
            DO 100 K = 2, KDIM
               RC(1,K,I,IDIR)      = 0.25E0*
     1                  ( R(1,K,I,IDIR)   + R(1,K,I-1,IDIR)
     2                  + R(1,K-1,I,IDIR) + R(1,K-1,I-1,IDIR) )
               RC(JDIM+1,K,I,IDIR) = 0.25E0*
     1                  ( R(JDIM,K,I,IDIR)   + R(JDIM,K,I-1,IDIR)
     2                  + R(JDIM,K-1,I,IDIR) + R(JDIM,K-1,I-1,IDIR) )
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
C
C      K = 1, K = KDIM+1 Planes
C
      DO 150 IDIR = 1, 3
         DO 140 I = 2, IDIM
            DO 130 J = 2, JDIM
               RC(J,1,I,IDIR)      = 0.25E0*
     1                  ( R(J,1,I,IDIR)   + R(J,1,I-1,IDIR)
     2                  + R(J-1,1,I,IDIR) + R(J-1,1,I-1,IDIR) )
               RC(J,KDIM+1,I,IDIR) = 0.25E0*
     1                  ( R(J,KDIM,I,IDIR)   + R(J,KDIM,I-1,IDIR)
     2                  + R(J-1,KDIM,I,IDIR) + R(J-1,KDIM,I-1,IDIR) )
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
C
C      I = 1, I = IDIM+1 Planes
C
      DO 180 IDIR = 1, 3
         DO 170 K = 2, KDIM
            DO 160 J = 2, JDIM
               RC(J,K,1,IDIR)      = 0.25E0*
     1                  ( R(J,K,1,IDIR)   + R(J-1,K,1,IDIR)
     2                  + R(J,K-1,1,IDIR) + R(J-1,K-1,1,IDIR) )
               RC(J,K,IDIM+1,IDIR) = 0.25E0*
     1                  ( R(J,K,IDIM,IDIR)   + R(J-1,K,IDIM,IDIR)
     2                  + R(J,K-1,IDIM,IDIR) + R(J-1,K-1,IDIM,IDIR) )
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
C
C 1-d. The interior of the domain
C
      DO 240 IDIR = 1, 3
         DO 230 I = 2, IDIM
            DO 220 K = 2, KDIM
               DO 210 J = 2, JDIM
                  RC(J,K,I,IDIR) = 0.125E0*
     1                     (  R(J,K,I,IDIR)     + R(J,K,I-1,IDIR)
     2                      + R(J-1,K,I,IDIR)   + R(J-1,K,I-1,IDIR)
     3                      + R(J,K-1,I,IDIR)   + R(J,K-1,I-1,IDIR)
     4                      + R(J-1,K-1,I,IDIR) + R(J-1,K-1,I-1,IDIR) )
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
C
C 2. Modify cell centers of wall function points to reflect proper location
C
      IF (IFWALF) THEN
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM
         IJKDIM(2) = JDIM
         IJKDIM(3) = KDIM
C
C Loop over all bcs
C
         DO 400 IBC = 1, NBCS
            IBCTYP   = IBCDAT(1,IBC)
            IF (IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
               IDIR     = IBCDAT(2,IBC)
               INROUT   = IBCDAT(3,IBC)
               ISTRT(1) = IBCDAT(4,IBC)
               IEND (1) = IBCDAT(5,IBC)
               ISTRT(2) = IBCDAT(6,IBC)
               IEND (2) = IBCDAT(7,IBC)
               ISTRT(3) = IBCDAT(8,IBC)
               IEND (3) = IBCDAT(9,IBC)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
               DO 300 L = 1, 3
                  IF (L .NE. IDIR) THEN
                     IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                     IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L)+1
                  END IF
  300          CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C
C
               IF (INROUT .EQ. 1) THEN
                  IBCPT  = 1
                  IBCN   = 2
                  IBCN2  = 3
                  SGNM   = 1.E0
               ELSE
                  IBCPT  = IJKDIM(IDIR) + 1
                  IBCN   = IJKDIM(IDIR) 
                  IBCN2  = IJKDIM(IDIR) - 1
                  SGNM   = - 1.E0
               ENDIF
C
               IC1 = ICY (IDIR, 1)
               IC2 = ICY (IDIR, 2)
               IC3 = ICY (IDIR, 3)
C
C Reset solution point location for wall functions
C
               DO 320 I3 = ISTRT(IC3), IEND(IC3)
                  DO 310 I2 = ISTRT(IC2), IEND(IC2)
C
                     IW   = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JW   = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KW   = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN2  = IBCN2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN2  = IBCN2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN2  = IBCN2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
C Calculate distance from 2 to wall, 3 to 2, 3 to interface
C
                     DS2TOW =SQRT((RC(JN ,KN ,IN ,1)-RC(JW,KW,IW,1))**2
     1                           +(RC(JN ,KN ,IN ,2)-RC(JW,KW,IW,2))**2
     2                           +(RC(JN ,KN ,IN ,3)-RC(JW,KW,IW,3))**2)
                     DS3TO2 =SQRT((RC(JN2,KN2,IN2,1)-RC(JN,KN,IN,1))**2
     1                           +(RC(JN2,KN2,IN2,2)-RC(JN,KN,IN,2))**2
     2                           +(RC(JN2,KN2,IN2,3)-RC(JN,KN,IN,3))**2)
                     DS3TOI = (DS2TOW + DS3TO2) - 2.E0 * DS2TOW
C
C Linearly interpolate cell center of first cell to new location
C
                     FACT   = 2.E0 * DS3TOI / DS3TO2
                     OMFACT = 1.E0 - FACT
                     XC     = RC(JN,KN,IN,1)
                     YC     = RC(JN,KN,IN,2)
                     ZC     = RC(JN,KN,IN,3)
                     RC(JN,KN,IN,1) = FACT   * XC
     1                              + OMFACT * RC(JN2,KN2,IN2,1)
                     RC(JN,KN,IN,2) = FACT   * YC
     1                              + OMFACT * RC(JN2,KN2,IN2,2)
                     RC(JN,KN,IN,3) = FACT   * ZC
     1                              + OMFACT * RC(JN2,KN2,IN2,3)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  400    CONTINUE
      ENDIF
C
C Finished with grid averaging
C
      RETURN
      END





      SUBROUTINE RCCOPY (IDIM, JDIM, KDIM, RCIN, RCOUT)
C
C Routine to copy the input cell centered grid (RCIN) to the output cell
C centered grid (RCOUT).
C
C IDIM,JDIM,KDIM : Dimensions of current block
C RCIN           : Input  position vector of grid vertices at cell centers
C RCOUT          : Output position vector of grid vertices at cell centers
C RC             : Calculated position of grid at cell centers
C
      include '../header/common.h'
C
      DIMENSION RCIN  (JDIM+1,KDIM+1,IDIM+1,3),
     1          RCOUT (JDIM+1,KDIM+1,IDIM+1,3)
C
C Loop and copy
C
      DO 130 K = 1, KDIM+1
         DO 120 I = 1, IDIM+1
            DO 110 J = 1, JDIM+1
               RCOUT(J,K,I,1) = RCIN (J,K,I,1)
               RCOUT(J,K,I,2) = RCIN (J,K,I,2)
               RCOUT(J,K,I,3) = RCIN (J,K,I,3)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C Finished copying
C
      RETURN
      END






      SUBROUTINE RDDATA (NXQ, NXCFS, NXLVLS, NXBLKS, NXBCS,NXCUTS,NXPRT,
     1                   NCFS, MGLVLS, NBLKS,  IDIM, JDIM, KDIM,
     2                   NBCS, IBCDAT, NCUTS, ICUTS, CUTNAM,
     3                   NPRINT, IPRINT, NPRMON, IPRMON, NITMON,
     4                   ITUPJ, ITSLVL,NITS,NSUBIT,NITFO,NITBBC,NITALG,
     5                   NITRST, NXSTAG, NSTAGE, RKALPH,
     6                   RESTRT, FRDRST, FWRRST, GRDFIL, IGRDTP, GRDBIN,
     7                   NXPROF, NPROFL, FPROFL, 
     8                   IFPLT3, BINPLT, FPLT3G, FPLT3Q, FPLT3F,FPLT3FN,
     9                   IORDER, LIMIT, RKAP, IFLUX, VISCOS,DIAG,IMODEL,
     A                   IMSTRT, IMEND, TOLER, VGNRN, SIGMA, IFDELQ,
     B                   CONTRB, IMODLX, QMIN, ITRSMX, IFFORC, IFRCPR, 
     C                   YNCALC, NXTRSG, NTRSEG, ITRSEG, RCNTR)
C
C Routine to set default values and read input data controlling run.
C
C This routine is organized as follows:
C     1. Set defaults
C     2. Read input deck
C     3. Set values based on input
C     4. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C NXQ    : Maximum number of Q variables
C NXCFS  :
C NXLVLS :
C NXBLKS :
C NXBCS  :
C NXCUTS :
C NXPRT  :
C NCFS   :
C MGLVLS :
C NBLKS  :
C IDIM   :
C JDIM   :
C KDIM   :
C NBCS   :
C IBCDAT :
C NCUTS  :
C ICUTS  :
C NPRINT :
C IPRINT :
C NPRMON :
C IPRMON :
C NITMON :
C ITUPJ  :
C ITSLVL :
C NITS   : Number of iterations to perform this run.
C NSUBIT : Number of sub-iterations per time step.
C NITFO  : Number of iterations first order.  First order in all directions
C          for time dependent calculations.  First order in the streamwise
C          (marching) direction for marching calculations.
C NITBBC : Number of iterations of the blunt body condition - NOT USED.
C NITALG : Number of iterations to update the eddy viscosity using the
C          algebraic model.  Used for the initial iterations of higher
C          order models to assist convergence.  The higher order models
C          continue to be updated - but the eddy viscosity used for the
C          flow equations is calculated from the algebraic model.
C NITRST : Number of iterations to save intermediate restart file.
C NXSTAG : Maximum number of stages
C NSTAGE : Number of stages for temporal integration scheme
C          RK : Multistage scheme   NSTAGE = 3 (default or typical value)
C          AF : Single stage scheme NSTAGE = 1
C RKALPH : Coefficients for multistage scheme
C RESTRT :
C FRDRST :
C FWRRST :
C GRDFIL :
C IGRDTP :
C GRDBIN :
C NXPROF : Maximum number of profile boundary condition segments
C NPROFL : Number of profile boundary condition segments
C FPROFL : Array of filenames for profile boundary condition segments
C IFPLT3 :
C BINPLT :
C FPLTSG : File name for PLOT3D grid     file
C FPLT3Q : File name for PLOT3D solution file
C FPLT3F : File name for PLOT3D function file
C FPLT3FN: File name for PLOT3D function name file
C IORDER :
C LIMIT  :
C RKAP   :
C IFLUX  :
C VISCOS :
C DIAG   : Control of diagonalized scheme
C IMODEL : Turbulence model array
C IMSTRT :
C IMEND  :
C TOLER  :
C VGNRN  :
C SIGMA  :
C IFDELQ : Array controlling calculation of delq
C CONTRB : Logical.  If true then convert from one turbulence model to
C          another when do restart.
C IMODLX : Turbulence model array defining model to convert FROM if 
C          CONTRB is true.
C QMIN   : Control array for minimum value of Q
C ITRSMX : Number of iterations between output of maximum residual
C IFFORC : Control of calculation on forces
C IFRCPR : Control of output of forces
C YNCALC : Logical to force recalculation of Y_normal
C NXTRSG : Maximum number of transition specification segments
C NTRSEG : Actual  number of transition specification segments
C ITRSEG : Data for transition specification segments
C RCNTR  : Location of point vortex for 2D farfield bc. Also location to 
C          calculate moments about.
C
      include '../header/common.h'
      include '../header/histry.h'
C
      CHARACTER*80 FRDRST, FWRRST, GRDFIL(NXBLKS), 
     1             FPLT3G, FPLT3Q, FPLT3F, FPLT3FN, FPROFL(NXPROF)
      CHARACTER*80 FNAME
      CHARACTER*25 VAR, VAR1, BCTYP, SOLVER
C
      LOGICAL RESTRT, IFPLT3, GRDBIN, BINPLT, VGNRN, CONTRB, YNCALC
C
      DIMENSION IDIM   (NXBLKS,NXLVLS), 
     1          JDIM   (NXBLKS,NXLVLS),
     2          KDIM   (NXBLKS,NXLVLS)
      DIMENSION IBCDAT (10,NXBCS,NXBLKS,NXLVLS), 
     1          NBCS   (NXBLKS)
      DIMENSION ICUTS  (21,NXCUTS,NXLVLS)
      INTEGER   NCUTS
      CHARACTER*20 CUTNAM(NXCUTS)
      DIMENSION IPRINT (10,NXPRT,NXBLKS), 
     1          NPRINT (NXBLKS),
     2          IPRMON (10,NXPRT,NXBLKS), 
     3          NPRMON (NXBLKS)
      DIMENSION NITS   (NXCFS), 
     1          MGLVLS (NXCFS), 
     2          ITSLVL (NXLVLS,NXCFS)
      DIMENSION RKALPH (NXSTAG)
      DIMENSION IORDER (3,NXBLKS,NXLVLS,NXCFS), 
     1          LIMIT  (3,NXBLKS), 
     2          RKAP   (3,NXBLKS), 
     3          IFLUX  (NXBLKS)
      DIMENSION IMODEL (NMDL),
     1          IMODLX (NMDL)
      LOGICAL   VISCOS (3,NXBLKS),
     1          DIAG   (3,NXBLKS)
      DIMENSION IFDELQ (NXQ)
      DIMENSION QMIN   (NXQ)
C
C     Transition specification
C
      DIMENSION ITRSEG (7,NXTRSG)
C
      DIMENSION RCNTR  (3)
C
C     Local variables
C
      DIMENSION ISTRT (3), IEND(3), IORD(3)
      DIMENSION IDUM  (3)
      DIMENSION ITORDR(3)
C
C     Temporary storage to read print data into to avoid overwriting memory
C
      DIMENSION IDUMPR(10)
C
C 1.  Set defaults
C     Initialize INPERR to be 0.  This is incremented if there is an error
C     in the input dek so that the entire input deck is scanned for potential
C     errors without running an incorrect case.  INPERR counts the number
C     of errors.
C
      INPERR  = 0
C
C     Freestream Quantities
C
      RHOINF = 1.E0
      AINF   = 1.E0
C
C     Initialize variables to default values
C
      NQ     = 5
      NF     = 5
      NP     = 5
      NRANK  = 5
      NTIME  = 1
      NTMTAU = 0
C
      RESTRT = .FALSE.
      INITTD = .FALSE.
      VGNRN  = .FALSE.
      SIGMA  = 0.95E0
C
      IFPLT3 = .FALSE.
      BINPLT = .TRUE.
C
      FOURTH = .FALSE.
      NITBBC = 0
      THNLYR = .FALSE.
      GLOBAL = .TRUE.
      THREED = .TRUE.
      AXISYM = .FALSE.
      TWSPE  = .FALSE.
      SIUNIT = .FALSE.
      YNCALC = .FALSE.
C
      AXIDTH = 1.E0
C
C     Initialize number of profile boundary condition segments to zero
C
      NPROFL = 0
C
C     Grid type
C
      IGRDTP = IGTSEP
      GRDBIN = .TRUE.
C
C     Flow Variables
C
      GAMMA  = 1.4E0
      FSMACH = 0.8E0
      ALPHA  = 0.E0
      BETA   = 0.E0
      RE     = 1.E06
      PR     = 0.7E0
      TINF   = 450.E0
      TWALL  = 1060.E0
      AREARF = 1.E0
C
C     Location of point vortex boundary condition and location for moments
C
      RCNTR(1) = 0.25E0
      RCNTR(2) = 0.E0
      RCNTR(3) = 0.E0
C
C     Subsonic inflow/outflow bc conditions.  Initialize to zero.  If input
C     does not reset to non-zero value, then set based on input freestream.
C     Input T_total as Kelvins or degrees R, 
C           P_back, P_total as p/p_inf 
C     and non-dimensionalize correctly at end of routine.
C
      PBAKBC = 0.E0
      PTOTBC = 0.E0
      TTOTBC = 0.E0
C
C     Jet conditions 
C     Velocities are input non-dimensionalized by ainf
C     Pressure is input as p(jet)/pinf 
C        - this is later converted to p(jet)/rhoinf*ainf**2
C     Temperature is input as T(jet)/Tinf
C
      UJET   = 0.E0
      VJET   = 0.E0
      WJET   = 0.E0
      PJET   = 1.E0
      TJET   = 1.E0
C
      TKEJET = 1.E-04
      EPSJET = 1.E-04
      OMEJET = 1.E-04
C
      TAUJET(1) = 2.E0 / 3.E0 * TKEJET
      TAUJET(2) = 2.E0 / 3.E0 * TKEJET
      TAUJET(3) = 2.E0 / 3.E0 * TKEJET
      TAUJET(4) = 0.E0
      TAUJET(5) = 0.E0
      TAUJET(6) = 0.E0
C
C     MUSCL Variables, diagonalized control
C
      DO 3 IBLK = 1, NXBLKS
         RKAP  (1,IBLK) = 0.E0
         RKAP  (2,IBLK) = 0.E0
         RKAP  (3,IBLK) = 0.E0
         LIMIT (1,IBLK) = ILMNMD
         LIMIT (2,IBLK) = ILMNMD
         LIMIT (3,IBLK) = ILMNMD
         IFLUX (IBLK)   = IFROE
         DIAG  (1,IBLK) = .TRUE.
         DIAG  (2,IBLK) = .TRUE.
         DIAG  (3,IBLK) = .TRUE.
    3 CONTINUE
C
C     Entropy Fix
C
      CNTRPY = 0.E0
C
C     Artificial dissipation for central difference scheme
C
      ADKP46 = 1.E0 / 64.E0
C
C     Constants for dQ/dt difference formula - default to first-order
C
      TDTHET = 1.E0
      TDPHI  = 0.E0
      TDPHIP = 0.E0
      ITMTYP = ITTS
C
C     Viscous/Inviscid
C
      DO 6 IBLK = 1, NXBLKS
         VISCOS(1,IBLK) = .FALSE.
         VISCOS(2,IBLK) = .FALSE.
         VISCOS(3,IBLK) = .FALSE.
    6 CONTINUE
C
C     Turn off all calculations of delq
C
      DO 7 L = 1, NXQ
         IFDELQ(L) = 0
    7 CONTINUE
C
C     Turbulence model - set to inviscid flow, no wall functions
C                        turbulence production limit = 100
C
      ITURB     = ITNVSD
      PRT       = 0.9E0
      CONTRB    = .FALSE.
      IFWALF    = .FALSE.
      IDAMP     = 0
      IEPSEQ    = 0
      IEPSLN    = 0
      ISRCW     = INO
      IADRM     = IADRNO
      IASM      = IASMBU
      IEPSC     = IECNO
      IPDIL     = IPDNO
      ICMUST    = 0
      PRDLIM    = 100.E0
C
C     Set default to no transition specification
C
      NTRSEG    = 0
C
C     Set default to limit the production on the epsilon/omega equation also
C
      PRDE      = 1.E0
      PRDEM1    = 1.E0 - PRDE
C
C     Set default to higher order advection on turbulence equations
C
      FOTURB    = .FALSE.
C
C     Set default to not enforce positivity of production
C
      POSPRD    = .FALSE.
C
C     Set default convert turbulence model to default turbulence model
C
      ITURBX    = ITURB
      IDAMPX    = IDAMP
      IASMX     = IASM
C
C     Storage locations for wall shear stress and heat transfer
C
      ITQTAU    = 3
      ITQDOT    = 4
C
C     Run Control Variables
C          NITS   <= MXHIST (maximum size of history array)
C          ITDFNL >= 2 to avoid divide by zero when evaluate (ITDFNL - 1)
C          NITRST  = MXHIST + 1 to turn off intermediate restart file output
C
C     Multigrid control
C
      NCFS      = 1
      DO 8 ICFS = 1, NXCFS
         MGLVLS(ICFS) = 1
         NITS  (ICFS) = MIN (200, MXHIST)
    8 CONTINUE
      CSMTH     = 0.25E0
      CMXCHG    = 0.25E0
C
      NITRST    = MXHIST + 1
      NSUBIT    = 1
      NSTAGE    = 1
      NITFO     = 0
      NITALG    = 0
      ISRCE     = .TRUE.
      ISOLVR    = IAF3F
      IMPSRC    = ISDIAG
      ITIMED    = ISTDY
      NBLKS     = 0
      NCUTS     = 0
      ITUPDT    = 1
      ITUPJ     = 1
      CFL       = 1.E0
      DT        = 0.E0
      CFLFNL    = 0.E0
      TIMEDT    = - 1.E0
      ITDFNL    = MAX (NITS(1), 2)
      ITRSMX    = MXHIST + 1
      IFFORC    = INO
      IFRCPR    = IYES
C
C     Control of FIX routines for positivity preservation
C
      IFIXQ     = IFIXMN
C
      QMIN( 1)  =   1.E-12
      QMIN( 2)  = - 1.E0
      QMIN( 3)  = - 1.E0
      QMIN( 4)  = - 1.E0
      QMIN( 5)  =   1.E-12
C
C     Output Control Variables
C
      NITMON    = 999999
      DO 9 IBLK = 1, NXBLKS
         NPRINT(IBLK) = 0
         NPRMON(IBLK) = 0
    9 CONTINUE
C
      FWRRST    = 'isaac.rst'
C
      DO 20 ICFS = 1, NXCFS
         DO 10 ILVL = 1, NXLVLS
            ITSLVL(ILVL,ICFS) = 1
   10    CONTINUE
   20 CONTINUE
C
C 2.  Read input deck
C     This is the beginning of the loop to read in and process the input
C
      WRITE (IOUT,1001)
C
C     Read card, set variable based on input
C
   50 CONTINUE
      READ (INPUT,*,END=990,ERR=995) VAR, VALUE
      WRITE (IOUT,1000) VAR, VALUE
      IF (VAR .EQ. 'END') THEN
         GO TO 900
C
C     Restart Control
C
      ELSE IF (VAR .EQ. 'RESTART IN') THEN
         RESTRT = .TRUE.
         READ (INPUT,*,END=990,ERR=995) FRDRST
         WRITE (IOUT,1025) FRDRST
      ELSE IF (VAR .EQ. 'RESTART OUT') THEN
         READ (INPUT,*,END=990,ERR=995) FWRRST
         WRITE (IOUT,1005) FWRRST
      ELSE IF (VAR .EQ. 'SAVE RESTART') THEN
         NITRST = INT(VALUE)
C
C     Plot Output
C
      ELSE IF (VAR .EQ. 'PLOT3D') THEN
         IFPLT3 = .TRUE.
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'BINARY') THEN
            BINPLT = .TRUE.
         ELSE IF (VAR .EQ. 'FORMATTED') THEN
            BINPLT = .FALSE.
         ELSE
            WRITE (IOUT,1010) VAR
            INPERR = INPERR + 1
         ENDIF
         READ (INPUT,*,END=990,ERR=995) FNAME
         FPLT3G  = FNAME
         FPLT3Q  = FNAME
         FPLT3F  = FNAME
         FPLT3FN = FNAME
         CALL CHAPND (FPLT3G,  '.g')
         CALL CHAPND (FPLT3Q,  '.q')
         CALL CHAPND (FPLT3F,  '.qt')
         CALL CHAPND (FPLT3FN, '.qt.nam')
         IF (BINPLT) THEN
            CALL CHAPND (FPLT3G, '.unf')
            CALL CHAPND (FPLT3Q, '.unf')
            CALL CHAPND (FPLT3F, '.unf')
         ELSE
            CALL CHAPND (FPLT3G, '.fmt')
            CALL CHAPND (FPLT3Q, '.fmt')
            CALL CHAPND (FPLT3F, '.fmt')
         ENDIF
         WRITE (IOUT,1015) FPLT3G, FPLT3Q, FPLT3F
C
C     Solution Procedure
C
      ELSE IF (VAR .EQ. 'SOLVER') THEN
         READ (INPUT,*,END=990,ERR=995) SOLVER
         IF (SOLVER .EQ. 'RK') THEN
            ISOLVR    = IRKN
            NSTAGE    = INT (VALUE)
            IF (NSTAGE .LE. 1) THEN
C
C     Default coefficients from Rick Gaffney optimized for kappa = 1/3 scheme
C     give maximum CFL = 1.76
C
               NSTAGE    = 3
               RKALPH(1) = 0.27E0
               RKALPH(2) = 0.50E0
               RKALPH(3) = 1.0E0
            ELSE
C
C     Read input coefficients
C
               DO 96 ISTAGE = 1, NSTAGE
                  READ (INPUT,*,END=990,ERR=995) RKALPH(ISTAGE)
   96          CONTINUE
            ENDIF
            IF (NSTAGE .GT. NXSTAG) THEN
               WRITE (IOUT,1020) NSTAGE, NXSTAG
               INPERR = INPERR + 1
            ENDIF
            WRITE (IOUT,1021) NSTAGE, (RKALPH(ISTAGE),ISTAGE=1,NSTAGE)
         ELSE IF (SOLVER .EQ. 'AF3F') THEN
            ISOLVR = IAF3F
            DO 98 IBLK = 1, NXBLKS
               DIAG(1,IBLK) = .FALSE.
               DIAG(2,IBLK) = .FALSE.
               DIAG(3,IBLK) = .FALSE.
   98       CONTINUE
         ELSE IF (SOLVER .EQ. 'DIAGONAL AF3F') THEN
            ISOLVR = IAF3F
            DO 99 IBLK = 1, NXBLKS
               DIAG(1,IBLK) = .TRUE.
               DIAG(2,IBLK) = .TRUE.
               DIAG(3,IBLK) = .TRUE.
   99       CONTINUE
         ELSE IF (SOLVER .EQ. 'LGS') THEN
            ISOLVR = ILGS
         ELSE IF (SOLVER .EQ. 'MARCH') THEN
            ISOLVR = IMARCH
            READ (INPUT,*,END=990,ERR=995) IMSTRT, IMEND
            WRITE (IOUT,1022) IMSTRT, IMEND
         ELSE
            WRITE (IOUT,1023) SOLVER
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'INITIALIZE TIME DEPENDENT') THEN
         INITTD = .TRUE.
      ELSE IF (VAR .EQ. 'IMPLICIT TIME DEPENDENT') THEN
         ITIMED = IUNSTD
         TDTHET = 1.0E0
         TDPHI  = 0.5E0
         TDPHIP = 0.E0
         NTIME  = 2
         NSUBIT = 2
         DO 53 IIMP = 1, INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'TIME STEP') THEN
               TIMEDT = VALUE
            ELSE IF (VAR .EQ. 'SUB ITERATIONS') THEN
               NSUBIT = INT(VALUE)
            ELSE IF (VAR .EQ. 'T-TS') THEN
               ITMTYP = ITTS
            ELSE IF (VAR .EQ. 'TAU-TS') THEN
               ITMTYP = ITAUTS
               TDPHIP = 0.E0
               NTMTAU = 0
            ELSE IF (VAR .EQ. 'TAU-TS 2ND ORDER') THEN
               ITMTYP = ITAUTS
               TDPHIP = 0.5E0
               NTMTAU = 1
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
   53    CONTINUE
      ELSE IF (VAR .EQ. 'CFL') THEN
         CFL = VALUE
      ELSE IF (VAR .EQ. 'DT') THEN
         DT  = VALUE
      ELSE IF (VAR .EQ. 'CFLFNL') THEN
         CFLFNL = VALUE
      ELSE IF (VAR .EQ. 'ITUPDT') THEN
         ITUPDT = INT(VALUE)
      ELSE IF (VAR .EQ. 'ITDFNL') THEN
         ITDFNL = MAX (INT(VALUE), 2)
      ELSE IF (VAR .EQ. 'ITUPJ') THEN
         ITUPJ = INT(VALUE)
      ELSE IF (VAR .EQ. 'NITS') THEN
         NITS(1) = INT(VALUE)
         IF (NITS(1) .GT. MXHIST) THEN
            WRITE (IOUT,1024) NITS(1), MXHIST
            NITS(1) = MXHIST
         ENDIF
      ELSE IF (VAR .EQ. 'ENTROPY FIX') THEN
         CNTRPY = VALUE
         IF (CNTRPY .LT. 0.E0 .OR. CNTRPY .GT. 1.E0) THEN
            WRITE (IOUT,1124) 
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'PRINT MAX RESIDUAL') THEN
         ITRSMX = INT(VALUE)
      ELSE IF (VAR .EQ. 'RECALCULATE YNORMAL') THEN
         ITEMP  = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            YNCALC = .FALSE.
         ELSE
            YNCALC = .TRUE.
         ENDIF
      ELSE IF (VAR .EQ. 'PRINT FORCES') THEN
         ITEMP  = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            IFRCPR = INO
         ELSE IF (ITEMP .EQ. 1) THEN
            IFRCPR = IYES
         ELSE
            WRITE (IOUT,1026) ITEMP
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'REFERENCE AREA') THEN
         AREARF = VALUE
      ELSE IF (VAR .EQ. 'MULTIGRID LEVELS') THEN
         DO 54 ICFS = 1, NXCFS
            MGLVLS(ICFS) = MAX ((INT(VALUE) + 1 - ICFS), 1)
   54    CONTINUE
      ELSE IF (VAR .EQ. 'IMPLICIT SMOOTHING') THEN
         CSMTH  = VALUE
      ELSE IF (VAR .EQ. 'MESH SEQUENCE') THEN
         NCFS    = INT(VALUE) + 1
         IF (NCFS .GT. NXCFS) THEN
            WRITE (IOUT,1210) NCFS, NXCFS, NCFS-1
            INPERR = INPERR + 1
C     Skip input lines if error
            DO 55 ICFS = NCFS, 2, -1
               READ (INPUT,*,END=990,ERR=995) 
   55       CONTINUE
         ELSE
            DO 60 ICFS = NCFS, 2, -1
               READ (INPUT,*,END=990,ERR=995) VALUE
               NITS  (ICFS) = INT(VALUE)
               WRITE (IOUT,1220) ICFS, NITS(ICFS)
               IF (NITS(ICFS) .GT. MXHIST) THEN
                  WRITE (IOUT,1024) NITS(ICFS), MXHIST
                  NITS(ICFS) = MXHIST
               ENDIF
   60       CONTINUE
         ENDIF
      ELSE IF (VAR .EQ. 'NITFO') THEN
         NITFO  = INT(VALUE)
      ELSE IF (VAR .EQ. 'FIRST ORDER TURBULENCE') THEN
         FOTURB = .TRUE.
      ELSE IF (VAR .EQ. 'ENFORCE POSITIVE PROD.') THEN
         POSPRD = .TRUE.
      ELSE IF (VAR .EQ. 'NITALG') THEN
         NITALG = INT(VALUE)
      ELSE IF (VAR .EQ. 'SOURCE') THEN
         ITEMP = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            ISRCE = .FALSE.
         ELSE IF (ITEMP .EQ. 1) THEN
            ISRCE = .TRUE.
         ELSE
            WRITE (IOUT,1027) ITEMP
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'IMPLICIT SOURCE') THEN
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1065) VAR
         IF (VAR .EQ. 'NONE') THEN
            IMPSRC = INO
         ELSE IF (VAR .EQ. 'DIAGONAL') THEN
            IMPSRC = ISDIAG
         ELSE IF (VAR .EQ. 'BLOCK') THEN
            IMPSRC = ISBLOC
         ELSE
            INPERR = INPERR + 1
            WRITE (IOUT,1028) VAR
         ENDIF
C
C     Control of FIX positivity routines
C
      ELSE IF (VAR .EQ. 'FIXQ') THEN
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'NONE') THEN
            IFIXQ  = INO
         ELSE IF (VAR .EQ. 'MINIMUM') THEN
            IFIXQ  = IFIXMN
         ELSE IF (VAR .EQ. 'AVERAGE') THEN
            IFIXQ  = IFIXAV
         ELSE 
            WRITE (IOUT,1229) VAR
            INPERR = INPERR + 1
         ENDIF
C
C     Input transition specification segments
C
      ELSE IF (VAR .EQ. 'ENFORCE LAMINAR FLOW') THEN
         NTRSEG = INT (VALUE)
         IF (NTRSEG .GT. NXTRSG) THEN
            WRITE (IOUT,1545) NTRSEG, NXTRSG
            INPERR = INPERR + 1
            GO TO 900
         ENDIF
C
C      Input each segment specification; Check that ending index is >=
C      beginning index; Switch from cell node to cell center.
C
         DO ISEG = 1, NTRSEG
            READ (INPUT,*,END=990,ERR=995) IBLOCK, (ISTRT(L),L=1,3), 
     1                                             (IEND (M),M=1,3)
            WRITE (IOUT,1580) IBLOCK, (ISTRT(L),L=1,3),
     1                                (IEND (M),M=1,3)
            DO L = 1, 3
               IF (IEND(L) .LT. ISTRT(L)) THEN
                  WRITE (IOUT,1590) L, IBLOCK, (ISTRT(LL),LL=1,3),
     1                                         (IEND (MM),MM=1,3)
                  INPERR = INPERR + 1
               ENDIF
            END DO
            ITRSEG(1,ISEG) = IBLOCK
            ITRSEG(2,ISEG) = ISTRT(1) + 1
            ITRSEG(3,ISEG) = ISTRT(2) + 1
            ITRSEG(4,ISEG) = ISTRT(3) + 1
            ITRSEG(5,ISEG) = IEND (1)
            ITRSEG(6,ISEG) = IEND (2)
            ITRSEG(7,ISEG) = IEND (3)
         END DO
C
C     Control of Marching Algorithm
C
      ELSE IF (VAR .EQ. 'VIGNERON') THEN
         VGNRN  = .TRUE.
      ELSE IF (VAR .EQ. 'TOLER') THEN
         TOLER = VALUE
      ELSE IF (VAR .EQ. 'SIGMA') THEN
         SIGMA = VALUE
c     ELSE IF (VAR .EQ. 'BLUNT BODY') THEN
c        BLUNT  = .TRUE.
c        NITBBC = INT(VALUE)
C
C     Viscous Terms
C
      ELSE IF (VAR .EQ. 'THIN LAYER') THEN
         THNLYR = .TRUE.
      ELSE IF (VAR .EQ. 'NAVIER STOKES') THEN
         THNLYR = .FALSE.
C
C     Turbulence Models
C
      ELSE IF (VAR .EQ. 'BALDWIN LOMAX') THEN
         ITURB  = ITBLMX
C
C     SIGK and SIGT2 are defined as the inverse of sigma(k) and sigma(epsilon)
C
      ELSE IF (VAR .EQ. 'K EPSILON') THEN
C
C     Set default k-epsilon model to be Zhang, So, Gatski, and Speziale
C
C     SAA : Set k-epsilon turbulence constants from Speziale, Abid 
C          and Anderson AIAA 90-1481 as
C             sigma k       = 1.36
C             sigma k       = 1.00 - CORRECTED TO THIS 
C                             (Abid had recommended 0.75 for supersonic flows)
C             sigma epsilon = 1.36
C             C2 epsilon    = 1.83
C        ICASE 90-9
C             alf1 kepsilon = 1.0 without pressure dilatation
C             alf1 kepsilon = 0.6    with pressure dilatation (Rick Gaffney)
C        ICASE 90-42
C             alf2 kepsilon = 0.4 
C             alf3 kepsilon = 0.2
C        CPDFRO,CPDFT = 0 as these models are not good
C     ZSGS: Zhang, So, Gatski, Speziale in Near-Wall Turbulent Flows (1993)
C     RNG : Yakhot, Orszag, Thangam, Gatski, Speziale in Phys. Fluids A
C
C     Default values
C
         ITURB  = ITKE
         IEPSEQ = IEEZSG
         IDAMP  = IFMZSG
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 7
         IFDL2K = 1
         NQ     = 7
         NP     = 5
         NF     = 7
         NRANK  = 7
         SIGK   = 1.0E0 / 1.00E0
         SIGT2  = 1.0E0 / 1.45E0
         SIGRHO = 1.0E0 / 0.95E0
         CMU    = 0.096E0
         CEPS1  = 1.50E0
         CEPS2  = 1.83E0
         A2KE   = 4.9E0
         ALF1   = 0.0E0
         ALF2   = 0.0E0
         ALF3   = 0.0E0
C--->    CPDFRO = SIGRHO
C--->    CPDFT  = 1.E0 / PRT
         CPDFRO = 0.E0
         CPDFT  = 0.E0
         CPVELK = 0.E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-04 
         RMUTNF = 1.0E0
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-8
         QMIN( 7)  =   1.E-8
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 101 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'SAA') THEN
               IEPSEQ = IEESAA
               IDAMP  = IFMSAA
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.36E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.09E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.6E0
               ALF2   = 0.4E0
               ALF3   = 0.2E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ZSGS') THEN
               IEPSEQ = IEEZSG
               IDAMP  = IFMZSG
               IFDL2K = 1
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.45E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.096E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ZSSL') THEN
               IEPSEQ = IEEZSL
               IDAMP  = IFMZSG
               IFDL2K = 1
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.45E0
               SIGRHO = 1.0E0 / 0.50E0
               CMU    = 0.096E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'HIGH RE') THEN
               IEPSEQ = IEEHR
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.30E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.09E0
               CEPS1  = 1.44E0
               CEPS2  = 1.92E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'RNG') THEN
               IEPSEQ = IEERNG
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.39E0
               SIGT2  = 1.39E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.085E0
               CEPS1  = 1.42E0
               CEPS2  = 1.68E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ASM') THEN
               IASM   = IASMGS
               ICMUST = 5
               IEPSEQ = IEEABD
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGRHO = 0.0E0
               CMU    = 0.081E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               A2KE   = 5.5E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
               SIGT2  = 1.0E0 / ( 0.41e0**2/((CEPS2-CEPS1)*SQRT(CMU)) )
            ELSE IF (VAR .EQ. 'ADRM') THEN
               IADRM  = IADRGS
               IEPSEQ = IEEHR
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.87E0
               SIGRHO = 0.0E0
               CMU    = 0.094E0
               CEPS1  = 1.20E0
               CEPS2  = 1.83E0
               A2KE   = 0.0E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMAE') THEN
               SIGT2  = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'CEPS1') THEN
               CEPS1  = VALUE
            ELSE IF (VAR .EQ. 'CEPS2') THEN
               CEPS2  = VALUE
            ELSE IF (VAR .EQ. 'A2KE') THEN
               A2KE   = VALUE
            ELSE IF (VAR .EQ. 'ALF1') THEN
               ALF1   = VALUE
            ELSE IF (VAR .EQ. 'ALF2') THEN
               ALF2   = VALUE
            ELSE IF (VAR .EQ. 'ALF3') THEN
               ALF3   = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL EPSILON PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN(6) = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM EPSILON') THEN
               QMIN(7) = VALUE
            ELSE IF (VAR .EQ. 'COMPRESSIBLE DISS') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IEPSC  = IECNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IEPSC  = IECSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IEPSC  = IECRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1068) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE DILATATION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IPDIL  = IPDNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IPDIL  = IPDSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IPDIL  = IPDRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1069) VAR
               ENDIF
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  101    CONTINUE
C
C     These variables are used only for output
C
         IF (IEPSEQ .EQ. IEESAA) THEN
            WRITE (IOUT,1070) 'Speziale-Abid-Anderson        '
         ELSE IF (IEPSEQ .EQ. IEEHR ) THEN
            WRITE (IOUT,1070) 'High Reynolds Number          '
         ELSE IF (IEPSEQ .EQ. IEEZSG) THEN
            WRITE (IOUT,1070) 'Zhang-So-Gatski-Speziale      '
         ELSE IF (IEPSEQ .EQ. IEEZSL) THEN
            WRITE (IOUT,1070) 'Zhang-So-Speziale-Lai         '
         ELSE IF (IEPSEQ .EQ. IEERNG) THEN
            WRITE (IOUT,1070) 'RNG (Yak-Ors-Than-Gat-Spez)   '
         ELSE IF (IEPSEQ .EQ. IEEABD) THEN
            WRITE (IOUT,1070) 'Abid ASM (Abid-Morr-Gat-Spez) '
         ENDIF
         SK     = 1.E0 / SIGK
         ST2    = 1.E0 / SIGT2
         WRITE (IOUT,1075) CEPS1, CEPS2, CMU, SK, ST2, SIGRHO, A2KE,
     1                     ALF1, ALF2, ALF3, CPDFRO, CPDFT, CPVELK,
     2                     TKEINF, RMUTNF
C
C     k-omega constants from Wilcox, AIAA Journal, Vol. 26, No. 11
C
      ELSE IF (VAR .EQ. 'K OMEGA') THEN
C
C     Default values
C
         ITURB  = ITKW
         IEPSEQ = 0
         IDAMP  = IFMHR
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 7
         IFDL2K = 0
         NQ     = 7
         NP     = 5
         NF     = 7
         NRANK  = 7
         SIGK   = 0.5E0
         SIGT2  = 0.5E0
         SIGRHO = 0.0E0
         BSTRKW = 0.09E0
         BKW    = 3.0E0 / 40.0E0
         CMU    = 1.0E0
         GKW    = 5.0E0 / 9.0E0
         ROUGHK = 0.0E0
         CPDFRO = 0.0E0
         CPDFT  = 0.0E0
         CPVELK = 0.0E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-06 
         RMUTNF = 1.0E-01
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-12
         QMIN( 7)  =   1.E-8
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 102 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'ASM') THEN
               IASM   = IASMGS
               ICMUST = 5
               BSTRKW = 1.0E0
               CMU    = 0.088E0
               SIGK   = 1.E0 / 1.4E0
               SIGT2  = 1.E0 / 2.0E0
               BKW    = 0.83E0
               GKW    = BKW - (0.41E0)**2 / (SQRT(CMU) / SIGT2)
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = VALUE
            ELSE IF (VAR .EQ. 'SIGMAW') THEN
               SIGT2  = VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'BSTRKW') THEN
               BSTRKW = VALUE
            ELSE IF (VAR .EQ. 'BKW') THEN
               BKW    = VALUE
            ELSE IF (VAR .EQ. 'GKW') THEN
               GKW    = VALUE
            ELSE IF (VAR .EQ. 'ROUGHK') THEN
               ROUGHK = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL OMEGA PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN(6) = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM OMEGA') THEN
               QMIN(7) = VALUE
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  102    CONTINUE
C
C     These variables are used only for output
C
         WRITE (IOUT,1080) SIGK, SIGT2, BSTRKW, BKW, CMU, GKW, ROUGHK,
     1                     SIGRHO, CPDFRO, CPDFT, CPVELK,
     2                     TKEINF, RMUTNF
C
C     Reynolds stress model
C
      ELSE IF (VAR .EQ. 'REYNOLDS STRESS') THEN
C
C     Set Reynolds Stress turbulence constants 
C     AIAA 90-1481 as
C          sigma k       =  1.36
C          sigma k       =  0.75 CORRECTED TO THIS
C          sigma epsilon =  1.36
C          C  epsilon    =  0.10
C          C1 epsilon    =  1.50
C          C2 epsilon    =  1.83
C          C3 epsilon    = -2.00
C          C4 epsilon    = -1.50
C     ICASE 90-9
C          alf1 kepsilon = 1.0 without pressure dilatation
C          alf1 kepsilon = 0.6    with pressure dilatation (Rick Gaffney)
C          alf1          = 0   for Lai-So
C     ICASE 90-42
C          alf2          = 0   for Lai-So
C          alf3          = 0   for Lai-So
C          alf2 kepsilon = 0.4 
C          alf3 kepsilon = 0.2
C     CPDFRO,CPDFT = 0 as these models are not good
C     Pressure-Strain Correlation
C       Lai-So
C          cpi1 = 3.0
C          cpi2 = 0.6
C
C     Default values
C
         ITURB  = ITRS
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 12
         IFDL2K = 1
         NQ     = 12
         NP     =  5
         NF     = 12
         NRANK  = 12
         SIGK   = 1.0E0 / 0.75E0
         SIGT2  = 1.0E0 / 1.36E0
         SIGRHO = 1.0E0 / 0.50E0
         CMU    = 0.09E0
         CEPSLN = 0.10E0
         CEPS1  = 1.50E0
         CEPS2  = 1.83E0
         CEPS3  = - 2.0E0
         CEPS4  = - 1.5E0
         CEPS5  = 0.0E0
         CSUBS  = 0.11E0
         A2KE   = 4.9E0
         ALF1   = 0.0E0
         ALF2   = 0.0E0
         ALF3   = 0.0E0
C--->    CPDFRO = SIGRHO
C--->    CPDFT  = 1.E0 / PRT
         CPDFRO = 0.E0
         CPDFT  = 0.E0
         CPVELK = 0.E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-04 * SQRT(10.E0)
         RMUTNF = 1.0E0 * 10.E0
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-8
         QMIN( 7)  =   1.E-8
         QMIN( 8)  =   1.E-8
         QMIN( 9)  = - 1.E0
         QMIN(10)  = - 1.E0
         QMIN(11)  = - 1.E0
         QMIN(12)  =   1.E-8
C
C     Pressure-Strain Correlation
C
C     Lai-So form (Zhang's) Model
C
         IDMMDL = IDMSO
         IDTMDL = IDTHL
         IPSMDL = IPLS
         CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 103 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'LRR') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.15E0
               CEPS1  = 1.44E0
               CEPS2  = 1.90E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPLRR1
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SSG') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.1842E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSSG
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SSG-ADRM') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.094E0
               CEPSLN = 0.217E0
               CEPS1  = 1.00E0
               CEPS2  = 1.83E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 0.0E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSSG
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
               IEPSIJ = IEADGS
            ELSE IF (VAR .EQ. 'FLT') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.15E0
               CEPS1  = 1.45E0
               CEPS2  = 1.90E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPFLT
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'ZSGS') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.10E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               CEPS3  = - 2.0E0
               CEPS4  = - 1.5E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPLS
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SO 95') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.12E0
               CEPS1  = 1.50E0
               CEPS2  = 1.90E0
               CEPS3  = 7.E0/(3.E0*SQRT(15.E0)) * 0.0066E0
               CEPS4  = 2.95E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSO95
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'GIBSON LAUNDER') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.16E0
               CEPS1  = 1.44E0
               CEPS2  = 1.92E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPGL
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMAE') THEN
               SIGT2  = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'CEPSLN') THEN
               CEPSLN = VALUE
            ELSE IF (VAR .EQ. 'CEPS1') THEN
               CEPS1  = VALUE
            ELSE IF (VAR .EQ. 'CEPS2') THEN
               CEPS2  = VALUE
            ELSE IF (VAR .EQ. 'CEPS3') THEN
               CEPS3  = VALUE
            ELSE IF (VAR .EQ. 'CEPS4') THEN
               CEPS4  = VALUE
            ELSE IF (VAR .EQ. 'CEPS5') THEN
               CEPS5  = VALUE
            ELSE IF (VAR .EQ. 'CSUBS') THEN
               CSUBS  = VALUE
            ELSE IF (VAR .EQ. 'A2KE') THEN
               A2KE   = VALUE
            ELSE IF (VAR .EQ. 'ALF1') THEN
               ALF1   = VALUE
            ELSE IF (VAR .EQ. 'ALF2') THEN
               ALF2   = VALUE
            ELSE IF (VAR .EQ. 'ALF3') THEN
               ALF3   = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL EPSILON PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'PSC1') THEN
               PSC1   = VALUE
            ELSE IF (VAR .EQ. 'PSC2') THEN
               PSC2   = VALUE
            ELSE IF (VAR .EQ. 'PSCW') THEN
               PSCW   = VALUE
            ELSE IF (VAR .EQ. 'PSALFA') THEN
               PSALFA = VALUE
            ELSE IF (VAR .EQ. 'PSBETA') THEN
               PSBETA = VALUE
            ELSE IF (VAR .EQ. 'PSGAMA') THEN
               PSGAMA = VALUE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN( 6) = 2.E0 / 3.E0 * VALUE
               QMIN( 7) = 2.E0 / 3.E0 * VALUE
               QMIN( 8) = 2.E0 / 3.E0 * VALUE
            ELSE IF (VAR .EQ. 'MINIMUM EPSILON') THEN
               QMIN(12) = VALUE
            ELSE IF (VAR .EQ. 'TURBULENT DIFFUSION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'DH') THEN
                  IDTMDL = IDTDH
               ELSE IF (VAR .EQ. 'HL') THEN
                  IDTMDL = IDTHL
               ELSE IF (VAR .EQ. 'MH') THEN
                  IDTMDL = IDTMH
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1066) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE STRAIN') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'LRR1') THEN
                  IPSMDL = IPLRR1
               ELSE IF (VAR .EQ. 'LRR2') THEN
                  IPSMDL = IPLRR2
               ELSE IF (VAR .EQ. 'SSG') THEN
                  IPSMDL = IPSSG
               ELSE IF (VAR .EQ. 'FLT') THEN
                  IPSMDL = IPFLT
               ELSE IF (VAR .EQ. 'LAI SO') THEN
                  IPSMDL = IPLS
               ELSE IF (VAR .EQ. 'GIBSON LAUNDER') THEN
                  IPSMDL = IPGL
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1067) VAR
               ENDIF
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'COMPRESSIBLE DISS') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IEPSC  = IECNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IEPSC  = IECSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IEPSC  = IECRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1068) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE DILATATION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1227) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IPDIL  = IPDNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IPDIL  = IPDSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IPDIL  = IPDRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1069) VAR
               ENDIF
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  103    CONTINUE
C
C     These variables are used only for output
C
         SK     = 1.E0 / SIGK
         ST2    = 1.E0 / SIGT2
         WRITE (IOUT,1090) CEPSLN, CEPS1, CEPS2, CEPS3, CEPS4, CEPS5, 
     1                     CSUBS, CMU, SK, ST2, SIGRHO, A2KE,
     2                     ALF1, ALF2, ALF3, CPDFRO, CPDFT, CPVELK,
     3                     TKEINF, RMUTNF
         IF (IPSMDL .EQ. IPLRR1) THEN
            WRITE (IOUT,1095) 'LRR (Launder-Reece-Rodi)      '
         ELSE IF (IPSMDL .EQ. IPLRR2) THEN
            WRITE (IOUT,1095) 'LRR2 (Launder-Reece-Rodi 2)   '
         ELSE IF (IPSMDL .EQ. IPSSG ) THEN
            WRITE (IOUT,1095) 'SSG (Speziale-Sarkar-Gatski)  '
         ELSE IF (IPSMDL .EQ. IPFLT ) THEN
            WRITE (IOUT,1095) 'FLT (Fu-Launder-Tselepidakis) '
         ELSE IF (IPSMDL .EQ. IPLS  ) THEN
            WRITE (IOUT,1095) 'LS (LAI-SO)                   '
         ELSE IF (IPSMDL .EQ. IPGL  ) THEN
            WRITE (IOUT,1095) 'GL (Gibson-Launder)           '
         ENDIF
         WRITE (IOUT,1097) PSC1, PSC2, PSCW, PSALFA, PSBETA, PSGAMA
C
C     Large-Eddy Simulation (LES)
C
      ELSE IF (VAR .EQ. 'LES') THEN
C
C     Set default LES model to Smagorinsky high Reynolds number
C
         ITURB = ITLES
         ISGS  = ISGSMG
         IDAMP = IFMHR
         NQ    = 5
         NP    = 5
         NF    = 5
         NRANK = 5
         CMU   = 0.2E0
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read additional input lines
C
         NCON   = INT (VALUE)
         DO 203 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'SMAGORINSKY') THEN
               ISGS  = ISGSMG
               IDAMP = IFMHR
               CMU   = 0.2E0
            ELSE IF (VAR .EQ. 'STRUCTURE FUNCTION MODEL') THEN
               ISGS  = ISGSFM
               IDAMP = IFMHR
               CKSFM = 1.4E0
               CMU   = 0.105E0 * CKSFM**(-3.E0/2.E0)
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  203    CONTINUE
C
         IF (ISGS .EQ. ISGSMG) THEN
            WRITE (IOUT,1098) 'Smagorinsky                   '
         ELSE IF (ISGS .EQ. ISGSFM) THEN
            WRITE (IOUT,1098) 'Structure Function Model      '
         ENDIF
         WRITE (IOUT,1099) CMU
C
C     Convert turbulence model on restart
C
      ELSE IF (VAR .EQ. 'CONVERT') THEN
         CONTRB = .TRUE.
         READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'LAMINAR') THEN
            ITURBX = ITLMNR
            WRITE (IOUT,1110) VAR
         ELSE IF (VAR .EQ. 'BALDWIN LOMAX') THEN
            ITURBX = ITBLMX
            WRITE (IOUT,1110) VAR
         ELSE IF (VAR .EQ. 'K EPSILON') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITKE
            IDAMPX = IFMZSG
            IASMX  = IASMBU
            NCON   = INT(VALUE)
            DO 104 ICON = 1, NCON
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'SAA') THEN
                  IDAMPX = IFMSAA
               ELSE IF (VAR .EQ. 'ZSGS') THEN
                  IDAMPX = IFMZSG
               ELSE IF (VAR .EQ. 'ZSSL') THEN
                  IDAMPX = IFMZSG
               ELSE IF (VAR .EQ. 'HIGH RE') THEN
                  IDAMPX = IFMHR
               ELSE IF (VAR .EQ. 'RNG') THEN
                  IDAMPX = IFMHR
               ELSE IF (VAR .EQ. 'ASM') THEN
                  ICMUST = 5
                  IASMX  = IASMGS
               ELSE
                  WRITE (IOUT,1105) VAR
                  INPERR = INPERR + 1
               END IF
  104       CONTINUE
         ELSE IF (VAR .EQ. 'K OMEGA') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITKW
            IDAMPX = IFMHR
            IASMX  = IASMBU
            NCON   = INT(VALUE)
            DO 105 ICON = 1, NCON
               READ (INPUT,*,END=990,ERR=995) VAR
               IF (VAR .EQ. 'ASM') THEN
                  ICMUST = 5
                  IASMX  = IASMGS
               ELSE
                  WRITE (IOUT,1105) VAR
                  INPERR = INPERR + 1
               END IF
  105       CONTINUE
         ELSE IF (VAR .EQ. 'REYNOLDS STRESS') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITRS
            IDAMPX = IFMZSG
            IASMX  = IASMBU
         ELSE
            WRITE (IOUT,1100) VAR
            INPERR = INPERR + 1
         ENDIF
C
C        Set up IMODLX to control conversion model 
C        Only require ITURBX, IDAMPX, and IASMX set all else to error code
C
         DO 106 IM = 1, NMDL
            IMODLX(IM) = - 1
  106    CONTINUE
         IMODLX( 1) = ITURBX
         IMODLX( 3) = IDAMPX
         IMODLX(10) = IASMX
C
C     End of Turbulence Models
C
C     Metric/English units
C
      ELSE IF (VAR .EQ. 'SIUNIT') THEN
         SIUNIT = .TRUE.
      ELSE IF (VAR .EQ. 'ENGLISH') THEN
         SIUNIT = .FALSE.
C
C     Two/Three-Dimensional
C
      ELSE IF (VAR .EQ. 'TWOD') THEN
         THREED = .FALSE.
         AXISYM = .FALSE.
      ELSE IF (VAR .EQ. 'THREED') THEN
         THREED = .TRUE.
         AXISYM = .FALSE.
      ELSE IF (VAR .EQ. 'AXISYMMETRIC') THEN
         THREED = .FALSE.
         AXISYM = .TRUE.
      ELSE IF (VAR .EQ. 'GRID') THEN
         READ (INPUT,*,END=990,ERR=995) GRDFIL(1)
         WRITE (IOUT,1040) GRDFIL(1)
         IGRDTP = IGTP3D
      ELSE IF (VAR .EQ. 'GRID BINARY') THEN
         GRDBIN = .TRUE.
      ELSE IF (VAR .EQ. 'GRID FORMATTED') THEN
         GRDBIN = .FALSE.
C
C     Flow Quantities
C
      ELSE IF (VAR .EQ. 'GAMMA') THEN
         GAMMA  = VALUE
      ELSE IF (VAR .EQ. 'MACH') THEN
         FSMACH = VALUE
      ELSE IF (VAR .EQ. 'ALPHA') THEN
         ALPHA = VALUE
      ELSE IF (VAR .EQ. 'BETA') THEN
         BETA = VALUE
      ELSE IF (VAR .EQ. 'RE') THEN
         RE = VALUE
      ELSE IF (VAR .EQ. 'PR') THEN
         PR = VALUE
      ELSE IF (VAR .EQ. 'PRT') THEN
         PRT = VALUE
      ELSE IF (VAR .EQ. 'TINF') THEN
         TINF = VALUE
      ELSE IF (VAR .EQ. 'TWALL') THEN
         TWSPE = .TRUE.
         TWALL = VALUE
C
C     MUSCL Control
C
      ELSE IF (VAR .EQ. 'KAPPA') THEN
         DO 108 IBLK = 1, NXBLKS
            RKAP(1,IBLK) = VALUE
            RKAP(2,IBLK) = VALUE
            RKAP(3,IBLK) = VALUE
  108    CONTINUE
      ELSE IF (VAR .EQ. 'LIMIT') THEN
         LIM  = INT(VALUE)
         LIMX = INT(VALUE)
         IF (LIMX .EQ. ILNONE) THEN
            LIM = ILNONE
         ELSE IF (LIMX .EQ. ILSMTH) THEN
            LIM = ILSMTH
         ELSE IF (LIMX .EQ. ILMNMD) THEN
            LIM = ILMNMD
         ELSE IF (LIMX .EQ. ILVNKT) THEN
            LIM = ILVNKT
         ELSE
            WRITE (IOUT,1058) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
         DO 110 IBLK = 1, NXBLKS
            LIMIT(1,IBLK) = LIM
            LIMIT(2,IBLK) = LIM
            LIMIT(3,IBLK) = LIM
  110    CONTINUE
C
C     Fourth order
C
      ELSE IF (VAR .EQ. 'CENTRAL DIFFERENCE') THEN
  114    READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'FOURTH ORDER') THEN
            FOURTH = .TRUE.
         ELSE IF (VAR .EQ. 'ARTIFICIAL DISSIPATION') THEN
            ADKP46 = VALUE
         ELSE IF (VAR .EQ. 'END CENTRAL DIFFERENCE') THEN
            GO TO 118
         ELSE
            WRITE (IOUT,1059) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Loop for next control of central difference
C
         GO TO 114
C
C     Finished with central difference control
C
  118    CONTINUE
C
C     Subsonic Inflow/Outflow BC conditions
C
      ELSE IF (VAR .EQ. 'BACK PRESSURE') THEN
         PBAKBC = VALUE
      ELSE IF (VAR .EQ. 'INFLOW PTOTAL') THEN
         PTOTBC = VALUE
      ELSE IF (VAR .EQ. 'INFLOW TTOTAL') THEN
         TTOTBC = VALUE
C
C     Jet Conditions
C
      ELSE IF (VAR .EQ. 'JET CONDITIONS') THEN
         JET  = INT(VALUE)
C
C     Read in Jet conditions
C
  120    READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'U') THEN
            UJET = VALUE
         ELSE IF (VAR .EQ. 'V') THEN
            VJET = VALUE
         ELSE IF (VAR .EQ. 'W') THEN
            WJET = VALUE
         ELSE IF (VAR .EQ. 'P') THEN
            PJET = VALUE
         ELSE IF (VAR .EQ. 'T') THEN
            TJET = VALUE
         ELSE IF (VAR .EQ. 'K') THEN
            TKEJET = VALUE
         ELSE IF (VAR .EQ. 'EPSILON') THEN
            EPSJET = VALUE
         ELSE IF (VAR .EQ. 'OMEGA') THEN
            OMEJET = VALUE
         ELSE IF (VAR .EQ. 'TAUXX') THEN
            TAUJET(1) = VALUE
         ELSE IF (VAR .EQ. 'TAUYY') THEN
            TAUJET(2) = VALUE
         ELSE IF (VAR .EQ. 'TAUZZ') THEN
            TAUJET(3) = VALUE
         ELSE IF (VAR .EQ. 'TAUXY') THEN
            TAUJET(4) = VALUE
         ELSE IF (VAR .EQ. 'TAUXZ') THEN
            TAUJET(5) = VALUE
         ELSE IF (VAR .EQ. 'TAUYZ') THEN
            TAUJET(6) = VALUE
         ELSE IF (VAR .EQ. 'END JET') THEN
            GO TO 129
         ELSE
            WRITE (IOUT,1060) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Loop for next Jet condition
C
         GO TO 120
C
C     Finished inputing Jet conditions
C
  129    CONTINUE
C
C     Printout control
C
      ELSE IF (VAR .EQ. 'PRINT') THEN
         NPRT = INT(VALUE)
         DO 300 IPRT = 1, NPRT
            READ (INPUT,*,END=990,ERR=995) IBLOCK, (IDUMPR(LL),LL=1,10)
            WRITE (IOUT,2005)              IBLOCK, (IDUMPR(LL),LL=1,10)
            NPRINT(IBLOCK) = NPRINT(IBLOCK) + 1
            IF (NPRINT(IBLOCK) .GT. NXPRT) THEN
               WRITE (IOUT,1540) NPRINT(IBLOCK), IBLOCK, NXPRT
               GO TO 900
            ENDIF
            DO 299 LL = 1, 10
               IPRINT(LL,NPRINT(IBLOCK),IBLOCK) = IDUMPR(LL)
 299        CONTINUE
 300     CONTINUE
      ELSE IF (VAR .EQ. 'MONITOR') THEN
         NPRT = INT(VALUE)
         READ (INPUT,*,END=990,ERR=995) NITMON
         DO 310 IPRT = 1, NPRT
            READ (INPUT,*,END=990,ERR=995) IBLOCK,
     1                (IPRMON(LL,NPRMON(IBLOCK)+1,IBLOCK),LL=1,10)
            WRITE (IOUT,2005) IBLOCK, 
     1                (IPRMON(LL,NPRMON(IBLOCK)+1,IBLOCK),LL=1,10)
            NPRMON(IBLOCK) = NPRMON(IBLOCK) + 1
 310     CONTINUE
C
C     CUT boundary condition data
C
      ELSE IF (VAR .EQ. 'CUT') THEN
         NCUTS = INT(VALUE)
         IF (NCUTS*2 .GT. NXCUTS) THEN
            WRITE (IOUT,1530) NCUTS*2, NXCUTS
            GO TO 900
         ENDIF
         DO 329 ICUT = 1, NCUTS
            IC1 = 2 * ICUT - 1
            IC2 = 2 * ICUT
c           READ (INPUT,*,END=990,ERR=995) VAR, IBLK, (ISTRT(L),L=1,3), 
c    1                                        (IEND(M),M=1,3)
            read (input,*,END=990,ERR=995) var, iblk,istrt(1),istrt(2),
     1           istrt(3), iend(1), iend(2), iend(3)
            WRITE (IOUT,1600) VAR, IBLK, (ISTRT(L),L=1,3), 
     1                                   (IEND(M),M=1,3)
C
C     Determine CUT face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a cut. In case an error occurs, set the
C     face direction to 1 to finish checking the CUT data.
C
            IDIR = 0
            DO 320 L = 1, 3
               IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
  320       CONTINUE
            IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
               WRITE (IOUT,1630) 
               INPERR = INPERR + 1
               IDIR   = 1
            ENDIF
C
C     Modify start or end appropriately to account for cut input as grid
C     points but cut to be cell-centers.  Do this only for the indices that
C     the cut loops over, not the idir direction.
C
            DO 321 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (IEND(L) .GT. ISTRT(L)) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ELSE
                     IEND(L)  = IEND(L)  + 1
                  ENDIF
               ENDIF
  321       CONTINUE
C
            CUTNAM(IC1) = VAR
            CUTNAM(IC2) = VAR
C
            ICUTS( 1,IC1,1) = IBLK
            ICUTS( 2,IC1,1) = IDIR
            ICUTS( 3,IC1,1) = ISTRT(IDIR)
            ICUTS( 4,IC1,1) = ISTRT(1)
            ICUTS( 5,IC1,1) = IEND (1)
            ICUTS( 6,IC1,1) = ISTRT(2)
            ICUTS( 7,IC1,1) = IEND (2)
            ICUTS( 8,IC1,1) = ISTRT(3)
            ICUTS( 9,IC1,1) = IEND (3)
C
            ICUTS(10,IC2,1) = IBLK
            ICUTS(11,IC2,1) = IDIR
            ICUTS(12,IC2,1) = ISTRT(IDIR)
            ICUTS(13,IC2,1) = ISTRT(1)
            ICUTS(14,IC2,1) = IEND (1)
            ICUTS(15,IC2,1) = ISTRT(2)
            ICUTS(16,IC2,1) = IEND (2)
            ICUTS(17,IC2,1) = ISTRT(3)
            ICUTS(18,IC2,1) = IEND (3)
C
            READ (INPUT,*,END=990,ERR=995) VAR1, IBLK, (ISTRT(L),L=1,3), 
     1                             (IEND(M),M=1,3), (IORD(N),N=1,3)
c           read (input,*,END=990,ERR=995) var1, iblk,istrt(1),istrt(2),
c    1           istrt(3), iend(1), iend(2), iend(3)
            WRITE (IOUT,1605) VAR1, IBLK, (ISTRT(L),L=1,3), 
     1                             (IEND(M),M=1,3), (IORD(N),N=1,3)
            IF (VAR .NE. VAR1) THEN
               INPERR = INPERR + 1
               WRITE (IOUT,1610) 
            ENDIF
C
            DO 324 L = 1, 3
               IF (IORD(L) .LT. 1 .OR. IORD(L) .GT. 3) THEN
                  INPERR = INPERR + 1
                  WRITE (IOUT,1620) L, IORD(L)
               ENDIF
  324       CONTINUE
C
C     Determine CUT face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a cut. In case an error occurs, set the
C     face direction to 1 to finish checking the CUT data.
C
            IDIR = 0
            DO 330 L = 1, 3
               IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
  330       CONTINUE
            IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
               WRITE (IOUT,1630) 
               INPERR = INPERR + 1
               IDIR   = 1
            ENDIF
C
C     Modify start or end appropriately to account for cut input as grid
C     points but cut to be cell-centers.  Do this only for the indices that
C     the cut loops over, not the idir direction.
C
            DO 331 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (IEND(L) .GT. ISTRT(L)) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ELSE
                     IEND(L)  = IEND(L)  + 1
                  ENDIF
               ENDIF
  331       CONTINUE
C
            ICUTS(10,IC1,1) = IBLK
            ICUTS(11,IC1,1) = IDIR
            ICUTS(12,IC1,1) = ISTRT(IDIR)
            ICUTS(13,IC1,1) = ISTRT(1)
            ICUTS(14,IC1,1) = IEND (1)
            ICUTS(15,IC1,1) = ISTRT(2)
            ICUTS(16,IC1,1) = IEND (2)
            ICUTS(17,IC1,1) = ISTRT(3)
            ICUTS(18,IC1,1) = IEND (3)
            ICUTS(19,IC1,1) = IORD (1)
            ICUTS(20,IC1,1) = IORD (2)
            ICUTS(21,IC1,1) = IORD (3)
C
            ICUTS( 1,IC2,1) = IBLK
            ICUTS( 2,IC2,1) = IDIR
            ICUTS( 3,IC2,1) = ISTRT(IDIR)
            ICUTS( 4,IC2,1) = ISTRT(1)
            ICUTS( 5,IC2,1) = IEND (1)
            ICUTS( 6,IC2,1) = ISTRT(2)
            ICUTS( 7,IC2,1) = IEND (2)
            ICUTS( 8,IC2,1) = ISTRT(3)
            ICUTS( 9,IC2,1) = IEND (3)
            ICUTS(19,IC2,1) = ICY(IORD(1),1)
            ICUTS(20,IC2,1) = ICY(IORD(2),1)
            ICUTS(21,IC2,1) = ICY(IORD(3),1)
  329    CONTINUE
         NCUTS = NCUTS * 2
C
C     Data for each grid block/zone 
C
      ELSE IF (VAR .EQ. 'BLOCK') THEN
         IBLOCK = INT(VALUE)
         NBLKS  = NBLKS + 1
         IF (NBLKS .GT. NXBLKS) THEN
            WRITE (IOUT,1525) NBLKS, NXBLKS
            GO TO 900
         ENDIF
  500    CONTINUE
         READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'DIMENSIONS') THEN
            READ (INPUT,*,END=990,ERR=995) IDIM(IBLOCK,1), 
     1                             JDIM(IBLOCK,1), KDIM(IBLOCK,1)
            WRITE (IOUT,1035) IDIM(IBLOCK,1), 
     1                             JDIM(IBLOCK,1), KDIM(IBLOCK,1)
         ELSE IF (VAR .EQ. 'KAPPA') THEN
            READ (INPUT,*,END=990,ERR=995) (RKAP(LL,IBLOCK),LL=1,3)
            WRITE (IOUT,2000) ' Kappa ', (RKAP(LL,IBLOCK),LL=1,3)
C
C N.B. this input needs modified to use the limiter mnemonics
C
         ELSE IF (VAR .EQ. 'LIMIT') THEN
            READ (INPUT,*,END=990,ERR=995) (LIMIT(LL,IBLOCK),LL=1,3)
            WRITE (IOUT,2002) ' Limit ', (LIMIT(LL,IBLOCK),LL=1,3)
C-       ELSE IF (VAR .EQ. 'ORDER') THEN
C-          READ (INPUT,*,END=990,ERR=995) (ITORDR(LL),LL=1,3)
C-          WRITE (IOUT,2002) ' Order ', (ITORDR(LL),LL=1,3)
C-          DO 509 ILVL = 1, NXLVLS
C-             IORDER(1,IBLOCK,ILVL,ILVL) = ITORDR(1)
C-             IORDER(2,IBLOCK,ILVL,ILVL) = ITORDR(2)
C-             IORDER(3,IBLOCK,ILVL,ILVL) = ITORDR(3)
C-509       CONTINUE
         ELSE IF (VAR .EQ. 'FLUX') THEN
            READ (INPUT,*,END=990,ERR=995) VAR
            WRITE (IOUT,2010) VAR
            IF (VAR .EQ. 'ROE') THEN
               IFLUX(IBLOCK) = IFROE
            ELSE
               WRITE (IOUT,1030)
               INPERR = INPERR + 1
            ENDIF
         ELSE IF (VAR .EQ. 'VISCOUS') THEN
            IV = INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) (IDUM(LL),LL=1,IV)
            WRITE (IOUT,1035) (IDUM(LL),LL=1,IV)
            DO 90 LL = 1, IV
               VISCOS(IDUM(LL),IBLOCK) = .TRUE.
 90         CONTINUE
            ITURB = MAX (ITURB, ITLMNR)
         ELSE IF (VAR .EQ. 'DIAGONAL') THEN
            IV = INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) (IDUM(LL),LL=1,IV)
            WRITE (IOUT,1035) (IDUM(LL),LL=1,IV)
            DO 91 LL = 1, IV
               DIAG  (IDUM(LL),IBLOCK) = .TRUE.
 91         CONTINUE
         ELSE IF (VAR .EQ. 'INVISCID') THEN
            DO 95 LL = 1, 3
               VISCOS(LL,IBLOCK) = .FALSE.
 95         CONTINUE
         ELSE IF (VAR .EQ. 'GRID') THEN
            READ (INPUT,*,END=990,ERR=995) GRDFIL(IBLOCK)
            WRITE (IOUT,1040) GRDFIL(IBLOCK)
            IGRDTP = IGTSEP
C
C     Boundary Condition Data
C
         ELSE IF (VAR .EQ. 'BC') THEN
            NBCS(IBLOCK) = INT (VALUE)
            IF (NBCS(IBLOCK) .GT. NXBCS) THEN
               INPERR = INPERR + 1
               WRITE (IOUT,1535) NBCS(IBLOCK), IBLOCK, NXBCS
               GO TO 900
            ENDIF
            DO 200 I = 1, NBCS(IBLOCK)
c              READ (INPUT,*,END=990,ERR=995) VAR, (ISTRT(L),L=1,3), 
c    1                                (IEND(M),M=1,3), BCTYP
               read (input,*,END=990,ERR=995) var, istrt(1), istrt(2),
     1              istrt(3), iend(1), iend(2), iend(3), bctyp
               WRITE (IOUT,1045) VAR, (ISTRT(L),L=1,3), 
     1                                (IEND(M),M=1,3), BCTYP
C
C     Determine BC face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a bc. In case an error occurs, set the
C     face direction to 1 to finish checking the BC data.
C
               IDIR = 0
               DO 150 L = 1, 3
                  IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
 150           CONTINUE
               IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
                  WRITE (IOUT,1049) 
                  INPERR = INPERR + 1
                  IDIR   = 1
               ENDIF
C
C     Set profile boundary segment identifier to zero for all other cases
C
               INPROF = 0
C
               IF (BCTYP .EQ. 'WALL') THEN
                  IBCDAT(1,I,IBLOCK,1) = IWALL
               ELSE IF (BCTYP .EQ. 'WALL FUNCTIONS') THEN
                  IBCDAT(1,I,IBLOCK,1) = IWALFN
                  IFWALF = .TRUE.
               ELSE IF (BCTYP .EQ. 'TANGENCY') THEN
                  IBCDAT(1,I,IBLOCK,1) = ITAN
               ELSE IF (BCTYP .EQ. 'FARFIELD') THEN
                  IBCDAT(1,I,IBLOCK,1) = IFAR
               ELSE IF (BCTYP .EQ. '2D FARFIELD VORTEX') THEN
                  IFFORC               = IYES
                  IBCDAT(1,I,IBLOCK,1) = IFAR2D
               ELSE IF (BCTYP .EQ. 'SYMMETRY XY') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMXY
               ELSE IF (BCTYP .EQ. 'SYMMETRY XZ') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMXZ
               ELSE IF (BCTYP .EQ. 'SYMMETRY YZ') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMYZ
               ELSE IF (BCTYP .EQ. 'EXTRAPOLATE') THEN
                  IBCDAT(1,I,IBLOCK,1) = IEXT
               ELSE IF (BCTYP .EQ. 'FIX') THEN
                  IBCDAT(1,I,IBLOCK,1) = IFIX
               ELSE IF (BCTYP .EQ. 'JET') THEN
                  IBCDAT(1,I,IBLOCK,1) = IJET
               ELSE IF (BCTYP .EQ. 'PERIODIC') THEN
                  IBCDAT(1,I,IBLOCK,1) = IPERD
               ELSE IF (BCTYP .EQ. 'HOMOGENEOUS SHEAR') THEN
                  IBCDAT(1,I,IBLOCK,1) = IHSHR
               ELSE IF (BCTYP .EQ. 'SUBSONIC INFLOW') THEN
                  IBCDAT(1,I,IBLOCK,1) = IINFLO
               ELSE IF (BCTYP .EQ. 'SUBSONIC OUTFLOW') THEN
                  IBCDAT(1,I,IBLOCK,1) = IOUTFL
               ELSE IF (BCTYP .EQ. 'AXISYMMETRIC') THEN
                  IBCDAT(1,I,IBLOCK,1) = IAXISM
               ELSE IF (BCTYP .EQ. 'AXISYMMETRIC CENTERLINE') THEN
                  IBCDAT(1,I,IBLOCK,1) = IAXICL
               ELSE IF (BCTYP .EQ. 'PROFILE CELL CENTER') THEN
                  IBCDAT(1,I,IBLOCK,1)  = IPROFL
                  NPROFL                = NPROFL + 1
                  INPROF                = NPROFL
                  READ (INPUT,*,END=990,ERR=995) FPROFL(NPROFL)
                  WRITE (IOUT,1047) FPROFL(NPROFL)
               ELSE IF (BCTYP .EQ. 'SUBSONIC PROFILE CC') THEN
                  IBCDAT(1,I,IBLOCK,1)  = IPROSB
                  NPROFL                = NPROFL + 1
                  INPROF                = NPROFL
                  READ (INPUT,*,END=990,ERR=995) FPROFL(NPROFL)
                  WRITE (IOUT,1048) FPROFL(NPROFL)
               ELSE
                  WRITE (IOUT,1050) BCTYP
                  INPERR = INPERR + 1
               ENDIF
C
C     Modify starting point appropriately to account for bc input as grid
C     points but bc to be cell-centers.  Do this only for the indices that
C     the bc loops over, not the idir direction.  BCs must be input with
C     indices in increasing order.
C     The starts are the max of the input starting point and 2 for the 
C     looping directions on a given boundary so that only loop over the 
C     interior boundary points to fit with the dimensioning of Q(+/-) etc 
C     as used in the MUSCL and FVSTL routines.
C
               DO 151 L = 1, 3
                  IF (L .NE. IDIR) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ENDIF
  151          CONTINUE
C
               IBCDAT( 2,I,IBLOCK,1) = IDIR
               IBCDAT( 3,I,IBLOCK,1) = ISTRT(IDIR)
               IBCDAT( 4,I,IBLOCK,1) = ISTRT(1)
               IBCDAT( 5,I,IBLOCK,1) = IEND (1)
               IBCDAT( 6,I,IBLOCK,1) = ISTRT(2)
               IBCDAT( 7,I,IBLOCK,1) = IEND (2)
               IBCDAT( 8,I,IBLOCK,1) = ISTRT(3)
               IBCDAT( 9,I,IBLOCK,1) = IEND (3)
               IBCDAT(10,I,IBLOCK,1) = INPROF
C
  200       CONTINUE
C
C     If read an END BLOCK, then have reached the end of the data for
C     the grid block/zone and are back in the main loop
C
         ELSE IF (VAR .EQ. 'END BLOCK') THEN
            GO TO 50
         ELSE
            WRITE (IOUT,1055) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Go back and get the next line of input from the block/zone data
C
         GO TO 500
C
C     End of data for grid block/zone
C
      ELSE
C
C     Error in main loop
C
         WRITE (IOUT,1060) VAR, VALUE
         INPERR = INPERR + 1
      ENDIF
C
C     Go back and read the next line
C
      GO TO 50
C
C     End of input file
C
  900 CONTINUE
      WRITE (IOUT,1002)
C
C 3.  Set values based on input
C
C     Iterative, implicit items
C     Total additional levels of time storage
C
      NTMLVL = NTIME + NTMTAU
C
C     Set calculation of forces: if output of forces is requested or the
C                                2D farfield vortex bc is used then calculate
C
      IF (IFRCPR .EQ. IYES) IFFORC = IYES
C
C     If unsteady, then verify using AF3F or DIAGONAL AF3F scheme.
C     Set TIMEDT to DT if not previously set; if DT not set then output error.
C     t-ts scheme must have DT = TIMEDT.
C
      IF (ITIMED .EQ. IUNSTD) THEN
         IF ( .NOT. (ISOLVR .EQ. IAF3F) ) THEN
            WRITE (IOUT,2015)
            INPERR = INPERR + 1
         ENDIF
C
         IF (TIMEDT .LE. 0.E0) THEN
            IF (DT .GT. 0.E0) THEN
               TIMEDT = DT
               WRITE (IOUT,2020) DT
            ELSE
               WRITE (IOUT,2030)
               INPERR = INPERR + 1
            ENDIF
         ENDIF
C
         IF (ITMTYP .EQ. ITTS) THEN
            IF (DT .NE. TIMEDT) THEN
               WRITE (IOUT,2040)
               INPERR = INPERR + 1
            ENDIF
         ENDIF
      ENDIF
C
C     Check for errors in the input deck
C
      IF (INPERR .GT. 0) THEN
         WRITE (IOUT,2050) INPERR
         STOP
      ENDIF
C
C     Set-up unit normal to use for collapsed cell faces
C     SHAT = (1,1,1) / sqrt(3) for 3D
C            (1,1,0) / sqrt(2) for 3D
C
      IF (THREED) THEN
         DENOM = SQRT (3.E0)
         SHATX = 1.E0 / DENOM
         SHATY = 1.E0 / DENOM
         SHATZ = 1.E0 / DENOM
      ELSE
         DENOM = SQRT (2.E0)
         SHATX = 1.E0 / DENOM
         SHATY = 1.E0 / DENOM
         SHATZ = 0.E0 
      ENDIF
C
C     Set-up IMODEL
C
      IMODEL( 1) = ITURB
      IMODEL( 2) = IEPSEQ
      IMODEL( 3) = IDAMP
      IMODEL( 4) = IADRM
      IMODEL( 5) = IFDL2K
      IF (ITURB .EQ. ITRS) THEN
         IMODEL( 6) = IPSMDL
         IMODEL( 7) = IEPSIJ
         IMODEL( 8) = IDMMDL
         IMODEL( 9) = IDTMDL
      ELSE IF (ITURB .EQ. ITLES) THEN
         IMODEL( 6) = ISGS
      ENDIF
      IMODEL(10) = IASM
      IMODEL(11) = IEPSC
      IMODEL(12) = IPDIL
      IMODEL(13) = ISRCW
C
C     Set values that depend on the entire range of input
C
      IF (ITURB .EQ. ITKE) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         EPSINF = CMU * RHOINF * TKEINF * TKEINF / RMUTNF * RE / FSMACH
      ELSE IF (ITURB .EQ. ITKW) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         OMEINF = CMU * RHOINF * TKEINF / RMUTNF * RE / FSMACH
      ELSE IF (ITURB .EQ. ITRS) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         EPSINF = CMU * RHOINF * TKEINF * TKEINF / RMUTNF * RE / FSMACH
C
C Isotropic turbulence
C
         TAUINF(1) = 2.E0 / 3.E0 * TKEINF
         TAUINF(2) = 2.E0 / 3.E0 * TKEINF
         TAUINF(3) = 2.E0 / 3.E0 * TKEINF
         TAUINF(4) = 0.E0
         TAUINF(5) = 0.E0
         TAUINF(6) = 0.E0
      ENDIF
C
C     Set PBAKBC, PTOTBC, TTOTBC if not input.
C     Convert PBAKBC, PTOTBC, TTOTBC from dimensional to non-dimensional.
C
      ARG    = 1.E0 + 0.5E0 * (GAMMA-1.E0) * FSMACH * FSMACH
      IF (PBAKBC .LE. 0.E0) PBAKBC = 1.E0
      IF (PTOTBC .LE. 0.E0) PTOTBC = ARG ** (GAMMA/(GAMMA-1.E0))
      IF (TTOTBC .LE. 0.E0) TTOTBC = TINF * ARG
C
      PBAKBC = PBAKBC / GAMMA
      PTOTBC = PTOTBC / GAMMA
      TTOTBC = TTOTBC / TINF
C
C     Convert PJET from P(jet) / Pinf to P(jet) / (rhoinf*ainf**2) 
C     Set the value of the Jet density
C
      PJET   = PJET / GAMMA
      RHOJET = GAMMA * PJET / TJET
C
C     Set TDPDTJ, TDPDEL for t-ts or tau-ts time stepping.
C
      IF (ITMTYP .EQ. ITTS) THEN
         TDPDTJ = TDPHI  + 1.E0
         TDPDEL = 0.E0
      ELSE
         TDPDTJ = TDPHIP + 1.E0
         TDPDEL = TDPHI  + 1.E0
      ENDIF
C
C     Sort the BCs into ascending order based on BC type
C
      WRITE (IOUT,2100)
      DO 980 IBLOCK = 1, NBLKS
         DO 970 I = 1, NBCS(IBLOCK)
            IBCTYP = IBCDAT(1,I,IBLOCK,1) 
            DO 960 J = I+1, NBCS(IBLOCK)
               IBCTST = IBCDAT(1,J,IBLOCK,1) 
               IF (IBCTST .LT. IBCTYP) THEN
                  IBCTYP = IBCTST
                  DO 950 L = 1, 10
                     IBCTMP               = IBCDAT(L,I,IBLOCK,1)
                     IBCDAT(L,I,IBLOCK,1) = IBCDAT(L,J,IBLOCK,1)
                     IBCDAT(L,J,IBLOCK,1) = IBCTMP
  950             CONTINUE
               ENDIF
  960       CONTINUE
  970    CONTINUE
  980 CONTINUE
C
C     Finished
C
      RETURN
C
C 4.   Error handling section
C 4-a. End-Of-File Handler
C
  990 CONTINUE
      WRITE (IOUT,9990)
      STOP
C
C 4-b. Error Handler
C
  995 CONTINUE
      WRITE (IOUT,9995)
      STOP
C
 1001 FORMAT (' ','RDDATA: Begin Input Echo')
 1002 FORMAT (' ','RDDATA: End Input Echo')
 1000 FORMAT (5X,A25,5X,E14.7)
 1005 FORMAT (' ',4X,'Restart being written to file:',
     1       /' ',9X,A70)
 1010 FORMAT (' ','Illegal PLOT3D format option = ',A25,
     1       /' ','Please specify either BINARY or FORMATTED ')
 1015 FORMAT (' ',4X,'Output PLOT3D Grid    : ',A50
     1       /' ',4X,'              Solution: ',A50,
     2       /' ',4X,'              Function: ',A50)
 1020 FORMAT (' ','Specified number of stages, NSTAGE ',I2,
     1        ' ','greater than allowed, MXSTAG ',I2,
     2       /' ','Reset MXSTAG in main and recompile.')
 1021 FORMAT (' ',4X,'Runge-Kutta: Number of Stages = ',I2,
     1       /' ',4X,'             Coefficients     = ',6F10.6)
 1022 FORMAT (5X,'IMSTRT = ',I5,' IMEND = ',I5)
 1023 FORMAT (' ','Illegal solution procedure specified ',A25)
 1024 FORMAT (' ','Specified number of iterations, NITS ',I10,
     1        ' ','greater than allowed, MXHIST ',I10,
     2       /' ','Reset to MXHIST and continue.')
 1124 FORMAT (' ','RDDATA: Error-> Invalid coefficient for entropy fix', 
     1       /' ','                CNTRPY must be 0 <= CNTRPY <= 1.')
 1025 FORMAT (' ',4X,'Restarting from file:',
     1       /' ',9X,A70)
 1026 FORMAT (' ','Illegal specification of force calculation = ',I10,
     1       /' ','     PRINT FORCES = 0 DO NOT calculate forces ',
     2       /' ','     PRINT FORCES = 1 Calculate forces ')
 1027 FORMAT (' ','Illegal specification of source term control = ',I10,
     1       /' ','     SOURCE = 0 DO NOT calculate source terms ',
     2       /' ','     SOURCE = 1 Calculate source terms ')
 1028 FORMAT (10X,'Invalid Implicit Source Term treatment ',A25)
 1210 FORMAT (' ','Required mesh levels ',I3,' exceed maximum ',I3,
     1       /' ','     Increase MXCFS in MAIN or decrease mesh ',
     2            'sequencing.'
     3       /' ','     Skipping next ',I3,' input lines.')
 1220 FORMAT (10X,'Mesh Sequence Level:',I4,' Number of Iterations:',I6)
 1227 FORMAT (10X,A25)
 1229 FORMAT (' ','Illegal specification of FIXQ  pos. control = ',I10,
     1       /' ','     FIXQ   = NONE    -> DO NOT call FIXQ ',
     2       /' ','     FIXQ   = MINIMUM -> Call FIXQ with Minimum ',
     3       /' ','     FIXQ   = AVERAGE -> Call FIXQ with Average ')
 1029 FORMAT (' ','Illegal turbulent constant specified ',A25)
 1030 FORMAT (///' ','***** ILLEGAL FLUX EVALUATION ROUTINE ****')
 1035 FORMAT (20X,3I5)
 1040 FORMAT (10X,'Grid File ',5X,A80)
 1045 FORMAT (10X,A25,6(1X,I4),2X,A25)
 1047 FORMAT (15X,'Profile File: ',5X,A80)
 1048 FORMAT (15X,'Subsonic Profile File: ',5X,A80)
 1049 FORMAT (' ','RDDATA: ERROR-> BC Face not specified. One pair of ',
     1            'BC indices must match',
     2       /' ','                to define the BC Face.')
 1050 FORMAT (' ','INVALID BOUNDARY CONDITION TYPE ',5X,A25)
 1055 FORMAT (10X,'Invalid Option in Block/Zone input ',A25,5X,E14.7)
 1058 FORMAT (' ','Invalid Limiter Specified ',5X,A25,5X,E14.7)
 1059 FORMAT (' ','Invalid Central Difference Option ',5X,A25,5X,E14.7)
 1060 FORMAT (' ','INVALID COMMAND',5X,A25,5X,E14.7)
 1061 FORMAT (10X,A25,5X,E14.7)
 1065 FORMAT (15X,A25)
 1066 FORMAT (10X,'Invalid Turbulent-DIffusion model ',A25)
 1067 FORMAT (10X,'Invalid Pressure-Strain model ',A25)
 1068 FORMAT (10X,'Invalid Compressible Dissipation model ',A25)
 1069 FORMAT (10X,'Invalid Pressure-Dilatation  model ',A25)
 1070 FORMAT (' ',5X,' K-Epsilon Model: ',A30)
 1075 FORMAT (' ',5X,' K-Epsilon Modelling Constants ',
     1       /' ',5X,'      CEPS1        = ',F12.8,
     2       /' ',5X,'      CEPS2        = ',F12.8,
     3       /' ',5X,'      CMU          = ',F12.8,
     4       /' ',5X,'      SIGMAK       = ',F12.8,
     5       /' ',5X,'      SIGMAE       = ',F12.8,
     6       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     7       /' ',5X,'      A2KE         = ',F12.8,
     8       /' ',5X,'      ALF1         = ',F12.8,
     9       /' ',5X,'      ALF2         = ',F12.8,
     A       /' ',5X,'      ALF3         = ',F12.8,
     B       /' ',5X,'      CPDFRO       = ',F12.8,
     C       /' ',5X,'      CPDFT        = ',F12.8,
     D       /' ',5X,'      CPVELK       = ',F12.8,
     E       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     F       /' ',5X,'      RMUTNF       = ',F12.8)
 1080 FORMAT (' ',5X,' K-Omega Modelling Constants ',
     1       /' ',5X,'      SIGMAK       = ',F12.8,
     2       /' ',5X,'      SIGMAW       = ',F12.8,
     3       /' ',5X,'      BSTRKW       = ',F12.8,
     4       /' ',5X,'      BKW          = ',F12.8,
     5       /' ',5X,'      CMU          = ',F12.8,
     6       /' ',5X,'      GKW          = ',F12.8,
     7       /' ',5X,'      ROUGHK       = ',F12.8,
     8       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     9       /' ',5X,'      CPDFRO       = ',F12.8,
     A       /' ',5X,'      CPDFT        = ',F12.8,
     B       /' ',5X,'      CPVELK       = ',F12.8,
     C       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     D       /' ',5X,'      RMUTNF       = ',F12.8)
 1090 FORMAT (' ',5X,' Reynolds Stress Modelling Constants ',
     1       /' ',5X,'      CEPSLN       = ',F12.8,
     2       /' ',5X,'      CEPS1        = ',F12.8,
     3       /' ',5X,'      CEPS2        = ',F12.8,
     4       /' ',5X,'      CEPS3        = ',F12.8,
     5       /' ',5X,'      CEPS4        = ',F12.8,
     6       /' ',5X,'      CEPS5        = ',F12.8,
     7       /' ',5X,'      CSUBS        = ',F12.8,
     8       /' ',5X,'      CMU          = ',F12.8,
     9       /' ',5X,'      SIGMAK       = ',F12.8,
     A       /' ',5X,'      SIGMAE       = ',F12.8,
     B       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     C       /' ',5X,'      A2KE         = ',F12.8,
     D       /' ',5X,'      ALF1         = ',F12.8,
     E       /' ',5X,'      ALF2         = ',F12.8,
     F       /' ',5X,'      ALF3         = ',F12.8,
     G       /' ',5X,'      CPDFRO       = ',F12.8,
     H       /' ',5X,'      CPDFT        = ',F12.8,
     I       /' ',5X,'      CPVELK       = ',F12.8,
     J       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     K       /' ',5X,'      RMUTNF       = ',F12.8)
 1095 FORMAT (' ',10X,' Pressure-Strain Correlation: ',A25)
 1097 FORMAT (' ',10X,'      PSC1         = ',F12.8,
     1       /' ',10X,'      PSC2         = ',F12.8,
     2       /' ',10X,'      PSCW         = ',F12.8,
     3       /' ',10X,'      PSALFA       = ',F12.8,
     4       /' ',10X,'      PSBETA       = ',F12.8,
     5       /' ',10X,'      PSGAMA       = ',F12.8)
 1098 FORMAT (' ',5X,' LES Subgrid-Scale Stress Model: ',A30)
 1099 FORMAT (' ',5X,' Large-Eddy Simulation Model Constants ',
     1       /' ',5X,'      C_s          = ',F12.8)
C
 1100 FORMAT (' ',5X,'Invalid conversion turbulence model specified',
     1               ' = ',A25)
 1105 FORMAT (' ',10X,'Invalid model specified for conversion ',A25)
 1110 FORMAT (' ',5X,'Convert on Restart From ',A25)
C
 1525 FORMAT (/' ','RDDATA: Number of blocks required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXBLKS in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1530 FORMAT (/' ','RDDATA: Number of cuts   required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXCUTS in MAIN and re-submit.',
     3        /' ', 8X,'   (Storage requires 2*number of cuts)',
     4        /' ', 8X,'   Unable to continue reading input.',
     5        /' ', 8X,'   RUN ABORTING!')
 1535 FORMAT (/' ','RDDATA: Number of BCs    required ',I5,
     1             ' for block ',I3,' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXBCS  in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1540 FORMAT (/' ','RDDATA: Number of print  required ',I5,
     1             ' for block ',I3,' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXPRT  in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1545 FORMAT (/' ','RDDATA: Number of transition segments required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXTRSG in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
C
 1580 FORMAT (10X,'Block = ',I5,5X,'Start = ',3I4,5X,'End = ',3I4)
 1590 FORMAT (' ','RDDATA: ERROR-> Starting index is greater than ',
     1            'ending index for transition at:'
     2       /' ','IDIR = ',I2,'Block = ',I5,5X,'Start = ',3I4,5X,
     3            'End = ',3I4)
C
 1600 FORMAT (10X,A25,I5,6(1X,I4))
 1605 FORMAT (10X,A25,I5,9(1X,I4))
 1610 FORMAT (' ','CUT names do not match')
 1620 FORMAT (' ','IORD(',I1,') =',I5,' does not equal 1, 2, or 3.')
 1630 FORMAT (' ','RDDATA: ERROR-> CUT Face not specified. One pair ',
     1            'of CUT indices must match',
     2       /' ','                to define the CUT Face.')
 2000 FORMAT (10X,A25,3E16.7)
 2002 FORMAT (10X,A25,3I5)
 2005 FORMAT (5X,'   Block = ',I4,'  IDIR = ',I4,'  Start = ',3I4,
     1           '  End = ',3I4,'  Increment = ',3I4)
 2010 FORMAT (10X,A25,5X,E14.7)
 2015 FORMAT (/' ','RDDATA: ERROR-> Iterative, implicit time ',
     1             'dependent currently only supports the',
     2        /' ','                AF3F or DIAGONAL AF3F solvers.',
     3        /' ','                Please switch to one of these and ',
     4             'resubmit.')
 2020 FORMAT (/' ','RDDATA: WARNING-> Iterative, implicit time ',
     1             'dependent time step set to DT ',1PE14.7)
 2030 FORMAT (/' ','RDDATA: ERROR-> Iterative, implicit time dependent',
     1             ' requires time step.',
     2        /' ',8X,'Input time step with TIME STEP or DT.')
 2040 FORMAT (/' ','RDDATA: ERROR-> t-ts Iterative, implicit time ',
     1             'dependent requires DT = TIME STEP',
     2        /' ',8X,'Check input TIME STEP and DT.')
 2050 FORMAT (/' ','RDDATA: ERROR-> There were ',I5,' errors in the ',
     1             'input deck.',
     2        /' ', 8X,'RUN ABORTING!')
 2100 FORMAT (//' ','RDDATA: Sorting boundary conditions into ',
     1              'ascending order based on type',//)
 9990 FORMAT (' ','RDDATA: ERROR-> Unexpected End-Of-File encountered.',
     1       /' ', 8X,'RUN ABORTING!')
 9995 FORMAT (' ','RDDATA: ERROR-> Syntax error in input line.'
     1       /' ', 8X,'RUN ABORTING!')
C
      END





      SUBROUTINE SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW)
C
C This routine sets the following constants for the various pressure-strain 
C correlations:
C      IEPSIJ, IEPSEQ, IDAMP, ISRCW, PSC1, PSC2, PSCW, PSALFA, PSBETA, PSGAMA
C
C IPSMDL         : Input pressure strain model to set constants for
C IEPSIJ         : Epsilon_ij model
C IEPSEQ         : Epsilon equation model
C IDAMP          : Damping used for mu_turbulent
C ISRCW          : Control of additional wall normal dependent source routines
C
      include '../header/common.h'
C
      IF (IPSMDL .EQ. IPLRR1) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 1.0E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPLRR2) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPSSG ) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPFLT ) THEN
         IEPSIJ = IEFLT
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPLS  ) THEN
         IEPSIJ = IELAIS
         IEPSEQ = IEEZSG
         IDAMP  = IFMZSG
         ISRCW  = IYES
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.01363540E0
         PSALFA = 0.45E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPSO95) THEN
         IEPSIJ = IESO95
         IEPSEQ = IEES95
         IDAMP  = IFMZSG
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.0E0
         PSALFA = - 0.36E0
         PSBETA =   0.0005E0
         PSGAMA =   0.075E0
      ELSE IF (IPSMDL .EQ. IPGL  ) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.6E0
         PSC2   = 0.6E0
         PSCW   = 0.E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE
         WRITE (IOUT,1000) IPSMDL
      ENDIF
C
C Finished
C
      RETURN
 1000 FORMAT (' ','SETPIJ: ERROR-> Invalid pressure-strain model: ',
     1            'IPSMDL = ',I5)
      END






      SUBROUTINE RDGRID (NBLKS, IDIM, JDIM, KDIM, IOFFS, GRDFIL, R, 
     1                   IDIMIN, JDIMIN, KDIMIN, IGRDTP, GRDBIN, 
     2                   IEOF, IERRCD)
C
C Routine to read in the complete grid data.  This routine will either read
C the grid for each block from a separate file (IGRDTP = IGTSEP) or
C will read the entire grid from a single PLOT3D file in single or multiblock
C format as appropriate (IGRDTP = IGTP3D).
C
C This routine is organized as follows:
C     1. Read Grid
C        a. Read Separate File (IGRDTP = IGTSEP) [formatted or binary]
C           For each block {
C           i.   Open file
C           ii.  Read header information
C           iii. Read block grid 
C           iv.  Close file
C           }
C        b. Read PLOT3D File   (IGRDTP = IGTP3D) [formatted or binary]
C           i.   Open file
C           ii.  Read header information for all blocks
C           iii. Verify dimensions
C           iv.  Loop over all blocks and read grids
C           v.   Close file
C        c. Error: Invalid type
C     2. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C NBLKS          : Number of blocks
C IDIM,JDIM,KDIM : Dimensions of all blocks on finest level
C IOFFS          : Offset pointer to each block
C GRDFIL         : Grid file name(s)
C R              : Position vector of grid vertices (x,y,z)
C IDIMIN,JDIMIN,
C KDIMIN         : Dimensions of all blocks read in the grid
C IGRDTP         : Grid type (separate files or Plot3D multiblock)
C GRDBIN         : Logical for binary (.TRUE.) or formatted (.FALSE.) file
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
C
      DIMENSION IDIM  (NBLKS),
     1          JDIM  (NBLKS),
     2          KDIM  (NBLKS),
     3          IOFFS (NBLKS)
C
      DIMENSION IDIMIN(NBLKS),
     1          JDIMIN(NBLKS),
     2          KDIMIN(NBLKS)
C
      CHARACTER*80 GRDFIL(NBLKS),
     1             FERROR
C
      DIMENSION R     (*)
C
      LOGICAL GRDBIN
C
      CHARACTER*80 FILTYP
C
C Initialize error handling
C
      IEOF   = INO
      IERRCD = INO
C
C Initialize file type
C
      IF (GRDBIN) THEN
         FILTYP = 'UNFORMATTED'
      ELSE
         FILTYP = 'FORMATTED'
      ENDIF
C
C 1.   Read Grid
C 1-a. Read Separate Files (IGRDTP = IGTSEP) [formatted or binary]
C
      IF (IGRDTP .EQ. IGTSEP) THEN
         DO 190 IBLK = 1, NBLKS
            FERROR = GRDFIL(IBLK)
            OPEN (IGRID,FILE=GRDFIL(IBLK),STATUS='OLD',FORM=FILTYP,
     1                  ERR=920)
            REWIND (IGRID)
            WRITE (IOUT,1000) IBLK, GRDFIL(IBLK)
            IR   = IOFFS(IBLK) * 3     + 1
C
C      Binary  files
C
            IF (GRDBIN) THEN
               IF (THREED) THEN
                  READ (IGRID,END=900,ERR=910) IDIM1, JDIM1, KDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK) .OR.
     2                KDIM1 .NE. KDIM(IBLK)) THEN
                     WRITE (IOUT,1100) IBLK, 
     1                                 IDIM(IBLK),JDIM(IBLK),KDIM(IBLK),
     2                                 IDIM1,     JDIM1,     KDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ELSE
                  READ (IGRID,END=900,ERR=910) IDIM1, JDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK)) THEN
                     WRITE (IOUT,1105) IBLK, IDIM(IBLK), JDIM(IBLK), 
     1                                       IDIM1,      JDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ENDIF
C
C      Formatted files
C
            ELSE
               IF (THREED) THEN
                  READ (IGRID,*,END=900,ERR=910) IDIM1, JDIM1, KDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK) .OR.
     2                KDIM1 .NE. KDIM(IBLK)) THEN
                     WRITE (IOUT,1100) IBLK, 
     1                                 IDIM(IBLK),JDIM(IBLK),KDIM(IBLK),
     2                                 IDIM1,     JDIM1,     KDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ELSE
                  READ (IGRID,*,END=900,ERR=910) IDIM1, JDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK)) THEN
                     WRITE (IOUT,1105) IBLK, IDIM(IBLK), JDIM(IBLK), 
     1                                       IDIM1,      JDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ENDIF
            ENDIF
            CLOSE (IGRID)
            IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
  190    CONTINUE
C
C 1-b. Read PLOT3D File   (IGRDTP = IGTP3D) [formatted or binary]
C 1-b-i.   Open file
C
      ELSE IF (IGRDTP .EQ. IGTP3D) THEN
         FERROR = GRDFIL(1)
         OPEN (IGRID,FILE=GRDFIL(1),STATUS='OLD',FORM=FILTYP,
     1               ERR=920)
         REWIND (IGRID)
C
C 1-b-ii.  Read grid header file
C
C      Binary
C
         IF (GRDBIN) THEN
C
C      Read multiblock data if more than 1 block
C
            IF (NBLKS .GT. 1) THEN
               READ (IGRID,END=900,ERR=910) NBLKIN
               IF (NBLKIN .NE. NBLKS) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1110) NBLKS, NBLKIN
                  RETURN
               ENDIF
            ENDIF
C
C      3D or 2D
C
            IF (THREED) THEN
               READ (IGRID,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                       KDIMIN(L), L=1,NBLKS)
            ELSE
               READ (IGRID,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                       L=1,NBLKS)
            ENDIF
C
C      Formatted
C
         ELSE
C
C      Read multiblock data if more than 1 block
C
            IF (NBLKS .GT. 1) THEN
               READ (IGRID,*,END=900,ERR=910) NBLKIN
               IF (NBLKIN .NE. NBLKS) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1110) NBLKS, NBLKIN
                  RETURN
               ENDIF
            ENDIF
C
C      3D or 2D
C
            IF (THREED) THEN
               READ (IGRID,*,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                         KDIMIN(L), L=1,NBLKS)
            ELSE
               READ (IGRID,*,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                         L=1,NBLKS)
            ENDIF
         ENDIF
C
C 1-b-iii. Verify dimensions
C
         IF (THREED) THEN
            DO 210 IBLK = 1, NBLKS
               IF (IDIM(IBLK) .NE. IDIMIN(IBLK) .OR.
     1             JDIM(IBLK) .NE. JDIMIN(IBLK) .OR.
     2             KDIM(IBLK) .NE. KDIMIN(IBLK) ) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1100) IBLK, 
     1                           IDIM  (IBLK),JDIM  (IBLK),KDIM  (IBLK),
     2                           IDIMIN(IBLK),JDIMIN(IBLK),KDIMIN(IBLK)
               ENDIF
  210       CONTINUE
         ELSE
            DO 220 IBLK = 1, NBLKS
               IF (IDIM(IBLK) .NE. IDIMIN(IBLK) .OR.
     1             JDIM(IBLK) .NE. JDIMIN(IBLK) ) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1105) IBLK, 
     1                              IDIM  (IBLK),JDIM  (IBLK),
     2                              IDIMIN(IBLK),JDIMIN(IBLK)
               ENDIF
  220       CONTINUE
         ENDIF
         IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
C
C 1-b-iv.  Loop over all blocks reading grid
C
         DO 290 IBLK = 1, NBLKS
            WRITE (IOUT,1000) IBLK, GRDFIL(1)
            IR   = IOFFS(IBLK) * 3     + 1
            CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                   KDIM(IBLK), R(IR), IEOF, IERRCD)
            IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
  290    CONTINUE
C
C 1-b-v.   Close file
C
         CLOSE (IGRID)
C
C 1-c. Error
C
      ELSE
         IERRCD = IYES
         WRITE (IOUT,1120) IGRDTP
         RETURN
      ENDIF
C
C     Succesful completion reading grid - return
C
      RETURN
C
C 2.   Error handling section
C 2-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      WRITE (IOUT,1200) 
      RETURN
C
C 2-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1210) 
      RETURN
C
C 2-c. Open Error Handler
C
  920 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1220) FERROR
      RETURN
C
 1000 FORMAT (' ','RDGRID: Reading grid for block ',I3,' from file ',
     1            A80)
 1100 FORMAT (' ','RDGRID: ERROR-> Reading input 3D grid for Block ',I3,
     1       /' ','                     IDIM, JDIM, KDIM = ',3I6,
     2       /' ','                Grid IDIM, JDIM, KDIM = ',3I6,
     3       /' ', 8X,'RUN ABORTING!')
 1105 FORMAT (' ','RDGRID: ERROR-> Reading input 2D grid for Block ',I3,
     1       /' ','                     IDIM, JDIM = ',2I6,
     2       /' ','                Grid IDIM, JDIM = ',2I6,
     3       /' ', 8X,'RUN ABORTING!')
 1110 FORMAT (' ','RDGRID: ERROR-> Reading input grid number of blocks',
     1       /' ','                     NBLKS = ',I6,
     2       /' ','                Grid NBLKS = ',I6,
     3       /' ', 8X,'RUN ABORTING!')
 1120 FORMAT (' ','RDGRID: ERROR-> Invalid grid type specified ',
     1            'IGRDTP = ',I6,
     2       /' ', 8X,'RUN ABORTING!')
C
 1200 FORMAT (' ','RDGRID: ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1210 FORMAT (' ','RDGRID: ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1220 FORMAT (' ','RDGRID: ERROR-> Error opening file ',A80,
     1       /' ', 8X,'RUN ABORTING!')
C
      END






      SUBROUTINE RDP3D  (GRDBIN, IDIM, JDIM, KDIM, R, IEOF, IERRCD)
C
C Routine to read in one block of the grid in PLOT3D format (either formatted
C or binary).  The loops run over kdim ( jdim ( idim ) ) ).  
C
C This routine is organized as follows:
C     1. Read grid
C        a. Binary    (3D or 2D grid)
C        b. Formatted (3D or 2D grid)
C     2. If two-dimensional grid then copy input plane to plane 2
C     3. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C GRDBIN         : Logical for binary (.TRUE.) or formatted (.FALSE.) file
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector of grid vertices (x,y,z)
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3)
C
      LOGICAL GRDBIN
C
C 1.   Read grid
C 1-a. Binary   PLOT3D grid
C
      IF (GRDBIN) THEN
         IF (THREED) THEN
            READ (IGRID,END=900,ERR=910) 
     1              ((((R(J,K,I,L),I=1,IDIM),J=1,JDIM),K=1,KDIM),L=1,3)
         ELSE
            READ (IGRID,END=900,ERR=910) 
     1               (((R(J,1,I,L),I=1,IDIM),J=1,JDIM),L=1,2)
         ENDIF
C
C 1-b. Formatted PLOT3D grid
C
      ELSE
         IF (THREED) THEN
            READ (IGRID,*,END=900,ERR=910) 
     1              ((((R(J,K,I,L),I=1,IDIM),J=1,JDIM),K=1,KDIM),L=1,3)
         ELSE
            READ (IGRID,*,END=900,ERR=910) 
     1               (((R(J,1,I,L),I=1,IDIM),J=1,JDIM),L=1,2)
         ENDIF
      ENDIF
C
C 2.   Axisymmetric or Two-Dimensional Grid Modifications
C 2-a. Axisymmetric: Rotate K=1, K=2 grid planes by AXIDTH degrees
C
      IF (AXISYM) THEN
         PI     = ACOS(-1.E0)
         STHETA = SIN (AXIDTH * PI / 180.E0)
         DO 220 I = 1, IDIM
            DO 210 J = 1, JDIM
               Z          =   R(J,1,I,2) * STHETA
               R(J,1,I,3) = - Z
               R(J,2,I,3) =   Z
               R(J,2,I,1) =   R(J,1,I,1)
               R(J,2,I,2) =   R(J,1,I,2)
 210        CONTINUE
 220     CONTINUE
C
C 2-b. Two-Dimensional: Copy input plane to plane 2 
C
      ELSE IF (.NOT. THREED) THEN
         DO 240 I = 1, IDIM
            DO 230 J = 1, JDIM
               R(J,1,I,3) = 0.E0
               R(J,2,I,3) = 1.E0
               R(J,2,I,1) = R(J,1,I,1)
               R(J,2,I,2) = R(J,1,I,2)
 230        CONTINUE
 240     CONTINUE
      ENDIF
C
C     Succesful completion reading block grid - return
C
      RETURN
C
C 3.   Error handling section
C 3-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      WRITE (IOUT,1000) 
      RETURN
C
C 3-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1005) 
      RETURN
C
 1000 FORMAT (' ','RDP3D : ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1005 FORMAT (' ','RDP3D : ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
C
      END






      SUBROUTINE RDPRFL (NPTS, QBC, IBCDAT, IERRCD)
C
C Routine to read in a profile and store it into the QBC array to
C be used in the boundary condition routine BCPROF.
C
C NPTS           : Number of points in the QBC array for current segment
C QBC            : Primitive variables read from a file
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT( 1,ibc) = Boundary Condition Type
C                          2      = Direction (1=i, 2=j, 3=k)
C                          3      = i(1 or idim)
C                          4      = istrt
C                          5      = iend
C                          6      = jstrt
C                          7      = jend
C                          8      = kstrt
C                          9      = kend
C                         10      = profile identifier
C IERRCD         : Error Code
C                  IERRCD = INO  if NO storage error or error reading file
C                  IERRCD = IYES if    storage error or error reading file
C
      include '../header/common.h'
C
      DIMENSION QBC   (NPTS,NQ)
C
      DIMENSION IBCDAT(10)
C
C     Store boundary condition data in mnemonic names
C
      IBCTYP = IBCDAT(1)
      IDIR   = IBCDAT(2)
      INROUT = IBCDAT(3)
      ISTRT  = IBCDAT(4)
      IEND   = IBCDAT(5)
      JSTRT  = IBCDAT(6)
      JEND   = IBCDAT(7)
      KSTRT  = IBCDAT(8)
      KEND   = IBCDAT(9)
C
C     I Direction
C
      IF (IDIR .EQ. 1) THEN
         WRITE (IOUT,1010) 'K', KSTRT, KEND
         WRITE (IOUT,1010) 'J', JSTRT, JEND
         DO 120 K = KSTRT, KEND
            DO 110 J = JSTRT, JEND
               JKBCPT = (J-JSTRT+1) + (K-KSTRT)*(JEND-JSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(JKBCPT,L),L=1,NQ)
  110       CONTINUE
  120    CONTINUE
C
C     J Direction
C
      ELSE IF (IDIR .EQ. 2) THEN
         WRITE (IOUT,1010) 'I', ISTRT, IEND
         WRITE (IOUT,1010) 'K', KSTRT, KEND
         DO 220 I = ISTRT, IEND
            DO 210 K = KSTRT, KEND
               IKBCPT = (K-KSTRT+1) + (I-ISTRT)*(KEND-KSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(IKBCPT,L),L=1,NQ)
  210       CONTINUE
  220    CONTINUE
C
C     K Direction
C
      ELSE IF (IDIR .EQ. 3) THEN
         WRITE (IOUT,1010) 'I', ISTRT, IEND
         WRITE (IOUT,1010) 'J', JSTRT, JEND
         DO 320 I = ISTRT, IEND
            DO 310 J = JSTRT, JEND
               IJBCPT = (J-JSTRT+1) + (I-ISTRT)*(JEND-JSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(IJBCPT,L),L=1,NQ)
  310       CONTINUE
  320    CONTINUE
      ELSE
         IERRCD = IYES
         WRITE (IOUT,1200) IDIR
         RETURN
      ENDIF
C
C     Finished reading profile boundary data
C
      RETURN
C
C 2.   Error handling section
C 2-a. End-Of-File Handler
C
  900 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1205)
      RETURN
C
C 2-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1210)
      RETURN
C
 1010 FORMAT (' ','RDPRFL: Looping over ',A1,' from ',I4,' to ',I4)
 1200 FORMAT (' ','RDPRFL: ERROR-> Invalid direction specified ',I4,
     1       /' ', 8X,'RUN ABORTING!')
 1205 FORMAT (' ','RDPRFL: ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1210 FORMAT (' ','RDPRFL: ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
C
      END






      SUBROUTINE RDREST (ITURB, IBLOCK, ILVL, IDIM, JDIM, KDIM, 
     1                   Q, QN, PROPS, IEOF,  IERRCD)
C
C Routine to read in the restart file for a block.
C
C Note: To read a restart file from ISAAC versions earlier than 3.5.15,
C       set RDPROP = FALSE.
C
C N.B.  QN is converted to PRIMITIVE variables for restart file.
C
C This routine is organized as follows:
C     If (IBLOCK = 1 AND ILVL = 1) 
C     {
C        1. Read Residual History Data 
C           a. Read input residual history data
C           b. Copy residual history data to output restart file
C           c. Echo residual history to listing
C        2. Read Restart File Header Information 
C           a. Number of variables
C           b. Number of properties
C           c. Number of additional time levels and time for unsteady
C     }
C     3. Read Flow Variables and Properties
C        a. Read Flow Variables 
C        b. Read Properties
C        c. Read flow variables at additional time levels for unsteady flow
C     4. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C ITURB          : Turbulence model
C IBLOCK         : Block number
C ILVL           : Grid level
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTIME-1),
     2          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Control of whether or not to read properties data and additional time levels
C
      LOGICAL   RDPROP, RDUNST
C
C Initialize error switches
C
      IEOF   = INO
      IERRCD = INO
C
      IF (ITURB  .GE. ITLMNR) THEN
         RDPROP = .TRUE.
      ELSE
         RDPROP = .FALSE.
      ENDIF
C
      IF (ITIMED .EQ. IUNSTD .AND. .NOT.(INITTD)) THEN
         RDUNST = .TRUE.
      ELSE
         RDUNST = .FALSE.
      ENDIF
C
C 1.  Read residual history and file header only for the first block
C
      IF (IBLOCK .EQ. 1 .AND. ILVL .EQ. 1) THEN
C
C     Initial residual output increment
C
         INC   = 1
C
C     Initialize variables in case restart with no residual history data
C
         ITOLD = 0
         R2OLD = 0.E0
         R2ONE = 1.E0
C
C 1-a. Read Residual History Data
C 1-b. Copy Residual History Data to output restart file
C
  100    CONTINUE
         READ (IRDRST,END=900,ERR=910)  IT, R2T
         IF (IT .LE. 0) GO TO 190
         ITOLD  = IT
         R2OLD  = R2T
         IF (IT .LT. MXHIST) THEN
            R2(IT) = R2T
         ELSE
            IERRCD = IYES
         ENDIF
         IF (IT/INC*INC .EQ. IT) THEN
            WRITE (IOUT,1000) IT, R2T
         ENDIF
         IF (IT .GT. 10)  INC = 10
         IF (IT .GT. 100) INC = 100
         GO TO 100
C
C 1-c. Echo the last residual to listing
C
  190    CONTINUE
         ITTOT = ITOLD 
         WRITE (IOUT,1000) ITOLD, R2OLD 
C
         IF (IERRCD .NE. INO) THEN
            WRITE (IOUT,1010) ITTOT, MXHIST
            RETURN
         ENDIF
C
         R2ONE = R2(1)
C
C 2.   Read Restart File Header Information
C 2-a. NQ
C
         READ (IRDRST,END=900,ERR=910) NQIN
C
         IF (NQIN .NE. NQ) THEN
            WRITE (IOUT,1100) NQIN, NQ
            STOP
         ENDIF
C
C 2-b. NP (not if inviscid)
C
         IF (RDPROP) THEN
            READ (IRDRST,END=900,ERR=910) NPIN
C
            IF (NPIN .NE. NP) THEN
               WRITE (IOUT,1105) NPIN, NP
               STOP
            ENDIF
         ENDIF
C
C 2-c. Number of additional time levels and time for unsteady
C
         IF (RDUNST) THEN
            READ (IRDRST,END=900,ERR=910) NTIMIN
            READ (IRDRST,END=900,ERR=910) TIME
C
            IF (NTIMIN .NE. (NTIME-1)) THEN
               WRITE (IOUT,1110) NTIMIN, NTIME-1
               STOP
            ENDIF
         ENDIF
C
C     End of history and file header section
C
      ENDIF
C
C 3.   Read Flow Variables and Properties
C 3-a. Flow Variables (Q)
C
      READ (IRDRST,END=900,ERR=910) Q
C
C 3-b. Properites (PROPS)
C
      IF (RDPROP) THEN
         READ (IRDRST,END=900,ERR=910) PROPS
      ENDIF
C
C 3-c. Read flow variables at additional time levels for unsteady flow
C
      IF (RDUNST) THEN
         READ (IRDRST,END=900,ERR=910) QN
      ENDIF
C
C     Succesful completion reading restart - return
C
      RETURN
C
C 4.   Error handling section
C 4-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      RETURN
C
C 4-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1900)
      RETURN
C
 1000 FORMAT (' ','   ITER = ',I5,'   R2 = ',E15.8)
 1010 FORMAT (' ','RDREST: ERROR-> Length of residual history in ',
     1            'restart file (ITTOT) ',I6,
     2       /' ',15X,' exceeds maximum allocated      (MXHIST ',
     3            'in hstry.h) ',I6,
     4       /' ',15X,' Increase MXHIST to be greater than ',
     5                'ITTOT + NITS and recompile.',
     6       /' ', 8X,'RUN ABORTING!')
 1100 FORMAT (' ','RDREST: ERROR-> Number of variables  (NQ) in ',
     1            'restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'Check input for incorrect model and resubmit.',
     4       /' ', 8X,'RUN ABORTING!')
 1105 FORMAT (' ','RDREST: ERROR-> Number of properties (NP) in ',
     1            'restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'Check input for incorrect model and resubmit.',
     4       /' ', 8X,'RUN ABORTING!')
 1110 FORMAT (' ','RDREST: ERROR-> Number of additional time levels ',
     1            '(NTIME-1) in restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'This may be caused by attempting to read a ',
     4                'steady restart file', 
     5       /' ', 8X,'without specifying the INITIALIZE TIME ',
     6                'DEPENDENT input option to ',
     7       /' ', 8X,'convert to an unsteady solution.',
     7      //' ', 8X,'Check input for incorrect scheme and resubmit.',
     8       /' ', 8X,'RUN ABORTING!')
 1900 FORMAT (' ','RDREST: ERROR-> Unclassifiable error reading ',
     1            'restart file.',
     2       /' ', 8X,'Check specified restart file and retry.',
     3       /' ', 8X,'RUN ABORTING!')
      END







      SUBROUTINE SFHT  (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                  Q, PROPS, S, VOL, IBCDAT, RC, IBLOCK)
C
C Routine to output the skin friction, heat transfer, and boundary layer
C quantities.
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C RC             : Cell center grid
C IBLOCK         : Block number to identify output
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C Set control for higher order evaluation of du/dn and dT/dn
C     IHIORD = 0 for first  order
C     IHIORD = 1 for second order
C
      IHIORD = 1
C
C Output headers
C
      IF (IHIORD .EQ. 0) THEN
         WRITE (IOUT,1000) IBLOCK, '1st'
      ELSE
         WRITE (IOUT,1000) IBLOCK, '2nd'
      ENDIF
      WRITE (IOUT,1005) 
C
C 1. Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C Store boundary condition data in mnemonic names
C
      IBCTYP   = IBCDAT(1)
      IDIR     = IBCDAT(2)
      INROUT   = IBCDAT(3)
      ISTRT(1) = MAX (IBCDAT(4), IPSTRT)
      IEND (1) = MIN (IBCDAT(5), IPEND)
      ISTRT(2) = IBCDAT(6)
      IEND (2) = IBCDAT(7)
      ISTRT(3) = IBCDAT(8)
      IEND (3) = IBCDAT(9)
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCN2  = 3
         IBCM   = 1
         IBCEND = IJKDIM(IDIR) + 1
         IBCINC = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
         IBCEND = 1
         IBCINC = - 1
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Loop over all points on the boundary surface and calculate skin friction,
C    heat transfer, and boundary layer quantities for output.
C
      DO 220 I2 = ISTRT(IC2), IEND(IC2)
         DO 210 I3 = ISTRT(IC3), IEND(IC3)
C
            IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IN2  = IBCN2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN2  = IBCN2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN2  = IBCN2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            RHOW = Q(JBC,KBC,IBC,1) 
            RHON = Q(JN ,KN ,IN ,1) 
            RHON2= Q(JN2,KN2,IN2,1) 
C
            PW   = Q(JBC,KBC,IBC,5) 
            PN   = Q(JN ,KN ,IN ,5) 
            PN2  = Q(JN2,KN2,IN2,5) 
C
C Calculate T = T / TINF  at the wall and cell center off the wall
C
            TBC  = GAMMA * PW  / RHOW
            TN   = GAMMA * PN  / RHON
            TN2  = GAMMA * PN2 / RHON2
C
C Calculate the velocity parallel to the wall at the two cell centers off 
C the wall (to allow for second order accurate expressions).
C UBAR = U(parallel) = U - UDOTN where U is the total velocity vector
C and UDOTN is the velocity at the cell center normal to the wall.
C (Velocity on the wall should be identically zero)
C
            UDOTN  = S(JM,KM,IM,1,IDIR) * Q(JN ,KN ,IN ,2)
     1             + S(JM,KM,IM,2,IDIR) * Q(JN ,KN ,IN ,3)
     2             + S(JM,KM,IM,3,IDIR) * Q(JN ,KN ,IN ,4)
            UBARN  = SQRT (
     1               (Q(JN ,KN ,IN ,2) - UDOTN *S(JM,KM,IM,1,IDIR))**2
     2             + (Q(JN ,KN ,IN ,3) - UDOTN *S(JM,KM,IM,2,IDIR))**2
     3             + (Q(JN ,KN ,IN ,4) - UDOTN *S(JM,KM,IM,3,IDIR))**2 )
C
            UDOTN2 = S(JM,KM,IM,1,IDIR) * Q(JN2,KN2,IN2,2)
     1             + S(JM,KM,IM,2,IDIR) * Q(JN2,KN2,IN2,3)
     2             + S(JM,KM,IM,3,IDIR) * Q(JN2,KN2,IN2,4)
            UBARN2 = SQRT (
     1               (Q(JN2,KN2,IN2,2) - UDOTN2*S(JM,KM,IM,1,IDIR))**2
     2             + (Q(JN2,KN2,IN2,3) - UDOTN2*S(JM,KM,IM,2,IDIR))**2
     3             + (Q(JN2,KN2,IN2,4) - UDOTN2*S(JM,KM,IM,3,IDIR))**2 )
C
C Viscosity (RMUW) at the wall and the cell center
C
            RMUW   = PROPS(JBC,KBC,IBC,1)
            RMUN   = PROPS(JN, KN ,IN ,1)
C
C Distances
C
            DN     = SQRT ( (RC(JN ,KN ,IN ,1) - RC(JBC,KBC,IBC,1))**2
     1                    + (RC(JN ,KN ,IN ,2) - RC(JBC,KBC,IBC,2))**2
     2                    + (RC(JN ,KN ,IN ,3) - RC(JBC,KBC,IBC,3))**2 )
C
            DN2    = SQRT ( (RC(JN2,KN2,IN2,1) - RC(JN ,KN ,IN ,1))**2
     1                    + (RC(JN2,KN2,IN2,2) - RC(JN ,KN ,IN ,2))**2
     2                    + (RC(JN2,KN2,IN2,3) - RC(JN ,KN ,IN ,3))**2 )
C
C Calculate first derivative terms using first or second order differencing
C
            DUDN   = (1-IHIORD) * (UBARN      / DN)
     1               +  IHIORD  * (- DN*DN*UBARN2 + (DN+DN2)**2*UBARN)
     2                            / (DN*DN2*(DN+DN2))
            DTDN   = (1-IHIORD) * ((TN - TBC) / DN)
     1               +  IHIORD  * (- DN*DN*TN2 + (DN+DN2)**2*TN 
     2                             - DN2*(2.E0*DN+DN2)*TBC)
     3                            / (DN*DN2*(DN+DN2))
C
C Pressure Coefficients Pw/Qinf =  P_w          / ( 1/2 * rhoinf * uinf ** 2 )
C                       Cp      = (P_w - P_inf) / ( 1/2 * rhoinf * uinf ** 2 )
C
            PWQINF =  PW         * 2.E0 / ( FSMACH * FSMACH )
            CP     = (PW - PREF) * 2.E0 / ( FSMACH * FSMACH )
C
C Find the maximum velocity magnitude in the profile
C
            UMAX   = 0.E0
CDIR$ IVDEP
            DO 100 I1 = IBCPT, IBCEND, IBCINC
C
               II   = I1    *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JJ   = I1    *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KK   = I1    *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               UMAG = SQRT ( Q(JJ,KK,II,2) * Q(JJ,KK,II,2) +
     1                       Q(JJ,KK,II,3) * Q(JJ,KK,II,3) +
     2                       Q(JJ,KK,II,4) * Q(JJ,KK,II,4) )
               UMAX = MAX  ( UMAX, UMAG )
  100       CONTINUE
C
C Calculate the edge velocity.  Based on Chris Rumsey's notes with the
C following assumptions: 
C      1. H   = H_inf 
C      2. P   = P_wall 
C      3. Isentropic from edge to wall to get density
C      4. TKE = TKE_inf
C If these assumptions result in imaginary UE, then set UE to UMAX
C
            HEDGE  = EINF + PREF / RHOINF
            PEDGE  = PW
            RHOEDG = RHOINF * (PEDGE / PREF)**(1.E0/GAMMA)
            TKEDGE = TKEINF
C
            ARGUE  = 2.E0*(HEDGE - GAMMA/GAMM1*PEDGE/RHOEDG - TKEDGE)
            IF (ARGUE .LE. 0.E0) THEN
               UE  = UMAX
            ELSE
               UE  = SQRT (ARGUE)
            ENDIF
C
C Make sure that UE is bounded by UMAX
C
            UE     = MIN (UE, UMAX)
C
C Skin Friction Cf      = Tau_wall / ( 1/2 * rho_inf  * u_inf  ** 2 )
C               Cf_edge = Tau_wall / ( 1/2 * rho_edge * u_edge ** 2 )
C
            IF ( IFWALF ) THEN
               TAUW   = PROPS (JBC2,KBC2,IBC2,ITQTAU)
            ELSE
               TAUW   = RMUW * DUDN
            ENDIF
            CF     = 2.E0 * TAUW / ( RE * FSMACH )
            CFEDGE = CF * (RHOINF/RHOEDG) * (FSMACH/UE)**2
C
C Heat Transfer Ch      = Qdot_wall / ( rho_inf  * u_inf  * (H_inf  - H_wall) )
C               Ch_edge = Qdot_wall / ( rho_edge * u_edge * (H_edge - H_wall) )
C
            IF ( IFWALF ) THEN
               QW     = PROPS (JBC2,KBC2,IBC2,ITQDOT)
            ELSE
               QW     = RMUW * DTDN / PR
            ENDIF
            DENOM = RE*(1.E0-TWALL/TINF+0.5E0*GAMM1*FSMACH*FSMACH)
C
C Check to see if can be non-dimensionalized this way
C
            IF (DENOM .EQ. 0.E0) THEN
               CH     = RMAX
               CHEDGE = RMAX
            ELSE
               CH     = QW / DENOM
               CHEDGE = CH * (RHOINF/RHOEDG) * (FSMACH/UE)**2
            ENDIF
C
C Ustar / ainf = sqrt (Tauwall / rho) * sqrt (Minf / Re)
C
            USTAR = SQRT (ABS(TAUW) / RHOW) * SQRT (FSMACH / RE)
C
C Yplus = y ustar / nu : ustar = sqrt (Tauwall / rho) 
C
            YPLUS = PROPS(JN ,KN ,IN ,3)
C
C Calculate boundary layer data:
C 1. Find the boundary layer edge defined as:
C    i.   location where U = .995 U_edge
C 2. Integrate using trapezoidal rule to find:
C    i.   boundary layer thickness (DELTA)
C    ii.  displacement thickness   (DSTAR)
C    iii. momentum thickness       (THETA)
C Formulas for DELTA, DSTAR, THETA from White, 'Viscous Fluid Flow', p.606
C
C JEDGE is the last cell-center before the boundary layer edge
C
            UEDGE  = .995E0 * UE
            DELTA  = 0.E0
            DSTAR  = 0.E0
            THETA  = 0.E0
CDIR$ IVDEP
            DO 110 I1 = IBCPT, IBCEND, IBCINC
C
               II    = I1    *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JJ    = I1    *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KK    = I1    *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IP    = II + IBCINC*IKD(IC1,1)
               JP    = JJ + IBCINC*IKD(IC1,2)
               KP    = KK + IBCINC*IKD(IC1,3)
C
               U     = SQRT ( Q(JJ,KK,II,2) * Q(JJ,KK,II,2) + 
     1                        Q(JJ,KK,II,3) * Q(JJ,KK,II,3) + 
     2                        Q(JJ,KK,II,4) * Q(JJ,KK,II,4) ) 
               UP    = SQRT ( Q(JP,KP,IP,2) * Q(JP,KP,IP,2) + 
     1                        Q(JP,KP,IP,3) * Q(JP,KP,IP,3) + 
     2                        Q(JP,KP,IP,4) * Q(JP,KP,IP,4) ) 
               IF (U  .LT. UEDGE .AND. UP .GE. UEDGE) THEN
                  I1EDGE = I1
                  GO TO 119
               ENDIF
C
               DY    = SQRT ( (RC(JP,KP,IP,1)-RC(JJ,KK,II,1))**2
     1                      + (RC(JP,KP,IP,2)-RC(JJ,KK,II,2))**2
     2                      + (RC(JP,KP,IP,3)-RC(JJ,KK,II,3))**2 )
C
               U     = U  / FSMACH
               UP    = UP / FSMACH
               RHOU  = Q(JJ,KK,II,1) * U
               RHOUP = Q(JP,KP,IP,1) * UP
C
               DELTA = DELTA + DY
               DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                 + (1.E0 - RHOUP) )*DY
               THETA = THETA + 0.5E0 * ( (1.E0 - U) * RHOU
     1                                 + (1.E0 - UP)* RHOUP )*DY
  110       CONTINUE
C
C Could not find bl edge.  Output error.
C
            DELTA  = - RMAX
            DSTAR  = - RMAX
            THETA  = - RMAX
C
  119       CONTINUE
C
C Output Cf, Ch, Cp, Delta, Delta*, Theta
C 
  209       CONTINUE
            WRITE (IOUT,1010) IBC, JBC, KBC, RC(JBC,KBC,IBC,1), 
     1                        RC(JBC,KBC,IBC,2), RC(JBC,KBC,IBC,3), 
     2                        PWQINF, CP, CF, CFEDGE, CH, CHEDGE, DN, 
     3                        YPLUS, USTAR, DELTA, DSTAR, THETA, UE
  210    CONTINUE
  220 CONTINUE
C
C Finished with Skin Friction and Heat Transfer
C
      RETURN
 1000 FORMAT (//' ',' Surface Data for BLOCK ',I5,10X,'(',A,
     1              ' Order accurate derivatives)')
 1005 FORMAT ( /' ',3X,'I',3X,'J',3X,'K',6X,'X',12X,'Y',12X,'Z',9X,
     1              'Pw/Qinf',8X,'Cp',11X,'Cf',9X,'Cf_edge',8X, 
     2              'Ch',9X,'Ch_edge',8X,'DN',10X,'Yplus',8X,'Ustar',8X,
     3              'Delta',8X,'Delta*',7X,'Theta',8X,'U_edge'/)
 1010 FORMAT (' ',3I4,16(1PE13.5))
      END
      SUBROUTINE xSFHTJ (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                  Q, PROPS, S, VOL, IBCDAT, RC, IBLOCK)
C
C Routine to calculate the skin friction and heat transfer using 
C the thin layer Navier-Stokes assumption for the J direction.
C
C This routine is not presently suited for the multiblock
C approach.  It is close to being ready but is not.  The problem
C lies that the metrics for one cell into the adjacent block
C are needed to be completely general.  These are currently not available.
C The routine metbc needs to be completed.
C
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C RC             : Cell center grid
C IBLOCK         : Block number to identify output
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
C Store boundary condition data in mnemonic names
C
      IBCTYP = IBCDAT(1)
      IDIR   = IBCDAT(2)
      INROUT = IBCDAT(3)
      ISTRT  = MAX (IBCDAT(4), IPSTRT)
      IEND   = MIN (IBCDAT(5), IPEND)
      JSTRT  = IBCDAT(6)
      JEND   = IBCDAT(7)
      KSTRT  = IBCDAT(8)
      KEND   = IBCDAT(9)
c
      WRITE (IOUT,1000) IBLOCK
      WRITE (IOUT,1005) 
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBC  is the boundary condition location
C      IBC2 is the boundary condition location outside the block from IBC
C      IN   is the point to the inside of the block from IBC
C      IM   is the location of the metrics
C
C
C J Direction
C
      IF (IDIR .EQ. 2) THEN
         IF (INROUT .EQ. 1) THEN
            JBC  = 1
            JBC2 = 0
            JN   = 2
            JM   = 1
            JMP1 = 2
            JBEG = 1
            JEND = JDIM + 1
            JINC = 1
         ELSE
            JBC  = JDIM + 1
            JBC2 = JDIM + 2
            JN   = JDIM
            JM   = JDIM
            JMP1 = JDIM - 1
            JBEG = JDIM + 1
            JEND = 1
            JINC = -1
         ENDIF
C
C Loop over the J = JBC plane
C
         DO 220 I = ISTRT, IEND
            DO 210 K = KSTRT, KEND
C
C Calculate T = T / TINF  at the wall and cell center off the wall
C
               TBC = GAMMA * Q(JBC,K,I,5) / Q(JBC,K,I,1)
               TN  = GAMMA * Q(JN,K,I,5) / Q(JN,K,I,1)
C
C Calculate the velocity parallel to the wall at 
C the cell center off the wall.
C U(parallel) = U - UBAR where U is the total velocity vector
C and UBAR is the velocity normal to the wall.
C (Velocity on the wall should be identically zero)
C
               UBARN = S(JM,K,I,1,2) * Q(JN,K,I,2)
     1               + S(JM,K,I,2,2) * Q(JN,K,I,3)
     2               + S(JM,K,I,3,2) * Q(JN,K,I,4)
               VELN  = SQRT ( (Q(JN,K,I,2) - UBARN*S(JM,K,I,1,2))**2
     1                      + (Q(JN,K,I,3) - UBARN*S(JM,K,I,2,2))**2
     2                      + (Q(JN,K,I,4) - UBARN*S(JM,K,I,3,2))**2 )
C
C Viscosity (RMUW) at the wall and the cell center
C
               RMUW = PROPS(JBC,K,I,1)
               RMUN = PROPS(JN, K,I,1)
C
C Average metrics and cell volume
C
               VOLAUX = 0.5E0 * VOL(JMP1,K,I)
               SPJ1   = 0.5E0 * ( S(JMP1,K,I,1,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,1,2)*S(JM,K,I,4,2) )
               SPJ2   = 0.5E0 * ( S(JMP1,K,I,2,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,2,2)*S(JM,K,I,4,2) )
               SPJ3   = 0.5E0 * ( S(JMP1,K,I,3,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,3,2)*S(JM,K,I,4,2) )
               SMJ1   = S(JM,K,I,1,2)*S(JM,K,I,4,2)
               SMJ2   = S(JM,K,I,2,2)*S(JM,K,I,4,2)
               SMJ3   = S(JM,K,I,3,2)*S(JM,K,I,4,2)
C
C Calculate first derivative terms using Divergence Theorem
C
               DVELDX = VELN*SPJ1 / VOLAUX
               DVELDY = VELN*SPJ2 / VOLAUX
               DVELDZ = VELN*SPJ3 / VOLAUX
               DTDX   = (TN*SPJ1 - TBC*SMJ1) / VOLAUX
               DTDY   = (TN*SPJ2 - TBC*SMJ2) / VOLAUX
               DTDZ   = (TN*SPJ3 - TBC*SMJ3) / VOLAUX
C
C Pressure Coefficients Pw/Qinf =  P_w          / ( 1/2 * rhoinf * uinf ** 2 )
C                       Cp      = (P_w - P_inf) / ( 1/2 * rhoinf * uinf ** 2 )
C
               PW     = Q(JBC,K,I,5) 
               PWQINF =  PW         * 2.E0 / ( FSMACH * FSMACH )
               CP     = (PW - PREF) * 2.E0 / ( FSMACH * FSMACH )
C
C Calculate the edge velocity.  Based on Chris Rumsey's notes with the
C following assumptions: 
C      1. H   = H_inf 
C      2. P   = P_wall 
C      3. Isentropic from edge to wall to get density
C      4. TKE = TKE_inf
C
               HEDGE  = EINF + PREF / RHOINF
               PEDGE  = PW
               RHOEDG = RHOINF * (PEDGE / PREF)**(1.E0/GAMMA)
               TKEDGE = TKEINF
C
               ARGUE  = 2.E0*(HEDGE - GAMMA/GAMM1*PEDGE/RHOEDG - TKEDGE)
               IF (ARGUE .LE. 0.E0) THEN
                  UE     = FSMACH
               ELSE
                  UE     = SQRT (ARGUE)
               ENDIF
C
C Skin Friction Cf      = Tau_wall / ( 1/2 * rho_inf  * u_inf  ** 2 )
C               Cf_edge = Tau_wall / ( 1/2 * rho_edge * u_edge ** 2 )
C
               DVELDN = DVELDX * S(JM,K,I,1,2) + DVELDY * S(JM,K,I,2,2)
     1                + DVELDZ * S(JM,K,I,3,2)
               IF ( IFWALF ) THEN
                  TAUW   = PROPS (JBC2,K,I,ITQTAU)
               ELSE
                  TAUW   = RMUW * DVELDN
               ENDIF
               CF     = 2.E0 * TAUW / ( RE * FSMACH )
               CFEDGE = CF * (RHOINF/RHOEDG) * (FSMACH/UE)**2
C
C Heat Transfer Ch      = Qdot_wall / ( rho_inf  * u_inf  * (H_inf  - H_wall) )
C               Ch_edge = Qdot_wall / ( rho_edge * u_edge * (H_edge - H_wall) )
C
               DTDN  = DTDX * S(JM,K,I,1,2) + DTDY * S(JM,K,I,2,2)
     1               + DTDZ * S(JM,K,I,3,2)
               IF ( IFWALF ) THEN
                  QW     = PROPS (JBC2,K,I,ITQDOT)
               ELSE
                  QW     = RMUW * DTDN / PR
               ENDIF
               DENOM = RE*(1.E0-TWALL/TINF+0.5E0*GAMM1*FSMACH*FSMACH)
C
C Check to see if can be non-dimensionalized this way
C
               IF (DENOM .EQ. 0.E0) THEN
                  CH     = RMAX
                  CHEDGE = RMAX
               ELSE
                  CH     = QW / DENOM
                  CHEDGE = CH * (RHOINF/RHOEDG) * (FSMACH/UE)**2
               ENDIF
C
C Ustar / ainf = sqrt (Tauwall / rho) * sqrt (Minf / Re)
C
               USTAR = SQRT (TAUW / Q(JBC,K,I,1)) *
     1                 SQRT (FSMACH / RE)
C
C Yplus = y ustar / nu : ustar = sqrt (Tauwall / rho) 
C
               DN    = VOLAUX / (0.5E0*(S(JMP1,K,I,4,2)+S(JM,K,I,4,2)))
               YPLUS = PROPS(JN,K,I,3)
C
C Find the boundary layer edge defined as:
C      1. Location where U = .995 U_edge
C JEDGE is the last cell-center before the boundary layer edge
C
               UE     = .995E0 * UE
               DO 110 J = JBEG, JEND-JINC, JINC
                  U  = SQRT ( Q(J     ,K,I,2)*Q(J     ,K,I,2) + 
     1                        Q(J     ,K,I,3)*Q(J     ,K,I,3) + 
     2                        Q(J     ,K,I,4)*Q(J     ,K,I,4) ) 
                  UP = SQRT ( Q(J+JINC,K,I,2)*Q(J+JINC,K,I,2) + 
     1                        Q(J+JINC,K,I,3)*Q(J+JINC,K,I,3) + 
     2                        Q(J+JINC,K,I,4)*Q(J+JINC,K,I,4) ) 
                  IF (U  .LT. UE .AND. UP .GE. UE) THEN
                     JEDGE = J
                     GO TO 119
                  ENDIF
  110          CONTINUE
C
C Could not find bl edge.  Output error.
C
               DELTA = - RMAX
               DSTAR = - RMAX
               THETA = - RMAX
               GO TO 209
C
C Integrate to find:
C    the boundary layer thickness (DELTA),
C    the displacement thickness (DSTAR),
C    the momentum thickness (THETA),
C using the trapezoidal rule.  Formulas for DELTA, DSTAR, THETA
C given on page 606 of White's 'Viscous Fluid Flow' for compressible fluids.
C
  119          CONTINUE
               DELTA = 0.E0
               DSTAR = 0.E0
               THETA = 0.E0
               DO 120 J = JBEG, JEDGE-JINC, JINC
                  DY    = SQRT ( (RC(J+JINC,K,I,1)-RC(J,K,I,1))**2
     1                         + (RC(J+JINC,K,I,2)-RC(J,K,I,2))**2
     2                         + (RC(J+JINC,K,I,3)-RC(J,K,I,3))**2 )
                  U     = SQRT (Q(J     ,K,I,2)*Q(J     ,K,I,2) + 
     1                          Q(J     ,K,I,3)*Q(J     ,K,I,3) + 
     2                          Q(J     ,K,I,4)*Q(J     ,K,I,4) )/FSMACH
                  RHOU  = Q(J,K,I,1) * U
                  UP    = SQRT (Q(J+JINC,K,I,2)*Q(J+JINC,K,I,2) + 
     1                          Q(J+JINC,K,I,3)*Q(J+JINC,K,I,3) + 
     2                          Q(J+JINC,K,I,4)*Q(J+JINC,K,I,4) )/FSMACH
                  RHOUP = Q(J+JINC,K,I,1) * UP
                  DELTA = DELTA + DY
                  DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                    + (1.E0 - RHOUP) )*DY
                  THETA = THETA + 0.5E0 * ( (1.E0 - U) * RHOU
     1                                    + (1.E0 - UP)* RHOUP )*DY
  120          CONTINUE
C
C Add contribution from JEDGE to DELTA to DELTA, DSTAR and THETA
C
               UJE   = SQRT ( Q(JEDGE     ,K,I,2)*Q(JEDGE     ,K,I,2) + 
     1                        Q(JEDGE     ,K,I,3)*Q(JEDGE     ,K,I,3) + 
     2                        Q(JEDGE     ,K,I,4)*Q(JEDGE     ,K,I,4) )
               UJEP  = SQRT ( Q(JEDGE+JINC,K,I,2)*Q(JEDGE+JINC,K,I,2) + 
     1                        Q(JEDGE+JINC,K,I,3)*Q(JEDGE+JINC,K,I,3) + 
     2                        Q(JEDGE+JINC,K,I,4)*Q(JEDGE+JINC,K,I,4) )
               FACT  = (UE - UJE) / (UJEP - UJE)
               DY    = FACT *
     1              SQRT ( (RC(JEDGE+JINC,K,I,1)-RC(JEDGE,K,I,1))**2
     2                   + (RC(JEDGE+JINC,K,I,2)-RC(JEDGE,K,I,2))**2
     3                   + (RC(JEDGE+JINC,K,I,3)-RC(JEDGE,K,I,3))**2 )
               RHOU  = Q(JEDGE,K,I,1) * UJE / FSMACH
               RHOE  = (1.E0 - FACT)*Q(JEDGE,K,I,1) 
     1                 + FACT*Q(JEDGE+JINC,K,I,1)
               RHOUP = RHOE * UE / FSMACH
               DELTA = DELTA + DY
               DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                 + (1.E0 - RHOUP) )*DY
               THETA = THETA + 0.5E0 * ( (1.E0 - RHOU) * RHOU
     1                                 + (1.E0 - RHOUP)* RHOUP )*DY
C
C Print data
C 
  209          CONTINUE
               WRITE (IOUT,1010) I, JBC, K, RC(JBC,K,I,1), 
     1                           RC(JBC,K,I,2), RC(JBC,K,I,3), 
     2                           PWQINF, CP, CF, CFEDGE, CH, CHEDGE, DN, 
     3                           YPLUS, USTAR, DELTA, DSTAR, THETA
  210       CONTINUE
  220    CONTINUE
      ENDIF
C
C Finished with Skin Friction and Heat Transfer
C
      RETURN
 1000 FORMAT (//' ',' Surface Data for BLOCK ',I5)
 1005 FORMAT ( /' ',3X,'I',3X,'J',3X,'K',6X,'X',12X,'Y',12X,'Z',9X,
     1              'Pw/Qinf',8X,'Cp',11X,'Cf',9X,'Cf_edge',8X,
     2              'Ch',9X,'Ch_edge',8X,'DN',10X,'Yplus',8X,
     3              'Ustar',8X,'Delta',8X,'Delta*',7X,'Theta'/)
 1010 FORMAT (' ',3I4,15(1PE13.5))
      END





      SUBROUTINE WRREST (NXBLKS, NXLVLS, NXPTS,  NXQNPT, NXPPTS, 
     1                   NXBCS,  NXSIZE, NXSECT,
     2                   IOFFS,  IOFFQ,  IOFFP,  IOFFQN,
     3                   NBLKS,  NLVLS,  IMODEL, IFDELQ, ITRUN,
     4                   IDIM,   JDIM,   KDIM, 
     5                   Q, QN, PROPS, S, VOL, NBCS, IBCDAT,
     6                   RES, QP, QM, DQP, DQM, F, FM, WORKX, STEMP,
     7                   FWRRST, IERRCD)
C
C Routine to write in the restart file for a block.
C
C N.B.  QN is converted to PRIMITIVE variables in restart file. 
C       QN is converted back to conserved variables before exiting this
C          routine to continue calculation.
C
C This routine is organized as follows:
C     0. Setup logicals to control additional output to restart file
C     1. Open restart file
C     2. Write (new) Residual History Data
C     3. Write Restart File Header Information 
C        a. Number of variables
C        b. Number of properties
C        c. Number of additional time levels and time for unsteady
C     4. For each level and block Write Flow Variables and Properties
C     5. Close restart file
C     6. Error handling section
C
C NXBLKS         : Maximum number of blocks dimensioned in main
C NXLVLS         : Maximum number of grid levels dimensioned in main
C NXPTS          : Maximum number of grid points dimensioned in main
C NXQNPT         : Maximum number of additional time level grid points 
C                  dimensioned in main
C NXPPTS         : Maximum number of property points dimensioned in main
C NXBCS          : Maximum number of boundary conditions dimensioned in main
C NXSIZE         : Maximum size of a single block as dimensioned in main
C NXSECT         : Maximum size of a cross section as dimensioned in main
C NBLKS          : Number of blocks
C NLVLS          : Number of grid levels
C IMODEL         : Data controlling turbulence model
C IFDELQ         : Array controlling calculation of delq
C ITRUN          : Number of iterations to write residual history
C IDIM,JDIM,KDIM : Dimensions of each block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C S              : Metrics
C VOL            : Cell volumes
C NBCS           : Number of boundary conditions for each block
C IBCDAT         : Data controlling boundary conditions
C RES,QP,QM,DQP,DQM,F,FM,WORKX,STEMP
C                : Used for temporary space
C FWRRST         : Restart filename
C IERRCD         : Error code
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION IDIM  (NXBLKS,NXLVLS),
     1          JDIM  (NXBLKS,NXLVLS),
     2          KDIM  (NXBLKS,NXLVLS)
C
      DIMENSION IOFFS (NXBLKS,NXLVLS),
     1          IOFFQ (NXBLKS,NXLVLS),
     2          IOFFP (NXBLKS,NXLVLS),
     3          IOFFQN(NXBLKS,NXLVLS)
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NXPTS *NQ),
     1          QN    (NXQNPT*NQ*NTMLVL),
     2          PROPS (NXPPTS*NP),
     3          S     (NXPTS *4*3),
     4          VOL   (NXPTS)
C
      DIMENSION IBCDAT(10,NXBCS,NXBLKS,NXLVLS),  
     1          NBCS  (NXBLKS)
C
      DIMENSION RES   (NXSIZE*NQ),
     1          QP    (NXSECT*NQ),
     2          QM    (NXSECT*NQ),
     3          DQP   (NXSECT*NQ),
     4          DQM   (NXSECT*NQ),
     5          F     (NXSECT*NQ),
     6          FM    (NXSECT*NQ),
     7          WORKX (NXSECT*NW),
     8          STEMP (NXSECT)
C
      DIMENSION IFDELQ(NQ)
C
      CHARACTER*80 FWRRST,
     1             FERROR
C
C Control of whether or not to write properties data and additional time levels
C
      LOGICAL   WRPROP, WRUNST
C
C 0.   Setup logicals to control additional output to restart file
C
      ITURB = IMODEL(1)
C
      IF (ITURB  .GE. ITLMNR) THEN
         WRPROP = .TRUE.
      ELSE
         WRPROP = .FALSE.
      ENDIF
C
      IF (ITIMED .EQ. IUNSTD) THEN
         WRUNST = .TRUE.
      ELSE
         WRUNST = .FALSE.
      ENDIF
C
C 1.   Open restart file
C
      IERRCD = INO
      FERROR = FWRRST
      OPEN (IWRRST,FILE=FWRRST,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1             ERR=900)
      REWIND (IWRRST)
C
C 2.   Write (new) Residual History Data
C
      DO 200 IT = 1, ITRUN + ITTOT
         WRITE (IWRRST) IT, R2(IT)
  200 CONTINUE
      ITF = -999
      R2F = -999.99
      WRITE (IWRRST) ITF, R2F
C
C 3.   Write Restart File Header Information
C 3-a. Number of flow variables
C
      WRITE (IWRRST) NQ
C
C 3-b. Number of properties
C
      IF (WRPROP) THEN
         WRITE (IWRRST) NP
      ENDIF
C
C 3-c. Number of additional time levels and time for unsteady
C
      IF (WRUNST) THEN
         WRITE (IWRRST) NTIME-1
         WRITE (IWRRST) TIME
      ENDIF
C
C 4.   For each level and block Write Flow Variables and Properties
C
      WRITE (IOUT,1000)
      DO 420 ILVL = 1, NLVLS
         DO 410 IBLK = 1, NBLKS
            IQ  = IOFFQ (IBLK,ILVL) * NQ          + 1
            IQN = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
            IP  = IOFFP (IBLK,ILVL) * NP          + 1
            IS  = IOFFS (IBLK,ILVL) * 4 * 3       + 1
            IV  = IOFFS (IBLK,ILVL)               + 1
            IF (ITURB .GE. ITLMNR) THEN
               IPSTRT = 1
               IPEND  = IDIM(IBLK,ILVL) + 1
               CALL PRPRTY (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), IPSTRT, IPEND, NXSECT,
     2                      Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                      IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     4                      IMODEL, IFDELQ,
     5                      STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
            ENDIF
            IF (ITIMED .EQ. IUNSTD) THEN
               NNPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1                 (KDIM(IBLK,ILVL)+3) * (NTIME-1)
               CALL CONC2P (ITURB, NNPTS, QN(IQN), RES)
            ENDIF
C
            WRITE (IOUT,1010) IBLK, ILVL
C
            CALL WRRSTB (ITURB, ITRUN, IBLK, ILVL, IDIM(IBLK,ILVL),
     1                  JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     2                  Q(IQ), QN(IQN), PROPS(IP), WRPROP, WRUNST)
C
            IF (ITIMED .EQ. IUNSTD) THEN
               CALL CONP2C (ITURB, NNPTS, QN(IQN), RES)
            ENDIF
  410    CONTINUE
  420 CONTINUE
C
C 5.   Close restart file
C
      CLOSE (IWRRST)
C
C     Finished writing restart file
C
      WRITE (IOUT,1000)
      RETURN
C
C 6.   Error handling section
C
  900 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1100) FERROR
      RETURN
C
 1000 FORMAT (////' ')
 1010 FORMAT (' ','WRREST: Writing Restart File for Block ',I4,
     1            '    Level ',I4)
 1100 FORMAT (/' ','WRREST: ERROR-> Error opening file ',A80,
     1        /' ', 8X,'RUN ABORTING!')
C
      END



      SUBROUTINE WRRSTB (ITURB, ITER, IBLOCK, ILVL, IDIM, JDIM, KDIM, 
     1                   Q, QN, PROPS, WRPROP, WRUNST)
C
C Routine to write the restart data for a block.
C
C N.B.  QN is converted to PRIMITIVE variables for restart file.
C
C This routine is organized as follows:
C     1. Write Flow Variables and Properties
C        a. Write Flow Variables
C        b. Write Properties
C        c. Write flow variables at additional time levels for unsteady flow
C
C ITURB          : Turbulence model
C ITER           : Number of iterations to write residual history
C IBLOCK         : Block number
C ILVL           : Grid level
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C WRPROP         : Logical controlling whether to write property data to file
C                  WRPROP = .true.  Write property data
C                         = .false. Do not write property data
C WRUNST         : Logical controlling whether to write unsteady data to file
C                  WRUNST = .true.  Write unsteady data
C                         = .false. Do not write unsteady data
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTIME-1),
     2          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Control of whether or not to write properties data and additional time levels
C
      LOGICAL   WRPROP, WRUNST
C
C 1.   Write Flow Variables and Properties
C 1-a. Write Flow Variables
C
      WRITE (IWRRST) Q
C
C 1-b. Write Properties
C
      IF (WRPROP) THEN
         WRITE (IWRRST) PROPS
      ENDIF
C
C 1-c. Write flow variables at additional time levels for unsteady flow
C
      IF (WRUNST) THEN
         WRITE (IWRRST) QN
      ENDIF
C
      RETURN
      END
