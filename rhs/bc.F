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
C     $Date: 2001/10/29 03:25:30 $
C     $State: Exp $
C     $Log: bc.F,v $
C     Revision 4.4  2001/10/29 03:25:30  jhmorr
C     Updated license information
C
C     Revision 4.3  2001/06/08 06:25:38  jhmorr
C     Fixed error in BCPRSB where S, PROPS are dimensioned but not passed
C
C     Revision 4.2  2001/06/08 04:56:15  jhmorr
C     Added notice
C
C     Revision 4.1.1.16  1998/09/22 14:16:17  jhmorr
C     Removed S from bchshr routine; it is not used.
C
c Revision 4.1.1.15  1998/06/20  18:55:58  jhmorr
c Modify shear condition to match Tavoularis and Corrsin.
c
c Revision 4.1.1.14  1998/06/11  00:29:32  jhmorr
c Removed debug output.
c
c Revision 4.1.1.13  1998/06/10  23:58:38  jhmorr
c Corrections to shear and domain constants in homogeneous shear BC.
c
c Revision 4.1.1.12  1998/06/10  19:19:40  jhmorr
c Added numbers for the homogeneous shear test case.
c
c Revision 4.1.1.11  1998/06/10  14:01:21  jhmorr
c Added check for INC being too large in homogeneous shear BC.
c
c Revision 4.1.1.10  1998/06/06  20:23:33  jhmorr
c Updates to corrected homogeneous shear BC.
c
c Revision 4.1.1.9  1998/06/06  18:55:11  jhmorr
c Corrections to calculate the initial interpolating point for the
c homogeneous shear BC.
c
c Revision 4.1.1.8  1998/06/06  14:38:51  jhmorr
c Corrected error for the second BC on the outer boundary for the homogeneous
c shear BC due to setting deltai = - deltai.
c
c Revision 4.1.1.7  1998/06/06  14:12:11  jhmorr
c Moved invariant calculation of interpolating coefficients outside
c of the loop.
c
c Revision 4.1.1.6  1998/06/06  13:54:44  jhmorr
c Modified calculation of DX0 to make it independent of the loop index.
c
c Revision 4.1.1.5  1998/06/06  13:42:46  jhmorr
c Correction in calculation of DX0 in homogeneous shear BC. DELTAI
c must have the sign of the difference.
c
c Revision 4.1.1.4  1998/06/04  21:24:38  jhmorr
c Corrections to homogeneous shear BC. Also correction to PERIODIC
c BC to loop over np for the properties array.
c
c Revision 4.1.1.3  1998/06/04  13:50:00  jhmorr
c Correction to indices and calculated dxs for homogeneous shear.
c
c Revision 4.1.1.2  1998/05/30  15:40:22  jhmorr
c Added call to homogeneous shear bc.
c
c Revision 4.1.1.1  1998/05/30  15:31:27  jhmorr
c Added initial coding of homogeneous shear boundary condition.
c
c Revision 4.1  1998/04/14  20:34:58  jhmorr
c Alpha 4.1.
c
C


      SUBROUTINE BC     (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   NBCS, IBCDAT, NXPROF, NPROFL, IOFQBC, QBC, 
     3                   VISCOS, IMODEL, WORK, MXSECT, CL, RCNTR)
C
C Routine to set boundary conditions. This routine calls the individual
C boundary condition routines to set each type of boundary condition.
C The boundary condition routines extend the BCs from I=2 to I=1 and
C I=IDIM to I=IDIM+1 and similarly for J and K faces to set the corner
C data. Therefore, calling this routine with IPSTRT=2 is identical to
C using IPSTRT=1 and IPEND=IDIM is identical to IPEND=IDIM+1.
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C NPPTS          : Number of points at which PROPS are stored.  This allows
C                  the storage of PROPS for viscous flows and the avoidance
C                  of the additional storage for inviscid flows.
C S              : Metrics
C VOL            : Cell volumes
C RC             : Position vector at cell centers
C NBCS           : Number of boundary conditions for current block
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
C NXPROF         : Maximum number of profile boundary condition segments
C NPROFL         : Number of profile boundary condition segments
C IOFQBC         : Offsets for QBC array
C QBC            : Boundary condition array - for PROFILE data
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence model array
C WORK           : Work array
C MXSECT         : Maximum size of the work array
C CL             : Lift coefficient. Used for 2D farfield vortex bc.
C RCNTR          : Center of airfoil. Used as location of point vortex.
C                  Also used as airfoil location at which to calculate moments.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (1:JDIM+1,1:KDIM+1,1:IDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IOFQBC(NXPROF)
      DIMENSION QBC   (*) 
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION WORK  (MXSECT,NQ)
C
      LOGICAL   VISCOS(3)
C
      DIMENSION RCNTR (3)
C
C Need to add Logic to avoid calling the IDIR = 1 (I) direction boundary
C     conditions when in the interior of the domain and using a planar
C     time integration scheme (e.g. LGS, MARCH).  Need to allow for the
C     call for other time integration schemes and when using a planar
C     scheme on the first plane for a conical solution procedure.
C
      ITURB = IMODEL(1)
C
      DO 900 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
C
C     Do not call K-direction boundary conditions for 2D flows
C
         IDIR   = IBCDAT(2,IBC)
         IF (.NOT.THREED .AND. .NOT.AXISYM .AND. IDIR .EQ. 3) GO TO 900
C
C     Wall boundary condition
C
         IF (IBCTYP .EQ. IWALL) THEN
            CALL BCWALL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Wall Functions boundary condition
C
         ELSE IF (IBCTYP .EQ. IWALFN) THEN
            CALL BCWALF (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, VOL, IBCDAT(1,IBC))
C
C     Tangency boundary condition
C
         ELSE IF (IBCTYP .EQ. ITAN) THEN
            CALL BCTAN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Farfield boundary condition
C
         ELSE IF (IBCTYP .EQ. IFAR) THEN
            CALL BCFAR  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     2D Farfield boundary condition with point vortex correction
C
         ELSE IF (IBCTYP .EQ. IFAR2D) THEN
            CALL BCFR2D (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, RC, IBCDAT(1,IBC), CL, RCNTR)
C
C     Extrapolation boundary condition
C     This is general for any number of equations
C
         ELSE IF (IBCTYP .EQ. IEXT) THEN
            CALL BCEXT  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Fix boundary condition
C
         ELSE IF (IBCTYP .EQ. IFIX) THEN
            CALL BCFIX  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Jet boundary condition
C
         ELSE IF (IBCTYP .EQ. IJET) THEN
            CALL BCJET  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Symmetry boundary condition
C
         ELSE IF (IBCTYP .EQ. ISYMXY .OR. IBCTYP .EQ. ISYMXZ
     1                               .OR. IBCTYP .EQ. ISYMYZ) THEN
            CALL BCSYM  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC), WORK)
C     
C     Periodic boundary condition
C
         ELSE IF (IBCTYP .EQ. IPERD) THEN
            CALL BCPERD (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Profile boundary condition
C
         ELSE IF (IBCTYP .EQ. IPROFL) THEN
            ISTRT  = IBCDAT( 4,IBC)
            IEND   = IBCDAT( 5,IBC)
            JSTRT  = IBCDAT( 6,IBC)
            JEND   = IBCDAT( 7,IBC)
            KSTRT  = IBCDAT( 8,IBC)
            KEND   = IBCDAT( 9,IBC)
            INPROF = IBCDAT(10,IBC)
            IQBC   = IOFQBC(INPROF) * NQ + 1
            NPTS   = (IEND-ISTRT+1) * (JEND-JSTRT+1) * (KEND-KSTRT+1)
            CALL BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC(IQBC), S, IBCDAT(1,IBC), IERROR)
C
C     Subsonic Profile boundary condition
C     - call profile bc first and then update the pressure
C
         ELSE IF (IBCTYP .EQ. IPROSB) THEN
            ISTRT  = IBCDAT( 4,IBC)
            IEND   = IBCDAT( 5,IBC)
            JSTRT  = IBCDAT( 6,IBC)
            JEND   = IBCDAT( 7,IBC)
            KSTRT  = IBCDAT( 8,IBC)
            KEND   = IBCDAT( 9,IBC)
            INPROF = IBCDAT(10,IBC)
            IQBC   = IOFQBC(INPROF) * NQ + 1
            NPTS   = (IEND-ISTRT+1) * (JEND-JSTRT+1) * (KEND-KSTRT+1)
            CALL BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC(IQBC), S, IBCDAT(1,IBC), IERROR)
            CALL BCPRSB (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, IBCDAT(1,IBC), IERROR)
C
C     Subsonic Inflow  boundary condition
C
         ELSE IF (IBCTYP .EQ. IINFLO) THEN
            CALL BCSIN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Subsonic Outflow boundary condition
C
         ELSE IF (IBCTYP .EQ. IOUTFL) THEN
            CALL BCSOUT (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Axisymmetric boundary condition
C
         ELSE IF (IBCTYP .EQ. IAXISM) THEN
            CALL BCAXI  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Axisymmetric centerline boundary condition
C
         ELSE IF (IBCTYP .EQ. IAXICL) THEN
            CALL BCAXCL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C     
C     Homogeneous Shear Periodic boundary condition
C
C N.B. Need to fix the WORK array here
C
         ELSE IF (IBCTYP .EQ. IHSHR) THEN
            CALL BCHSHR (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, IBCDAT(1,IBC), 
     2                   WORK(1,1), WORK(1,2))
C
C     Error
C
         ELSE
            WRITE (IOUT,1000) IBC, IBCTYP
         ENDIF
 900  CONTINUE
C
C     Finished with boundary condition routine
C
      RETURN
 1000 FORMAT (' ','BC    : ERROR-> An illegal boundary condition ',
     1            'type occurs in the ',I3,
     2            ' boundary condition and is = ',I3)
      END




      SUBROUTINE BCAXI  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to set the axisymmetric boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C 0. Constants for rotation (2 theta) as the grid is rotated theta and the
C    bc data point is a ghost point.
C
      PI        = ACOS(-1.E0)
      STHETA    = SIN (2.E0 * AXIDTH * PI / 180.E0)
      CTHETA    = COS (2.E0 * AXIDTH * PI / 180.E0)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C      SGNM   is sign for the rotation (negative for K = 1    boundary, 
C                                       positive for K = KDIM boundary)
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCM   = 1
         SGNM   = - 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = 1.E0
      ENDIF
C
      STHETA = STHETA * SGNM 
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            U    = Q(JN, KN, IN, 2)
            V    = Q(JN, KN, IN, 3)
            W    = Q(JN, KN, IN, 4)
C
            UAXI = U
            VAXI = V * CTHETA 
            WAXI = V * STHETA
C
            Q(JBC ,KBC ,IBC ,1) = Q(JN  ,KN  ,IN  ,1)
            Q(JBC ,KBC ,IBC ,2) = UAXI
            Q(JBC ,KBC ,IBC ,3) = VAXI
            Q(JBC ,KBC ,IBC ,4) = WAXI
            Q(JBC ,KBC ,IBC ,5) = Q(JN  ,KN  ,IN  ,5)
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon, k-omega
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
               Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
C Stresses in cylindrical coordinates are invariant in theta
C => Tau_ij_bc = Tau_ij_n (cylindrical coordinates) 
C    Note also that Tau_r,theta = Tau_theta,z = 0
C
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               TAUXXN = Q(JN  ,KN  ,IN  , 6)
               TAUYYN = Q(JN  ,KN  ,IN  , 7)
               TAUZZN = Q(JN  ,KN  ,IN  , 8)
               TAUXYN = Q(JN  ,KN  ,IN  , 9)
               TAUXZN = Q(JN  ,KN  ,IN  ,10)
               TAUYZN = Q(JN  ,KN  ,IN  ,11)
C
C Cylindrical coordinates (r,theta,z correspond to R,T(theta),S(streamwise))
C
               TAURR  = TAUYYN
               TAUTT  = TAUZZN
               TAUSS  = TAUXXN
               TAURS  = TAUXYN
C
C Back to Cartesian
C
               TAUXX  =   TAUSS
               TAUYY  =   CTHETA*CTHETA*TAURR + STHETA*STHETA*TAUTT 
               TAUZZ  =   STHETA*STHETA*TAURR + CTHETA*CTHETA*TAUTT 
               TAUXY  =   CTHETA * TAURS
               TAUXZ  = - STHETA * TAURS
               TAUYZ  =   CTHETA*STHETA*(TAUTT - TAURR)
C
               Q(JBC ,KBC ,IBC , 6) = TAUXX
               Q(JBC ,KBC ,IBC , 7) = TAUYY
               Q(JBC ,KBC ,IBC , 8) = TAUZZ
               Q(JBC ,KBC ,IBC , 9) = TAUXY
               Q(JBC ,KBC ,IBC ,10) = TAUXZ
               Q(JBC ,KBC ,IBC ,11) = TAUYZ
               Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with FIX boundary condition
C
      RETURN
      END




      SUBROUTINE BCAXCL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to set the axisymmetric centerline boundary condition for y = 0.
C N.B. ONLY FOR Y = 0.
C Set scalars equal to interior value, reflect v velocity.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C      SGNM   is sign for the rotation (negative for K = 1    boundary, 
C                                       positive for K = KDIM boundary)
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCN2  = 3
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C    Set scalars equal to interior value, reflect v velocity about y = 0.
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            Q(JBC ,KBC ,IBC ,1) =   Q(JN  ,KN  ,IN  ,1)
            Q(JBC ,KBC ,IBC ,2) =   Q(JN  ,KN  ,IN  ,2)
            Q(JBC ,KBC ,IBC ,3) = - Q(JN  ,KN  ,IN  ,3)
            Q(JBC ,KBC ,IBC ,4) =   Q(JN  ,KN  ,IN  ,4)
            Q(JBC ,KBC ,IBC ,5) =   Q(JN  ,KN  ,IN  ,5)
C
            Q(JBC2,KBC2,IBC2,1) =   Q(JN2 ,KN2 ,IN2 ,1)
            Q(JBC2,KBC2,IBC2,2) =   Q(JN2 ,KN2 ,IN2 ,2)
            Q(JBC2,KBC2,IBC2,3) = - Q(JN2 ,KN2 ,IN2 ,3)
            Q(JBC2,KBC2,IBC2,4) =   Q(JN2 ,KN2 ,IN2 ,4)
            Q(JBC2,KBC2,IBC2,5) =   Q(JN2 ,KN2 ,IN2 ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon, k-omega
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
               Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
               Q(JBC2,KBC2,IBC2,6) = Q(JN2 ,KN2 ,IN2 ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JN2 ,KN2 ,IN2 ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
C Infinity conditions
C
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               TAUXX = Q(JN  ,KN  ,IN  , 6)
               TAUYY = Q(JN  ,KN  ,IN  , 7)
               TAUZZ = Q(JN  ,KN  ,IN  , 8)
               TAUXY = Q(JN  ,KN  ,IN  , 9)
               TAUXZ = Q(JN  ,KN  ,IN  ,10)
               TAUYZ = Q(JN  ,KN  ,IN  ,11)
C
               Q(JBC ,KBC ,IBC , 6) = TAUXX
               Q(JBC ,KBC ,IBC , 7) = TAUYY
               Q(JBC ,KBC ,IBC , 8) = TAUZZ
               Q(JBC ,KBC ,IBC , 9) = TAUXY
               Q(JBC ,KBC ,IBC ,10) = TAUXZ
               Q(JBC ,KBC ,IBC ,11) = TAUYZ
               Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JN2 ,KN2 ,IN2 , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JN2 ,KN2 ,IN2 , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JN2 ,KN2 ,IN2 , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JN2 ,KN2 ,IN2 , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JN2 ,KN2 ,IN2 ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JN2 ,KN2 ,IN2 ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JN2 ,KN2 ,IN2 ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with FIX boundary condition
C
      RETURN
      END





      SUBROUTINE BCEXT  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do First Order Extrapolation boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for all variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for ALL Variables
C
      DO 130 L = 1, NQ
         DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,L) = Q(JN  ,KN  ,IN  ,L)
               Q(JBC2,KBC2,IBC2,L) = Q(JN  ,KN  ,IN  ,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C Finished with First Order Extrapolation boundary condition
C
      RETURN
      END




      SUBROUTINE BCFAR  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do FARFIELD (1-dimensional Riemann invariants) boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C Conditions at infinity
C
      SINF  = PREF / RHOINF**GAMMA
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C       U (dot) n > 0  for inflow  to   domain
C                 < 0  for outflow from domain
C The use of the SGNM variable sets up one set of logic
C to handle both the inner and outer boundary cases.
C
            SX    = S(JM,KM,IM,1,IDIR) * SGNM
            SY    = S(JM,KM,IM,2,IDIR) * SGNM
            SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
            QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1            + Q(JN,KN,IN,4) * SZ
            QINF  = UINF * SX + VINF * SY + WINF * SZ
            RHON  = Q(JN,KN,IN,1)
            PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
            AEXT  = SQRT (GAMMA * PN / RHON)
            REXT  = QEXT - 2.E0 * AEXT / GAMM1
            RINF  = QINF + 2.E0 * AINF / GAMM1
            QDOTN = 0.5E0 * (REXT + RINF)
            A     = 0.25E0 * GAMM1 * (RINF - REXT)
            RMACH = QDOTN / A
            SEXT  = PN   / RHON  **GAMMA
C
C Inflow
C
            IF (QDOTN .GT. 0.E0) THEN
C
C     Supersonic Inflow - Data Specified
C
               IF (ABS(RMACH) .GE. 1.0E0) THEN
                  Q(JBC ,KBC ,IBC ,1) = RHOINF
                  Q(JBC ,KBC ,IBC ,2) = UINF
                  Q(JBC ,KBC ,IBC ,3) = VINF
                  Q(JBC ,KBC ,IBC ,4) = WINF
                  Q(JBC ,KBC ,IBC ,5) = PREF
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ELSE
C
C     Subsonic Inflow - One Dimensional Riemann Problem
C     From Jameson and Baker, AIAA 83-1929 
C
                  RHOBC  = (A*A / (GAMMA*SINF))**(1.E0/GAMM1)
                  Q(JBC ,KBC ,IBC ,1) = RHOBC
                  Q(JBC ,KBC ,IBC ,2) = UINF + (QDOTN-QINF)*SX
                  Q(JBC ,KBC ,IBC ,3) = VINF + (QDOTN-QINF)*SY
                  Q(JBC ,KBC ,IBC ,4) = WINF + (QDOTN-QINF)*SZ
                  Q(JBC ,KBC ,IBC ,5) = RHOBC * A * A / GAMMA
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ENDIF
C
C Outflow
C
            ELSE
C
C     Supersonic Outflow - Data Extrapolated
C
               IF (ABS(RMACH) .GE. 1.0E0) THEN
                  Q(JBC ,KBC ,IBC ,1) = Q(JN  ,KN  ,IN  ,1)
                  Q(JBC ,KBC ,IBC ,2) = Q(JN  ,KN  ,IN  ,2)
                  Q(JBC ,KBC ,IBC ,3) = Q(JN  ,KN  ,IN  ,3)
                  Q(JBC ,KBC ,IBC ,4) = Q(JN  ,KN  ,IN  ,4)
                  Q(JBC ,KBC ,IBC ,5) = Q(JN  ,KN  ,IN  ,5)
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ELSE
C
C     Subsonic Outflow - One Dimensional Riemann Problem
C     From Jameson and Baker, AIAA 83-1929
C
                  RHOBC  = (A*A / (GAMMA*SEXT))**(1.E0/GAMM1)
                  Q(JBC ,KBC ,IBC ,1) = RHOBC
                  Q(JBC ,KBC ,IBC ,2) = Q(JN,KN,IN,2) + (QDOTN-QEXT)*SX
                  Q(JBC ,KBC ,IBC ,3) = Q(JN,KN,IN,3) + (QDOTN-QEXT)*SY
                  Q(JBC ,KBC ,IBC ,4) = Q(JN,KN,IN,4) + (QDOTN-QEXT)*SZ
                  Q(JBC ,KBC ,IBC ,5) = RHOBC * A * A / GAMMA
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ENDIF
            ENDIF
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC ,6) = TKEINF
c-                Q(JBC ,KBC ,IBC ,7) = EPSINF
C
c-                Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
c-                Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
                  Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
                  Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
                  Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
c-             ENDIF
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC ,6) = TKEINF
c-                Q(JBC ,KBC ,IBC ,7) = OMEINF
C
c-                Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
c-                Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
                  Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
                  Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
                  Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
c-             ENDIF
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC , 6) = TAUINF(1)
c-                Q(JBC ,KBC ,IBC , 7) = TAUINF(2)
c-                Q(JBC ,KBC ,IBC , 8) = TAUINF(3)
c-                Q(JBC ,KBC ,IBC , 9) = TAUINF(4)
c-                Q(JBC ,KBC ,IBC ,10) = TAUINF(5)
c-                Q(JBC ,KBC ,IBC ,11) = TAUINF(6)
c-                Q(JBC ,KBC ,IBC ,12) = EPSINF
C
c-                Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
c-                Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
c-                Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
c-                Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
c-                Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
c-                Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
c-                Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC , 6) = Q(JN  ,KN  ,IN  , 6)
                  Q(JBC ,KBC ,IBC , 7) = Q(JN  ,KN  ,IN  , 7)
                  Q(JBC ,KBC ,IBC , 8) = Q(JN  ,KN  ,IN  , 8)
                  Q(JBC ,KBC ,IBC , 9) = Q(JN  ,KN  ,IN  , 9)
                  Q(JBC ,KBC ,IBC ,10) = Q(JN  ,KN  ,IN  ,10)
                  Q(JBC ,KBC ,IBC ,11) = Q(JN  ,KN  ,IN  ,11)
                  Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
C
                  Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
                  Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
                  Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
                  Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
                  Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
                  Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
                  Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
c-             ENDIF
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with FARFIELD boundary condition
C
      RETURN
      END




      SUBROUTINE BCFR2D (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, RC, IBCDAT, CL, RCNTR)
C
C Routine to do 2D FARFIELD (1-dimensional Riemann invariants) boundary 
C condition with a point vortex correction. The theory of the farfield bc
C with a point vortex correction is available in AIAA Paper 85-0020 by
C J.L. Thomas and M. D. Salas.
C
C N.B. This BC is ONLY valid for a 2D airfoil flow.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
C RC             : Position vector at cell centers
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
C CL             : Lift coefficient.  Used to calculate circulation.
C RCNTR          : Center of airfoil. Used as location of point vortex.
C                  Also used as airfoil location at which to calculate moments.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          RC    (1:JDIM+1,1:KDIM+1,1:IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
      DIMENSION RCNTR (3)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C Conditions at infinity
C
      SINF   = PREF / RHOINF**GAMMA
C
C Calculate circulation, compressibility factor, and constants for the
C point vortex correction.
C     AREARF is the airfoil chord in 2D.
C
      PI     = ACOS(-1.E0)
      ALPHAR = ALPHA * PI / 180.E0
C
      CBETA  = SQRT ( ABS (1.E0 - FSMACH * FSMACH) )
      FCONST = 0.25E0 * FSMACH * CL * AREARF * CBETA / PI
c-----write (iout,'(a,2e15.7)') ' cbeta, fconst: ', cbeta, fconst
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C       U (dot) n > 0  for inflow  to   domain
C                 < 0  for outflow from domain
C The use of the SGNM variable sets up one set of logic
C to handle both the inner and outer boundary cases.
C
            SX    = S(JM,KM,IM,1,IDIR) * SGNM
            SY    = S(JM,KM,IM,2,IDIR) * SGNM
            SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
C Calculate the farfield velocity as the freestream velocity plus a correction
C due to a point vortex on the airfoil. 2D ONLY => WFAR = 0.
C      atan2 returns -pi < theta < pi
C
            XFAR  = RC(JBC,KBC,IBC,1) - RCNTR(1)
            YFAR  = RC(JBC,KBC,IBC,2) - RCNTR(2)
            RFAR  = SQRT (XFAR * XFAR + YFAR * YFAR)
            ANGLE = ATAN2 (YFAR, XFAR)
C
            FVORT = FCONST / 
     1              (RFAR * (1.E0 - (FSMACH*SIN(ANGLE-ALPHAR))**2))
            UFAR  = UINF + FVORT * SIN (ANGLE)
            VFAR  = VINF - FVORT * COS (ANGLE)
c-----write (iout,'(a,3i5,4e13.5)') ' i,j,k,x,y,r,a ', ibc, jbc, kbc,
c----1          xfar, yfar, rfar, angle
c-----write (iout,'(a,15x,3e13.5)') ' fvort,ufar,vf ', fvort, ufar, vfar
C
            QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
            QINF  = UFAR * SX + VFAR * SY 
            RHON  = Q(JN,KN,IN,1)
            PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
            AEXT  = SQRT (GAMMA * PN / RHON)
            REXT  = QEXT - 2.E0 * AEXT / GAMM1
            RINF  = QINF + 2.E0 * AINF / GAMM1
            QDOTN = 0.5E0 * (REXT + RINF)
            A     = 0.25E0 * GAMM1 * (RINF - REXT)
            RMACH = QDOTN / A
            SEXT  = PN   / RHON  **GAMMA
C
C Inflow
C
            IF (QDOTN .GT. 0.E0) THEN
C
C     Supersonic Inflow - Data Specified
C
               IF (ABS(RMACH) .GE. 1.0E0) THEN
                  Q(JBC ,KBC ,IBC ,1) = RHOINF
                  Q(JBC ,KBC ,IBC ,2) = UINF
                  Q(JBC ,KBC ,IBC ,3) = VINF
                  Q(JBC ,KBC ,IBC ,4) = WINF
                  Q(JBC ,KBC ,IBC ,5) = PREF
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ELSE
C
C     Subsonic Inflow - One Dimensional Riemann Problem
C     From Jameson and Baker, AIAA 83-1929 
C
                  RHOBC  = (A*A / (GAMMA*SINF))**(1.E0/GAMM1)
                  Q(JBC ,KBC ,IBC ,1) = RHOBC
                  Q(JBC ,KBC ,IBC ,2) = UFAR + (QDOTN-QINF)*SX
                  Q(JBC ,KBC ,IBC ,3) = VFAR + (QDOTN-QINF)*SY
                  Q(JBC ,KBC ,IBC ,4) = 0.E0
                  Q(JBC ,KBC ,IBC ,5) = RHOBC * A * A / GAMMA
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ENDIF
C
C Outflow
C
            ELSE
C
C     Supersonic Outflow - Data Extrapolated
C
               IF (ABS(RMACH) .GE. 1.0E0) THEN
                  Q(JBC ,KBC ,IBC ,1) = Q(JN  ,KN  ,IN  ,1)
                  Q(JBC ,KBC ,IBC ,2) = Q(JN  ,KN  ,IN  ,2)
                  Q(JBC ,KBC ,IBC ,3) = Q(JN  ,KN  ,IN  ,3)
                  Q(JBC ,KBC ,IBC ,4) = Q(JN  ,KN  ,IN  ,4)
                  Q(JBC ,KBC ,IBC ,5) = Q(JN  ,KN  ,IN  ,5)
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ELSE
C
C     Subsonic Outflow - One Dimensional Riemann Problem
C     From Jameson and Baker, AIAA 83-1929
C
                  RHOBC  = (A*A / (GAMMA*SEXT))**(1.E0/GAMM1)
                  Q(JBC ,KBC ,IBC ,1) = RHOBC
                  Q(JBC ,KBC ,IBC ,2) = Q(JN,KN,IN,2) + (QDOTN-QEXT)*SX
                  Q(JBC ,KBC ,IBC ,3) = Q(JN,KN,IN,3) + (QDOTN-QEXT)*SY
                  Q(JBC ,KBC ,IBC ,4) = Q(JN,KN,IN,4) + (QDOTN-QEXT)*SZ
                  Q(JBC ,KBC ,IBC ,5) = RHOBC * A * A / GAMMA
C
                  Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
                  Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
                  Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
                  Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
                  Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
               ENDIF
            ENDIF
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC ,6) = TKEINF
c-                Q(JBC ,KBC ,IBC ,7) = EPSINF
C
c-                Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
c-                Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
                  Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
                  Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
                  Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
c-             ENDIF
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC ,6) = TKEINF
c-                Q(JBC ,KBC ,IBC ,7) = OMEINF
C
c-                Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
c-                Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
                  Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
                  Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
                  Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
c-             ENDIF
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Check U (dot) n for inflow/outflow determination
C
               SX    = S(JM,KM,IM,1,IDIR) * SGNM
               SY    = S(JM,KM,IM,2,IDIR) * SGNM
               SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1               + Q(JN,KN,IN,4) * SZ
               QINF  = UINF * SX + VINF * SY + WINF * SZ
               RHON  = Q(JN,KN,IN,1)
               PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariants
C - the sign of U dot n handles the sign change in the Riemann invariants
C
               AEXT  = SQRT (GAMMA * PN / RHON)
               REXT  = QEXT - 2.E0 * AEXT / GAMM1
               RINF  = QINF + 2.E0 * AINF / GAMM1
               QDOTN = 0.5E0 * (REXT + RINF)
C
C Inflow - Data specified from freestream
C
c-             IF (QDOTN .GT. 0.E0) THEN
C
c-                Q(JBC ,KBC ,IBC , 6) = TAUINF(1)
c-                Q(JBC ,KBC ,IBC , 7) = TAUINF(2)
c-                Q(JBC ,KBC ,IBC , 8) = TAUINF(3)
c-                Q(JBC ,KBC ,IBC , 9) = TAUINF(4)
c-                Q(JBC ,KBC ,IBC ,10) = TAUINF(5)
c-                Q(JBC ,KBC ,IBC ,11) = TAUINF(6)
c-                Q(JBC ,KBC ,IBC ,12) = EPSINF
C
c-                Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
c-                Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
c-                Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
c-                Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
c-                Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
c-                Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
c-                Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
C
C Outflow - Data extrapolated from interior
C
c-             ELSE
C
                  Q(JBC ,KBC ,IBC , 6) = Q(JN  ,KN  ,IN  , 6)
                  Q(JBC ,KBC ,IBC , 7) = Q(JN  ,KN  ,IN  , 7)
                  Q(JBC ,KBC ,IBC , 8) = Q(JN  ,KN  ,IN  , 8)
                  Q(JBC ,KBC ,IBC , 9) = Q(JN  ,KN  ,IN  , 9)
                  Q(JBC ,KBC ,IBC ,10) = Q(JN  ,KN  ,IN  ,10)
                  Q(JBC ,KBC ,IBC ,11) = Q(JN  ,KN  ,IN  ,11)
                  Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
C
                  Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
                  Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
                  Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
                  Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
                  Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
                  Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
                  Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
c-             ENDIF
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with FARFIELD boundary condition
C
      RETURN
      END




      SUBROUTINE BCFIX  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to FIX the boundary condition to the freestream condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
            IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            Q(JBC ,KBC ,IBC ,1) = RHOINF
            Q(JBC ,KBC ,IBC ,2) = UINF
            Q(JBC ,KBC ,IBC ,3) = VINF
            Q(JBC ,KBC ,IBC ,4) = WINF
            Q(JBC ,KBC ,IBC ,5) = PREF
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = EPSINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = OMEINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
C Infinity conditions
C
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC , 6) = TAUINF(1)
               Q(JBC ,KBC ,IBC , 7) = TAUINF(2)
               Q(JBC ,KBC ,IBC , 8) = TAUINF(3)
               Q(JBC ,KBC ,IBC , 9) = TAUINF(4)
               Q(JBC ,KBC ,IBC ,10) = TAUINF(5)
               Q(JBC ,KBC ,IBC ,11) = TAUINF(6)
               Q(JBC ,KBC ,IBC ,12) = EPSINF
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with FIX boundary condition
C
      RETURN
      END




      SUBROUTINE BCJET  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do JET boundary condition.
C Presently, this is only useful where the jet conditions are prescribed, 
C i.e. supersonic jet. This can be used for single/multiple jets of the 
C same inflow conditions or for a shear layer simulation.
C
C The jet conditions are passed in through a common in non-dimensional form as:
C      UJET   = u(jet) / ainf
C      VJET   = v(jet) / ainf
C      WJET   = w(jet) / ainf
C      TJET   = T(jet) / Tinf
C      PJET   = p(jet) / (rhoinf * ainf**2)
C The value of the Jet density is then derived as:
C      RHOJET = GAMMA * PJET / TJET
C
C               gamma * p(jet) / (rhoinf * ainf**2)
C             = -----------------------------------
C                       T(jet) / Tinf
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
            IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            Q(JBC ,KBC ,IBC ,1) = RHOJET
            Q(JBC ,KBC ,IBC ,2) = UJET
            Q(JBC ,KBC ,IBC ,3) = VJET
            Q(JBC ,KBC ,IBC ,4) = WJET
            Q(JBC ,KBC ,IBC ,5) = PJET
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEJET
               Q(JBC ,KBC ,IBC ,7) = EPSJET
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEJET
               Q(JBC ,KBC ,IBC ,7) = OMEJET
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC , 6) = TAUJET(1)
               Q(JBC ,KBC ,IBC , 7) = TAUJET(2)
               Q(JBC ,KBC ,IBC , 8) = TAUJET(3)
               Q(JBC ,KBC ,IBC , 9) = TAUJET(4)
               Q(JBC ,KBC ,IBC ,10) = TAUJET(5)
               Q(JBC ,KBC ,IBC ,11) = TAUJET(6)
               Q(JBC ,KBC ,IBC ,12) = EPSJET
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with JET boundary condition
C
      RETURN
      END





      SUBROUTINE BCPERD (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do PERIODIC boundary condition in a coordinate direction.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for all variables
C     3. Set boundary conditions for all properties
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCPR  is the periodic location corresponding to IBCPT
C      IBCPR2 is the periodic location corresponding to IBCPT2
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCPR  = IJKDIM(IDIR)
         IBCPR2 = IJKDIM(IDIR) - 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCPR  = 2
         IBCPR2 = 3
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for ALL Variables
C
      DO 130 L = 1, NQ
         DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IPR  = IBCPR *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JPR  = IBCPR *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KPR  = IBCPR *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IPR2 = IBCPR2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JPR2 = IBCPR2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KPR2 = IBCPR2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,L) = Q(JPR ,KPR ,IPR ,L)
               Q(JBC2,KBC2,IBC2,L) = Q(JPR2,KPR2,IPR2,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C 3. Boundary conditions for ALL properties
C
      DO 330 L = 1, NP
         DO 320 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 310 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IPR  = IBCPR *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JPR  = IBCPR *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KPR  = IBCPR *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IPR2 = IBCPR2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JPR2 = IBCPR2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KPR2 = IBCPR2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               PROPS(JBC ,KBC ,IBC ,L) = PROPS(JPR ,KPR ,IPR ,L)
               PROPS(JBC2,KBC2,IBC2,L) = PROPS(JPR2,KPR2,IPR2,L)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
C
C Finished with PERIODIC boundary condition
C
      RETURN
      END




      SUBROUTINE BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC, S, IBCDAT, IERROR)
C
C Routine to set boundary condition from input profile
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C NPTS           : Number of points in QBC array for this bc segment
C Q              : Primitive variables at cell centers
C QBC            : Boundary data from the PROFILE data
C S              : Metrics
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
C IERROR         : Error code
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QBC   (NPTS,NQ),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C     Store boundary condition data in mnemonic names
C     N.B. Do not modify ISTRT or IEND as they are used to calculate indices
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
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C          IBC  is the boundary condition location
C          IBC2 is the boundary condition location outside the block from IBC
C          IN   is the point to the inside of the block from IBC
C
C
C     I Direction
C
      IF (IDIR .EQ. 1) THEN
         IF (INROUT .EQ. 1) THEN
            IBC  = 1
            IBC2 = 0
            IN   = 2
         ELSE
            IBC  = IDIM + 1
            IBC2 = IDIM + 2
            IN   = IDIM
         ENDIF
         DO 130 L = 1, NQ
            DO 120 K = KSTRT, KEND
               DO 110 J = JSTRT, JEND
                  JKBCPT = (J-JSTRT+1) + (K-KSTRT)*(JEND-JSTRT+1) 
                  Q(J,K,IBC,L)  = QBC(JKBCPT,L)
                  Q(J,K,IBC2,L) = Q(J,K,IBC,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C     J Direction
C
      ELSE IF (IDIR .EQ. 2) THEN
         IF (INROUT .EQ. 1) THEN
            JBC  = 1
            JBC2 = 0
            JN   = 2
         ELSE
            JBC  = JDIM + 1
            JBC2 = JDIM + 2
            JN   = JDIM
         ENDIF
         DO 230 L = 1, NQ
            DO 220 I = MAX(ISTRT, IPSTRT), MIN(IEND, IPEND)
               DO 210 K = KSTRT, KEND
                  IKBCPT = (K-KSTRT+1) + (I-ISTRT)*(KEND-KSTRT+1) 
                  Q(JBC,K,I,L)  = QBC(IKBCPT,L)
                  Q(JBC2,K,I,L) = Q(JBC,K,I,L)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     K Direction
C
      ELSE IF (IDIR .EQ. 3) THEN
         IF (INROUT .EQ. 1) THEN
            KBC  = 1
            KBC2 = 0
            KN   = 2
         ELSE
            KBC  = KDIM + 1
            KBC2 = KDIM + 2
            KN   = KDIM
         ENDIF
         DO 330 L = 1, NQ
            DO 320 I = MAX(ISTRT, IPSTRT), MIN(IEND, IPEND)
               DO 310 J = JSTRT, JEND
                  JIBCPT = (J-JSTRT+1) + (I-ISTRT)*(JEND-JSTRT+1) 
                  Q(J,KBC,I,L)  = QBC(JIBCPT,L)
                  Q(J,KBC2,I,L) = Q(J,KBC,I,L)
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
      ELSE
         WRITE (IOUT,1000) IDIR
         IERROR = IYES
      ENDIF
C
C     Finished with PROFILE boundary condition
C
      RETURN
 1000 FORMAT (' ','BC    : ERROR-> Illegal direction = ',I3,
     1            ' in BCPRFL.')
      END




      SUBROUTINE BCPRSB (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, IBCDAT, IERROR)
C
C Routine to modify the pressure for the SUBSONIC PROFILE CC boundary 
C condition. The routine bcprfl must first be called to set the profile
C data. This routine then applies a subsonic condition on the pressure.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for pressure
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
c--  1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
c--  2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for pressure
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            Q(JBC ,KBC ,IBC ,5) = Q(JN  ,KN  ,IN  ,5)
            Q(JBC2,KBC2,IBC2,5) = Q(JN  ,KN  ,IN  ,5)
  110    CONTINUE
  120 CONTINUE
C
C Finished with SUBSONIC PROFILE CC boundary condition
C
      RETURN
      END




      SUBROUTINE BCSIN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do SUBSONIC INFLOW boundary condition.  This routine is based
C on Jeff White's recommendation to use the following: specify total pressure
C and total Temperature, specify flow direction, extrapolate velocities from
C the interior, solve for the speed of sound, get M, calculate static
C conditions using isentropic relations.  Then specify turbulence quantities.
C
C This routine should be used for nozzles, tunnels, etc.
C
C NOTE: This routine cannot currently be used on a boundary that adjoins
C       a WALL or WALL FUNCTION boundary.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCM   = 1
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C Conditions at infinity
C
      SINF  = PREF / RHOINF**GAMMA
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Set the flow direction to use.  Jeff recommends three options:
C     i. Normal to boundary face - default
C    ii. Aligned with local grid - better for highly curved nozzles
C   iii. User input              - to handle other cases
C Current implementation is for i.
C
C The use of the SGNM variable sets up one set of logic
C to handle both the inner and outer boundary cases.
C
            SX    = S(JM,KM,IM,1,IDIR) * SGNM
            SY    = S(JM,KM,IM,2,IDIR) * SGNM
            SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
            QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1            + Q(JN,KN,IN,4) * SZ
            RHON  = Q(JN,KN,IN,1)
            PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariant
C - the sign of U dot n handles the sign change in the Riemann invariants
C
            AEXT  = SQRT (GAMMA * PN / RHON)
            REXT  = QEXT - 2.E0 * AEXT / GAMM1
            HTEXT = (AEXT*AEXT/GAMM1 + 0.5E0*QEXT*QEXT)
            A     = ( - REXT + SQRT(REXT*REXT 
     1                          +(GAMMA+1.E0)*(HTEXT-0.5E0*REXT*REXT)) )
     2              * GAMM1 / (GAMMA + 1.E0)
            QDOTN = REXT + 2.E0 * A / GAMM1
            RMACH = QDOTN / A
C
C Isentropic relationships
C
            ARG   = 1.E0 + 0.5E0 * GAMM1 * RMACH * RMACH
            P     = PTOTBC / ARG**(GAMMA/GAMM1)
            T     = TTOTBC / ARG
            RHO   = GAMMA * P / T
C
C Mean Variables
C
            Q(JBC ,KBC ,IBC ,1) = RHO
            Q(JBC ,KBC ,IBC ,2) = QDOTN * SX
            Q(JBC ,KBC ,IBC ,3) = QDOTN * SY
            Q(JBC ,KBC ,IBC ,4) = QDOTN * SZ
            Q(JBC ,KBC ,IBC ,5) = P
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = EPSINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = OMEINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC , 6) = TAUINF(1)
               Q(JBC ,KBC ,IBC , 7) = TAUINF(2)
               Q(JBC ,KBC ,IBC , 8) = TAUINF(3)
               Q(JBC ,KBC ,IBC , 9) = TAUINF(4)
               Q(JBC ,KBC ,IBC ,10) = TAUINF(5)
               Q(JBC ,KBC ,IBC ,11) = TAUINF(6)
               Q(JBC ,KBC ,IBC ,12) = EPSINF
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with SUBSONIC INFLOW boundary condition
C
      RETURN
      END




      SUBROUTINE BCSOUT (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do SUBSONIC OUTFLOW boundary condition.  This is based on Jeff
C White's recommendation for the following: specify static pressure and
C extrapolate everything else.
C
C This routine should be used for nozzles, tunnels, etc.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            Q(JBC ,KBC ,IBC ,1) = Q(JN  ,KN  ,IN  ,1)
            Q(JBC ,KBC ,IBC ,2) = Q(JN  ,KN  ,IN  ,2)
            Q(JBC ,KBC ,IBC ,3) = Q(JN  ,KN  ,IN  ,3)
            Q(JBC ,KBC ,IBC ,4) = Q(JN  ,KN  ,IN  ,4)
            Q(JBC ,KBC ,IBC ,5) = PBAKBC
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
               Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
               Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC , 6) = Q(JN  ,KN  ,IN  , 6)
               Q(JBC ,KBC ,IBC , 7) = Q(JN  ,KN  ,IN  , 7)
               Q(JBC ,KBC ,IBC , 8) = Q(JN  ,KN  ,IN  , 8)
               Q(JBC ,KBC ,IBC , 9) = Q(JN  ,KN  ,IN  , 9)
               Q(JBC ,KBC ,IBC ,10) = Q(JN  ,KN  ,IN  ,10)
               Q(JBC ,KBC ,IBC ,11) = Q(JN  ,KN  ,IN  ,11)
               Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with SUBSONIC OUTFLOW boundary condition
C
      RETURN
      END




      SUBROUTINE BCSYM  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT, SYM)
C
C Routine to do SYMMETRY boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Setup variables for even/odd symmetry
C     3. Set boundary conditions for all variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C SYM            : Variable for even/odd symmetry
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
      DIMENSION SYM   (NQ)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCN2  is the second point inside of the block from IBCPT
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCN2  = 3
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Set up the symmetry array such that all variables are
C    reflected across the symmetry plane.
C
      IF ( IBCTYP .EQ. ISYMXY ) THEN
         ISYMDR = 3
         IASYM1 = 10
         IASYM2 = 11
      ELSE IF ( IBCTYP .EQ. ISYMXZ ) THEN
         ISYMDR = 2
         IASYM1 = 9
         IASYM2 = 11
      ELSE IF ( IBCTYP .EQ. ISYMYZ ) THEN
         ISYMDR = 1
         IASYM1 = 9
         IASYM2 = 10
      ENDIF
      DO 100 L = 1, NQ
         SYM(L) = 1.E0
 100  CONTINUE
C
      SYM(ISYMDR+1) = -1.E0
      IF (ITURB .EQ. ITRS) THEN
         SYM(IASYM1) = -1.E0
         SYM(IASYM2) = -1.E0
      ENDIF
C
C 3. Boundary conditions for all Variables
C
      DO 130 L = 1, NQ
         DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
               Q(JBC ,KBC ,IBC ,L) = SYM(L) * Q(JN  ,KN  ,IN  ,L)
               Q(JBC2,KBC2,IBC2,L) = SYM(L) * Q(JN2 ,KN2 ,IN2 ,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C
C Finished with SYMMETRY boundary condition
C
      RETURN
      END





      SUBROUTINE BCTAN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do TANGENCY (inviscid wall) boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            RHON  = Q(JN ,KN ,IN ,1)
            PN    = Q(JN ,KN ,IN ,5)
            PN2   = Q(JN2,KN2,IN2,5)
C
c-1st order PW    = PN
            PW    = (9.E0 * PN - PN2) / 8.E0
            RHOW  = ( PW / (PN/RHON**GAMMA) )**(1.E0/GAMMA)
C
            UBAR  = Q(JN,KN,IN,2) * S(JM,KM,IM,1,IDIR)
     1            + Q(JN,KN,IN,3) * S(JM,KM,IM,2,IDIR)
     2            + Q(JN,KN,IN,4) * S(JM,KM,IM,3,IDIR)
C
            Q(JBC ,KBC ,IBC ,1) = RHOW
            Q(JBC ,KBC ,IBC ,2) = Q(JN,KN,IN,2)-UBAR*S(JM,KM,IM,1,IDIR)
            Q(JBC ,KBC ,IBC ,3) = Q(JN,KN,IN,3)-UBAR*S(JM,KM,IM,2,IDIR)
            Q(JBC ,KBC ,IBC ,4) = Q(JN,KN,IN,4)-UBAR*S(JM,KM,IM,3,IDIR)
            Q(JBC ,KBC ,IBC ,5) = PW
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               TKEN = Q(JN,KN,IN,6)
               EPSN = Q(JN,KN,IN,7)
               TKEW = TKEN
               EPSW = EPSN
C
               Q(JBC ,KBC ,IBC ,6) = TKEW
               Q(JBC ,KBC ,IBC ,7) = EPSW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
               TKEN = Q(JN,KN,IN,6)
               OMEN = Q(JN,KN,IN,7)
               TKEW = TKEN
               OMEW = OMEN
C
               Q(JBC ,KBC ,IBC ,6) = TKEW
               Q(JBC ,KBC ,IBC ,7) = OMEW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               TXXN = Q(JN,KN,IN, 6)
               TYYN = Q(JN,KN,IN, 7)
               TZZN = Q(JN,KN,IN, 8)
               TXYN = Q(JN,KN,IN, 9)
               TXZN = Q(JN,KN,IN,10)
               TYZN = Q(JN,KN,IN,11)
               EPSN = Q(JN,KN,IN,12)
               TXXW = TXXN
               TYYW = TYYN
               TZZW = TZZN
               TXYW = TXYN
               TXZW = TXZN
               TYZW = TYZN
               EPSW = EPSN
C
               Q(JBC ,KBC ,IBC , 6) = TXXW
               Q(JBC ,KBC ,IBC , 7) = TYYW
               Q(JBC ,KBC ,IBC , 8) = TZZW
               Q(JBC ,KBC ,IBC , 9) = TXYW
               Q(JBC ,KBC ,IBC ,10) = TXZW
               Q(JBC ,KBC ,IBC ,11) = TYZW
               Q(JBC ,KBC ,IBC ,12) = EPSW
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with TANGENCY boundary condition
C
      RETURN
      END




      SUBROUTINE BCWALL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do VISCOUS WALL boundary condition.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C Set control for higher order evaluation of PW
C     IHIORD = 0 for first  order
C     IHIORD = 1 for higher order
C
      IF (FOURTH) THEN
         IHIORD = 1
      ELSE
         IHIORD = 0
      ENDIF
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
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT(3) = 2
         IEND (3) = 2
      ENDIF
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
         IBCN3  = 4
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCN3  = IJKDIM(IDIR) - 2
         IBCM   = IJKDIM(IDIR)
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
            IN3  = IBCN3 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN3  = IBCN3 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN3  = IBCN3 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
c-1st order PW    = Q(JN,KN,IN,5)
            PN    = Q(JN ,KN ,IN ,5)
            PN2   = Q(JN2,KN2,IN2,5)
            PN3   = Q(JN3,KN3,IN3,5)
c-2nd order PW    = (9.E0 * PN - PN2) / 8.E0
c-3rd order PW    = (225.E0 * PN - 50.E0 * PN2 + 9.E0 * PN3) / 184.E0
            PW    = (1 - IHIORD) * PN +
     1              IHIORD * (225.E0*PN - 50.E0*PN2 + 9.E0*PN3) / 184.E0
            IF ( TWSPE ) THEN
               RHOW = GAMMA * PW / ( TWALL / TINF )
            ELSE
               RHOW = Q(JN,KN,IN,1)
            ENDIF
            Q(JBC ,KBC ,IBC ,1) = RHOW
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PW
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
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
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               TKEN   = Q    (JN ,KN ,IN ,6)
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
C Wall boundary condition for epsilon from Speziale
C epsilon+(wall) = 2 NU K+(1) / Y+(1)**2
C
               EPSLNW = 2.E0 * RMUW / RHOW * TKEN
     1                  / (YNORML * YNORML + RSMASQ) * FSMACH / RE
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = EPSLNW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
C-----Begin wall BC based on Wilcox's rough wall approach
C     Calculate the velocity parallel to the wall at the cell center off 
C     the wall.  U(parallel) = U - UBAR where U is the total velocity vector
C     and UBAR is the velocity normal to the wall.  (Velocity on the wall 
C     should be identically zero)
C     VELN is magnitude of the parallel velocity at point N
C
c              UBARN = S(JM,KM,IM,1,IDIR) * Q(JN,KN,IN,2)
c    1               + S(JM,KM,IM,2,IDIR) * Q(JN,KN,IN,3)
c    2               + S(JM,KM,IM,3,IDIR) * Q(JN,KN,IN,4)
c              VELN  = SQRT((Q(JN,KN,IN,2)-UBARN*S(JM,KM,IM,1,IDIR))**2
c    1                     +(Q(JN,KN,IN,3)-UBARN*S(JM,KM,IM,2,IDIR))**2
c    2                     +(Q(JN,KN,IN,4)-UBARN*S(JM,KM,IM,3,IDIR))**2)
C
C     Tau(wall) = mu * du/dn
C
c              DVELDN = VELN / (YNORML + RSMALL)
c              TAUW   = RMUW * DVELDN
C
C     Kr+         = utau *  rho(wall) * Kr / mu(wall)
C     Omega(wall) = utau**2 rho(wall) * Sr / mu(wall)
C     If Kr = 0, set Kr+ = 5 -> hydraulically smooth
C
c              IF (ROUGHK .GT. 0.E0) THEN
c                 RKRPLS = SQRT (RHOW * TAUW) * ROUGHK / RMUW
c              ELSE
c                 RKRPLS = 5.E0
c              ENDIF
C
c              IF (RKRPLS .GE. 25.0E0) THEN
c                 SR =  100.E0 / RKRPLS
c              ELSE
c                 SR = 2500.E0 / (RKRPLS * RKRPLS)
c              ENDIF
C
c              OMEGAW = DVELDN * SR
C-----End   wall BC based on Wilcox's rough wall approach
C-----Begin wall BC based on Menter (NASA TM 103975)
               OMEGAW = 60.E0*RMUW / (BKW*RHOW*YNORML*YNORML+RSMASQ)
     1                  *FSMACH/RE
C-----End   wall BC based on Menter (NASA TM 103975)
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = OMEGAW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
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
               RHOW  = Q    (JBC,KBC,IBC,1)
               RMUW  = PROPS(JBC,KBC,IBC,1)
C
               TKEN  = 0.5E0*(Q(JN,KN,IN,6)+Q(JN,KN,IN,7)+Q(JN,KN,IN,8))
               YNORML= ABS (PROPS(JN ,KN ,IN ,4))
C
C Wall boundary condition for epsilon from Speziale
C epsilon+(wall) = 2 NU K+(1) / Y+(1)**2
C
               EPSLNW = 2.E0 * RMUW / RHOW * TKEN
     1                  / (YNORML * YNORML + RSMASQ) * FSMACH / RE
C
               Q(JBC ,KBC ,IBC , 6) = 0.E0
               Q(JBC ,KBC ,IBC , 7) = 0.E0
               Q(JBC ,KBC ,IBC , 8) = 0.E0
               Q(JBC ,KBC ,IBC , 9) = 0.E0
               Q(JBC ,KBC ,IBC ,10) = 0.E0
               Q(JBC ,KBC ,IBC ,11) = 0.E0
               Q(JBC ,KBC ,IBC ,12) = EPSLNW
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with WALL boundary condition
C
      RETURN
      END





      SUBROUTINE BCWALF (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, VOL, IBCDAT)
C
C Routine to do TURBULENT WALL FUNCTION boundary condition.
C
C This routine is organized as follows:
C     1. Setup constants for turbulent wall functions
C     2. Setup variables for looping
C     3. Set boundary conditions for the mean variables
C     4. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3),
     3          ISTRTI(3),  
     4          IENDI (3)
C
C 1. Control of iterative process
C
      MAXIT  = 50
C
C 32-Bit machine can't converge 1x10-7
C
      TOLER  = 1.E-06
C
C No underrelaxation results in more planes not converging,
C underrelaxation with RELAX = .9 converges more rapidly (11 to 20 cycles)
C than none but some planes still do not converge (get in limit cycle)
C RELAX = .5 requires 14 iterations instead of 11 but is more robust
C
      RELAX  = 0.5E0
C
C Constants for wall functions
C
      RKAPVK = 0.41E0
      ESTAR  = 9.0E0
      YV     = 11.6E0
      SCALE  = SQRT (RE / FSMACH)
C
C 2. Setup memory for looping.
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
C ISTRTI, IENDI correspond to the Interior of the boundary condition face
C i.e., this does not extend to the edge/corner 
C
      DO 5 L = 1, 3
         ISTRTI(L) = ISTRT(L)
         IENDI (L) = IEND (L)
    5 CONTINUE
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT (3) = 2
         IEND  (3) = 2
         ISTRTI(3) = 2
         IENDI (3) = 2
      ENDIF
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
         IBCM2  = 2
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
         IBCM2  = IJKDIM(IDIR) - 1
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2-b. Set the boundary conditions on the edges/corners 
C
      IF (ISTRT(IC2) .EQ. 1) THEN
         I2 = ISTRT(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  210    CONTINUE
      ENDIF
C
      IF (IEND (IC2) .EQ. (IJKDIM(IC2)+1)) THEN
         I2 = IEND (IC2)
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  220    CONTINUE
      ENDIF
C
      IF (ISTRT(IC3) .EQ. 1) THEN
         I3 = ISTRT(IC3)
         DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  230    CONTINUE
      ENDIF
C
      IF (IEND (IC3) .EQ. (IJKDIM(IC3)+1)) THEN
         I3 = IEND (IC3)
         DO 240 I2 = ISTRT(IC2), IEND(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  240    CONTINUE
      ENDIF
C
C 3. Boundary conditions for Mean Variables using Wall Functions
C    Compute using the wall functions only in the INTERIOR of the boundary 
C
      DO 320 I3 = ISTRTI(IC3), IENDI(IC3)
CDIR$ IVDEP
         DO 310 I2 = ISTRTI(IC2), IENDI(IC2)
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
            IM2  = IBCM2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM2  = IBCM2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM2  = IBCM2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
c---------- RHOW   = Q(JBC,KBC,IBC,1)
            RHON   = Q(JN ,KN ,IN ,1)
            PN     = Q(JN ,KN ,IN ,5)
C
C Calculate velocity magnitude parallel to the wall
C
            UDOTN  = ( Q(JN,KN,IN,2) * S(JM,KM,IM,1,IDIR) +
     1                 Q(JN,KN,IN,3) * S(JM,KM,IM,2,IDIR) +
     2                 Q(JN,KN,IN,4) * S(JM,KM,IM,3,IDIR) ) * SGNM
            UPN    = Q(JN,KN,IN,2) - UDOTN * SGNM * S(JM,KM,IM,1,IDIR)
            VPN    = Q(JN,KN,IN,3) - UDOTN * SGNM * S(JM,KM,IM,2,IDIR)
            WPN    = Q(JN,KN,IN,4) - UDOTN * SGNM * S(JM,KM,IM,3,IDIR)
            UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
C
            UDOTN2 = ( Q(JN2,KN2,IN2,2) * S(JM2,KM2,IM2,1,IDIR) +
     1                 Q(JN2,KN2,IN2,3) * S(JM2,KM2,IM2,2,IDIR) +
     2                 Q(JN2,KN2,IN2,4) * S(JM2,KM2,IM2,3,IDIR) ) * SGNM
            UPN2   = Q(JN2,KN2,IN2,2)-UDOTN2*SGNM*S(JM2,KM2,IM2,1,IDIR)
            VPN2   = Q(JN2,KN2,IN2,3)-UDOTN2*SGNM*S(JM2,KM2,IM2,2,IDIR)
            WPN2   = Q(JN2,KN2,IN2,4)-UDOTN2*SGNM*S(JM2,KM2,IM2,3,IDIR)
            UPARN2 = SQRT( UPN2*UPN2 + VPN2*VPN2 + WPN2*WPN2 )
C
            TN     = GAMMA * PN / RHON
C
            RMUW   = PROPS(JBC,KBC,IBC,1)
            RMUN   = PROPS(JN ,KN ,IN ,1)
            YNORML = ABS (PROPS(JN ,KN ,IN ,4))
            YNORML = MAX (YNORML, RSMALL)
C
C Calculate normal Reynolds stress
C
            IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
               TKEN   = Q(JN,KN,IN,6)
               TAUNN  = 2.E0 / 3.E0 * TKEN
            ELSE IF (ITURB .EQ. ITRS) THEN
C
               RNX    = S(JM,KM,IM,1,IDIR) * SGNM
               RNY    = S(JM,KM,IM,2,IDIR) * SGNM
               RNZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               TAUXX  = Q(JN,KN,IN, 6)
               TAUYY  = Q(JN,KN,IN, 7)
               TAUZZ  = Q(JN,KN,IN, 8)
               TAUXY  = Q(JN,KN,IN, 9)
               TAUXZ  = Q(JN,KN,IN,10)
               TAUYZ  = Q(JN,KN,IN,11)
C
C Transform stresses to local streamwise, normal coordinates
C
               TAUNN  =   RNX*RNX*TAUXX + RNY*RNY*TAUYY + RNZ*RNZ*TAUZZ 
     1            + 2.E0*(RNX*RNY*TAUXY + RNX*RNZ*TAUXZ + RNY*RNZ*TAUYZ)
C
            ELSE
               TAUNN  = 0.E0
            ENDIF
c+++++++
c+++        IF (S(JM,KM,IM,4,IDIR) .LE. 0.E0 .AND.
c+++ 1          VOL(JM,KM,IM)        .LE. 0.E0) THEN
c+++           WRITE (iout,'(a,3i5)') ' i, j, k ', im, jm, km
c+++        END IF
c+++++++
C
C Wall pressure boundary condition and wall density
C
            PW     = PN + RHON * TAUNN
            IF ( TWSPE ) THEN
               TW    = TWALL / TINF
            ELSE
               TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
            ENDIF
            RHOW = GAMMA * PW / TW
C
C Initialize Tau_w for iterative procedure to previous result
C
            TAUW   = PROPS (JBC2,KBC2,IBC2,ITQTAU)
C
C Iterate on tau_w
C
            DO 300 ITER = 1, MAXIT
               YPLUSN = SCALE * YNORML * SQRT (ABS(TAUW) * RHOW) / RMUW
C
C Case 1: Y+_N > YV ==> u+ = Ln(y+) / kappa + B
C
               IF (YPLUSN .GT. YV) THEN
                  RMUTN = SCALE * RKAPVK * YNORML * SQRT(ABS(TAUW)*RHON)
                  IF ( TWSPE ) THEN
                     TW    = TWALL / TINF
                     QW    = TAUW * ( (TW - TN) / (PRT * UPARN + RSMALL) 
     1                              - 0.5E0 * GAMM1 * UPARN ) 
                  ELSE
                     QW    = 0.E0
                     TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
                  ENDIF
                  RHOW = GAMMA * PW / TW
C
C Calculate U_C_N
C
                  DEN1  = SQRT( QW*QW + 2.E0*TW/(GAMM1*PRT)*TAUW*TAUW )
                  DEN1  = MAX (DEN1, 1.e-08)
                  ARG11 = (QW + UPARN*TAUW) / DEN1
                  ARG12 = (QW) / DEN1
                  ARG11 = MIN (  1.E0, ARG11 )
                  ARG11 = MAX ( -1.E0, ARG11 )
                  ARG12 = MIN (  1.E0, ARG12 )
                  ARG12 = MAX ( -1.E0, ARG12 )
                  UCN   = SQRT( 2.E0*TW/(GAMM1*PRT) ) 
     1                    * ( ASIN(ARG11) - ASIN(ARG12) )
C
C Calculate Tau_w
C
                  TAUWN = RHOW*(RKAPVK*UCN/LOG(ESTAR*YPLUSN))**2
     1                    *RE/FSMACH
C
C Case 2: Y+_N < YV ==> u+ = y+
C
               ELSE
                  IF ( TWSPE ) THEN
                     TW    = TWALL / TINF
                     QW    = TAUW * ( (TW - TN) / (PRT * UPARN + RSMALL) 
     1                              - 0.5E0 * GAMM1 * UPARN ) 
                  ELSE
                     QW    = 0.E0
                     TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
                  ENDIF
                  RHOW = GAMMA * PW / TW
C
C Calculate Tau_w
C
                  TAUWN = (UPARN * RMUN / YNORML)
               ENDIF
C
C Update Tauw
C
               ERROR = TAUWN - TAUW
               TAUW  = (1.E0 - RELAX) * TAUW + RELAX * TAUWN
C
C Converged?
C
               RELERR = ABS(ERROR) / (ABS(TAUW) + 1.E-10)
               IF (RELERR .LT. TOLER) GO TO 305
  300       CONTINUE
            WRITE (IOUT,1010) IBC, JBC, KBC, MAXIT, RELERR, TOLER,YPLUSN
  305       CONTINUE
C
C Pass back the wall shear stress and heat transfer
C Note: this will always develop TAUW >= 0.  The direction of TAUW is
C       assumed to align with the local velocity direction.  SGNM is
C       required here to handle DIM end of coordinate direction.
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = TAUW
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = QW
C
C Calculate wall velocity so that du/dy at wall-function cell is from log-law.
C Velocity at the interface of the first cell (wall-function cell) is
C calculated as the average of the values on either side of the interface
C since the solution point in the wall-function cell is equidistant from
C the interface.  This calculates the magnitude of the wall velocity.  Its
C direction is taken as aligned with the velocity parallel to the wall.
C Therefore, UINT, DUDY2, DYINT are all positive quantities.
C
            UINT   = 0.5E0 * ( UPARN + UPARN2 )
            UTAUW  = SQRT (ABS(TAUW) / RHOW * FSMACH / RE)
            DUDY2  = UTAUW / (RKAPVK * YNORML)
            DYINT  = 0.5E0 * (YNORML + PROPS(JN2,KN2,IN2,4))
            UPWALL = UINT - DUDY2 * DYINT
C
C Set wall mean quantities
C
            UPARN  = MAX ( UPARN, RSMALL )
            Q(JBC ,KBC ,IBC ,1) = RHOW
            Q(JBC ,KBC ,IBC ,2) = UPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,3) = VPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,4) = WPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,5) = PW
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  310    CONTINUE
  320 CONTINUE
C
C 4-a. Boundary conditions for Two-Equation variables
C
C 4-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 420 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 410 I2 = ISTRT(IC2), IEND(IC2)
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
               RHOW   = Q    (JBC ,KBC ,IBC ,1)
C
               TAUW   = PROPS(JBC2,KBC2,IBC2,ITQTAU)
               YNORML = ABS (PROPS(JN  ,KN  ,IN  ,4))
               YNORML = MAX ( YNORML, RSMALL )
C
C Boundary condition for tke     = abs (tau_w) / rho / sqrt(C_mu)
C Boundary condition for epsilon = (sqrt(C_mu) k)^(3/2) / kappa y
C
               CMURT  = SQRT (CMU)
               TKEN   = ABS (TAUW) / RHOW / CMURT * FSMACH / RE
               TKEN   = MAX (TKEN, RSMALL)
C
               TERM   = CMURT * TKEN
               EPSLNN = TERM * SQRT (TERM) / (RKAPVK * YNORML)
               EPSLNN = MAX (EPSLNN, RSMALL)
               EPSLNW = EPSLNN
C
               Q(JN  ,KN  ,IN  ,6) = TKEN
               Q(JN  ,KN  ,IN  ,7) = EPSLNW
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = EPSLNW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  410       CONTINUE
  420    CONTINUE
C
C 4-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
      write (iout,'(a)') ' Wall functions not implemented for k-w'
      stop
         DO 440 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 430 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
  430       CONTINUE
  440    CONTINUE
C
C 4-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 460 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 450 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC   = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC   = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC   = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2  = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2  = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2  = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IN    = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN    = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN    = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IN2   = IBCN2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN2   = IBCN2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN2   = IBCN2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IM    = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM    = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM    = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Set up normals for local coordinate system (s=SP, n=RN, t) where
C n is normal to the surface, s is parallel to the surface aligned with
C the flow, and t is perpendicular to s, n (and is not needed)
C
               RNX    = S(JM,KM,IM,1,IDIR) * SGNM
               RNY    = S(JM,KM,IM,2,IDIR) * SGNM
               RNZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               UDOTN  = Q(JN,KN,IN,2) * RNX + Q(JN,KN,IN,3) * RNY +
     1                  Q(JN,KN,IN,4) * RNZ
               UPN    = Q(JN,KN,IN,2) - UDOTN * RNX
               VPN    = Q(JN,KN,IN,3) - UDOTN * RNY
               WPN    = Q(JN,KN,IN,4) - UDOTN * RNZ
               UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
               IF (UPARN .EQ. 0.E0) THEN
                  SPX    = SHATX
                  SPY    = SHATY
                  SPZ    = SHATZ
               ELSE
                  SPX    = UPN / UPARN
                  SPY    = VPN / UPARN
                  SPZ    = WPN / UPARN
               ENDIF
C
               RHOW   = Q    (JBC ,KBC ,IBC ,1)
C
               TAUW   = PROPS(JBC2,KBC2,IBC2,ITQTAU)
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
               YNORML = MAX ( YNORML, RSMALL )
C
               RHO    = Q(JN,KN,IN, 1)
               TAUXX  = Q(JN,KN,IN, 6)
               TAUYY  = Q(JN,KN,IN, 7)
               TAUZZ  = Q(JN,KN,IN, 8)
               TAUXY  = Q(JN,KN,IN, 9)
               TAUXZ  = Q(JN,KN,IN,10)
               TAUYZ  = Q(JN,KN,IN,11)
               TKEN   = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
               UTAUW  = SQRT (ABS(TAUW) / RHOW * FSMACH / RE)
               DUPDN  = UTAUW / (RKAPVK * YNORML)
C
C Transform stresses to local streamwise, normal coordinates
C
               TAUNN  =   RNX*RNX*TAUXX + RNY*RNY*TAUYY + RNZ*RNZ*TAUZZ 
     1            + 2.E0*(RNX*RNY*TAUXY + RNX*RNZ*TAUXZ + RNY*RNZ*TAUYZ)
               TAUSN  =   SPX*RNX*TAUXX + SPY*RNY*TAUYY + SPZ*RNZ*TAUZZ 
     1              + (SPX*RNY+SPY*RNX)*TAUXY + (SPX*RNZ+SPZ*RNX)*TAUXZ 
     2              + (SPY*RNZ+SPZ*RNY)*TAUYZ
C
C Set sign of du_parallel/dy_normal as opposite of sign of tausn
C
               DUPDN  = SIGN (DUPDN, (-TAUSN))
C
C Calculate production = - rho * du_i/dx_j * Tau_ij
C
               PROD   = - RHO * TAUSN * DUPDN
C
C Boundary condition for epsilon from Younis' thesis as restriction of
C the epsilon transport equation
C
               CRAT   = CEPS1 / CEPS2 * PROD / RHO
               TAUNN  = ABS (TAUNN)
               EPSLNN = 0.5E0 * (CRAT + SQRT (CRAT*CRAT + 
     1              4.E0*CEPSLN/CEPS2*TKEN*TKEN*TAUNN/(YNORML*YNORML)) )
               EPSLNW = EPSLNN
C
               Q(JN  ,KN  ,IN  ,12) = EPSLNW
C
               Q(JBC ,KBC ,IBC , 6) = 0.E0
               Q(JBC ,KBC ,IBC , 7) = 0.E0
               Q(JBC ,KBC ,IBC , 8) = 0.E0
               Q(JBC ,KBC ,IBC , 9) = 0.E0
               Q(JBC ,KBC ,IBC ,10) = 0.E0
               Q(JBC ,KBC ,IBC ,11) = 0.E0
               Q(JBC ,KBC ,IBC ,12) = EPSLNW
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  450       CONTINUE
  460    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with TURBULENT WALL FUNCTIONS boundary condition
C
      RETURN
 1010 FORMAT (' ','BCWALF: Tau_w not converged at (I:J:K) ', 3I5, 
     1            ' after ', I5, ' iterations.',
     2       /' ','  Rel. Error ', E15.7, 5X, 'Tolerance ', E15.7,
     3            ' Y+ = ',E15.7)
      END





      SUBROUTINE CUT    (IDIM1, JDIM1, KDIM1, Q1,
     1                   IDIM2, JDIM2, KDIM2, Q2, ICUTS, NV)
C
C Routine to do cut boundary condition.
C Flow data stored in Q1 are updated from data in the interior of Q2.
C This routine works for Q or PROPS and needs the size of the array.
C
C IDIM1,JDIM1,KDIM1 : Block 1 dimensions 
C Q1                : Block 1 data (Q or PROPS)
C IDIM2,JDIM2,KDIM2 : Block 2 dimensions 
C Q2                : Block 2 data (Q or PROPS)
C ICUTS             : Array controlling beginning/ending of cut
C NV                : Number of variables (NV = NQ for Q, NP for PROPS)
C
      include '../header/common.h'
      DIMENSION Q1   (0:JDIM1+2,0:KDIM1+2,0:IDIM1+2,NV),
     1          Q2   (0:JDIM2+2,0:KDIM2+2,0:IDIM2+2,NV)
C
      DIMENSION ICUTS(21)
C
      DIMENSION IJKDM1(3), IJKDM2(3)
      DIMENSION ISTR1(3),  IEND1(3), 
     1          ISTR2(3),  IEND2(3), 
     2          ISGN1(3),  ISGN2(3),
     3          IORD(3),   LEN(3)
C
C     Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDM1(1) = IDIM1
      IJKDM1(2) = JDIM1
      IJKDM1(3) = KDIM1
      IJKDM2(1) = IDIM2
      IJKDM2(2) = JDIM2
      IJKDM2(3) = KDIM2
C
C     Store boundary condition data in mnemonic names
C
      IDIR1    = ICUTS( 2)
      INOUT1   = ICUTS( 3)
      ISTR1(1) = ICUTS( 4)
      IEND1(1) = ICUTS( 5)
      ISTR1(2) = ICUTS( 6)
      IEND1(2) = ICUTS( 7)
      ISTR1(3) = ICUTS( 8)
      IEND1(3) = ICUTS( 9)
C
      IBLK2    = ICUTS(10)
      IDIR2    = ICUTS(11)
      INOUT2   = ICUTS(12)
      ISTR2(1) = ICUTS(13)
      IEND2(1) = ICUTS(14)
      ISTR2(2) = ICUTS(15)
      IEND2(2) = ICUTS(16)
      ISTR2(3) = ICUTS(17)
      IEND2(3) = ICUTS(18)
C
      IORD(1)  = ICUTS(19)
      IORD(2)  = ICUTS(20)
      IORD(3)  = ICUTS(21)
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INOUT1 = 1) or the outer boundary (INOUT1 > 1)
C          IBCPT  = boundary point of block 1
C          IBCINC = increment to second boundary point of block 1
C                   IBCPT2 = IBCPT + IBCINC
C          INR    = interior point of block 2.
C          INRINC = increment to second interior point of block 2
C                   INR2 = INR + INRINC
C
      IF (INOUT1 .EQ. 1) THEN
         IBCPT  =   1
         IBCINC = - 1
      ELSE
         IBCPT  = IJKDM1(IDIR1) + 1
         IBCINC = 1
      ENDIF
C
      IF (INOUT2 .EQ. 1) THEN
         INR    =   2
         INRINC =   1
      ELSE
         INR    =   IJKDM2(IDIR2)
         INRINC = - 1
      ENDIF
C
C     Find the length of the two outer loops and loop over these using
C     offsets and delta function to set the two cut data points to the
C     two interior data points of block 2.
C
      DO 100 L = 1, 3
C
         LEN(L) = ABS ( IEND1(L) - ISTR1(L) )
C
C     Increment/Decrement 
C
         IF ( IEND1(L) .GT. ISTR1(L) ) THEN
            ISGN1(L) =   1
         ELSE
            ISGN1(L) = - 1
         ENDIF
C
C     Increment/Decrement 
C
         IF ( IEND2(L) .GT. ISTR2(L) ) THEN
            ISGN2(L) =   1
         ELSE
            ISGN2(L) = - 1
         ENDIF
  100 CONTINUE
C
C     II1 first  boundary point of block 1
C     II2 second boundary point of block 1
C     IN1 first  interior point of block 2
C     IN2 second interior point of block 2
C
      IC1 = ICY (IDIR1, 1)
      IC2 = ICY (IDIR1, 2)
      IC3 = ICY (IDIR1, 3)
C
      JC1 = IORD (IC1)
      JC2 = IORD (IC2)
      JC3 = IORD (IC3)
C
      DO 800 I3 = 0, LEN(IC3)
         DO 700 I2 = 0, LEN(IC2)
            II1 = IBCPT                        * IKD (IC1, 1) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 1) 
            JJ1 = IBCPT                        * IKD (IC1, 2) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 2) 
            KK1 = IBCPT                        * IKD (IC1, 3) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 3) 
C
            II2 = II1 + IBCINC * IKD (IC1, 1)
            JJ2 = JJ1 + IBCINC * IKD (IC1, 2)
            KK2 = KK1 + IBCINC * IKD (IC1, 3)
C
            IN1 = INR                          * IKD (JC1, 1) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 1) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 1) 
            JN1 = INR                          * IKD (JC1, 2) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 2) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 2) 
            KN1 = INR                          * IKD (JC1, 3) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 3) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 3) 
C
            IN2 = IN1 + INRINC * IKD (JC1, 1)
            JN2 = JN1 + INRINC * IKD (JC1, 2)
            KN2 = KN1 + INRINC * IKD (JC1, 3)
C
CDIR$ SHORTLOOP
            DO 600 L = 1, NV
               Q1(JJ1,KK1,II1,L) = Q2(JN1,KN1,IN1,L)
               Q1(JJ2,KK2,II2,L) = Q2(JN2,KN2,IN2,L)
 600        CONTINUE
 700     CONTINUE
 800  CONTINUE
C
C     Finished with Cut boundary condition
C
      RETURN
      END





      SUBROUTINE BCHSHR (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, IBCDAT, QBC, PBC)
C
C Routine to do PERIODIC boundary condition for HOMOGENEOUS SHEAR in the
C J-direction ONLY.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for all variables
C     3. Set boundary conditions for all properties
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
C QBC,PBC        : Temporary storage for Q, PROPS for interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
      DIMENSION IBCDAT(10)
C
C Temporary storage for Q and PROPS for interpolation
C
      DIMENSION QBC   (3-IDIM:2*IDIM-1,NQ),
     1          PBC   (3-IDIM:2*IDIM-1,NP)
C
C 1. Setup memory for looping.
C
C
C Store boundary condition data in mnemonic names
C
      IBCTYP   = IBCDAT(1)
      IDIR     = IBCDAT(2)
      INROUT   = IBCDAT(3)
C
C This routine is set-up for the J-Direction ONLY!
C
      IF (IDIR .NE. 2) THEN
         WRITE (IOUT,1000) IDIR
         STOP
      END IF
C
C      Calculate the increment in the I coordinate for interpolation
C
      DLENGTH = 1.E0
c--   DUDY    = 0.
      DUDY    = 0.6847688E-01
      DELTAU  = DUDY * DLENGTH
C
      DX      = DLENGTH / FLOAT(IDIM-1)
      DELTA   = DELTAU * TIMEDT
      DELTAI  = DELTA / DX
C
C      Check to see that DELTAI < IDIM; Stop if an error
C
      IF (INT(DELTAI) .GT. IDIM) THEN
         WRITE (IOUT,1005) INT(DELTAI), IDIM
         STOP
      ENDIF
C
C 2.   Boundary conditions for ALL Variables and ALL Properties
C
C      Loop over all K planes
C
      DO K = 1, KDIM+1
C
C 2-a. Set up the arrays to interpolate in
C
         DO JLOOP = 0, 1
C
C      Set needed variables depending on whether the boundary is
C      the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C
            IF (INROUT .EQ. 1) THEN
               DINC   = DELTAI
               INC    = INT(DINC)
               DU     = - DELTAU
               IF (JLOOP .EQ. 0) THEN
                  JBC = 0
                  JPR = JDIM - 1
               ELSE
                  JBC = 1
                  JPR = JDIM
               END IF
            ELSE
               DINC   = - DELTAI
               INC    = INT(DINC)
               DU     = DELTAU
               IF (JLOOP .EQ. 0) THEN
                  JBC = JDIM + 1
                  JPR = 2
               ELSE
                  JBC = JDIM + 2
                  JPR = 3
               END IF
            END IF
C
C      Calculate the initial point in the interpolating stencil as the third
C      grid index to the left of the interpolating point.
C
            IF (DINC .LT. 0) INC = INC - 1
C
C      Calculate the interpolation coefficients. Due to the design of the
C      homogeneous shear case, the interpolation is independent of the
C      location. The interpolation coefficients are then invariant.
C      5th order Lagrange interpolating polynomial. Carnahan, Luther, and
C      Wilkes, p. 27.
C
            DX0 = DINC - FLOAT(INC) + 2.E0
            DX1 = DX0    - 1.E0
            DX2 = DX0    - 2.E0
            DX3 = DX0    - 3.E0
            DX4 = DX0    - 4.E0
            DX5 = DX0    - 5.E0
C
            A0  = -       DX1 * DX2 * DX3 * DX4 * DX5 / 120.E0
            A1  =   DX0       * DX2 * DX3 * DX4 * DX5 /  24.E0
            A2  = - DX0 * DX1       * DX3 * DX4 * DX5 /  12.E0
            A3  =   DX0 * DX1 * DX2       * DX4 * DX5 /  12.E0
            A4  = - DX0 * DX1 * DX2 * DX3       * DX5 /  24.E0
            A5  =   DX0 * DX1 * DX2 * DX3 * DX4       / 120.E0
C
C      Store Q and PROPS in a periodic one-dimensional array to allow
C      for a simple form of the interpolation.
C
            DO I = 2, IDIM
               DO L = 1, NQ
                  QBC(I       ,L) = Q    (JPR,K,I,L)
                  QBC(I+IDIM-1,L) = Q    (JPR,K,I,L)
                  QBC(I-IDIM+1,L) = Q    (JPR,K,I,L)
               END DO
               DO L = 1, NP
                  PBC(I       ,L) = PROPS(JPR,K,I,L)
                  PBC(I+IDIM-1,L) = PROPS(JPR,K,I,L)
                  PBC(I-IDIM+1,L) = PROPS(JPR,K,I,L)
               END DO
            END DO
C
C 2-b. Interpolate
C      5th order Lagrange interpolating polynomial. Carnahan, Luther, and
C      Wilkes, p. 27. The II index is the cell center location of the
C      first point in the interpolating polynomial corresponding to the
C      A0 coefficient.
C
            DO I = 1, IDIM+1
               II = I + INC - 2
C
               DO L = 1, NQ
                  Q    (JBC,K,I,L) = A0 * QBC(II  ,L) + A1 * QBC(II+1,L)
     1                             + A2 * QBC(II+2,L) + A3 * QBC(II+3,L)
     2                             + A4 * QBC(II+4,L) + A5 * QBC(II+5,L)
               END DO
               DO L = 1, NP
                  PROPS(JBC,K,I,L) = A0 * PBC(II  ,L) + A1 * PBC(II+1,L)
     1                             + A2 * PBC(II+2,L) + A3 * PBC(II+3,L)
     2                             + A4 * PBC(II+4,L) + A5 * PBC(II+5,L)
               END DO
            END DO
C
C     Add/Subtract the mean velocity gradient from U
C
            DO I = 1, IDIM+1
               Q(JBC,K,I,2) = Q(JBC,K,I,2) + DU
            END DO
C
         END DO
      END DO
C
C Finished with PERIODIC HOMOGENEOUS SHEAR boundary condition
C
      RETURN
 1000 FORMAT (' ','ERROR-> BCHSHR: This routine is for the J-direction',
     1            ' only.',
     2       /' ',15x,' An invalid direction ',I2,' has been specified',
     3       /' ',15x,' Run Aborting.')
 1005 FORMAT (' ','ERROR-> BCHSHR: The calculated offset for the ',
     1                             'J-direction ',I5,
     2       /' ','                exceeds the allowed (IDIM) ',I5,
     3       /' ',15x,' Run Aborting.')
      END
