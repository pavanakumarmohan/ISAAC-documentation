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
C     $Revision: 4.3 $
C     $Author: jhmorr $
C     $Date: 2001/10/29 03:25:34 $
C     $State: Exp $
C     $Log: resid.F,v $
C     Revision 4.3  2001/10/29 03:25:34  jhmorr
C     Updated license information
C
C     Revision 4.2  2001/06/08 04:56:18  jhmorr
C     Added notice
C
C     Revision 4.1  1998/04/14 20:34:58  jhmorr
C     Alpha 4.1.
C
C




      SUBROUTINE DQDT   (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, QN, VOL, DTJ, RES, QCP, TKE, IMODEL)
C
C Routine to calculate the dQ/dt term to add to residual.
C     t-ts: Physical Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C
C     tau-ts: Pseudo Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C                 + [phi' (Q^(m) - Q^(m-1))] / dtau
C
C     TDTHET : Theta for Time Dependent calculation (physical time term)
C     TDPHI  : Phi   for Time Dependent calculation (physical time term)
C     TDPHIP : Phi'  for Time Dependent calculation (pseudo   time term)
C
C +-------+-----------------+-------------+-------------------+---------------+
C |Scheme | Order           | Theta  Phi  |  Phi' store dQ^(m)| TDPDTJ TDPDEL |
C +-------+-----------------+-------------+-------------------+---------------+
C |t-ts   | 1st             |   1     0   |   -         No    | 1+Phi     0   |
C |       | 2nd             |   1    1/2  |   -         No    | 1+Phi     0   |
C |tau-ts | 2nd, 1st pseudo |   1    1/2  |   0         No    | 1+Phi'  1+Phi |
C |       | 2nd, 2nd pseudo |   1    1/2  |  1/2       Yes    | 1+Phi'  1+Phi |
C +-------+-----------------+-------------+-------------------+---------------+
C
C N.B. This routine assumes that the previous time level solutions are
C      passed as conserved variables rather than primitive variables.
C
C This routine is organized as follows:
C     1. Copy Q to work array
C     2. Convert Q from primitive to conserved variables
C     3. Calculate dQ/dt and add to residual array
C     4. Q^(m) - Q^(m-1) term for tau-ts 
C        a. Add Q^(m) - Q^(m-1) term for tau-ts with 2nd order psuedo 
C           time stepping
C        b. Shift Q^(m) to Q^(m-1) for next subiteration
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space)
C Q              : Primitive variables at iterative time level p
C QN             : Conserved variables at time levels n, n-1, etc.
C VOL            : Cell volumes
C DTJ            : Time step divided by cell volume
C RES            : Residual
C QCP            : Work array for conserved Q at level p
C TKE            : Work array for turbulent kinetic energy
C IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      INTEGER   IDIM, JDIM, KDIM
C
      DIMENSION Q   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN  (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTMLVL),
     2          VOL (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     3          DTJ (JDIM,KDIM,IDIM),
     4          RES (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
      DIMENSION QCP (2:JDIM,2:KDIM,NQ),
     1          TKE (2:JDIM,2:KDIM)
C
      DIMENSION IMODEL(NMDL)
C
C Local variables
C
      INTEGER   I, J, K, L, NPTS
C
C     Turbulence model
C
      ITURB = IMODEL(1)
C
C     Constants for temporal accuracy 
C     (see e.g. Pulliam AIAA 93-3360, Rai AIAA 95-0583, 
C               Rumsey et al. Computers & Fluids 1996)
C          CONP constant for level p, CONN constant for time n, 
C          CONNM1 constant for time n-1
C
      CONP   = 1.E0 + TDPHI
      CONN   = 1.E0 + 2.E0 * TDPHI
      CONNM1 = TDPHI
C
C      Loop over each cross section
C
      DO 900 I = IPSTRT, IPEND
C
C 1.   Copy Q to work array
C
         DO 130 L = 1, NQ
            DO 120 K = 2, KDIM
               DO 110 J = 2, JDIM
                  QCP(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C 2.   Convert Q from primitive to conserved
C
         NPTS = (JDIM-1) * (KDIM-1)
         CALL CONP2C (ITURB, NPTS, QCP, TKE)
C
C 3.   Calculate dQ/dt
C
         DO 330 L = 1, NQ
            DO 320 K = 2, KDIM
               DO 310 J = 2, JDIM
                  DQDTP = ( CONP   * QCP(J,K,L) - CONN * QN(J,K,I,L,1) 
     1                    + CONNM1 * QN(J,K,I,L,2) ) * VOL(J,K,I)/TIMEDT
                  RES(J,K,I,L) = RES(J,K,I,L) + DQDTP
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C
C 4.   Q^(m) - Q^(m-1) term for tau-ts 
C
         IF (ITMTYP .EQ. ITAUTS .AND. TDPHIP .NE. 0.E0) THEN
            IT = NTIME + 1
C
C 4-a. Add Q^(m) - Q^(m-1) term for tau-ts with 2nd order psuedo time stepping
C
            DO 430 L = 1, NQ
               DO 420 K = 2, KDIM
                  DO 410 J = 2, JDIM
                     DQMDT = TDPHIP * (QCP(J,K,L) - QN(J,K,I,L,IT))
     1                       / DTJ(J,K,I)
                     RES(J,K,I,L) = RES(J,K,I,L) + DQMDT
  410             CONTINUE
  420          CONTINUE
  430       CONTINUE
C
C 4-b. Shift Q^(m) to Q^(m-1) for next subiteration
C
            DO 460 L = 1, NQ
               DO 450 K = 2, KDIM
                  DO 440 J = 2, JDIM
                     QN(J,K,I,L,IT) = QCP(J,K,L)
  440             CONTINUE
  450          CONTINUE
  460       CONTINUE
         ENDIF
C
C      Next cross section
C
  900 CONTINUE
C
C      Finished calculating dQ/dt
C
      RETURN
      END





      SUBROUTINE FLUXI (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, VOL, 
     2                  RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                  IBCDAT, NBCS, 
     4                  IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
C     Routine to calculate the contriubtion to the residual from the
C     i direction.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
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
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RES            : Residual
C     DQP,DQM        : Storage for delta q = q(i) - q(i-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage for Reynolds stresses at interfaces
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
C     IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (MXSECT,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (IPSTRT-1:IPEND,2:JDIM,NF), 
     7          FM    (IPSTRT-1:IPEND,2:JDIM,NF), 
     8          STEMP (IPSTRT-1:IPEND,2:JDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over K calculating Q(+,-), F(Q(+,-)) and updating the Residual.
C
      DO 500 K = 2, KDIM
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
            CALL FDS (IFLUX, ITURB, NPTS, QP, QM, STEMP, DQP, DQM, F)
C
            DO 180 L = 1, NF
               DO 160 I = IPSTRT, IPEND
                  DO 140 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) + 
     1                              ( F(I,J,L)   * S(J,K,I,4,1) -
     2                                F(I-1,J,L) * S(J,K,I-1,4,1) )
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
         ENDIF
C
C     Viscous Residuals
C
         IF ( VISCOS ) THEN
C
C     Calculate the derivatives and averages of stored variables at the
C     interfaces.
C     Temporary variables passed to the viscous flux evaluation routine:
C         QP    for Q      (variables at interface)
C         FM    for PROPSI (properties at interface)
C               NOTE: This works properly only as long as NP <= NQ
C         QM    for DQDX   (derivatives of Q in x direction at interface)
C         DQP   for DQDY   (derivatives of Q in y direction at interface)
C         DQM   for DQDZ   (derivatives of Q in z direction at interface)
C     Temporary variables used internally in DELQ:
C         F(2-NQ) for SI (STEMP is dimensioned incorrectly to be used for
C                     the averaged metrics in the J direction.  F is large
C                     enough to provide the necessary space.  Do not use
C                     F for anything else in this call).
C         F(1)    for VOLAUX 
C
            CALL DELQI (IDIM, JDIM, KDIM, IPSTRT, IPEND, K,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  F(1,2,2), F(1,2,1), QP, FM, QM, DQP, DQM)
C
C     Viscous Terms 
C
            CALL FVIS(IMODEL, NPTS, QP, FM, QM, DQP, DQM, STEMP, TAU, F)
C
C Wall Functions - if not integrating to the wall
C
            IF ( IFWALF ) THEN
C
C Loop over WALL FUNCTION boundary conditions for I faces - Calculate 
C flux boundary conditions for Wall Functions 
C
               DO 300 IBC = 1, NBCS
                  IBCTYP = IBCDAT(1,IBC)
                  IDIR   = IBCDAT(2,IBC)
                  ISTRT  = IBCDAT(4,IBC)
                  IEND   = IBCDAT(5,IBC)
                  IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 1) THEN
C
                     INROUT = IBCDAT(3,IBC)
                     JSTRT  = MAX (IBCDAT(6,IBC), 2)
                     JEND   = MIN (IBCDAT(7,IBC), JDIM)
C
                     IF (INROUT .EQ. 1) THEN
                        IN     = 2
                        IM     = 1
                        IW     = 1
                        IBC2   = 0
                        SGNM   = 1.E0
                     ELSE
                        IN     = IDIM
                        IM     = IDIM
                        IW     = IDIM
                        IBC2   = IDIM + 2
                        SGNM   = - 1.E0
                     ENDIF
C
C Modify viscous wall fluxes - set all viscous fluxes to zero then correct
C
                     DO 270 L = 1, NQ
                        DO 260 J = JSTRT, JEND
                           F(IW,J,L) = 0.E0
  260                   CONTINUE
  270                CONTINUE
C
                     DO 280 J = JSTRT, JEND
C
C Non-dimensionalize variables for flux terms, SGNM to account for inner
C or outer end of coordinate 
C
                        TAUW   = PROPS(J,K,IBC2,ITQTAU)*SGNM*FSMACH/RE
                        QW     = PROPS(J,K,IBC2,ITQDOT)*SGNM*FSMACH/RE
C
C Tangency vector => aligned with parallel velocity at adjoining cell center
C Note: SGNM used in specification of unit vector parallel to wall in BCWALF
C       cancels and is removed here for computational efficiency.
C
                        UDOTN  = Q(J,K,IN,2) * S(J,K,IM,1,IDIR) +
     1                           Q(J,K,IN,3) * S(J,K,IM,2,IDIR) +
     2                           Q(J,K,IN,4) * S(J,K,IM,3,IDIR)
                        UPN    = Q(J,K,IN,2) - UDOTN * S(J,K,IM,1,IDIR)
                        VPN    = Q(J,K,IN,3) - UDOTN * S(J,K,IM,2,IDIR)
                        WPN    = Q(J,K,IN,4) - UDOTN * S(J,K,IM,3,IDIR)
                        UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
                        UPARN  = MAX ( UPARN, RSMALL )
C
C Calculate the viscous wall flux terms
C
                        F(IW,J,2) = TAUW * UPN / UPARN
                        F(IW,J,3) = TAUW * VPN / UPARN
                        F(IW,J,4) = TAUW * WPN / UPARN
                        F(IW,J,5) = - QW
  280                CONTINUE
                  ENDIF
  300          CONTINUE
            ENDIF
C
C     Add fluxes to residual
C
            DO 380 L = 1, NF
               DO 360 I = IPSTRT, IPEND
                  DO 340 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) - 
     1                              ( F(I,J,L)   * S(J,K,I,4,1) - 
     2                                F(I-1,J,L) * S(J,K,I-1,4,1) )
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
C
         ENDIF
  500 CONTINUE
C
C     Finished with i direction residual computation.
C
      RETURN
      END






      SUBROUTINE FLUXJ (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, VOL, 
     2                  RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                  IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX, 
     4                  VISCOS, IMODEL)
C
C     Routine to calculate the contriubtion to the residual from the
C     j direction.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
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
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RES            : Residual
C     DQP,DQM        : Storage for delta q = q(j) - q(j-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage of Reynolds stresses at interfaces
C     F,FM           : Fluxes at cell interfaces -
C                      if Flux Difference Splitting then F = Flux
C                      if Flux Vector Splitting then F = Fplus, FM = Fminus
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for J direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C     RKAP           : Kappa to use in MUSCL interpolation
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (MXSECT,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (JDIM,2:KDIM,NF), 
     7          FM    (JDIM,2:KDIM,NF), 
     8          STEMP (JDIM,2:KDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over I calculating Q(+,-), F(Q(+,-)) and updating the Residual.
C
      DO 500 I = IPSTRT, IPEND
C
C     Inviscid Residuals
C
         CALL MUSCLJ (IDIM, JDIM, KDIM, I,
     1                Q, DQP, DQM, QP, QM,
     2                IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C     Calculate the fluxes from q(plus) and q(minus).  Pass only
C     the appropriate metrics from the current direction.
C
         DO 120 K = 2, KDIM
            DO 110 J = 1, JDIM
               DO 100 L = 1, 3
                  STEMP(J,K,L) = S(J,K,I,L,2)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         NPTS = JDIM * (KDIM - 1)
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
            CALL FDS (IFLUX, ITURB, NPTS, QP, QM, STEMP, DQP, DQM, F)
C
            DO 180 L = 1, NF
               DO 160 K = 2, KDIM
                  DO 140 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) + 
     1                              ( F(J,K,L)   * S(J,K,I,4,2) -
     2                                F(J-1,K,L) * S(J-1,K,I,4,2) )
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
         ENDIF
C
C     Viscous Residuals
C
         IF ( VISCOS ) THEN
C
C     Calculate the derivatives and averages of stored variables at the
C     interfaces.
C     Temporary variables passed to the viscous flux evaluation routine:
C         QP    for Q      (variables at interface)
C         FM    for PROPSI (properties at interface)
C               NOTE: This works properly only as long as NP <= NQ
C         QM    for DQDX   (derivatives of Q in x direction at interface)
C         DQP   for DQDY   (derivatives of Q in y direction at interface)
C         DQM   for DQDZ   (derivatives of Q in z direction at interface)
C     Temporary variables used internally in DELQ:
C         F(2-NQ) for SJ (STEMP is dimensioned incorrectly to be used for
C                     the averaged metrics in the J direction.  F is large
C                     enough to provide the necessary space.  Do not use
C                     F for anything else in this call).
C         F(1)    for VOLAUX 
C
            CALL DELQJ (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  F(1,2,2), F(1,2,1), QP, FM, QM, DQP, DQM)
C
C     Viscous Terms 
C
            CALL FVIS(IMODEL, NPTS, QP, FM, QM, DQP, DQM, STEMP, TAU, F)
C
C Wall Functions - if not integrating to the wall
C
            IF ( IFWALF ) THEN
C
C Loop over WALL FUNCTION boundary conditions for J faces - Calculate 
C flux boundary conditions for Wall Functions 
C
               DO 300 IBC = 1, NBCS
                  IBCTYP = IBCDAT(1,IBC)
                  IDIR   = IBCDAT(2,IBC)
                  ISTRT  = IBCDAT(4,IBC)
                  IEND   = IBCDAT(5,IBC)
                  IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 2 .AND.
     1                     I .GE. ISTRT .AND.    I .LE. IEND) THEN
C
                     INROUT = IBCDAT(3,IBC)
                     KSTRT  = MAX (IBCDAT(8,IBC), 2)
                     KEND   = MIN (IBCDAT(9,IBC), KDIM)
C
                     IF (INROUT .EQ. 1) THEN
                        JN     = 2
                        JM     = 1
                        JW     = 1
                        JBC2   = 0
                        SGNM   = 1.E0
                     ELSE
                        JN     = JDIM
                        JM     = JDIM
                        JW     = JDIM
                        JBC2   = JDIM + 2
                        SGNM   = - 1.E0
                     ENDIF
C
C Modify viscous wall fluxes - set all viscous fluxes to zero then correct
C
                     DO 270 L = 1, NQ
                        DO 260 K = KSTRT, KEND
                           F(JW,K,L) = 0.E0
  260                   CONTINUE
  270                CONTINUE
C
                     DO 280 K = KSTRT, KEND
C
C Non-dimensionalize variables for flux terms, SGNM to account for inner
C or outer end of coordinate 
C
                        TAUW   = PROPS(JBC2,K,I,ITQTAU)*SGNM*FSMACH/RE
                        QW     = PROPS(JBC2,K,I,ITQDOT)*SGNM*FSMACH/RE
C
C Tangency vector => aligned with parallel velocity at adjoining cell center
C Note: SGNM used in specification of unit vector parallel to wall in BCWALF
C       cancels and is removed here for computational efficiency.
C
                        UDOTN  = Q(JN,K,I,2) * S(JM,K,I,1,IDIR) +
     1                           Q(JN,K,I,3) * S(JM,K,I,2,IDIR) +
     2                           Q(JN,K,I,4) * S(JM,K,I,3,IDIR)
                        UPN    = Q(JN,K,I,2) - UDOTN * S(JM,K,I,1,IDIR)
                        VPN    = Q(JN,K,I,3) - UDOTN * S(JM,K,I,2,IDIR)
                        WPN    = Q(JN,K,I,4) - UDOTN * S(JM,K,I,3,IDIR)
                        UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
                        UPARN  = MAX ( UPARN, RSMALL )
C
C Calculate the viscous wall flux terms
C
                        F(JW,K,2) = TAUW * UPN / UPARN
                        F(JW,K,3) = TAUW * VPN / UPARN
                        F(JW,K,4) = TAUW * WPN / UPARN
                        F(JW,K,5) = - QW
  280                CONTINUE
                  ENDIF
  300          CONTINUE
            ENDIF
C
C     Add fluxes to residual
C
            DO 380 L = 1, NF
               DO 360 K = 2, KDIM
                  DO 340 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) - 
     1                              ( F(J,K,L)   * S(J,K,I,4,2) - 
     2                                F(J-1,K,L) * S(J-1,K,I,4,2) )
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
C
         ENDIF
C
  500 CONTINUE
C
C     Finished with j direction residual computation.
C
      RETURN
      END






      SUBROUTINE FLUXK (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, VOL, 
     2                  RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                  IBCDAT, NBCS, 
     4                  IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
C     Routine to calculate the contriubtion to the residual from the
C     k direction.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
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
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RES            : Residual
C     DQP,DQM        : Storage for delta q = q(k) - q(k-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage for Reynolds stresses at interfaces
C     F,FM           : Fluxes at cell interfaces -
C                      if Flux Difference Splitting then F = Flux
C                      if Flux Vector Splitting then F = Fplus, FM = Fminus
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for K direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C     RKAP           : Kappa to use in MUSCL interpolation
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     IMODEL         : Turbulence model array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (MXSECT,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (KDIM,2:JDIM,NF), 
     7          FM    (KDIM,2:JDIM,NF), 
     8          STEMP (KDIM,2:JDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over I calculating Q(+,-), F(Q(+,-)) and updating the Residual.
C
      DO 500 I = IPSTRT, IPEND
C
C     Inviscid Residuals
C
         CALL MUSCLK (IDIM, JDIM, KDIM, I,
     1                Q, DQP, DQM, QP, QM,
     2                IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C     Calculate the fluxes from q(plus) and q(minus).  Pass only
C     the appropriate metrics from the current direction.
C
         DO 120 J = 2, JDIM
            DO 110 K = 1, KDIM
               DO 100 L = 1, 3
                  STEMP(K,J,L) = S(J,K,I,L,3)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         NPTS = KDIM * (JDIM - 1)
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
            CALL FDS (IFLUX, ITURB, NPTS, QP, QM, STEMP, DQP, DQM, F)
C
            DO 180 L = 1, NF
               DO 160 J = 2, JDIM
                  DO 140 K = 2, KDIM
                     RES(J,K,I,L) = RES(J,K,I,L) + 
     1                              ( F(K,J,L)   * S(J,K,I,4,3) -
     2                                F(K-1,J,L) * S(J,K-1,I,4,3) )
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
         ENDIF
C
C     Viscous Residuals
C
         IF ( VISCOS ) THEN
C
C     Calculate the derivatives and averages of stored variables at the
C     interfaces.
C     Temporary variables passed to the viscous flux evaluation routine:
C         QP    for Q      (variables at interface)
C         FM    for PROPSI (properties at interface)
C               NOTE: This works properly only as long as NP <= NQ
C         QM    for DQDX   (derivatives of Q in x direction at interface)
C         DQP   for DQDY   (derivatives of Q in y direction at interface)
C         DQM   for DQDZ   (derivatives of Q in z direction at interface)
C     Temporary variables used internally in DELQ:
C         F(2-NQ) for SK (STEMP is dimensioned incorrectly to be used for
C                     the averaged metrics in the J direction.  F is large
C                     enough to provide the necessary space.  Do not use
C                     F for anything else in this call).
C         F(1)    for VOLAUX 
C
            CALL DELQK (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  F(1,2,2), F(1,2,1), QP, FM, QM, DQP, DQM)
C
C     Viscous Terms 
C
            CALL FVIS(IMODEL, NPTS, QP, FM, QM, DQP, DQM, STEMP, TAU, F)
C
C Wall Functions - if not integrating to the wall
C
            IF ( IFWALF ) THEN
C
C Loop over WALL FUNCTION boundary conditions for K faces - Calculate 
C flux boundary conditions for Wall Functions 
C
               DO 300 IBC = 1, NBCS
                  IBCTYP = IBCDAT(1,IBC)
                  IDIR   = IBCDAT(2,IBC)
                  ISTRT  = IBCDAT(4,IBC)
                  IEND   = IBCDAT(5,IBC)
                  IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 3 .AND.
     1                     I .GE. ISTRT .AND.    I .LE. IEND) THEN
C
                     INROUT = IBCDAT(3,IBC)
                     JSTRT  = MAX (IBCDAT(6,IBC), 2)
                     JEND   = MIN (IBCDAT(7,IBC), JDIM)
C
                     IF (INROUT .EQ. 1) THEN
                        KN     = 2
                        KM     = 1
                        KW     = 1
                        KBC2   = 0
                        SGNM   = 1.E0
                     ELSE
                        KN     = KDIM
                        KM     = KDIM
                        KW     = KDIM
                        KBC2   = KDIM + 2
                        SGNM   = - 1.E0
                     ENDIF
C
C Modify viscous wall fluxes - set all viscous fluxes to zero then correct
C
                     DO 270 L = 1, NQ
                        DO 260 J = JSTRT, JEND
                           F(KW,J,L) = 0.E0
  260                   CONTINUE
  270                CONTINUE
C
                     DO 280 J = JSTRT, JEND
C
C Non-dimensionalize variables for flux terms, SGNM to account for inner
C or outer end of coordinate 
C
                        TAUW   = PROPS(J,KBC2,I,ITQTAU)*SGNM*FSMACH/RE
                        QW     = PROPS(J,KBC2,I,ITQDOT)*SGNM*FSMACH/RE
C
C Tangency vector => aligned with parallel velocity at adjoining cell center
C Note: SGNM used in specification of unit vector parallel to wall in BCWALF
C       cancels and is removed here for computational efficiency.
C
                        UDOTN  = Q(J,KN,I,2) * S(J,KM,I,1,IDIR) +
     1                           Q(J,KN,I,3) * S(J,KM,I,2,IDIR) +
     2                           Q(J,KN,I,4) * S(J,KM,I,3,IDIR)
                        UPN    = Q(J,KN,I,2) - UDOTN * S(J,KM,I,1,IDIR)
                        VPN    = Q(J,KN,I,3) - UDOTN * S(J,KM,I,2,IDIR)
                        WPN    = Q(J,KN,I,4) - UDOTN * S(J,KM,I,3,IDIR)
                        UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
                        UPARN  = MAX ( UPARN, RSMALL )
C
C Calculate the viscous wall flux terms
C
                        F(KW,J,2) = TAUW * UPN / UPARN
                        F(KW,J,3) = TAUW * VPN / UPARN
                        F(KW,J,4) = TAUW * WPN / UPARN
                        F(KW,J,5) = - QW
  280                CONTINUE
                  ENDIF
  300          CONTINUE
            ENDIF
C
C     Add fluxes to residual
C
            DO 380 L = 1, NF
               DO 360 J = 2, JDIM
                  DO 340 K = 2, KDIM
                     RES(J,K,I,L) = RES(J,K,I,L) - 
     1                              ( F(K,J,L)   * S(J,K,I,4,3) -
     2                                F(K-1,J,L) * S(J,K-1,I,4,3) )
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
C
         ENDIF
C
  500 CONTINUE
C
C     Finished with k direction residual computation.
C
      RETURN
      END




      SUBROUTINE FLUXKA (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, 
     4                   IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
C Routine to calculate the axisymmetric flux contriubtion to the residual 
C from the k direction.  This routine is used for either the time dependent 
C calculation or the marching calculation.  
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
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
C RES            : Residual
C DQP,DQM        : Storage for delta q = q(k) - q(k-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses at interfaces
C F,FM           : Fluxes at cell interfaces -
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Storage for metrics passed to flux evaluation routine.
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for K direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence model array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (KDIM,2:JDIM,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (KDIM,2:JDIM,NF), 
     7          FM    (KDIM,2:JDIM,NF), 
     8          STEMP (KDIM,2:JDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over I calculating axisymmetric fluxes and updating the residual.
C
      DO 500 I = IPSTRT, IPEND
C
C 1.  Inviscid fluxes: Pressure only on the v- and w-momentum equations
C
         DO 110 J = 2, JDIM
            P            = Q(J,2,I,5)
C
            F(1,J,3)     = P * S(J,1,I,2,3)
            F(2,J,3)     = P * S(J,2,I,2,3)
            RES(J,2,I,3) = RES(J,2,I,3) + 
     1                     ( F(2,J,3) * S(J,2,I,4,3) -
     2                       F(1,J,3) * S(J,1,I,4,3) )
C
            F(1,J,4)     = P * S(J,1,I,3,3)
            F(2,J,4)     = P * S(J,2,I,3,3)
            RES(J,2,I,4) = RES(J,2,I,4) + 
     1                     ( F(2,J,4) * S(J,2,I,4,3) -
     2                       F(1,J,4) * S(J,1,I,4,3) )
  110    CONTINUE
C
C 2.  Viscous Residuals
C
         IF ( VISCOS ) THEN
C
C     Calculate the derivatives and averages of stored variables at the
C     interfaces.
C     Temporary variables passed to the viscous flux evaluation routine:
C         QP    for Q      (variables at interface)
C         FM    for PROPSI (properties at interface)
C               NOTE: This works properly only as long as NP <= NQ
C         QM    for DQDX   (derivatives of Q in x direction at interface)
C         DQP   for DQDY   (derivatives of Q in y direction at interface)
C         DQM   for DQDZ   (derivatives of Q in z direction at interface)
C     Temporary variables used internally in DELQ:
C         F(2-NQ) for SK (STEMP is dimensioned incorrectly to be used for
C                     the averaged metrics in the J direction.  F is large
C                     enough to provide the necessary space.  Do not use
C                     F for anything else in this call).
C         F(1)    for VOLAUX 
C
            CALL DELQK (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  F(1,2,2), F(1,2,1), QP, FM, QM, DQP, DQM)
C
C     Metrics from the current direction.
C
            DO 320 K = 1, KDIM
               DO 310 J = 2, JDIM
                  STEMP(K,J,1) = S(J,K,I,1,3)
                  STEMP(K,J,2) = S(J,K,I,2,3)
                  STEMP(K,J,3) = S(J,K,I,3,3)
  310          CONTINUE
  320       CONTINUE
            NPTS = KDIM * (JDIM - 1)
C
C     Viscous Terms 
C
            CALL FVIS(IMODEL, NPTS, QP, FM, QM, DQP, DQM, STEMP, TAU, F)
C
C     Add fluxes to residual
C
            DO 380 L = 1, NF
               DO 360 K = 2, KDIM
                  DO 340 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) - 
     1                              ( F(K,J,L)   * S(J,K,I,4,3) -
     2                                F(K-1,J,L) * S(J,K-1,I,4,3) )
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
C
         ENDIF
C
  500 CONTINUE
C
C     Finished with k direction axisymmetric residual computation.
C
      RETURN
      END




      SUBROUTINE FLUXM (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, VOL, RES, 
     2                  DQP, DQM, QP, QM, FP, FM, STEMP,
     3                  IBCDAT, NBCS, 
     4                  IORDER, LIMIT, VGNRN, SIGMA, IMODEL)
C
C     Routine to calculate the contribution to the residual from the
C     i direction for marching.  Uses fully upwind either first or
C     second order accurate in i.
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
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RES            : Residual
C     DQP,DQM        : Used as temporary storage for FLUX presently
C     QP,QM          : Q at a cell interface for MUSCL scheme
C                      QP = Q(i+1/2), QM = Q(i-1/2)
C     FP,FM          : Fluxes at cell interfaces
C                      FP = F(i+1/2), FM = F(i-1/2)
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for I direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C                      If LIMIT > 0 then do catastrophic limiting in the
C                      I direction for marching
C     VGNRN          : Logical for Vigneron approach in marching flux 
C                      evaluation
C     SIGMA          : Safety factor for Vigneron technique when marching
C     IMODEL         : Turbulence model array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (2:JDIM,2:KDIM,NQ), 
     2          DQM   (2:JDIM,2:KDIM,NQ), 
     3          QP    (2:JDIM,2:KDIM,NQ), 
     4          QM    (2:JDIM,2:KDIM,NQ), 
     5          FP    (2:JDIM,2:KDIM,NF), 
     6          FM    (2:JDIM,2:KDIM,NF), 
     7          STEMP (2:JDIM,2:KDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VGNRN
C
      ITURB = IMODEL(1)
C
C     ICAT counts the number of times that the catastrophic limiter is on
C
      ICAT = 0
C
C     Calculate everything at the cross section I = IPSTRT = IPEND
C
      IF (IPSTRT .NE. IPEND) THEN
         WRITE (IOUT,1000) IPSTRT, IPEND
         STOP
      ENDIF
C
      I = IPSTRT
C
C     Inviscid Residuals
C
      IF (IORDER .EQ. 1) THEN
         DO 60 L = 1, NQ
            DO 50 K = 2, KDIM
               DO 40 J = 2, JDIM
                  QP(J,K,L) = Q(J,K,I,L) 
                  QM(J,K,L) = Q(J,K,I-1,L) 
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
      ELSE
         DO 90 L = 1, NQ
            DO 80 K = 2, KDIM
               DO 70 J = 2, JDIM
                  QP(J,K,L) = 1.5E0*Q(J,K,I,L)   - 0.5E0*Q(J,K,I-1,L)
                  QM(J,K,L) = 1.5E0*Q(J,K,I-1,L) - 0.5E0*Q(J,K,I-2,L)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
C
C     Catastrophic limiting for higher order - check density and pressure
C
         IF ( LIMIT .GT. 0 ) THEN
            EPS  = 0.E0
            ICAT = 0
            DO 95 K = 2, KDIM
               DO 94 J = 2, JDIM
                  IF (QP(J,K,1) .LE. EPS .OR. QP(J,K,5) .LE. EPS) THEN
                     ICAT = ICAT + 1
                     DO 91 L = 1, NQ
                        QP(J,K,L) = Q(J,K,I,L)
   91                CONTINUE
                  ENDIF
                  IF (QM(J,K,1) .LE. EPS .OR. QM(J,K,5) .LE. EPS) THEN
                     ICAT = ICAT + 1
                     DO 92 L = 1, NQ
                        QM(J,K,L) = Q(J,K,I-1,L)
   92                CONTINUE
                  ENDIF
   94          CONTINUE
   95       CONTINUE
            IF (ICAT .GT. 0) THEN
               WRITE (IOUT,1200) ICAT
            ENDIF
         ENDIF
      ENDIF
C
C     Calculate the fluxes F(i+1/2) from q(plus). Pass only
C     the appropriate metrics from the current direction.
C
      DO 103 L = 1, 3
         DO 102 K = 2, KDIM
            DO 101 J = 2, JDIM
               STEMP(J,K,L) = S(J,K,I,L,1)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
C
C     Call Flux evaluation routine.  For temporary storage use:
C        DQP => TKE
C        DQM => H
C
      NPTS = (JDIM - 1) * (KDIM - 1)
      CALL FLUX (ITURB, NPTS, QP, STEMP, DQP, DQM, FP)
C
C     Calculate the fluxes F(i-1/2) from q(minus). Pass only
C     the appropriate metrics from the current direction.
C
      DO 113 L = 1, 3
         DO 112 K = 2, KDIM
            DO 111 J = 2, JDIM
               STEMP(J,K,L) = S(J,K,I-1,L,1)
  111       CONTINUE
  112    CONTINUE
  113 CONTINUE
C
C     Call Flux evaluation routine.  For temporary storage use:
C        DQP => TKE
C        DQM => H
C
      NPTS = (JDIM - 1) * (KDIM - 1)
      CALL FLUX (ITURB, NPTS, QM, STEMP, DQP, DQM, FM)
C
C     Add the flux contributions to the residual.
C
      DO 180 L = 1, NF
         DO 160 K = 2, KDIM
            DO 140 J = 2, JDIM
               RES(J,K,I,L) = RES(J,K,I,L) + 
     1                        ( FP(J,K,L) * S(J,K,I,4,1) -
     2                          FM(J,K,L) * S(J,K,I-1,4,1) )
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
C
C     Source term for the Vigneron technique
C
      IF (VGNRN) THEN
         DO 290 K = 2, KDIM
            DO 280 J = 2, JDIM
               UBAR  = 0.5E0*(S(J,K,I,1,1)+S(J,K,I-1,1,1))*Q(J,K,I,2)
     1               + 0.5E0*(S(J,K,I,2,1)+S(J,K,I-1,2,1))*Q(J,K,I,3)
     2               + 0.5E0*(S(J,K,I,3,1)+S(J,K,I-1,3,1))*Q(J,K,I,4)
               AMSQ  = UBAR*UBAR / (GAMMA * Q(J,K,I,5) / Q(J,K,I,1))
               OMEGA = MIN(1.E0, SIGMA*GAMMA*AMSQ/(1.E0+GAMM1*AMSQ) )
C
               PPI   = QP(J,K,5)
               PMI   = QM(J,K,5)
               PPJ   = 0.5E0*(Q(J,K,I,5) + Q(J+1,K,I,5))
               PMJ   = 0.5E0*(Q(J,K,I,5) + Q(J-1,K,I,5))
               PPK   = 0.5E0*(Q(J,K,I,5) + Q(J,K+1,I,5))
               PMK   = 0.5E0*(Q(J,K,I,5) + Q(J,K-1,I,5))
C
               UPI   = S(J,K,I,1,1)   * QP(J,K,2) +
     1                 S(J,K,I,2,1)   * QP(J,K,3) +
     2                 S(J,K,I,3,1)   * QP(J,K,4)
               UMI   = S(J,K,I-1,1,1) * QM(J,K,2) + 
     1                 S(J,K,I-1,2,1) * QM(J,K,3) +
     2                 S(J,K,I-1,3,1) * QM(J,K,4)
C
               RES(J,K,I,2) = RES(J,K,I,2) + (OMEGA - 1.E0) *
     1                        ( PPI*S(J,K,I,1,1)*S(J,K,I,4,1)
     2                        - PMI*S(J,K,I-1,1,1)*S(J,K,I-1,4,1)
     3                        + PPJ*S(J,K,I,1,2)*S(J,K,I,4,2)
     4                        - PMJ*S(J-1,K,I,1,2)*S(J-1,K,I,4,2)
     5                        + PPK*S(J,K,I,1,3)*S(J,K,I,4,3)
     6                        - PMK*S(J,K-1,I,1,3)*S(J,K-1,I,4,3) )
c              RES(J,K,I,3) = RES(J,K,I,3) + (OMEGA - 1.E0) *
c    1                        ( PPI*S(J,K,I,2,1)*S(J,K,I,4,1)
c    2                        - PMI*S(J,K,I-1,2,1)*S(J,K,I-1,4,1)
c    3                        + PPJ*S(J,K,I,2,2)*S(J,K,I,4,2)
c    4                        - PMJ*S(J-1,K,I,2,2)*S(J-1,K,I,4,2)
c    5                        + PPK*S(J,K,I,2,3)*S(J,K,I,4,3)
c    6                        - PMK*S(J,K-1,I,2,3)*S(J,K-1,I,4,3) )
c              RES(J,K,I,4) = RES(J,K,I,4) + (OMEGA - 1.E0) *
c    1                        ( PPI*S(J,K,I,3,1)*S(J,K,I,4,1)
c    2                        - PMI*S(J,K,I-1,3,1)*S(J,K,I-1,4,1)
c    3                        + PPJ*S(J,K,I,3,2)*S(J,K,I,4,2)
c    4                        - PMJ*S(J-1,K,I,3,2)*S(J-1,K,I,4,2)
c    5                        + PPK*S(J,K,I,3,3)*S(J,K,I,4,3)
c    6                        - PMK*S(J,K-1,I,3,3)*S(J,K-1,I,4,3) )
  280       CONTINUE
  290    CONTINUE
      ENDIF
C
C     Finished with i direction residual computation for marching.
C
      RETURN
C
 1000 FORMAT (' ','Routine FLUXM - invalid IPSTRT not equal IPEND ',2I6,
     1       /' ','RUN STOPPING')
 1200 FORMAT (' ',' Catastrophic limiting for I direction ',I8,' times')
C
      END







      SUBROUTINE PRPRTY (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS, IMODEL, IFDELQ,
     2                   T, QC, PROPSC, VOLC, DQDX, DQDY, DQDZ, WORKX)
C
C Routine to calculate the PROPERTIES of the flowfield.
C This routine is used for either the time dependent calculation 
C or the marching calculation.  
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM+1)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IMODEL         : Turbulence model array
C IFDELQ         : Array to control calculation of delq
C T              : Temporary storage for the Temperature
C QC,PROPSC,VOLC : Q, PROPS and VOL in one-dimensional array
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C WORKX          : Work array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (0:JDIM+2,0:KDIM+2),
     1          QC    (1:JDIM+1,1:KDIM+1,NQ),
     2          PROPSC(1:JDIM+1,1:KDIM+1,NP),
     3          VOLC  (1:JDIM+1,1:KDIM+1),
     4          DQDX  (1:JDIM+1,1:KDIM+1,NQ),
     5          DQDY  (1:JDIM+1,1:KDIM+1,NQ),
     6          DQDZ  (1:JDIM+1,1:KDIM+1,NQ)
C
      DIMENSION WORKX (MXSECT,NW)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      ITURB = IMODEL(1)
      ISGS  = IMODEL(6)
C
C     Loop over I 
C
      DO 300 I = IPSTRT, IPEND
C
C     Calculating the molecular viscosity
C
         DO 220 K = 0, KDIM + 2
            DO 210 J = 0, JDIM + 2
               T(J,K) = GAMMA * Q(J,K,I,5) / Q(J,K,I,1)
  210       CONTINUE
  220    CONTINUE
C
C     Calculate the molecular viscosity and store it into the PROPS array
C
         NPTS = (JDIM + 3) * (KDIM + 3)
         CALL VSCSTY (NPTS, T, PROPS(0,0,I,1))
  300 CONTINUE
C
C     Calculate the values of Y+
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW .OR. 
     1    ITURB .EQ. ITRS .OR. ITURB .EQ. ITLES) THEN
         CALL YPLUSS (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                Q, PROPS, S, VOL, IBCDAT, NBCS)
      ENDIF
C
C     Calculate the Turbulent Eddy Viscosity
C
C     Loop over I
C
      DO 600 I = IPSTRT, IPEND
C
C     Baldwin-Lomax Turbulence Model
C     Update turbulent eddy viscosity with algebraic model if algebraic model 
C     is specified or if NITALG < ITER such that the algebraic model is used
C     in the initial transient to provide a better convergence for the
C     higher order models.
C
C     Only call Baldwin-Lomax for 2 < i < idim-1
C     PROPS(2) is RMUT
C     Temporary storage: Use WORKX(1)     for VORT
C                            WORKX(2-5)   for TVIS
C                            WORKX(6-9)   for DIST
C                            WORKX(10)    for RHO
C                            WORKX(11)    for VORTJ
C                            WORKX(12)    for RMUTJ
C                            WORKX(13)    for Y
C                            WORKX(14)    for FUN
C
         IF (MUTALG .OR. ITURB .EQ. ITBLMX) THEN
            IF (I .GE. 2 .AND. I .LE. IDIM) THEN
               MAXDIM = MAX (JDIM, KDIM)
               CALL TURB (IDIM, JDIM, KDIM, MAXDIM, I, ITURB,
     1                    Q, PROPS(0,0,0,2), S, VOL, IBCDAT, NBCS,
     2                    WORKX(1, 1), WORKX(1, 2), WORKX(1, 6),
     3                    WORKX(1,10), WORKX(1,11), WORKX(1,12), 
     4                    WORKX(1,13), WORKX(1,14))
            ENDIF
C
C     k-epsilon, Reynolds stress turbulence model
C     Temporary storage: Use WORKX(1)  for TKE
C                            WORKX(2)  for FDAMP
C
         ELSE IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITRS) THEN
C
            DO 330 L = 1, NQ
               DO 320 K = 1, KDIM+1
                  DO 310 J = 1, JDIM+1
                     QC(J,K,L) = Q(J,K,I,L)
  310             CONTINUE
  320          CONTINUE
  330       CONTINUE
C
            DO 360 L = 1, NP
               DO 350 K = 1, KDIM+1
                  DO 340 J = 1, JDIM+1
                     PROPSC(J,K,L) = PROPS(J,K,I,L)
  340             CONTINUE
  350          CONTINUE
  360       CONTINUE
C
            NPTS  = (JDIM + 1) * (KDIM + 1)
            CALL GETTKE (ITURB,  NPTS, QC, WORKX(1,1))
            CALL RMUTKE (IMODEL, NPTS, QC, PROPSC,WORKX(1,1),WORKX(1,2))
C
            DO 380 K = 1, KDIM + 1
               DO 370 J = 1, JDIM + 1
                  PROPS(J,K,I,2) = PROPSC(J,K,2)
  370          CONTINUE
  380       CONTINUE
C
C     k-omega turbulence model
C     PROPS(1) is RMU and PROPS(2) is RMUT
C
         ELSE IF (ITURB .EQ. ITKW) THEN
            CALL RMUTKW (IDIM, JDIM, KDIM, I, Q, PROPS(0,0,0,2)) 
C
C     Large-Eddy Simulation
C     Temporary storage: WORKX    for FDAMP
C
         ELSE IF (ITURB .EQ. ITLES) THEN
C
            DO 470 K = 1, KDIM+1
               DO 465 J = 1, JDIM+1
                  VOLC(J,K) = VOL(J,K,I)
  465          CONTINUE
  470       CONTINUE
C
            NPTS  = (JDIM + 1) * (KDIM + 1)
C
C     Temporary variables used internally in DELQC:
C         WORKX(1-2) for QJ
C         WORKX(3-4) for QK
C         WORKX(5-6) for QPI
C         WORKX(7-8) for QMI
C     QJ is dimensioned one larger than normal routines therefore assign it
C     2 planes of WORKX. Similarly for QK. QPI and QMI require only on plane
C     as they are dimensioned normally - but they are allocated one additional
C     plane as a precaution.
C
C     Calculate derivatives at interior cells and boundaries.
C
            NBOUND = 1
            CALL DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  WORKX(1,1), WORKX(1,3), WORKX(1,5), WORKX(1,7),
     3                  QC, PROPSC, DQDX, DQDY, DQDZ)
C
C     Calculate Structure function for SF model
C     Use WORKX(1) to store SF
C
            IF (ISGS .EQ. ISGSFM) THEN
               CALL GETSF2 (IMODEL, IDIM, JDIM, KDIM, I, Q, PROPS, S, 
     1                      VOL, IBCDAT, NBCS, WORKX(1,1), IERROR)
            ENDIF
C
C     Calculate SGS coefficient
C
            CALL RMUSGS (IMODEL, NPTS, QC, PROPSC, VOLC, DQDX, DQDY, 
     1                   DQDZ, WORKX(1,1), WORKX(1,2), IERROR)
            IF (IERROR .NE. INO) THEN
               STOP
            ENDIF
C
            DO 490 K = 1, KDIM + 1
               DO 480 J = 1, JDIM + 1
                  PROPS(J,K,I,2) = PROPSC(J,K,2)
  480          CONTINUE
  490       CONTINUE
C
C     No turbulence - set RMUT to zero
C
         ELSE
            DO 520 K = 0, KDIM + 2
               DO 510 J = 0, JDIM + 2
                  PROPS(J,K,I,2) = 0.E0
  510          CONTINUE
  520       CONTINUE
         ENDIF
C
  600 CONTINUE
C
C     Finished with the PROPERTIES
C
      RETURN
      END






      SUBROUTINE RESID (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, QN, NQNPTS, PROPS, NPPTS, S, VOL, RC,DTJ,RES, 
     2                  FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,WORKX,
     3                  IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX,VISCOS, 
     4                  IMODEL, INITRS, IFDELQ, MRCHNG, 
     5                  VGNRN, SIGMA, NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
C
C Routine to calculate the residual throughout the flowfield for either 
C the time dependent calculation or the marching calculation.  
C The input logical INITRS determines whether to initialize the residual
C to zero or to not initialize the residual. INITRS = false is used when
C the residual is initialized in the multigrid cycle to the coarse grid
C forcing function.
C
C Time dependent calculation: IPSTRT = 2, IPEND = IDIM
C Marching calculation      : IPSTRT = IPEND = I plane presently at
C
C This routine is organized as follows:
C     1. If INITRS then zero residuals
C     2. If iterative time dependent then add dQ/dt term
C     3. Residual in I direction - Add dF/dx term
C     4. Residual in J direction - Add dG/dy term
C     5. If 3D or axisymmetric then residual in K direction - add dH/dz term
C     6. Calculate Source Terms
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
C Q              : Primitive variables at cell centers
C QN             : Conserved variables at time levels n, n-1, etc.
C NQNPTS         : Number of QN points. This allows avoiding storing QN when
C                  running time accurate.
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
C DTJ            : Time step divided by cell volume
C RES            : Residual
C FHAT           : Storage for fluxes for 4th order
C DQP,DQM        : Storage for delta q = q(i) - q(i-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses to improve modularity
C F,FM           : Fluxes at cell interfaces -
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Temporary storage for metrics
C WORKX          : Work array for cross-section (used for damping, etc.)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for current direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence modeling scheme array
C INITRS         : Logical controlling initializing RES to zero for MG
C IFDELQ         : Array controlling calculation of delq
C MRCHNG         : Logical for Marching
C VGNRN          : Logical for Vigneron approach in marching flux evaluation
C SIGMA          : Safety factor for Vigneron technique when marching
C NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C IFL2DQ         : Logical to calculate L2 norm of dQ/dt term for unsteady
C R2DQDT         : L2 norm of dQ/dt term for unsteady
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (NQNPTS,NQ,NTMLVL),
     2          PROPS (NPPTS,NP),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     5          RC    (JDIM+1,KDIM+1,IDIM+1,3),
     6          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          FHAT  (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     2          DQP   (MXSECT,NQ),
     3          DQM   (MXSECT,NQ),
     4          QP    (MXSECT,NQ), 
     5          QM    (MXSECT,NQ), 
     6          TAU   (MXSECT,6), 
     7          F     (MXSECT,NF),
     8          FM    (MXSECT,NF),
     9          STEMP (MXSECT,3),
     A          WORKX (MXSECT,NW)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IORDER(3), LIMIT(3), RKAP(3)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      DIMENSION R2DQDT(NF+1)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
      LOGICAL   VISCOS(3), INITRS, MRCHNG, VGNRN, IFL2DQ
C
      ITURB = IMODEL(1)
C
C 1.  If INITRS then zero residuals
C
      IF (INITRS) THEN
         DO 140 L = 1, NF
            DO 130 I = IPSTRT - 1, IPEND + 1
               DO 120 K = 1, KDIM + 1
                  DO 110 J = 1, JDIM + 1
                     RES(J,K,I,L) = 0.E0
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
      ENDIF
C
C 2.  If iterative time dependent then add dQ/dt term and
C     calculate L2 Norm of dQ/dt residual
C
      IF (ITIMED .EQ. IUNSTD) THEN
         CALL DQDT   (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, QN, VOL, DTJ, RES, QP, QM, IMODEL)
         IF (IFL2DQ) THEN
            CALL L2NORM (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                   RES, VOL, R2DQDT)

         ENDIF
      ENDIF
C
C 3.  Residual in I direction - Add dF/dx term
C         FLUXM if marching calculation
C         FLUXI if time-dependent calculation
C
      IF ( MRCHNG ) THEN
         CALL FLUXM (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1               Q, PROPS, NPPTS, S, VOL, 
     2               RES, DQP, DQM, QP, QM, F, FM, STEMP,
     3               IBCDAT, NBCS, IORDER(1), LIMIT(1), 
     4               VGNRN, SIGMA, IMODEL)
      ELSE
         IF (FOURTH) THEN
            CALL FLUXI4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(1), LIMIT(1), RKAP(1), 
     4                   IFLUX, VISCOS(1), IMODEL)
         ELSE
            CALL FLUXI  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(1), LIMIT(1), RKAP(1), 
     4                   IFLUX, VISCOS(1), IMODEL)
         ENDIF
      ENDIF
C
C 4.  Residual in J direction - Add dG/dy term
C
      IF (FOURTH) THEN
         CALL FLUXJ4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(2), LIMIT(2), RKAP(2), 
     4                IFLUX, VISCOS(2), IMODEL)
      ELSE
         CALL FLUXJ  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(2), LIMIT(2), RKAP(2), 
     4                IFLUX, VISCOS(2), IMODEL)
      ENDIF
C
C 5.  If 3D or axisymmetric then Residual in K direction - add dH/dz term
C
      IF ( THREED ) THEN
         IF (FOURTH) THEN
            CALL FLUXK4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                   IFLUX, VISCOS(3), IMODEL)
         ELSE
            CALL FLUXK  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                   IFLUX, VISCOS(3), IMODEL)
         ENDIF
      ELSE IF ( AXISYM ) THEN
         CALL FLUXKA (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                IFLUX, VISCOS(3), IMODEL)
      ENDIF
C
C 6.  Calculate Source Terms
C     Temporary variables passed to the source term evaluation routine:
C         QP     for QC (variables in one-dimensional array)
C         FM(2-) for PROPSC (properties in one-dimensional array - this is 
C                            passed in using FM.  Must always be sure that 
C                            NP+1 <= NF or else will need to fix this)
C         STEMP        for SN   (Wall normals - needed in Reynolds stress model)
C         QM           for DQDX (derivatives of Q wrt x at cell centers)
C         DQP          for DQDY (derivatives of Q wrt y at cell centers)
C         DQM          for DQDZ (derivatives of Q wrt z at cell centers)
C         FM(1)        for DEL2K (Laplacian of turbulent kinetic energy)
C         WORKX(1)     for TKE
C         WORKX(2)     for ETILDE
C         WORKX(3)     for XSIW
C         WORKX(4-5)   for FDAMP
C         WORKX(6)     for EPSC
C         WORKX(7)     for PDIL
C         WORKX( 8-13) for SIJ
C         WORKX(14-19) for WIJ
C
      IF ( ISRCE ) THEN
         IF (ITURB.EQ.ITKE .OR. ITURB.EQ.ITKW .OR. ITURB.EQ.ITRS) THEN
            IF (NF .LT. (NP + 1)) THEN
               WRITE (IOUT,1000) NP+1, NF
               STOP
            ENDIF
C
            CALL SRCRES (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                   Q, PROPS, NPPTS, S, VOL, RC, RES,
     2                   QP, FM(1,2), STEMP, QM, DQP, DQM, FM(1,1), 
     3                   WORKX(1,1), WORKX(1,2), WORKX(1,3),WORKX(1,4),
     4                   WORKX(1,6), WORKX(1,7), WORKX(1,8),WORKX(1,14),
     5                   TAU, F, IBCDAT, NBCS, IFLUX, IMODEL, IFDELQ,
     6                   NPRLIM, NPRNEG)
         ENDIF
      ENDIF
C
C Wall Functions: reset residuals to zero on the turbulence quantities that
C are set from the wall functions (epsilon is set from wall functions, all
C other turbulence quantities, tau_ij or k, are calculated at cell center).
C The assumption in the wall functions is that we replace the calculated
C value of the turbulence quantities with an 'exact' value.  This
C then implies that the residual at these points is zero.  Necessary as
C the implicit time integration couples the mean and turbulent variables.
C
      IF ( IFWALF ) THEN
C
C Set up equation indices to reset residual to zero
C
         IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
            LSTRT = 6
            LEND  = IEPSLN
         ELSE IF (ITURB .EQ. ITRS) THEN
            LSTRT = IEPSLN
            LEND  = IEPSLN
         ELSE
            LSTRT = 1
            LEND  = 0
         ENDIF
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM
         IJKDIM(2) = JDIM
         IJKDIM(3) = KDIM
C
C Loop over all boundary conditions
C
         DO 490 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IF (IBCTYP .EQ. IWALFN) THEN
               IDIR      = IBCDAT(2,IBC)
               INROUT    = IBCDAT(3,IBC)
               ISTRT(1)  = MAX (IBCDAT(4,IBC), IPSTRT)
               IEND (1)  = MIN (IBCDAT(5,IBC), IPEND)
               ISTRT(2)  = IBCDAT(6,IBC)
               IEND (2)  = IBCDAT(7,IBC)
               ISTRT(3)  = IBCDAT(8,IBC)
               IEND (3)  = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  IBCN = 2
               ELSE
                  IBCN = IJKDIM(IDIR)
               ENDIF
C
               IC1 = ICY (IDIR, 1)
               IC2 = ICY (IDIR, 2)
               IC3 = ICY (IDIR, 3)
C
               DO 430 L = LSTRT, LEND
                  DO 420 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                     DO 410 I2 = ISTRT(IC2), IEND(IC2)
C
                        IN = IBCN*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                        JN = IBCN*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                        KN = IBCN*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                        RES(JN,KN,IN,L) = 0.E0
C
  410                CONTINUE
  420             CONTINUE
  430          CONTINUE
            ENDIF
  490    CONTINUE
      END IF
C
C     Finished with residual computation.
C
      RETURN
 1000 FORMAT (' ','***** FATAL ERROR - STOP in routine RESID ',
     1       /' ','Temporary storage required to evaluate source terms',
     2       /' ','is insufficient.  Contact ISAAC developer to fix.',
     3       /' ','NP+1 = ', I5, ' NF = ', I5)
      END






      SUBROUTINE SRCRES (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, NPPTS, S, VOL, RC, RES,
     2                   QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, NPRLIM, NPRNEG)
C
C     Routine to calculate the contriubtion to the residual from the
C     source terms.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RC             : Position vector at cell centers
C     RES            : Residual
C     QC             : Storage for Q in one-dimensional array
C     PROPSC         : Storage for PROPS array in one-dimensional array
C     SN             : Storage for wall normals in one-dimensional array
C     DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C     DEL2K          : Laplacian of turbulent kinetic energy
C     TKE            :
C     ETILDE         :
C     XSIW           :
C     FDAMP          :
C     EPSC           :
C     PDIL           :
C     SIJ            :
C     WIJ            :
C     TAU            : Storage for Reynolds stress terms
C     SRC            : Source terms for cell (J,K)
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     IMODEL         : Turbulence model array
C     IFDELQ         : Array controlling calculation of delq
C     NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          QC    (2:JDIM,2:KDIM,NQ), 
     2          PROPSC(2:JDIM,2:KDIM,NP),
     3          SN    (2:JDIM,2:KDIM,3),
     4          DQDX  (2:JDIM,2:KDIM,NQ),
     5          DQDY  (2:JDIM,2:KDIM,NQ), 
     6          DQDZ  (2:JDIM,2:KDIM,NQ),
     7          DEL2K (2:JDIM,2:KDIM),
     8          TKE   (2:JDIM,2:KDIM),
     9          ETILDE(2:JDIM,2:KDIM),
     A          XSIW  (2:JDIM,2:KDIM),
     B          FDAMP (2:JDIM,2:KDIM,2),
     C          EPSC  (2:JDIM,2:KDIM),
     D          PDIL  (2:JDIM,2:KDIM),
     E          SIJ   (2:JDIM,2:KDIM,6),
     F          WIJ   (2:JDIM,2:KDIM,6),
     G          TAU   (2:JDIM,2:KDIM,6),
     H          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
C     Local variable:
C
      LOGICAL   ISAVPR
C
      ITURB  = IMODEL( 1)
      IFDL2K = IMODEL( 5)
      IASM   = IMODEL(10)
C
C     Set up control to store CMUSTAR back into global array if 
C     Algebraic Reynolds stress model
C
c-    Current EASM model for GS does not include variable C_mu for diffusion
c-    IF (IASM .EQ. IASMGS) THEN
c-       ISAVPR = .TRUE.
c-    ELSE
         ISAVPR = .FALSE.
c-    ENDIF
C
C     Number of points in a cross-section to calculate source terms at
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C     Loop over I calculating the source terms
C
      DO 300 I = IPSTRT, IPEND
C
         CALL SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                Q, PROPS, NPPTS, S, VOL, RC,
     2                QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     4                TAU, SRC, IBCDAT, NBCS, 
     5                IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C     Add the source term to the residual - include volume
C
         DO 250 L = 1, NF
            DO 240 K = 2, KDIM
               DO 230 J = 2, JDIM
                  RES(J,K,I,L) = RES(J,K,I,L) - SRC(J,K,L)*VOL(J,K,I)
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
  300 CONTINUE
C
C     Finished adding source terms to residual
C
      RETURN
      END





      SUBROUTINE SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C Routine to calculate the contriubtion the source terms.  This routine is 
C used for either the time dependent calculation or the marching calculation.  
C This routine is used for explicit and implicit source term treatment.
C
C This routine is organized as follows:
C      1. Calculate the Laplacian of the turbulent kinetic energy
C      2. Calculate derivatives of Q at cell centers
C      3. Calculate turbulent kinetic energy in 1D array
C      4. Calculate the Reynold's stresses at the cell centers.  
C      5. Store CMU* from the local array back into the global storage
C      6. Calculate S_ij, W_ij
C      7. Call routines to calculate turbulence model source terms
C         a. K-epsilon
C         b. K-omega
C         c. Reynolds stress
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane (cell center) to calculate source terms at
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
C RES            : Residual
C QC             : Storage for Q in one-dimensional array
C PROPSC         : Storage for PROPS array in one-dimensional array
C SN             : Storage for wall normals in one-dimensional array
C DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C DEL2K          : Laplacian of turbulent kinetic energy
C TKE            :
C ETILDE         :
C XSIW           :
C FDAMP          :
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C TAU            : Storage for Reynolds stress terms
C SRC            : Source terms for cell (J,K)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IFLUX          : Flux evaluation scheme
C IMODEL         : Turbulence model array
C IFDELQ         : Array controlling calculation of delq
C ISAVPR         : Logical controlling storage of CMUSTAR from local
C                  array where it is calculated in TAUIJ back into
C                  global array.  This is done only when called from RESID.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION QC    (2:JDIM,2:KDIM,NQ), 
     1          PROPSC(2:JDIM,2:KDIM,NP),
     2          SN    (2:JDIM,2:KDIM,3),
     3          DQDX  (2:JDIM,2:KDIM,NQ),
     4          DQDY  (2:JDIM,2:KDIM,NQ), 
     5          DQDZ  (2:JDIM,2:KDIM,NQ),
     6          DEL2K (2:JDIM,2:KDIM),
     7          TKE   (2:JDIM,2:KDIM),
     8          ETILDE(2:JDIM,2:KDIM),
     9          XSIW  (2:JDIM,2:KDIM),
     A          FDAMP (2:JDIM,2:KDIM,2),
     B          EPSC  (2:JDIM,2:KDIM),
     C          PDIL  (2:JDIM,2:KDIM),
     D          SIJ   (2:JDIM,2:KDIM,6),
     E          WIJ   (2:JDIM,2:KDIM,6),
     F          TAU   (2:JDIM,2:KDIM,6),
     G          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      LOGICAL   ISAVPR
C
      ITURB  = IMODEL( 1)
      IFDL2K = IMODEL( 5)
      IASM   = IMODEL(10)
      IEPSC  = IMODEL(11)
      IPDIL  = IMODEL(12)
      ISRCW  = IMODEL(13)
C
C     Number of points in a cross-section to calculate source terms at
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C 1.  Calculate the Laplacian of the turbulent kinetic energy
C     Temporary variables used internally in DEL2KC:
C         TAU(1-3) for QJ  
C         SRC(1-3) for QK  
C         TAU(5->) for QPI 
C         SRC(4-6) for QMI 
C     QJ is dimensioned as JDIM,2:KDIM therefore need more room for it 
C     than 1 plane of TAU. Similarly for QK.  QPI and QMI require only 
C     one plane as they are dimensioned 2:JDIM,2:KDIM.
C     Use SRC for some of the temporary space to fix a memory allocation
C     problem.  Use DQDX,DQDY,DQDZ as temporary space here as it isn't
C     calculated until after this call.
C
      IF (IFDL2K .EQ. 1) THEN
         CALL DEL2KC (IDIM, JDIM, KDIM, I, ITURB,
     1                Q, S, VOL, IBCDAT, NBCS,
     2                TAU(2,2,1), SRC(2,2,1), TAU(2,2,5), SRC(2,2,4),
     3                DQDX, DQDY, DQDZ, DEL2K)
      ENDIF
C
C 2.  Calculate derivatives of Q at cell centers.
C     Temporary variables used internally in DELQC:
C         TAU(1-3) for QJ  
C         SRC(1-3) for QK  
C         TAU(5->) for QPI 
C         SRC(4-6) for QMI 
C     QJ is dimensioned as JDIM,2:KDIM therefore need more room for it 
C     than 1 plane of TAU. Similarly for QK.  QPI and QMI require only 
C     one plane as they are dimensioned 2:JDIM,2:KDIM.
C     Use SRC for some of the temporary space to fix a memory allocation
C     problem.
C
C     Calculate derivatives at interior cell centers only.
C
      NBOUND = 0
      CALL DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1            Q, PROPS, S, VOL, IBCDAT, NBCS,
     2            TAU(2,2,1), SRC(2,2,1), TAU(2,2,5), SRC(2,2,4),
     3            QC, PROPSC, DQDX, DQDY, DQDZ)
C
C 3.  Calculate turbulent kinetic energy in 1D array
C
      CALL GETTKE (ITURB, NPTS, QC, TKE)
C
C 4.  Calculate the Reynold's stresses at the cell centers.  
C                                      ~
C     The array TAU contains {- rho (u" u")}.
C                                     i  j
C     For laminar flow:            TAU_ij => zero.
C         algebraic models:        TAU_ij => modeled with Boussinesq approx.
C         two-equation models:     TAU_ij => modeled with Boussinesq approx.
C         Reynold's stress models: TAU_ij => computed from PDE.
C     Pass in RMUT and DUDX,DUDY,DUDZ for Boussinesq approximation.
C     Pass in QC for Q data at same location as DQDX,... derivative data.
C     The Reynold's stresses are stored as:
C         TAU(1) = TAUXX, TAU(2) = TAUYY, TAU(3) = TAUZZ,
C         TAU(4) = TAUXY, TAU(5) = TAUXZ, TAU(6) = TAUYZ
C
      CALL TAUIJ (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU)
C
C 5.  Store the value of CMU* from the local array back into the global
C     array if using Algebraic Reynolds stress model and ISAVPR is true
C
      IF (IASM .EQ. IASMGS  .AND.  ISAVPR) THEN
         CALL SAVEPR (IDIM, JDIM, KDIM, I, PROPSC, PROPS)
      ENDIF
C
C 6.  Calculate S_ij, W_ij
C
      CALL GETSW  (NPTS, QC, DQDX, DQDY, DQDZ, SIJ, WIJ)
C
C 7.  Turbulence model source terms
C
C 7-a. k-epsilon turbulence model
C
      IF (ITURB .EQ. ITKE) THEN
         CALL SRCKE  (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, 
     1                DEL2K, TKE, ETILDE, XSIW, FDAMP, 
     2                EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C 7-b. k-omega turbulence model
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         CALL SRCKW  (NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, SRC, 
     1                NPRLIM, NPRNEG)
C
C 7-c. Reynolds Stress turbulence model
C      SRCRS  calculates the source terms
C      SRCRSW calculates additional wall normal dependent source terms
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         CALL SRCRS  (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, 
     1                TAU, DEL2K, TKE, ETILDE, XSIW, FDAMP, 
     2                EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C     Calculate wall normal dependent source terms for models that include them
C     Temporary storage: TAU    for SRCW
C                        ETILDE for YNI
C                        XSIW   for YNITOT
C
         IF (ISRCW .EQ. IYES) THEN
            CALL SRCRSW (IMODEL, QC, PROPSC, SN, DQDX, DQDY, DQDZ, 
     1                   TKE, EPSC, SRC, TAU, 
     2                   IDIM, JDIM, KDIM, I, NBCS, IBCDAT,
     3                   S, RC, ETILDE, XSIW)
         ENDIF
      ENDIF
C
C     Finished with the source terms at I
C
      RETURN
      END




      SUBROUTINE SAVEPR (IDIM, JDIM, KDIM, I, PROPSC, PROPS)
C
C Routine to store CMUSTAR calculated in TAUIJ from local array in TAUIJ
C back into global array for the rest of the program.
C
C IDIM,JDIM,KDIM : Dimensions of block
C I              : Current I = const. plane
C PROPSC         : Local array for Properties (stored on an I = const. plane)
C PROPS          : Global array for Properties
C
      include '../header/common.h'
C
      DIMENSION PROPSC(2:JDIM,2:KDIM,NP),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Copy data from temporary storage into permanent PROPS array
C
      DO 120 K = 2, KDIM
         DO 110 J = 2, JDIM
            PROPS(J,K,I,ICMUST) = PROPSC(J,K,ICMUST)
  110    CONTINUE
  120 CONTINUE
C
C Do boundary conditions so that Cmu* is consistently defined on boundaries.
C
      DO 210 K = 2, KDIM
         PROPS(1     ,K,I,ICMUST) = PROPS(2   ,K,I,ICMUST)
         PROPS(JDIM+1,K,I,ICMUST) = PROPS(JDIM,K,I,ICMUST)
  210 CONTINUE
C
      DO 310 J = 1, JDIM+1
         PROPS(J,1     ,I,ICMUST) = PROPS(J,2   ,I,ICMUST)
         PROPS(J,KDIM+1,I,ICMUST) = PROPS(J,KDIM,I,ICMUST)
  310 CONTINUE
C
C If I = 2   , then do I = 1      boundary
C If I = IDIM, then do I = IDIM+1 boundary
C
      IF (I .EQ. 2) THEN
         DO 420 K = 1, KDIM+1
            DO 410 J = 1, JDIM+1
               PROPS(J,K,1     ,ICMUST) = PROPS(J,K,2   ,ICMUST)
  410       CONTINUE
  420    CONTINUE
C
      ELSE IF (I .EQ. IDIM) THEN
         DO 440 K = 1, KDIM+1
            DO 430 J = 1, JDIM+1
               PROPS(J,K,IDIM+1,ICMUST) = PROPS(J,K,IDIM,ICMUST)
  430       CONTINUE
  440    CONTINUE
      ENDIF
C
      RETURN
      END





      SUBROUTINE VSCSTY (NPTS, T, RMU)
C
C     Routine to calculate the molecular viscosity.
C
C
C     NPTS           : Number of points at which to calculate the molecular
C                      viscosity
C     T              : Temperature (T/Tinf)
C     RMU            : Molecular Viscosity
C
      include '../header/common.h'
C
      DIMENSION T     (NPTS),
     1          RMU   (NPTS)
C
C     Sutherland's Law
C
C     Constants for Sutherland's Law
C     Sref = 198.6 degrees Rankine.  If using SI units, correct to Kelvins.
C
      SREF  = 198.6E0
      IF ( SIUNIT ) SREF = SREF / 1.8E0
      STINF = SREF / TINF
      SNUM  = 1.E0 + STINF
C
C     Calculate molecular viscosity
C
      DO 100 I = 1, NPTS
         RMU(I) = T(I) * SQRT( T(I) ) * SNUM / ( T(I) + STINF )
  100 CONTINUE
C
C     Finished calculating molecular viscosity
C
      RETURN
      END





      SUBROUTINE YPLUSS (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS)
C
C Routine to calculate Y+ in the flowfield.
C This routine is used for either the time dependent calculation 
C or the marching calculation.  
C
C This routine is organized as follows:
C     1. Initialize Y+ = infinity 
C     2. Loop over all boundary conditions - for wall/wall functions do
C        a. Calculate Y+ = min (Y+, Y+ calculated from current wall)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT,IPEND   : I plane to begin,end calculations at
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 1. Initialize Y+ to <infinity> in case there are no walls
C
      DO 130 I = IPSTRT, IPEND
         DO 120 K = 1, KDIM + 1
            DO 110 J = 1, JDIM + 1
               PROPS(J,K,I,3) = RMAX
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C 2. Loop over WALL boundary conditions
C    Calculate Y+ = MIN (y+, Y+) to allow for multiple walls in
C    the same coordinate direction.
C
      DO 500 LBC = 1, NBCS
         IBCTYP = IBCDAT(1,LBC)
         IF (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN) THEN
            IDIR      = IBCDAT(2,LBC)
            INROUT    = IBCDAT(3,LBC)
            ISTRT(1)  = MAX (IBCDAT(4,LBC), IPSTRT)
            IEND (1)  = MIN (IBCDAT(5,LBC), IPEND)
            ISTRT(2)  = IBCDAT(6,LBC)
            IEND (2)  = IBCDAT(7,LBC)
            ISTRT(3)  = IBCDAT(8,LBC)
            IEND (3)  = IBCDAT(9,LBC)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 210 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
  210       CONTINUE
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
            DO 320 I3 = ISTRT(IC3), IEND(IC3)
               DO 310 I2 = ISTRT(IC2), IEND(IC2)
C
                  IBC  = IBCPT *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC  = IBCPT *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC  = IBCPT *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IBC2 = IBCPT2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC2 = IBCPT2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC2 = IBCPT2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IN   = IBCN  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN   = IBCN  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN   = IBCN  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IM   = IBCM  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JM   = IBCM  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KM   = IBCM  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  RHOW   = Q    (JBC,KBC,IBC,1)
                  RMUW   = PROPS(JBC,KBC,IBC,1)
                  YNORMW = ABS (PROPS(JBC,KBC,IBC,4))
C
                  YNORMN = ABS (PROPS(JN ,KN ,IN ,4))
C
C Calculate wall shear or use value from wall functions
C
                  IF ( IFWALF ) THEN
                     TAUW   = ABS (PROPS(JBC2,KBC2,IBC2,ITQTAU))
                  ELSE
C
C Calculate the wall value of du/dn
C
                     UBARN = S(JM,KM,IM,1,IDIR) * Q(JN,KN,IN,2)
     1                     + S(JM,KM,IM,2,IDIR) * Q(JN,KN,IN,3)
     2                     + S(JM,KM,IM,3,IDIR) * Q(JN,KN,IN,4)
                     VELN = SQRT ( 
     1                     (Q(JN,KN,IN,2)-UBARN*S(JM,KM,IM,1,IDIR))**2
     2                   + (Q(JN,KN,IN,3)-UBARN*S(JM,KM,IM,2,IDIR))**2
     3                   + (Q(JN,KN,IN,4)-UBARN*S(JM,KM,IM,3,IDIR))**2 )
C
C Calculate first derivative term as du/dn = delta(u_parallel) / delta(normal)
C
                     DVELDN = VELN / (ABS (YNORMN - YNORMW) + RSMALL)
C
C TAUW = mu * du/dy |(wall)
C
                     TAUW   = RMUW * DVELDN
                  ENDIF
C
C CTAUW = constants for non-dimensionalization * ustar
C         ustar = sqrt (Tau_wall / rho_wall)
C
                  CTAUW  = SQRT ( RE / FSMACH * TAUW / RHOW )
C
C Loop normal to the wall setting Y+
C
CDIR$ IVDEP
                  DO 300 I1 = 1, IJKDIM(IDIR) + 1
C
                     II = I1*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
                     JJ = I1*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
                     KK = I1*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Yplus = y ustar / nu 
C
                     RHO               = Q    (JJ,KK,II,1)
                     RMU               = PROPS(JJ,KK,II,1)
                     YNORML            = ABS (PROPS(JJ,KK,II,4))
                     YPL               = CTAUW * YNORML * RHO / RMU
                     PROPS(JJ,KK,II,3) = MIN (YPL, PROPS(JJ,KK,II,3))
  300             CONTINUE
  310          CONTINUE
  320       CONTINUE
C
         ENDIF
C
C Next boundary condition
C
  500 CONTINUE
C
C Finished calculating Y+
C
      RETURN
      END
