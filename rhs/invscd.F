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
C     $Date: 2001/10/29 03:25:31 $
C     $State: Exp $
C     $Log: invscd.F,v $
C     Revision 4.3  2001/10/29 03:25:31  jhmorr
C     Updated license information
C
C     Revision 4.2  2001/06/08 04:56:16  jhmorr
C     Added notice
C
C     Revision 4.1  1998/04/14 20:34:58  jhmorr
C     Alpha 4.1.
C
C





      SUBROUTINE FDS (IFLUX, ITURB, NPTS, QP, QM, S, DQP, DQM, F)
C
C     Routine to calculate interface fluxes by calling the 
C     appropriate Flux Difference Splitting scheme.
C
C
C     IFLUX          : Flux evaluation scheme
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     S              : Metrics
C     DQP,DQM        : Used for temporary storage in FROE
C     F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION QP    (NPTS,NQ),
     1          QM    (NPTS,NQ),
     2          S     (NPTS,3),
     3          DQP   (NPTS,NQ),
     4          DQM   (NPTS,NQ),
     5          F     (NPTS,NF)
C
C     Call the appropriate Flux Difference Splitting scheme
C     Temporary storage in FROE:
C        DQP(1)  => TKEP
C        DQP(2)  => TKEM
C        DQP(3)  => HP
C        DQP(4)  => HM
C        DQM     => WORK(1-4)
C
      IF (IFLUX .EQ. IFROE) THEN
         CALL FROE (ITURB, NPTS, QP, QM, S, 
     1              DQP(1,1), DQP(1,2), DQP(1,3), DQP(1,4), DQM, F)
      ENDIF
C
C     Finished calculating flux at the interface 
C
      RETURN
      END






      SUBROUTINE FLUX (ITURB, NPTS, Q, S, TKE, H, F)
C
C     Routine to calculate interface fluxes as the full advection flux.
C     This routine calculates the fluxes for the Favre averaged 
C     Navier-Stokes equations for inviscid, laminar, algebraic turbulence 
C     model, two-equation turbulence model or Reynolds stress turbulence 
C     model.  It computes all of the inviscid fluxes for the above models.
C
C     The two-equation and Reynolds stress models need the turbulent kinetic 
C     energy (K) in the mean flow equations.  Therefore, this routine is
C     organized into the following sections:
C          1. Calculate Turbulent Kinetic Energy (TKE)
C             a. Inviscid, laminar, 
C                algebraic model       => K = 0
C             b. Two-Equation model    => K = K (Q(6))
C             c. Reynolds stress model => K = (Tau_xx+Tau_yy+Tau_zz)/2
C          2. Calculate Total Enthalpy (H)
C             a. Perfect gas    => H = gamma P / (RHO (GAMMA-1)) + (q^2)/2 + K
C          3. Calculate Advection Fluxes for Mean Flow Equations
C          4. Calculate additional Advection Fluxes
C             a. Inviscid, laminar, 
C                algebraic model       => No additional fluxes
C             b. Two-Equation model    => Fluxes for K and T2 equations
C             c. Reynolds Stress model => Fluxes for Tau_ij and T2 equations
C
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     Q              : Q at a cell interface to calculate fluxes at
C     S              : Metrics
C     TKE            : Turbulent kinetic energy 
C     H              : Total enthalpy 
C     F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION Q   (NPTS,NQ),
     1          S   (NPTS,3),
     2          TKE (NPTS),
     3          H   (NPTS),
     4          F   (NPTS,NF)
C
C     1. Calculate the Turbulence Kinetic Energy (TKE)
C
C     1-b. Turbulence Kinetic Energy for Two Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKE(I) = Q(I,6)
  100    CONTINUE
C
C     1-c. Turbulence Kinetic Energy for Reynolds Stress Turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
C
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TKE(I) = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
  110    CONTINUE
C
C     1-a. Turbulence Kinetic Energy = 0 for Inviscid, Laminar and
C          Algebraic Turbulence models
C
      ELSE
         DO 120 I = 1, NPTS
            TKE(I) = 0.E0
  120    CONTINUE
      ENDIF
C
C     2. Calculate Total Enthalpy
C
C     2-a. Perfect Gas
C
      DO 200 I = 1, NPTS
C
         RHO   = Q(I,1)
         U     = Q(I,2) 
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
         H(I)  = GAMMA * P / (RHO * GAMM1) 
     1            + 0.5E0 * (U*U + V*V + W*W) + TKE(I)
  200 CONTINUE
C
C     3. Calculate F for Mean Flow Equations at the cell interfaces 
C
      DO 300 I = 1, NPTS
C
         RHO   = Q(I,1)
         U     = Q(I,2) 
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         UBAR  = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
C     Calculate F 
C     Area term is included in the RESID routine.
C
         F(I,1)  = RHO * UBAR 
         F(I,2)  = RHO * UBAR * U + S(I,1) * P
         F(I,3)  = RHO * UBAR * V + S(I,2) * P
         F(I,4)  = RHO * UBAR * W + S(I,3) * P
         F(I,5)  = RHO * UBAR * H(I)
 300  CONTINUE
C
C     4. Calculate Additional Advection Fluxes
C
C     4-a. Inviscid, Laminar, Algebraic Turbulence models => no additional flux
C
C     4-b. Two-Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 400 I = 1, NPTS
C
            RHO    = Q(I,1)
            U      = Q(I,2) 
            V      = Q(I,3)
            W      = Q(I,4)
C
            T2     = Q(I,7)
C
            UBAR   = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            F(I,6)  = RHO * UBAR * TKE(I)
            F(I,7)  = RHO * UBAR * T2 
  400    CONTINUE
C
C     4-c. Reynolds Stress Turbulence models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 410 I = 1, NPTS
C
            RHO    = Q(I,1)
            U      = Q(I,2) 
            V      = Q(I,3)
            W      = Q(I,4)
C
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            T2     = Q(I,12)
C
            UBAR   = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            F(I,6)  = RHO * UBAR * TAUXX
            F(I,7)  = RHO * UBAR * TAUYY
            F(I,8)  = RHO * UBAR * TAUZZ
            F(I,9)  = RHO * UBAR * TAUXY
            F(I,10) = RHO * UBAR * TAUXZ
            F(I,11) = RHO * UBAR * TAUYZ
            F(I,12) = RHO * UBAR * T2
  410    CONTINUE
      ENDIF
C
C     Finished calculating full flux at the interface 
C
      RETURN
      END




      SUBROUTINE FDISS (ITURB, NPTS, Q, DQC, S, ALAM, F)
C
C Routine to calculate interface fluxes as the full advection flux.
C This routine calculates the fluxes for the Favre averaged 
C Navier-Stokes equations for inviscid, laminar, algebraic turbulence 
C model, two-equation turbulence model or Reynolds stress turbulence 
C model.  It computes all of the inviscid fluxes for the above models.
C
C The two-equation and Reynolds stress models need the turbulent kinetic 
C energy (K) in the mean flow equations.  Therefore, this routine is
C organized into the following sections:
C      1. Calculate artificial dissipation for 2nd or 4th order
C         a. Scalar dissipation
C            i.   Mean flow equations: Inviscid, laminar, algebraic model
C            ii.  Two-equation model
C            iii. Reynolds stress model
C
C ITURB          : Turbulence model
C NPTS           : Size of array to calculate fluxes on.
C Q              : Q at a cell interface to calculate fluxes at
C DQC            : Delta Q conserved at interfaces for dissipation
C                  = 2nd order -> third undivided difference
C                    4th order -> fifth undivided difference
C S              : Metrics
C ALAM           : Eigenvalues (Lambda)
C F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION Q   (NPTS,NQ),
     1          DQC (NPTS,NQ),
     2          S   (NPTS,3),
     3          ALAM(NPTS),
     4          F   (NPTS,NF)
C
C 1.   Calculate artificial dissipation for 2nd or 4th order
C 1-a. Scalar dissipation
C 1-a-i.   Mean flow equations: Inviscid, laminar, algebraic model
C
c     IF () THEN
         DO 100 I = 1, NPTS
C
         RHO     = Q(I,1)
         U       = Q(I,2) 
         V       = Q(I,3)
         W       = Q(I,4)
         P       = Q(I,5)
C
         UBAR    = S(I,1)*U + S(I,2)*V + S(I,3)*W
         A       = SQRT (GAMMA * P / RHO)
C
         ALAM(I) = ABS(UBAR) + A
C
C      Calculate F = F - F_diss
C      Area term is included in the RESID routine.
C
            F(I,1)  = F(I,1) - ADKP46 * ALAM(I) * DQC(I,1)
            F(I,2)  = F(I,2) - ADKP46 * ALAM(I) * DQC(I,2)
            F(I,3)  = F(I,3) - ADKP46 * ALAM(I) * DQC(I,3)
            F(I,4)  = F(I,4) - ADKP46 * ALAM(I) * DQC(I,4)
            F(I,5)  = F(I,5) - ADKP46 * ALAM(I) * DQC(I,5)
 100     CONTINUE
C
C 1-a-ii.  Two-equation model
C
         IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
            DO 110 I = 1, NPTS
               F(I,6)  = F(I,6) - ADKP46 * ALAM(I) * DQC(I,6)
               F(I,7)  = F(I,7) - ADKP46 * ALAM(I) * DQC(I,7)
 110        CONTINUE
C
C 1-a-iii. Reynolds stress model
C
         ELSE IF (ITURB .EQ. ITRS) THEN
            DO 120 I = 1, NPTS
               F(I, 6)  = F(I, 6) - ADKP46 * ALAM(I) * DQC(I, 6)
               F(I, 7)  = F(I, 7) - ADKP46 * ALAM(I) * DQC(I, 7)
               F(I, 8)  = F(I, 8) - ADKP46 * ALAM(I) * DQC(I, 8)
               F(I, 9)  = F(I, 9) - ADKP46 * ALAM(I) * DQC(I, 9)
               F(I,10)  = F(I,10) - ADKP46 * ALAM(I) * DQC(I,10)
               F(I,11)  = F(I,11) - ADKP46 * ALAM(I) * DQC(I,11)
               F(I,12)  = F(I,12) - ADKP46 * ALAM(I) * DQC(I,12)
 120        CONTINUE
         ENDIF
c     ELSE
c     ENDIF
C
C     Finished calculating artificial dissipation flux at interface
C
      RETURN
      END





      SUBROUTINE FROE (ITURB, NPTS, QP, QM, S, TKEP, TKEM, HP, HM, 
     1                 WORK, F)
C
C     Routine to calculate interface fluxes using Roe's
C     approximate Riemann solver.  This routine calculates the
C     fluxes for the Favre averaged Navier-Stokes equations for inviscid,
C     laminar, algebraic turbulence model, two-equation turbulence model
C     or Reynolds stress turbulence model.  It computes all of the inviscid
C     fluxes for the above models.
C
C     The two-equation and Reynolds stress models need the turbulent kinetic 
C     energy (K) in the mean flow equations.  Therefore, this routine is
C     organized into the following sections:
C          1. Calculate Turbulent Kinetic Energy (TKEP,TKEM)
C             a. Inviscid, laminar, 
C                algebraic model       => K = 0
C             b. Two-Equation model    => K = K (Q(6))
C             c. Reynolds stress model => K = (Tau_xx+Tau_yy+Tau_zz)/2
C          2. Calculate Total Enthalpy (HP,HM)
C             a. Perfect gas    => H = gamma P / (RHO (GAMMA-1)) + (q^2)/2 + K
C          3. Calculate Advection Fluxes for Mean Flow Equations
C          4. Calculate additional Advection Fluxes
C             a. Inviscid, laminar, 
C                algebraic model       => No additional fluxes
C             b. Two-Equation model    => Fluxes for K and T2 equations
C             c. Reynolds Stress model => Fluxes for Tau_ij and T2 equations
C
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     S              : Metrics
C     TKEP,TKEM      : Turbulent kinetic energy for QP,QM
C     HP,HM          : Total enthalpy for QP,QM
C     WORK           : Temporary storage for the terms common to the 
C                      T |Lambda| T[-1] DQ term
C     F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION QP   (NPTS,NQ),
     1          QM   (NPTS,NQ),
     2          S    (NPTS,3),
     3          TKEP (NPTS),
     4          TKEM (NPTS),
     5          HP   (NPTS),
     6          HM   (NPTS),
     7          WORK (NPTS,4),
     8          F    (NPTS,NF)
C
C     1. Calculate the Turbulence Kinetic Energy (TKEP,TKEM)
C
C     1-b. Turbulence Kinetic Energy for Two Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKEP(I) = QP(I,6)
            TKEM(I) = QM(I,6)
  100    CONTINUE
C
C     1-c. Turbulence Kinetic Energy for Reynolds Stress Turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
C
            TAUXXP  = QP(I,6)
            TAUYYP  = QP(I,7)
            TAUZZP  = QP(I,8)
            TKEP(I) = 0.5E0 * (TAUXXP + TAUYYP + TAUZZP)
C
            TAUXXM  = QM(I,6)
            TAUYYM  = QM(I,7)
            TAUZZM  = QM(I,8)
            TKEM(I) = 0.5E0 * (TAUXXM + TAUYYM + TAUZZM)
  110    CONTINUE
C
C     1-a. Turbulence Kinetic Energy = 0 for Inviscid, Laminar and
C          Algebraic Turbulence models
C
      ELSE
         DO 120 I = 1, NPTS
            TKEP(I) = 0.E0
            TKEM(I) = 0.E0
  120    CONTINUE
      ENDIF
C
C     2. Calculate Total Enthalpy
C
C     2-a. Perfect Gas
C
      DO 200 I = 1, NPTS
C
         RHOP   = QP(I,1)
         UP     = QP(I,2) 
         VP     = QP(I,3)
         WP     = QP(I,4)
         PP     = QP(I,5)
         HP(I)  = GAMMA * PP / (RHOP * GAMM1) 
     1            + 0.5E0 * (UP*UP + VP*VP + WP*WP) + TKEP(I)
C
         RHOM   = QM(I,1)
         UM     = QM(I,2) 
         VM     = QM(I,3)
         WM     = QM(I,4)
         PM     = QM(I,5)
         HM(I)  = GAMMA * PM / (RHOM * GAMM1)
     1            + 0.5E0 * (UM*UM + VM*VM + WM*WM) + TKEM(I)
  200 CONTINUE
C
C     3. Calculate F for Mean Flow Equations at the cell interfaces 
C        using Roe's approximate Riemann Solver
C
      DO 300 I = 1, NPTS
C
C     Calculate F+
C
         RHOP   = QP(I,1)
         UP     = QP(I,2) 
         VP     = QP(I,3)
         WP     = QP(I,4)
         PP     = QP(I,5)
C
         UBARP  = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
C
C     Calculate F-
C
         RHOM   = QM(I,1)
         UM     = QM(I,2) 
         VM     = QM(I,3)
         WM     = QM(I,4)
         PM     = QM(I,5)
C
         UBARM  = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
C
C     Calculate T * ABS(LAMBDA) * T(-1) * DQ
C
         DRHO   = RHOP   - RHOM
         DU     = UP     - UM
         DV     = VP     - VM
         DW     = WP     - WM
         DUBAR  = UBARP  - UBARM
         DP     = PP     - PM
C
         DK     = TKEP(I) - TKEM(I)
C
         RHOMSR = SQRT (RHOM) 
         RHOPSR = SQRT (RHOP)
         RHOAV  = RHOMSR + RHOPSR
         UR     = ( RHOMSR*UM + RHOPSR*UP ) / RHOAV
         VR     = ( RHOMSR*VM + RHOPSR*VP ) / RHOAV
         WR     = ( RHOMSR*WM + RHOPSR*WP ) / RHOAV
         UBARR  = S(I,1)*UR + S(I,2)*VR + S(I,3)*WR
         HR     = ( RHOMSR*HM(I)   + RHOPSR*HP(I) )   / RHOAV
C
         TKER   = ( RHOMSR*TKEM(I) + RHOPSR*TKEP(I) ) / RHOAV
C
         QSQR   = UR*UR + VR*VR + WR*WR
         ASQR   = GAMM1 * (HR - 0.5E0*QSQR - TKER)
         AR     = SQRT (ASQR)
         RHOR   = RHOMSR*RHOPSR
C
C     Calculate the eigenvalues.
C
         ALAM1  = ABS (UBARR) 
         ALAM4  = ABS (UBARR + AR) 
         ALAM5  = ABS (UBARR - AR) 
C
C     Entropy Fix: Based on Yee, Warming, and Harten (1985) and Yee ()
C       Apply fix only to acoustic eigenvalues to preserve boundary layer
C       calculations. Avoidance of the carbuncle phenomena requires fix
C       applied to convective eigenvalue also. Details based on conversations
C       with Jeff White.
C       Begin {N.B. Comment out this section to remove the entropy fix
C
         TSTLAM = CNTRPY * (SQRT(QSQR) + AR)
         TSTLAM = MAX (TSTLAM, RSMALL)
         IF (ALAM4 .LT. TSTLAM) ALAM4=0.5E0*(ALAM4*ALAM4/TSTLAM+TSTLAM)
         IF (ALAM5 .LT. TSTLAM) ALAM5=0.5E0*(ALAM5*ALAM5/TSTLAM+TSTLAM)
C
C       End} Entropy Fix
C
C     Temporary constants used in T * ABS(LAMBDA) * T(-1) * DQ term
C
         WORK(I,1) = ALAM1 * (DRHO - DP/ASQR)
         WORK(I,2) = ALAM1 * RHOR
         WORK(I,3) = ALAM4 * 0.5E0 * (DP/ASQR + RHOR*DUBAR/AR)
         WORK(I,4) = ALAM5 * 0.5E0 * (DP/ASQR - RHOR*DUBAR/AR)
C
C     Calculate F = 0.5 * (F+ + F-) - 0.5 * T * ABS(LAMBDA) * T(-1) * DQ
C     Area term is included in the RESID routine.
C
         F(I,1)  = 0.5E0 * 
     1             ( RHOP*UBARP
     2             + RHOM*UBARM
     3             - ( WORK(I,1)
     4               + WORK(I,3)
     5               + WORK(I,4) ) )
         F(I,2)  = 0.5E0 *
     1             ( (RHOP*UBARP*UP + S(I,1)*PP)
     2             + (RHOM*UBARM*UM + S(I,1)*PM)
     3             - ( UR * WORK(I,1)
     4               + (DU - S(I,1)*DUBAR) * WORK(I,2)
     5               + (UR + AR*S(I,1))    * WORK(I,3)
     6               + (UR - AR*S(I,1))    * WORK(I,4) ) )
         F(I,3)  = 0.5E0 * 
     1             ( (RHOP*UBARP*VP + S(I,2)*PP)
     2             + (RHOM*UBARM*VM + S(I,2)*PM)
     3             - ( VR * WORK(I,1)
     4               + (DV - S(I,2)*DUBAR) * WORK(I,2)
     5               + (VR + AR*S(I,2))    * WORK(I,3)
     6               + (VR - AR*S(I,2))    * WORK(I,4) ) )
         F(I,4)  = 0.5E0 * 
     1             ( (RHOP*UBARP*WP + S(I,3)*PP)
     2             + (RHOM*UBARM*WM + S(I,3)*PM)
     3             - ( WR * WORK(I,1)
     4               + (DW - S(I,3)*DUBAR) * WORK(I,2)
     5               + (WR + AR*S(I,3))    * WORK(I,3)
     6               + (WR - AR*S(I,3))    * WORK(I,4) ) )
         F(I,5)  = 0.5E0 *
     1             ( RHOP*UBARP*HP(I)
     2             + RHOM*UBARM*HM(I)
     3             - ( (0.5E0*QSQR + TKER) * WORK(I,1)
     4               + (UR*DU+VR*DV+WR*DW-UBARR*DUBAR+DK)*WORK(I,2)
     5               + (HR + AR*UBARR)    * WORK(I,3)
     6               + (HR - AR*UBARR)    * WORK(I,4) ) )
 300  CONTINUE
C
C     4. Calculate Additional Advection Fluxes
C
C     4-a. Inviscid, Laminar, Algebraic Turbulence models => no additional flux
C
C     4-b. Two-Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 400 I = 1, NPTS
C
            RHOP    = QP(I,1)
            UP      = QP(I,2) 
            VP      = QP(I,3)
            WP      = QP(I,4)
C
            T2P     = QP(I,7)
C
            UBARP   = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
C
            RHOM    = QM(I,1)
            UM      = QM(I,2) 
            VM      = QM(I,3)
            WM      = QM(I,4)
C
            T2M     = QM(I,7)
C
            UBARM   = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
C
            DK      = TKEP(I) - TKEM(I)
            DT2     = T2P     - T2M
C
            RHOMSR  = SQRT (RHOM) 
            RHOPSR  = SQRT (RHOP)
            RHOAV   = RHOMSR + RHOPSR
            TKER    = ( RHOMSR*TKEM(I) + RHOPSR*TKEP(I) ) / RHOAV
            T2R     = ( RHOMSR*T2M     + RHOPSR*T2P )     / RHOAV
C
            F(I,6)  = 0.5E0 *
     1                ( RHOP*UBARP*TKEP(I)
     2                + RHOM*UBARM*TKEM(I)
     3                - ( TKER * WORK(I,1)
     4                  + DK   * WORK(I,2)
     5                  + TKER * WORK(I,3)
     6                  + TKER * WORK(I,4) ) )
            F(I,7)  = 0.5E0 *
     1                ( RHOP*UBARP*T2P
     2                + RHOM*UBARM*T2M
     3                - ( T2R * WORK(I,1)
     4                  + DT2 * WORK(I,2)
     5                  + T2R * WORK(I,3)
     6                  + T2R * WORK(I,4) ) )
  400    CONTINUE
C
C     4-c. Reynolds Stress Turbulence models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 410 I = 1, NPTS
C
            RHOP    = QP(I,1)
            UP      = QP(I,2) 
            VP      = QP(I,3)
            WP      = QP(I,4)
C
            TAUXXP  = QP(I,6)
            TAUYYP  = QP(I,7)
            TAUZZP  = QP(I,8)
            TAUXYP  = QP(I,9)
            TAUXZP  = QP(I,10)
            TAUYZP  = QP(I,11)
            T2P     = QP(I,12)
C
            UBARP   = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
C
            RHOM    = QM(I,1)
            UM      = QM(I,2) 
            VM      = QM(I,3)
            WM      = QM(I,4)
C
            TAUXXM  = QM(I,6)
            TAUYYM  = QM(I,7)
            TAUZZM  = QM(I,8)
            TAUXYM  = QM(I,9)
            TAUXZM  = QM(I,10)
            TAUYZM  = QM(I,11)
            T2M     = QM(I,12)
C
            UBARM   = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
C
            DTAUXX  = TAUXXP - TAUXXM
            DTAUYY  = TAUYYP - TAUYYM
            DTAUZZ  = TAUZZP - TAUZZM
            DTAUXY  = TAUXYP - TAUXYM
            DTAUXZ  = TAUXZP - TAUXZM
            DTAUYZ  = TAUYZP - TAUYZM
            DT2     = T2P    - T2M
C
            RHOMSR  = SQRT (RHOM) 
            RHOPSR  = SQRT (RHOP)
            RHOAV   = RHOMSR + RHOPSR
            TAUXXR  = ( RHOMSR*TAUXXM + RHOPSR*TAUXXP ) / RHOAV
            TAUYYR  = ( RHOMSR*TAUYYM + RHOPSR*TAUYYP ) / RHOAV
            TAUZZR  = ( RHOMSR*TAUZZM + RHOPSR*TAUZZP ) / RHOAV
            TAUXYR  = ( RHOMSR*TAUXYM + RHOPSR*TAUXYP ) / RHOAV
            TAUXZR  = ( RHOMSR*TAUXZM + RHOPSR*TAUXZP ) / RHOAV
            TAUYZR  = ( RHOMSR*TAUYZM + RHOPSR*TAUYZP ) / RHOAV
            T2R     = ( RHOMSR*T2M    + RHOPSR*T2P )    / RHOAV
C
            F(I,6)  = 0.5E0 *
     1                ( RHOP*UBARP*TAUXXP
     2                + RHOM*UBARM*TAUXXM
     3                - ( TAUXXR * WORK(I,1)
     4                  + DTAUXX * WORK(I,2)
     5                  + TAUXXR * WORK(I,3)
     6                  + TAUXXR * WORK(I,4) ) )
            F(I,7)  = 0.5E0 *
     1                ( RHOP*UBARP*TAUYYP
     2                + RHOM*UBARM*TAUYYM
     3                - ( TAUYYR * WORK(I,1)
     4                  + DTAUYY * WORK(I,2)
     5                  + TAUYYR * WORK(I,3)
     6                  + TAUYYR * WORK(I,4) ) )
            F(I,8)  = 0.5E0 *
     1                ( RHOP*UBARP*TAUZZP
     2                + RHOM*UBARM*TAUZZM
     3                - ( TAUZZR * WORK(I,1)
     4                  + DTAUZZ * WORK(I,2)
     5                  + TAUZZR * WORK(I,3)
     6                  + TAUZZR * WORK(I,4) ) )
            F(I,9)  = 0.5E0 *
     1                ( RHOP*UBARP*TAUXYP
     2                + RHOM*UBARM*TAUXYM
     3                - ( TAUXYR * WORK(I,1)
     4                  + DTAUXY * WORK(I,2)
     5                  + TAUXYR * WORK(I,3)
     6                  + TAUXYR * WORK(I,4) ) )
            F(I,10) = 0.5E0 *
     1                ( RHOP*UBARP*TAUXZP
     2                + RHOM*UBARM*TAUXZM
     3                - ( TAUXZR * WORK(I,1)
     4                  + DTAUXZ * WORK(I,2)
     5                  + TAUXZR * WORK(I,3)
     6                  + TAUXZR * WORK(I,4) ) )
            F(I,11) = 0.5E0 *
     1                ( RHOP*UBARP*TAUYZP
     2                + RHOM*UBARM*TAUYZM
     3                - ( TAUYZR * WORK(I,1)
     4                  + DTAUYZ * WORK(I,2)
     5                  + TAUYZR * WORK(I,3)
     6                  + TAUYZR * WORK(I,4) ) )
            F(I,12) = 0.5E0 *
     1                ( RHOP*UBARP*T2P
     2                + RHOM*UBARM*T2M
     3                - ( T2R * WORK(I,1)
     4                  + DT2 * WORK(I,2)
     5                  + T2R * WORK(I,3)
     6                  + T2R * WORK(I,4) ) )
  410    CONTINUE
C
      ENDIF
C
C     Finished calculating flux at the interface with Roe's scheme
C
      RETURN
      END





      SUBROUTINE LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LBEG, LEND, 
     1                   IERROR)
C
C The second order limited correction to the first order interpolation
C is calculated by this routine and stored in DQP for the q(+) term and 
C DQM for the q(-) term.  Delta(q(i)) = q(i+1) - q(i) is stored in DQP(i) 
C on entry.
C
C N.B. The values of DQP(NPTS) and DQM(1) are not defined properly on exit 
C as this leads to better vectorization and they are unnecessary.
C
C This routine is organized as follows:
C     0. Set constants for the interpolation
C     1. Venkatakrishnan's limiter
C     2. MinMod limiter
C     3. Smooth limiter
C     4. Unlimited scheme
C     5. Error
C
C NPTS           : Length of array to limit successive values of delta(q)
C DQP,DQM        : Storage for delta q = q(i) - q(i-1) on entry.
C                  Storage of the higher order correction to the 
C                  first order scheme on exit
C LIMIT          : Limiter to use for MUSCL interpolation
C                  ILNONE [0] => Unlimited Scheme
C                  ILSMTH [1] => Smooth Limiter
C                  ILMNMD [2] => Minmod Limiter
C                  ILVNKT [3] => Smooth limiter (fixed at kappa = 1/3) 
C                           from Venkatakrishnan, AIAA J. Vol. 29 No. 7
C                           also Thomas, J. L., 5th Copper Mtn. Multigrid
C RKAP           : Kappa to use in MUSCL interpolation
C LBEG,LEND      : Starting,stopping equation index for limiting
C IERROR         : Error code
C
      include '../header/common.h'
C
      DIMENSION DQP   (NPTS,NQ),
     1          DQM   (NPTS,NQ)
C
C     Functions defining the limiters.
C
      SMOOTH(D1, D2, EP, EPSQ) = ( 2.E0*D1*D2 + EP ) 
     1                           / ( D1*D1 + D2*D2 + EPSQ )
      RMINMD(D1, D2, B) = 0.5E0 * ( SIGN(1.E0,D1) + SIGN(1.E0,D2) )
     1                          * MIN ( ABS(D1), ABS(B*D2) )
c     RMINMD(D1, D2, B) = SIGN(1.E0, D1) * 
c    1      MAX( 0.E0, MIN( D1*SIGN(1.E0,D2), B*D2*SIGN(1.E0,D1) ) )
      SMI(D1, D2, EPSQ) = ( D1*(D2*D2+2.E0*EPSQ) 
     1                    + D2*(2.E0*D1*D1+EPSQ) )
     2                / ( 2.E0*D1*D1 - D1*D2 + 2.E0*D2*D2 + 3.E0*EPSQ )
C
C     Initialize error
C
      IERROR = INO
C
C 0.  Set constants for the interpolation.
C
      OMKAP = 1.E0 - RKAP
      OPKAP = 1.E0 + RKAP
C
C 1.  Venkat's Smooth Limiter
C
      IF (LIMIT .EQ. ILVNKT) THEN
c        EPSQ = 1.E-12
         EPSQ = 1.E-8
         DO L = LBEG, LEND
            DO I = 2, NPTS
               XP = SMI( DQP(I-1,L), DQP(I,L),   EPSQ )
               XM = SMI( DQP(I,L),   DQP(I-1,L), EPSQ ) 
               DQP(I-1,L) = - 0.5E0 * XP
               DQM(I,L)   =   0.5E0 * XM
            END DO
         END DO
C
C 2.  MinMod Limiter
C
      ELSE IF (LIMIT .EQ. ILMNMD) THEN
         B = (3.E0 - RKAP) / (1.E0 - RKAP)
         DO L = LBEG, LEND
            DO I = 2, NPTS
               XM = RMINMD( DQP(I-1,L), DQP(I,L),   B )
               XP = RMINMD( DQP(I,L),   DQP(I-1,L), B ) 
               DQP(I-1,L) = - 0.25E0 * ( OMKAP*XP + OPKAP*XM )
               DQM(I,L)   =   0.25E0 * ( OMKAP*XM + OPKAP*XP )
            END DO
         END DO
C
C 3.  Smooth Limiter
C
      ELSE IF (LIMIT .EQ. ILSMTH) THEN
      write (iout,'(a)') ' WARNING - THE SMOOTH LIMITER IS UNTESTED '
         EPS   = 1.E-06
         EPSSQ = EPS*EPS
         DO L = LBEG, LEND
            DO I = 2, NPTS
               XP         = DQP(I,L)
               XM         = DQP(I-1,L)
               SM         = SMOOTH (XM, XP, EPS, EPSSQ)
               DQP(I-1,L) = - 0.25E0 * SM * 
     1               ( (1.E0 - SM*RKAP)*XP
     2               + (1.E0 + SM*RKAP)*XM )
               DQM(I,L) =     0.25E0 * SM * 
     1               ( (1.E0 - SM*RKAP)*XM
     2               + (1.E0 + SM*RKAP)*XP )
            END DO
         END DO
C
C 4.  Unlimited Scheme
C
      ELSE IF (LIMIT .EQ. ILNONE) THEN
         DO L = LBEG, LEND
            DO I = 2, NPTS
               XP         = DQP(I,L)
               XM         = DQP(I-1,L)
               DQP(I-1,L) = - 0.25E0*(OMKAP*XP + OPKAP*XM)
               DQM(I,L)   =   0.25E0*(OMKAP*XM + OPKAP*XP)
            END DO
         END DO
C
C 5.  Error
C
      ELSE
         WRITE (IOUT,1000) LIMIT
         IERROR = IYES
         STOP
      ENDIF
C
C Finished with limiting routine.
C
      RETURN
 1000 FORMAT(' ','LIMITR  : Invalid value of LIMIT = ',I10,
     1      /' ','Run Aborting.')
      END






      SUBROUTINE MUSCLI (IDIM, JDIM, KDIM, K, IPSTRT, IPEND,
     1                   Q, DQP, DQM, QP, QM,
     2                   IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C Routine to do MUSCL interpolation in the I direction. This routine 
C interpolates for Q(+,-) on the I faces for i = ipstrt-1, ipend by 
C j = 2, jdim on the current K plane.  (IPSTRT and IPEND refer to cell centers, 
C therefore the interfaces are from IPSTRT-1 to IPEND.) This allows for greater 
C vectorization and also minimizes intermediate storage to a cross section.  
C
C Note the beginning and ending indices of the arrays.
C
C This routine is organized as follows:
C     0. Set up beginning/ending indices for first and second order sections
C     1. Calculate Q(+,-)
C        a. First order
C        b. Second order
C           i.   Calculate delta(Q)
C           ii.  Limit delta(Q)
C           iii. Construct Q(+,-)
C     2. Correct wall values of Q(+,-)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C K              : K Plane to do interpolation on
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C Q              : Primitive variables at cell centers
C DQP,DQM        : Storage for delta q = q(i) - q(i-1) and
C                  for the higher order interpolant to be added
C                  to the first order scheme to create QP,QM
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
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
C IORDER         : Order of MUSCL interpolation for I direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          DQP   (IPSTRT-2:IPEND+1,2:JDIM,NQ),
     2          DQM   (IPSTRT-2:IPEND+1,2:JDIM,NQ),
     3          QP    (IPSTRT-1:IPEND,2:JDIM,NQ), 
     4          QM    (IPSTRT-1:IPEND,2:JDIM,NQ) 
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local Storage
C
      LOGICAL   DOFRST,
     1          DOSCND
C
C 0.  Set up beginning/ending indices for first and second order sections
C
      LFOBEG = 0
      LFOEND = -1
      LSOBEG = 1
      LSOEND = NQ
      DOFRST = .FALSE.
      DOSCND = .TRUE.
C
      IF (IORDER .EQ. 1) THEN
         LFOBEG = 1
         LFOEND = NQ
         LSOBEG = 0
         LSOEND = -1
         DOFRST = .TRUE.
         DOSCND = .FALSE.
      ELSE IF (FOTURB) THEN
         LFOBEG = 6
         LFOEND = NQ
         LSOBEG = 1
         LSOEND = 5
         DOFRST = .TRUE.
         DOSCND = .TRUE.
      END IF
C
C 1.  Calculate Q(+,-)
C     Loop over j,i formulating q(extrapolated) - inner loop is over j for 
C     maximum vectorization for 2d and planar oriented algorithms. This also 
C     results in unit stride. I is stored first so that successive values of 
C     i (i = 1,2,3...) are stored together for limiting.
C
C 1-a. First Order
C
      IF (DOFRST) THEN
         DO 130 L = LFOBEG, LFOEND
            DO 120 I = IPSTRT - 1, IPEND
               DO 110 J = 2, JDIM
                  QP(I,J,L) = Q(J,K,I+1,L) 
                  QM(I,J,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
      END IF
C
C 1-b. Second Order or Third Order
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(i)) = q(i+1) - q(i) to formulate q(extrapolated) 
C        using limiting.
C
      IF (DOSCND) THEN
         DO 230 L = LSOBEG, LSOEND
            DO 220 I = IPSTRT - 2, IPEND + 1
               DO 200 J = 2, JDIM
                  DQP(I,J,L) = Q(J,K,I+1,L) - Q(J,K,I,L)
  200          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 1 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               KSTRT  = IBCDAT(8,IBC)
               KEND   = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  IBCPT = 1
               ELSE
                  IBCPT = IDIM
               ENDIF
               FACTR  = 2.E0
               IF (IBCTYP .EQ. IWALFN) FACTR = 0.E0
               IF (IBCPT .GE. IPSTRT-1 .AND. IBCPT .LE. IPEND .AND.
     1                 K .GE. KSTRT    .AND.     K .LE. KEND      ) THEN
                  DO 290 L = LSOBEG, LSOEND
                     DO 280 J = JSTRT, JEND
                        DQP(IBCPT,J,L) = FACTR * DQP(IBCPT,J,L)
  280                CONTINUE
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from i = ipstrt-2 to i = ipend+1 by j = 2 to j = jdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each i line, but are
C     calculated here to maximize vectorization on a ixj plane.
C
         NPTS = ( IPEND - IPSTRT + 4 ) * ( JDIM - 1 )
         CALL LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LSOBEG, LSOEND, 
     1                IERROR)
C
C 1-b-iii. Construct Q(+,-)
C     Calculate q(extrapolated,i+1/2) from q and the limited values of delta(q)
C
         DO 380 L = LSOBEG, LSOEND
            DO 370 I = IPSTRT - 1, IPEND
               DO 360 J = 2, JDIM
                  QP(I,J,L) = Q(J,K,I+1,L) + DQP(I,J,L)
                  QM(I,J,L) = Q(J,K,I,L)   + DQM(I,J,L)
  360          CONTINUE
  370       CONTINUE
  380    CONTINUE
      ENDIF
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
      DO 800 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( IDIR .EQ. 1 .AND.
     1        ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2          IBCTYP .EQ. IWALFN ) ) THEN
            INROUT = IBCDAT(3,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               IQPM  = 1
               IBCPT = 1
            ELSE
               IQPM  = IDIM
               IBCPT = IDIM + 1
            ENDIF
            IF (IBCPT .GE. IPSTRT-1 .AND. IBCPT .LE. IPEND .AND.
     1              K .GE. KSTRT    .AND.     K .LE. KEND        ) THEN
               DO 790 L = 1, NQ
                  DO 780 J = JSTRT, JEND
                     QP(IQPM,J,L) = Q(J,K,IBCPT,L)
                     QM(IQPM,J,L) = Q(J,K,IBCPT,L)
  780             CONTINUE
  790          CONTINUE
            ENDIF
         ENDIF
  800 CONTINUE
C
C Finished with Q(+/-) calculation for MUSCL scheme in i.
C
      RETURN
      END






      SUBROUTINE MUSCLJ (IDIM, JDIM, KDIM, I,
     1                   Q, DQP, DQM, QP, QM,
     2                   IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C Routine to do MUSCL interpolation in the J direction. This routine 
C interpolates for Q(+,-) on the J faces for j = 1, jdim by k = 2, kdim 
C on the current I plane. This allows for greater vectorization and also 
C minimizes intermediate storage to a cross section.  
C
C This routine is organized as follows:
C     0. Set up beginning/ending indices for first and second order sections
C     1. Calculate Q(+,-)
C        a. First order
C        b. Second order
C           i.   Calculate delta(Q)
C           ii.  Limit delta(Q)
C           iii. Construct Q(+,-)
C     2. Correct wall values of Q(+,-)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane to do interpolation on
C Q              : Primitive variables at cell centers
C DQP,DQM        : Storage for delta q = q(j) - q(j-1) and
C                  for the higher order interpolant to be added
C                  to the first order scheme to create QP,QM
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = j(1 or jdim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for I direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION DQP   (0:JDIM+1,2:KDIM,NQ),
     1          DQM   (0:JDIM+1,2:KDIM,NQ),
     2          QP    (JDIM,2:KDIM,NQ), 
     3          QM    (JDIM,2:KDIM,NQ) 
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local Storage
C
      LOGICAL   DOFRST,
     1          DOSCND
C
C 0.  Set up beginning/ending indices for first and second order sections
C
      LFOBEG = 0
      LFOEND = -1
      LSOBEG = 1
      LSOEND = NQ
      DOFRST = .FALSE.
      DOSCND = .TRUE.
C
      IF (IORDER .EQ. 1) THEN
         LFOBEG = 1
         LFOEND = NQ
         LSOBEG = 0
         LSOEND = -1
         DOFRST = .TRUE.
         DOSCND = .FALSE.
      ELSE IF (FOTURB) THEN
         LFOBEG = 6
         LFOEND = NQ
         LSOBEG = 1
         LSOEND = 5
         DOFRST = .TRUE.
         DOSCND = .TRUE.
      END IF
C
C 1.  Calculate Q(+,-)
C     Loop over j,k formulating q(extrapolated). J is stored first so that 
C     successive values of j (j=1,2,3...) are stored together for limiting.
C
C 1-a. First Order
C
      IF (DOFRST) THEN
         DO 130 L = LFOBEG, LFOEND
            DO 120 K = 2, KDIM
               DO 110 J = 1, JDIM
                  QP(J,K,L) = Q(J+1,K,I,L) 
                  QM(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
      END IF
C
C 1-b. Second Order or Third Order
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(j)) = q(j+1) - q(j) to formulate q(extrapolated) 
C        using limiting.
C
      IF (DOSCND) THEN
         DO 230 L = LSOBEG, LSOEND
            DO 220 K = 2, KDIM
               DO 200 J = 0, JDIM + 1
                  DQP(J,K,L) = Q(J+1,K,I,L) - Q(J,K,I,L)
  200          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 2 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               KSTRT  = IBCDAT(8,IBC)
               KEND   = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  JBCPT = 1
               ELSE
                  JBCPT = JDIM
               ENDIF
               FACTR  = 2.E0
               IF (IBCTYP .EQ. IWALFN) FACTR = 0.E0
               IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
                  DO 290 L = LSOBEG, LSOEND
                     DO 280 K = KSTRT, KEND
                        DQP(JBCPT,K,L) = FACTR * DQP(JBCPT,K,L)
  280                CONTINUE
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from j = 0 to j = jdim+1 by k = 2 to k = kdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each j line, but are
C     calculated here to maximize vectorization on a jxk plane.
C
         NPTS = ( JDIM + 2 ) * ( KDIM - 1 )
         CALL LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LSOBEG, LSOEND, 
     1                IERROR)
C
C     Calculate q(extrapolated,j+1/2) from q and the 
C     limited values of delta(q).
C
         DO 380 L = LSOBEG, LSOEND
            DO 370 K = 2, KDIM
               DO 360 J = 1, JDIM
                  QP(J,K,L) = Q(J+1,K,I,L) + DQP(J,K,L)
                  QM(J,K,L) = Q(J,K,I,L)   + DQM(J,K,L)
  360          CONTINUE
  370       CONTINUE
  380    CONTINUE
      ENDIF
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
C
      DO 800 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( IDIR .EQ. 2 .AND.
     1        ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2          IBCTYP .EQ. IWALFN ) ) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               JQPM  = 1
               JBCPT = 1
            ELSE
               JQPM  = JDIM
               JBCPT = JDIM + 1
            ENDIF
            IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
               DO 790 L = 1, NQ
                  DO 780 K = KSTRT, KEND
                     QP(JQPM,K,L) = Q(JBCPT,K,I,L)
                     QM(JQPM,K,L) = Q(JBCPT,K,I,L)
  780             CONTINUE
  790          CONTINUE
            ENDIF
         ENDIF
  800 CONTINUE
C
C Finished with Q(+/-) calculation for MUSCL scheme in j.
C
      RETURN
      END






      SUBROUTINE MUSCLK (IDIM, JDIM, KDIM, I,
     1                   Q, DQP, DQM, QP, QM,
     2                   IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C Routine to do MUSCL interpolation in the K direction. This routine 
C interpolates for Q(+,-) on the K faces for k = 1, kdim by j = 2, jdim 
C on the current I plane. This allows for greater vectorization and also 
C minimizes intermediate storage to a cross section.  
C
C This routine is organized as follows:
C     0. Set up beginning/ending indices for first and second order sections
C     1. Calculate Q(+,-)
C        a. First order
C        b. Second order
C           i.   Calculate delta(Q)
C           ii.  Limit delta(Q)
C           iii. Construct Q(+,-)
C     2. Correct wall values of Q(+,-)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane to do interpolation on
C Q              : Primitive variables at cell centers
C DQP,DQM        : Storage for delta q = q(k) - q(k-1) and
C                  for the higher order interpolant to be added
C                  to the first order scheme to create QP,QM
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = k(1 or kdim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for I direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION DQP   (0:KDIM+1,2:JDIM,NQ),
     1          DQM   (0:KDIM+1,2:JDIM,NQ),
     2          QP    (KDIM,2:JDIM,NQ), 
     3          QM    (KDIM,2:JDIM,NQ) 
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local Storage
C
      LOGICAL   DOFRST,
     1          DOSCND
C
C 0.  Set up beginning/ending indices for first and second order sections
C
      LFOBEG = 0
      LFOEND = -1
      LSOBEG = 1
      LSOEND = NQ
      DOFRST = .FALSE.
      DOSCND = .TRUE.
C
      IF (IORDER .EQ. 1) THEN
         LFOBEG = 1
         LFOEND = NQ
         LSOBEG = 0
         LSOEND = -1
         DOFRST = .TRUE.
         DOSCND = .FALSE.
      ELSE IF (FOTURB) THEN
         LFOBEG = 6
         LFOEND = NQ
         LSOBEG = 1
         LSOEND = 5
         DOFRST = .TRUE.
         DOSCND = .TRUE.
      END IF
C
C 1.  Calculate Q(+,-)
C     Loop over j,k formulating q(extrapolated). K is stored first so that 
C     successive values of k (k = 1,2,3...) are stored together for limiting.
C
C 1-a. First Order
C
      IF (DOFRST) THEN
         DO 130 L = LFOBEG, LFOEND
            DO 120 J = 2, JDIM
               DO 110 K = 1, KDIM
                  QP(K,J,L) = Q(J,K+1,I,L) 
                  QM(K,J,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
      END IF
C
C 1-b. Second Order or Third Order
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(k)) = q(k+1) - q(k) to formulate q(extrapolated) 
C        using limiting.
C
      IF (DOSCND) THEN
         DO 230 L = LSOBEG, LSOEND
            DO 220 J = 2, JDIM
               DO 200 K = 0, KDIM + 1
                  DQP(K,J,L) = Q(J,K+1,I,L) - Q(J,K,I,L)
  200          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 3 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               IF (INROUT .EQ. 1) THEN
                  KBCPT = 1
               ELSE
                  KBCPT = KDIM
               ENDIF
               FACTR  = 2.E0
               IF (IBCTYP .EQ. IWALFN) FACTR = 0.E0
               IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
                  DO 290 L = LSOBEG, LSOEND
                     DO 280 J = JSTRT, JEND
                        DQP(KBCPT,J,L) = FACTR * DQP(KBCPT,J,L)
  280                CONTINUE
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from k = 0 to k = kdim+1 by j = 2 to j = jdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each k line, but are
C     calculated here to maximize vectorization on a kxj plane.
C
         NPTS = ( KDIM + 2 ) * ( JDIM - 1 )
         CALL LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LSOBEG, LSOEND, 
     1                IERROR)
C
C 1-b-iii. Construct Q(+,-)
C     Calculate q(extrapolated,k+1/2) from q and the limited values of delta(q)
C
         DO 380 L = LSOBEG, LSOEND
            DO 370 J = 2, JDIM
               DO 360 K = 1, KDIM
                  QP(K,J,L) = Q(J,K+1,I,L) + DQP(K,J,L)
                  QM(K,J,L) = Q(J,K,I,L)   + DQM(K,J,L)
  360          CONTINUE
  370       CONTINUE
  380    CONTINUE
      ENDIF
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
C
      DO 800 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( IDIR .EQ. 3 .AND.
     1        ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2          IBCTYP .EQ. IWALFN ) ) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            IF (INROUT .EQ. 1) THEN
               KQPM  = 1
               KBCPT = 1
            ELSE
               KQPM  = KDIM
               KBCPT = KDIM + 1
            ENDIF
            IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
               DO 790 L = 1, NQ
                  DO 780 J = JSTRT, JEND
                     QP(KQPM,J,L) = Q(J,KBCPT,I,L)
                     QM(KQPM,J,L) = Q(J,KBCPT,I,L)
  780             CONTINUE
  790          CONTINUE
            ENDIF
         ENDIF
  800 CONTINUE
C
C Finished with Q(+/-) calculation for MUSCL scheme in k.
C
      RETURN
      END
