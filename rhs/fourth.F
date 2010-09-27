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
C     $Revision: 4.1 $
C     $Author: jhmorr $
C     $Date: 2001/11/13 04:45:54 $
C     $State: Exp $
C     $Log: fourth.stub.FOR,v $
C     Revision 4.1  2001/11/13 04:45:54  jhmorr
C     Renamed fourth.stub.F to fourth.stub.FOR to fix distribution.
C
C     Revision 4.3  2001/10/29 03:25:31  jhmorr
C     Updated license information
C
C     Revision 4.2  2001/06/08 04:56:16  jhmorr
C     Added notice
C
C     Revision 4.1  2001/06/08 04:17:19  jhmorr
C     Added fourth.stub.F
C
C





      SUBROUTINE DQCC   (IDIM, JDIM, KDIM, IPSTRT, IPEND, RES, DQ)
C
      RETURN
      END





      SUBROUTINE DELQI4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, K,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                   STEMP, VOLAUX, QI, PROPSI,
     3                   DQDX, DQDY, DQDZ)
C
      RETURN
      END




      SUBROUTINE DELQJ4 (IDIM, JDIM, KDIM, I,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                   STEMP, VOLAUX, QI, PROPSI,
     3                   DQDX, DQDY, DQDZ)
C
      RETURN
      END




      SUBROUTINE DELQK4 (IDIM, JDIM, KDIM, I,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                   STEMP, VOLAUX, QI, PROPSI,
     3                   DQDX, DQDY, DQDZ)
C
      RETURN
      END




      SUBROUTINE FLUXI4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, 
     4                   IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
      RETURN
      END




      SUBROUTINE FLUXJ4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, 
     4                   IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
      RETURN
      END




      SUBROUTINE FLUXK4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, 
     4                   IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
      RETURN
      END




      SUBROUTINE QICD   (IDIM, JDIM, KDIM, K, IPSTRT, IPEND, Q, QI, 
     1                   QC, DQC, IBCDAT, NBCS, IORDER)
C
      RETURN
      END





      SUBROUTINE QJCD   (IDIM, JDIM, KDIM, I, Q, QI, QC, DQC,
     1                   IBCDAT, NBCS, ITURB, IORDER)
C
      RETURN
      END






      SUBROUTINE QKCD   (IDIM, JDIM, KDIM, I, Q, QI, QC, DQC,
     1                   IBCDAT, NBCS, ITURB, IORDER)
C
      RETURN
      END
