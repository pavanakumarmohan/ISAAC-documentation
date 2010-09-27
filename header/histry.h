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
C     $Date: 2001/10/29 03:25:26 $
C     $State: Exp $
C




C
C     Common that holds the residual history for this run.
C     This common restricts the run to have 10000 iterations.
C
C     R2ONE  : The residual on the very first iteration of the case
C
      PARAMETER (MXHIST = 40000)
      COMMON /IHSTRY/ ITTOT
      COMMON /HSTRY / R2(MXHIST), R2ONE
