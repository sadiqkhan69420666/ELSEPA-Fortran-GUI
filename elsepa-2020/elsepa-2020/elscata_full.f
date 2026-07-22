CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             RRRRR     AA    DDDDD   IIII    AA    L                  C
C             R    R   A  A   D    D   II    A  A   L                  C
C             R    R  A    A  D    D   II   A    A  L                  C
C             RRRRR   AAAAAA  D    D   II   AAAAAA  L                  C
C             R  R    A    A  D    D   II   A    A  L                  C
C             R   R   A    A  DDDDB   IIII  A    A  LLLLLL             C
C                                                                      C
C                                                    (version 2018).   C
C                                                                      C
C  Numerical solution of the Schrodinger (S) and Dirac (D) radial      C
C  wave equations. Cubic spline interpolation + power series method.   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  References:
C  [1] F. Salvat and R. Mayol,
C      'Accurate numerical solution of Schrodinger and Dirac wave
C      equations for central fields'.
C      Comput. Phys. Commun. 62 (1991) 65-79.
C  [2] F. Salvat, J. M. Fernandez-Varea and W. Williamson, Jr.,
C      'Accurate numerical solution of the radial Schrodinger and Dirac
C      wave equations'.
C      Comput. Phys. Commun. 90 (1995) 151-158.
C  [3] F. Salvat and J.M. Fernandez-Varea,
C      'RADIAL: a FORTRAN subroutine package for the solution of the
C      radial Schrodinger and Dirac wave equations'.
C      Internal report, University of Barcelona, 2018.
C      This document describes the RADIAL subroutine package and its
C      operation. The PDF file is included in the distribution package.
C
C
C  It is assumed that the (central) potential energy V(R) is such that
C  the function R*V(R) is finite for all R and tends to a constant value
C  when R tends to infinity.
C
C****   All quantities are in atomic Hartree units.
C  For electrons and positrons
C     unit of length = A0 = 5.2917721092D-11 m (= Bohr radius),
C     unit of energy = E0 = 27.21138505 eV (= Hartree energy).
C  For particles of mass 'M' (in units of the electron mass) the atomic
C  units of length and energy are
C     unit of length = A0/M,
C     unit of energy = M*E0.
C
C
C  The calling sequence from the main program is:
C
C****   CALL VINT(R,RV,NV)
C
C  This is an initialization subroutine. It determines the natural cubic
C  spline that interpolates the table of values of the function R*V(R)
C  provided by the user.
C   Input arguments:
C     R(I) ..... input potential grid points. They must be in non-
C                decreasing order, i.e. R(I+1).GE.R(I). (Repeated
C                values are interpreted as discontinuities).
C     RV(I) .... R(I) times the potential energy at R=R(I).
C     NV ....... number of points in the table (.LE.NDIM), must be
C                greater than or equal to 4.
C  The R(I) grid _must_ include the origin (R=0), and extend up to
C  radial distances for which the function R*V(R) reaches its (constant)
C  asymptotic value.
C
C  The function RVSPL(R) gives the interpolated value of R*V(R) at R,
C  i.e., the potential function effectively used in the numerical
C  solution. The function VRANGE() gives the range of the 'inner
C  component' of the potential', which is the radius where the Coulomb
C  tail starts.
C
C****   CALL SBOUND(E,EPS,N,L) or DBOUND(E,EPS,N,K)
C
C  These subroutines solve the radial wave equations for bound states.
C   Input arguments:
C     E ........ estimated binding energy (a good initial estimate
C                speeds up the calculation).
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series. The EPS value
C                must be greater than 1.0D-15.
C     N ........ principal quantum number.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: 0.LE.L.LE.N-1, -N.LE.K.LE.N-1, K.NE.0)
C   Output argument:
C     E ........ binding energy.
C
C****   CALL SFREE(E,EPS,PHASE,L,IRWF) or DFREE(E,EPS,PHASE,K,IRWF)
C
C  These subroutines solve the radial wave equations for free states.
C   Input arguments:
C     E ........ kinetic energy.
C     EPS ...... global tolerance, i.e. allowed relative error in the
C                summation of the radial function series. The EPS value
C                must be greater than 1.0D-15.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: L.GE.0, K.NE.0)
C     IRWF ..... when =0 the radial function is not returned. Serves
C                to avoid unnecessary calculations when only the phase
C                shift is required.
C   Output arguments:
C     PHASE .... inner phase shift (in radians), caused by the short
C                range component of the potential. For modified Coulomb
C                potentials, we have
C                       total phase shift = PHASE + DELTA
C                where DELTA is the Coulomb phase shift (which is
C                delivered through the common block
C                       COMMON/OCOUL/RK,ETA,DELTA  ).
C
C  The values of the radial functions are delivered through the common
C  block
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C  with NDIM=25000 (if a larger number of grid points is needed, edit
C  the present source file and change the value of the parameter NDIM in
C  module CONSTANTS). The grid of radii RAD(I), where the radial wave
C  functions are tabulated, can be arbitrarily selected by the user. The
C  quantities
C     RAD(I) ... user's radial grid,
C     NGP ...... number of grid points (.LE.NDIM),
C  must be defined before calling the solution subroutines. Although
C  it is advisable to have the radial grid points sorted in increasing
C  order, this is not strictly necessary; the solution subroutines do
C  not alter the ordering of the input radii. The output quantities are:
C     P(I) ..... value of the radial function P(R) at the I-th grid
C                point, R=RAD(I).
C     Q(I) ..... value of the radial function Q(R) at the I-th grid
C                point (= P'(R) for Schrodinger particles).
C     ILAST .... *** Bound states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are set equal to 0.0D0.
C                *** Free states: for R.GT.RAD(ILAST), P(R) and Q(R)
C                are obtained in terms of the regular and irregular
C                asymptotic Coulomb functions as
C                  P(R)=COS(PHASE)*FU(R)+SIN(PHASE)*GU(R)
C                  Q(R)=COS(PHASE)*FL(R)+SIN(PHASE)*GL(R)
C                where FU, GU and FL, GL are calculated by subroutines
C                SCOULF and DCOULF with Z=RV(NV). When the absolute
C                value of RV(NV) is less than EPS, Z is set equal to
C                zero, so that the functions FU, GU, FL and GL then
C                reduce to spherical Bessel functions of integer order
C                (which are calculated by function SBESJN).
C     IER ...... error code. A value larger than zero indicates that
C                some fatal error has been found during the calculation.
C
C
C  Bound state wave functions are normalized to unity. The adopted
C  normalization for free states is such that P(R) oscillates with unit
C  amplitude in the asymptotic region (r --> infinity).
C
C
C****   Error codes (and tentative solutions...):
C    0 .... everything is OK.
C    1 .... EMIN.GE.0 in 'BOUND' (Use a denser grid. If the error
C           persists then probably such a bound state does not exist).
C    2 .... E=0 in 'BOUND' (Probably this bound state does not exist).
C    3 .... RAD(NGP) is too small in 'BOUND' (Extend the grid to larger
C           radii).
C    4 .... several zeros of P(R) in a single interval in 'BOUND' (Use
C           a denser grid).
C    5 .... E out of range in 'BOUND' (Accumulated round-off errors?).
C    6 .... RV(NGP)<<0 OR E>0 in 'BOUND' (Check the input potential
C           values).
C    7 .... E.LT.0.0001 in 'FREE'.
C    8 .... RAD(NGP) is too small in 'FREE' (Extend the grid to larger
C           radii). The subroutine tries to find a value of the outer
C           radius RAD(NGP) where it can match the inner and asymptotic
C           solutions.
C
C  The program stops when the input quantum numbers are out of range.
C
C  Some of the RADIAL subroutines generate output files through UNIT 33,
C  which is opened and closed within each individual subroutine. It is
C  advisable not to use this UNIT in the calling program.
C
C  NOTE: The present source file implements the theory described in the
C  accompanying manual, ref. [3] (see above). Numbers in parenthesis in
C  comment lines, with the format (ME-s.nn), indicate the relevant
C  equations in that manual.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  COMPLEX OPTICAL POTENTIALS:
C
C  Radial wave functions of free states for a complex optical potential
C  of the type V(R)+SQRT(-1)*W(R) can be calculated by using the
C  subroutines ZVINT and ZSFREE or ZDFREE.
C
C  It is assumed that V(R) is a modified Coulomb potential, i.e., such
C  that R*V(R) is finite for all R and tends to a constant value when R
C  tends to infinity. W(R) is a finite-range negative (i.e., absorptive)
C  potential such that R*W(R) is finite everywhere.
C
c  The calling sequence from the main program is:
C
C****   CALL ZVINT(R,RV,RW,NV)
C
C  This is the initialization subroutine for complex potentials. It
C  determines the natural cubic splines that interpolate the tables of
C  values of the functions R*V(R) and R*W(R) provided by the user.
C   Input arguments:
C     R(I) ..... input potential grid points. They must be in non-
C                decreasing order, i.e. R(I+1).GE.R(I). (Repeated
C                values are interpreted as discontinuities).
C     RV(I) .... R(I) times the real potential V(R) at R=R(I).
C     RW(I) .... R(I) times the imaginary potential W(R) at R=R(I).
C     NV ....... number of points in the table (.LE.NDIM), must be
C                greater than or equal to 4.
C  The R(I) grid _must_ include the origin (R=0), and extend up to
C  radial distances at which the function R*V(R) reaches its constant
C  asymptotic value and R*W(R) vanishes.
C
C  The subroutine ZRVSPL(R,RVS,RWS) gives the interpolated values of
C  R*V(R) and R*W(R) at R, i.e., the potential functions effectively
C  used in the numerical solution.
C
C****   CALL ZSFREE(E,EPS,PHASER,PHASEI,L,IRWF)    (Schrodinger)
C  or
C       CALL ZDFREE(E,EPS,PHASER,PHASEI,K,IRWF)    (Dirac)
C
C  These subroutines solve the radial wave equations for free states.
C   Input arguments:
C     E ........ kinetic energy.
C     EPS ...... global tolerance. Must be larger than 1.0D-15.
C     L ........ orbital angular momentum quantum number.
C     K ........ relativistic angular momentum quantum number, kappa.
C                (note: L.GE.0, K.NE.0)
C     IRWF ..... when =0 the radial functions are not returned. Serves
C                to avoid unnecessary calculations when only the phase
C                shift is required.
C   Output arguments:
C     PHASER and PHASEI .... real and imaginary parts of the inner phase
C                shift (in rad) caused by the short range component of
C                of the potential. Notice that
C                  total phase shift = PHASER + SQRT(-1)*PHASEI + DELTA
C                where DELTA is the Coulomb phase shift (which is
C                delivered through the common block
C                  COMMON/OCOUL/RK,ETA,DELTA  ).
C
C  The values of the radial functions are delivered through the common
C  blocks
C     COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C     COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
C  The input quantities
C     RAD(I) ... user's radial grid,
C     NGP ...... number of grid points (.LE.NDIM),
C  must be defined before calling the solution subroutines.
C  The output quantities are:
C     P(I),PIM(I) .... real and imaginary parts of the radial function
C                P(R) at the I-th grid point, R=RAD(I).
C     Q(I), QIM(I) .... real and imaginary parts of the radial function
C                function Q(R) at the I-th grid point [= P'(R) for
C                Schrodinger particles].
C     ILAST .... for R.GT.RAD(ILAST), P(R) and Q(R) are obtained in
C                in terms of the regular and irregular asymptotic
C                Coulomb functions as described in the manual.
C     IER ...... error code (the same as for the SFREE and DFREE
C                subroutines).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE CONSTANTS  ! Array dimensions and physical constants.
      SAVE  ! Saves all items in the module.
C  ****  Maximum radial grid dimension.
      INTEGER*4, PARAMETER :: NDIM=25000
C  ****  Maximum number of terms in asymptotic series.
      INTEGER*4, PARAMETER :: MNT=50
C  ----  Speed of light (1/alpha).
      DOUBLE PRECISION, PARAMETER :: SL=137.035999139D0
C  ----  Bohr radius (cm).
      DOUBLE PRECISION, PARAMETER :: A0B=5.2917721067D-9
C  ----  Hartree energy (eV).
      DOUBLE PRECISION, PARAMETER :: HREV=27.21138602D0
C  ----  Electron rest energy (eV).
      DOUBLE PRECISION, PARAMETER :: REV=510.9989461D3
      END MODULE CONSTANTS
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  *********************************************************************
C                       SUBROUTINE VINT
C  *********************************************************************
      SUBROUTINE VINT(R,RV,NV)
C
C     Natural cubic spline interpolation of R*V(R) from the input radii
C  and potential values (ME-4.3).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/RGRID/X(NPTG),RT(NPTG),VT(NPTG),IND(NPTG),NRT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION R(NV),RV(NV)
C
      IF(NV.GT.NDIM) THEN
        WRITE(6,2101) NV,NDIM
 2101   FORMAT(1X,'*** Error in VINT: input potential grid with NV = ',
     1    I5,' data points.',/5X,'NV must be less than NDIM = ',I5,'.')
        STOP
      ENDIF
      IF(NV.LT.4) THEN
        WRITE(6,2102) NV
 2102   FORMAT(1X,'*** Error in VINT: the input potential grid must ',
     1    /5X,'have more than 4 data points. NV =',I5,'.')
        STOP
      ENDIF
      IF(R(1).LT.0.0D0) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in VINT: R(1).LT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
      IF(R(1).GT.1.0D-15) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in VINT: R(1).GT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
C
      RT(1)=0.0D0
      VT(1)=RV(1)
      DO I=2,NV
        RT(I)=R(I)
        VT(I)=RV(I)
        IF(RT(I-1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I-1)))) THEN
          WRITE(6,2105)
 2105     FORMAT(1X,'*** Error in VINT: R values in',
     1      'decreasing order.',
     2      /5X,'Details in file ''VINT-error.dat''.')
          OPEN(33,FILE='VINT-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,NV
              WRITE(33,'(I5,1P,2E18.10)') J,R(J),RV(J)
            ENDDO
          CLOSE(33)
          STOP 'VINT: R values in decreasing order.'
        ENDIF
      ENDDO
C
C  ****  Coulomb tail.
C
      ZINF=RV(NV)
      TOL=MAX(ABS(ZINF)*1.0D-10,1.0D-10)
      NVI=NV
      DO I=NV,4,-1
        IF(ABS(VT(I-1)-ZINF).GT.TOL) THEN
          NVE=I
          IF(RT(I)-RT(I-1).LT.EPS*MAX(ABS(RT(I-1)),ABS(RT(I)))) THEN
            RT(I)=RT(I-1)
            VT(I)=ZINF
            NVE=I-1
          ELSE
            NVI=NVI+1  ! Add a discontinuity.
            DO J=NVI,NVE+1,-1
              RT(J)=RT(J-1)
            ENDDO
            VT(NVE+1)=ZINF
          ENDIF
          GO TO 10
        ELSE
          VT(I)=ZINF
        ENDIF
      ENDDO
      NVE=4
 10   CONTINUE
C
C  ****  Natural cubic spline interpolation, piecewise.
C
      IO=0
      I=0
      K=0
 1    I=I+1
      K=K+1
      X(K)=RT(I)
      Y(K)=VT(I)
      IF(I.EQ.NVE) GO TO 2
C  ****  Duplicated points are considered as discontinuities.
      IF(RT(I+1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I+1)))) GO TO 1
 2    CONTINUE
C
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      DO J=1,K-1
        IO=IO+1
        RG(IO)=X(J)
        RVG(IO)=Y(J)
        VA(IO)=A(J)
        VB(IO)=B(J)
        VC(IO)=C(J)
        VD(IO)=D(J)
      ENDDO
      IF(I.LT.NVE) THEN
        K=0
        GO TO 1
      ENDIF
C  ****  The last set of coefficients of the spline is replaced by those
C        of the potential tail.
      IO=IO+1
      NVT=IO
      RG(IO)=X(K)
      RVG(IO)=ZINF
      VA(IO)=ZINF
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RVSPL
C  *********************************************************************
      FUNCTION RVSPL(R)
C
C     This function gives the (natural cubic spline) interpolated value
C  of R*V(R) at R, i.e., the potential energy function effectively used
C  in the numerical solution.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C
      IF(R.LT.0.0D0) THEN
        RVSPL=0.0D0
      ELSE
        RVSPL=SPLVAL(R,RG,VA,VB,VC,VD,NVT)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION VRANGE
C  *********************************************************************
      FUNCTION VRANGE()
C
C     This function gives the range of the 'inner component' of the
C  potential, i.e., the radius where the Coulomb tail starts.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C
      VRANGE=RG(NVT)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SBOUND
C  *********************************************************************
      SUBROUTINE SBOUND(E,EPS,N,L)
C
C     This subroutine solves the Schrodinger radial equation for bound
C  states.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C  ****  Set IWR=1 to print partial results.
      DATA IWR/0/
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in SBOUND: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SBOUND: N.LT.1.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in SBOUND: L.LT.0.')
        STOP
      ENDIF
C
      IF(L.GE.N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in SBOUND: L.GE.N.')
        STOP
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=10.0D0*EPS
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2104) NDIM
 2104   FORMAT(1X,'*** Error in SBOUND: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (ME-5.1)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 in SBOUND: EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(RV(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
      ISUM=0
C
C  ************  New shot.
C
 2    CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 in SBOUND: E=0.',
     1    /5X,'(Probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (ME-5.2)
      DO K=2,NRT
        IOTP=NRT+2-K
        RN=R(IOTP)
        J=IND(IOTP)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
 3    CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 in SBOUND: RAD(NGP) is too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 in SBOUND: Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
        IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
        EMAX=E
        ICMAX=1
        E=0.5D0*(EMIN+E)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,L
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      IF(ISUM.EQ.0) THEN
        ISUM=1
        CALL SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
      ENDIF
      PO=PT(IOTP)
      QO=QT(IOTP)
C
C  ****  Inward solution.
C
      CALL SINW(E,EPS,SUMIN,L,IOTP,ISUM)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/PT(IOTP)
      DO I=IOTP,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
      ENDDO
      SUMIN=SUMIN*FACT**2
      QI=QT(IOTP)
      RLAST=R(ILAST)
C  ****  Normalization.
      SUM=SUMIN+SUMOUT
C  ****  Eigenvalue correction. (ME-5.11)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=PO*(QO-QI)/(SUM+SUM)
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,L
 2000   FORMAT(/2X,'Subroutine SBOUND.   N =',I3,'   L =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GE.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 in SBOUND: E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalization
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN  ! Prevents -0.0D0 values.
        DO I=ILAST+1,NRT
          PT(I)=0.0D0
          QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(PT(NRT)).GT.1.0D-5*ABS(PT(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBOUND
C  *********************************************************************
      SUBROUTINE DBOUND(E,EPS,N,K)
C
C     This subroutine solves the Dirac radial equation for bound states.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
C  ****  Set IWR=1 to print partial results.
      DATA IWR/0/
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in DBOUND: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(N.LT.1) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DBOUND: N.LT.1.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DBOUND: K.EQ.0.')
        STOP
      ENDIF
C
      IF(K.LT.-N) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in DBOUND: K.LT.-N.')
        STOP
      ENDIF
C
      IF(K.GE.N) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in DBOUND: K.GE.N.')
        STOP
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        ELSE
        L=K
      ENDIF
C  ****  Radial quantum number.
      NR=N-L-1
C
      DELL=10.0D0*EPS
      IF(E.GT.-1.0D-1) E=-1.0D-1
      FL1=0.5D0*L*(L+1)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2105) NDIM
 2105   FORMAT(1X,'*** Error in DBOUND: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Minimum of the effective radial potential. (ME-5.1)
C
      EMIN=1.0D0
      DO I=2,NRT
        RN=R(I)
        J=IND(I)
        RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
        EMIN=MIN(EMIN,(RVN+FL1/RN)/RN)
      ENDDO
      IF(EMIN.GT.-1.0D-35) THEN
        IER=1
        WRITE(6,1001)
 1001   FORMAT(1X,'*** Error 1 in DBOUND: EMIN.GE.0.'/5X,
     1    '(Use a denser grid. If the error persists then probably',
     2    /6X,'such a bound state does not exist).')
        RETURN
      ENDIF
      RN=R(NRT)
      J=IND(NRT)
      RVN=VA(J)+RN*(VB(J)+RN*(VC(J)+RN*VD(J)))
      ESUP=(RVN+FL1/RN)/RN
      IF(ESUP.GT.0.0D0) ESUP=0.0D0
      IF(L.EQ.0) THEN
        TEHYDR=-(RV(1)/N)**2
        EMIN=MIN(E,TEHYDR,-10.0D0)
      ENDIF
      IF(EMIN.LT.-SL*SL) EMIN=-SL*SL
      IF(E.GT.ESUP.OR.E.LT.EMIN) E=0.5D0*(ESUP+EMIN)
      EMAX=ESUP
      ICMIN=0
      ICMAX=0
      ISUM=0
C
C  ************  New shot.
C
 2    CONTINUE
      IF(E.GT.-1.0D-16) THEN
        IER=2
        WRITE(6,1002)
 1002   FORMAT(1X,'*** Error 2 in DBOUND: E=0.',
     1    /5X,'(probably this bound state does not exist).')
        RETURN
      ENDIF
C  ****  Outer turning point. (ME-5.2)
      DO J=2,NRT
        IOTP=NRT+2-J
        RN=R(IOTP)
        JJ=IND(IOTP)
        RVN=VA(JJ)+RN*(VB(JJ)+RN*(VC(JJ)+RN*VD(JJ)))
        EKIN=E-(RVN+FL1/RN)/RN
        IF(EKIN.GT.0.0D0) GO TO 3
      ENDDO
 3    CONTINUE
      IOTP=IOTP+1
      IF(IOTP.GT.NRT-1) THEN
        IER=3
        WRITE(6,1003)
 1003   FORMAT(1X,'*** Error 3 in DBOUND: RAD(NGP) is too small.'
     1    /5X,'(Extend the grid to larger radii).')
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      IF(NZMAX.GT.1.AND.NZERO.LE.NR) THEN
        IER=4
        WRITE(6,1004)
 1004   FORMAT(1X,'*** Error 4 in DBOUND: Several zeros of P(R)',
     1    /5X,'in a single interval (Use a denser grid).')
        RETURN
      ENDIF
C
C  ****  Too many nodes.
      IF(NZERO.GT.NR) THEN
        IF(ICMIN.EQ.0) EMIN=EMIN-2.0D0*(EMAX-EMIN)
        EMAX=E
        ICMAX=1
        E=0.5D0*(EMIN+E)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,K
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  Too few nodes.
      IF(NZERO.LT.NR) THEN
        ICMIN=1
        EMIN=E
        E=0.5D0*(E+EMAX)
        IF(IWR.NE.0) THEN
          WRITE(6,2000) N,K
          WRITE(6,2001) NR,NZERO,IOTP,NRT
          WRITE(6,2002) E
          WRITE(6,2004) EMIN,EMAX
        ENDIF
C
        IF(EMAX-EMIN.LT.DELL*ABS(EMIN)) THEN
          IER=5
          WRITE(6,1005)
          RETURN
        ENDIF
C
        GO TO 2
      ENDIF
C  ****  The correct number of nodes has been found.
      IF(ISUM.EQ.0) THEN
        ISUM=1
        CALL DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
      ENDIF
      PO=PT(IOTP)
      QO=QT(IOTP)
C
C  ****  Inward solution.
C
      CALL DINW(E,EPS,SUMIN,K,IOTP,ISUM)
      IF(IER.GT.0) RETURN
C  ****  Matching of the outward and inward solutions.
      FACT=PO/PT(IOTP)
      DO I=IOTP,NRT
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
      ENDDO
      SUMIN=SUMIN*FACT**2
      QI=QT(IOTP)
      RLAST=R(ILAST)
C  ****  Normalization integral.
      SUM=SUMIN+SUMOUT
C  ****  Eigenvalue correction. (ME-5.15)
      IF(SUM.LT.1.0D-15) SUM=1.0D0
      DE=SL*PO*(QO-QI)/SUM
      EP=E+DE
C
      IF(DE.LT.0.0D0) THEN
        ICMAX=1
        EMAX=E
      ENDIF
C
      IF(DE.GT.0.0D0) THEN
        ICMIN=1
        EMIN=E
      ENDIF
C
      IF(ICMIN.EQ.0.AND.EP.LT.EMIN) THEN
        EMIN=1.1D0*EMIN
        IF(EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      ENDIF
      IF(ICMIN.EQ.1.AND.EP.LT.EMIN) EP=0.5D0*(E+EMIN)
      IF(ICMAX.EQ.1.AND.EP.GT.EMAX) EP=0.5D0*(E+EMAX)
      IF(EP.GT.ESUP) EP=0.5D0*(ESUP+E)
C
      IF(IWR.NE.0) THEN
        WRITE(6,2000) N,K
 2000   FORMAT(/2X,'Subroutine DBOUND.   N =',I3,'   K =',I3)
        WRITE(6,2001) NR,NZERO,IOTP,NRT
 2001   FORMAT(2X,'NR =',I3,'   NZERO =',I3,'   IOTP = ',I5,
     1    '   NGP =',I5)
        WRITE(6,2002) EP
 2002   FORMAT(2X,'E new = ',1P,D22.15)
        WRITE(6,2003) E,DE
 2003   FORMAT(2X,'E old = ',1P,D22.15,'   DE = ',D11.4)
        WRITE(6,2004) EMIN,EMAX
 2004   FORMAT(2X,'EMIN = ',1P,D12.5,'   EMAX = ',D12.5)
      ENDIF
C
      IF(EP.GT.ESUP.AND.ABS(E-ESUP).LT.DELL*ABS(ESUP)) THEN
        IER=5
        WRITE(6,1005)
 1005   FORMAT(1X,'*** Error 5 in DBOUND: E out of range.'/5X,
     1    '(Accumulated round-off errors?).')
        RETURN
      ENDIF
      EO=E
      E=EP
      IF(MIN(ABS(DE),ABS(E-EO)).GT.ABS(E*DELL)) GO TO 2
C  ****  Normalization.
      FACT=1.0D0/SQRT(SUM)
      DO I=1,ILAST
        PT(I)=PT(I)*FACT
        QT(I)=QT(I)*FACT
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN  ! Prevents -0.0D0 values.
        DO I=ILAST+1,NRT
          PT(I)=0.0D0
          QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      IF(ABS(PT(NRT)).GT.1.0D-5*ABS(PT(IOTP))) THEN
        IER=3
        WRITE(6,1003)
      ENDIF
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFREE
C  *********************************************************************
      SUBROUTINE SFREE(E,EPS,PHASE,L,IRWF)
C
C     This subroutine solves the Schrodinger radial equation for free
C  states.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in SFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SFREE: L.LT.0.')
        STOP
      ENDIF
      FL1=0.5D0*L*(L+1)
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in SFREE: E.LT.0.0001')
        RETURN
      ENDIF
      RK=SQRT(E+E)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in SFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      ZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(ZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 3
          BNL1=SBESJN(2,L+1,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL=SBESJN(2,L,X)
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+1,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)
          QB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL SCOULF(ZINF,E,L,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (SFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in SFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      ISUM=0
      CALL SOUTW(E,EPS,SUMOUT,L,1,NZERO,ILAST,ISUM)
C
C  ****  Phase shift. (ME-6.8, ME-6.9)
C
      PO=PT(ILAST)
      POP=QT(ILAST)
      PIA=PA(ILAST)
      PIAP=QA(ILAST)
      PIB=PB(ILAST)
      PIBP=QB(ILAST)
C
      PHASE=ATAN2(POP*PIA-PO*PIAP,PO*PIBP-POP*PIB)
      CD=COS(PHASE)
      SD=SIN(PHASE)
      IF(ABS(PO).GT.EPS) THEN
        RNORM=(CD*PIA+SD*PIB)/PO
      ELSE
        RNORM=(CD*PIAP+SD*PIBP)/POP
      ENDIF
C  ****  The normalization factor RNORM is required to be positive.
      IF(RNORM.LT.0.0D0) THEN
        RNORM=-RNORM
        CD=-CD
        SD=-SD
        PHASE=PHASE+PI
      ENDIF
C  ****  The phase shift is reduced to the interval (-PI,PI).
      IF(PHASE.GT.PI-EPS) THEN
        PHASE=PHASE-TPI
      ELSE IF(PHASE.LT.-PI+EPS) THEN
        PHASE=PHASE+TPI
      ENDIF
C
      IF(IRWF.EQ.0) RETURN
C  ****  Normalized wave function. (ME-6.10)
      DO I=1,ILAST
        PT(I)=RNORM*PT(I)
        QT(I)=RNORM*QT(I)
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          PT(I)=CD*PA(I)+SD*PB(I)
          QT(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
          IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFREE
C  *********************************************************************
      SUBROUTINE DFREE(E,EPS,PHASE,K,IRWF)
C
C     This subroutine solves the Dirac radial equation for free states.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
      COMMON/NZT/NZMAX
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in DFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DFREE: K.EQ.0.')
        STOP
      ENDIF
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in DFREE: E.LT.0.0001')
        RETURN
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in DFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      ZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(ZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        FACTOR=SQRT(E/(E+2.0D0*SL*SL))
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*RN*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(RVN).GT.T) GO TO 3
          BNL=SBESJN(2,L,X)
          IF(ABS(BNL).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL1=SBESJN(2,L+KSIGN,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+KSIGN,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=-FACTOR*KSIGN*X*BJL1
          QB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(ZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(ABS(RVN-ZINF).GT.TAS) GO TO 3
          CALL DCOULF(ZINF,E,K,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (DFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in DFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      ISUM=0
      CALL DOUTW(E,EPS,SUMOUT,K,1,NZERO,ILAST,ISUM)
C
C  ****  Phase shift. (ME-6.8, ME-6.9)
C
      RM=R(ILAST)
      IL=IND(ILAST-1)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PO=PT(ILAST)
      POP=-K*PO/RM+FG*QT(ILAST)
      IL=IND(ILAST)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PIA=PA(ILAST)
      PIAP=-K*PIA/RM+FG*QA(ILAST)
      PIB=PB(ILAST)
      PIBP=-K*PIB/RM+FG*QB(ILAST)
C
      PHASE=ATAN2(POP*PIA-PO*PIAP,PO*PIBP-POP*PIB)
C  ****  The phase shift is reduced to the interval (-PI/2,PI/2).
      TT=ABS(PHASE)
      IF(TT.GT.PIH) PHASE=PHASE*(1.0D0-PI/TT)
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.10)
C
      CD=COS(PHASE)
      SD=SIN(PHASE)
      IF(ABS(PO).GT.EPS) THEN
        RNORM=(CD*PIA+SD*PIB)/PO
      ELSE
        RNORM=(CD*PIAP+SD*PIBP)/POP
      ENDIF
C
      DO I=1,ILAST
        PT(I)=RNORM*PT(I)
        QT(I)=RNORM*QT(I)
        IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
        IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          PT(I)=CD*PA(I)+SD*PB(I)
          QT(I)=CD*QA(I)+SD*QB(I)
          IF(ABS(PT(I)).LT.1.0D-99) PT(I)=0.0D0
          IF(ABS(QT(I)).LT.1.0D-99) QT(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=PT(J)
          Q(I)=QT(J)
        ELSE
          P(I)=PT(J+1)
          Q(I)=QT(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE SOUTW
C  *********************************************************************
      SUBROUTINE SOUTW(E,EPS,SUMOUT,L,NR,NZERO,IOTP,ISUM)
C
C     Outward solution of the Schrodinger radial equation for a  piece-
C  wise cubic potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AL=L
      N1=IOTP-1
C
      PT(1)=0.0D0
      QT(1)=0.0D0
      SUMOUT=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I)
        QQI=QT(I)
        CALL SCH(E,AL,EPS,ISUM)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        PT(I+1)=PF
        QT(I+1)=QF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=1,I
          PT(K)=PT(K)*FACT
          QT(K)=QT(K)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMOUT=SUMOUT*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SINW
C  *********************************************************************
      SUBROUTINE SINW(E,EPS,SUMIN,L,IOTP,ISUM)
C
C     Inward solution of the Schrodinger radial equation for a piece-
C  wise cubic potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (TRINF=22500.0D0)  ! TRINF=150.**2
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      AL=L
C  ****  WKB solution at the outer grid point. (ME-5.4, ME-5.5)
      N=NRT
 1    N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=2.0D0*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 in SBOUND: RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (ME-5.6)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=1.0D0-RN*(RVN+RN*RVNP-2*E*RN)/CMU
        PT(N)=1.0D0
        QT(N)=(0.5D0/RN)*CRAT-SQRT(CMU)/RN
        ILAST=N
      ELSE
        PT(N)=0.0D0
        QT(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      SUMIN=0.0D0
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I1)
        QQI=QT(I1)
        CALL SCH(E,AL,EPS,ISUM)
        PT(I)=PF
        QT(I)=QF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=I1,N
            PT(K)=PT(K)*FACT
            QT(K)=QT(K)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMIN=SUMIN*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMIN=SUMIN+RSUM
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DOUTW
C  *********************************************************************
      SUBROUTINE DOUTW(E,EPS,SUMOUT,K,NR,NZERO,IOTP,ISUM)
C
C     Outward solution of the Dirac radial equation for a piecewise
C  cubic potential. Power series method.
C
       USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AK=K
      IF(E.LT.0.0D0) THEN
        N1=NRT
      ELSE
        N1=IOTP-1
      ENDIF
C
      PT(1)=0.0D0
      QT(1)=0.0D0
      SUMOUT=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I)
        QQI=QT(I)
        CALL DIR(E,AK,EPS,ISUM)
        NZERO=NZERO+NCHS
        IF(NCHS.GT.NZMAX) NZMAX=NCHS
        IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
        PT(I+1)=PF
        QT(I+1)=QF
        IF(E.LT.0.0D0) THEN
C  ****  TCONV is the product of P and its second derivative at the
C        I-th grid point (positive if P is convex).
          TCONV=2.0D0*CA(3)*PPI
          IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
            IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
            IOTP=I+1
            RETURN
          ENDIF
        ENDIF
        IF(I.EQ.1) THEN
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
          GO TO 1
        ENDIF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO J=1,I
            PT(J)=PT(J)*FACT
            QT(J)=QT(J)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMOUT=SUMOUT*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMOUT=SUMOUT+RSUM
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DINW
C  *********************************************************************
      SUBROUTINE DINW(E,EPS,SUMIN,K,IOTP,ISUM)
C
C     Inward solution of the Dirac radial equation for a piecewise cubic
C  potential. Power series method.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (TRINF=22500.0D0)  ! TRINF=150.**2
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),PT(NPTG),QT(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
      ELSE
        L=K
      ENDIF
      AK=K
      AL=L
C  ****  WKB solution at the outer grid point. (ME-5.13, ME-5.14)
      N=NRT
      FACT=(E+2.0D0*SL*SL)/(SL*SL)
 1    N1=IND(N-1)
      RN=R(N)
      RVN=VA(N1)+RN*(VB(N1)+RN*(VC(N1)+RN*VD(N1)))
      RVNP=VB(N1)+RN*(2.0D0*VC(N1)+RN*3.0D0*VD(N1))
      CMU=FACT*RN*(RVN-E*RN)+AL*(AL+1)
      IF(CMU.LE.0.0D0) THEN
        IER=6
        WRITE(6,1006)
 1006   FORMAT(1X,'*** Error 6 in DBOUND: RV(NGP)<<0 OR E>0.',
     1    /5X,'(Check the input potential values).')
        RETURN
      ENDIF
C  ****  Practical infinity. (ME-5.2)
      IF(CMU.LT.TRINF.OR.N.EQ.IOTP+1) THEN
        CRAT=(0.5D0-SQRT(CMU))/RN-0.25D0*FACT*(RVN+RN*RVNP
     1      -2.0D0*RN*E)/CMU
        PT(N)=1.0D0
        QT(N)=SL*(CRAT+AK/RN)/(E+2.0D0*SL*SL)
        ILAST=N
      ELSE
        PT(N)=0.0D0
        QT(N)=0.0D0
        N=N-1
        GO TO 1
      ENDIF
C
      SUMIN=0.0D0
      N1=N-IOTP
      DO J=1,N1
        I=N-J
        I1=I+1
        RA=R(I1)
        RB=R(I)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        PPI=PT(I1)
        QQI=QT(I1)
        CALL DIR(E,AK,EPS,ISUM)
        PT(I)=PF
        QT(I)=QF
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO M=I1,N
            PT(M)=PT(M)*FACT
            QT(M)=QT(M)*FACT
          ENDDO
          IF(ISUM.EQ.1) SUMIN=SUMIN*FACT**2+RSUM
        ELSE
          IF(ISUM.EQ.1) SUMIN=SUMIN+RSUM
        ENDIF
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH
C  *********************************************************************
      SUBROUTINE SCH(E,AL,EPS,ISUM)
C
C     This subroutine solves the Schrodinger radial equation for a
C  central potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the radial function
C  and its derivative) at RA, the solution in the interval between RA
C  and RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AL .................... orbital angular momentum quantum number.
C     ISUM .................. normalization flag.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common SINOUT):
C     RA, RB ................ interval end points (input),
C     PPI, QQI .............. values of the radial function and its
C                             derivative at RA (input),
C     PF, QF ................ values of the radial function and its
C                             derivative at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     RSUM .................. normalization integral,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,S,NSUM
      NCHS=0
      RLN=0.0D0
      RSUM=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
        DIRECT=-1.0D0
      ELSE
        DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PPI
      Q1=QQI
 1    CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL SCH0(E,AL,EPS)
      HP=H
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C  ****  Normalization integral (ME-5.17)
      IF(ISUM.EQ.1) THEN
        RSUMI=0.0D0
        CDEMS=2.0D0*S-1.0D0
        DO I1=1,NSUM
          DO I2=1,NSUM
            RSUMI=RSUMI+CA(I1)*CA(I2)/(CDEMS+DBLE(I1+I2))
          ENDDO
        ENDDO
        RSUM=RSUM+DIRECT*HP*RSUMI
      ENDIF
C
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
        RSUM=RSUM/TST**2
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCH0
C  *********************************************************************
      SUBROUTINE SCH0(E,AL,EPS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/SSAVE/P0,Q0,P1,Q1,CA(60),R0,R1,S,NSUM
C
      RVE=RV1-E
      IF(R0.GT.1.0D-10) GO TO 2
C
C  ****  First interval. (ME-4.15 to ME-4.18)
C
      S=AL+1
      U0=AL*S
      U1=2*RV0*R1
      U2=2*RVE*R1**2
      U3=2*RV2*R1**3
      U4=2*RV3*R1**4
      UT=U0+U1+U2+U3+U4
C
      CA(1)=1.0D0
      CA(2)=U1*CA(1)/((S+1)*S-U0)
      CA(3)=(U1*CA(2)+U2*CA(1))/((S+2)*(S+1)-U0)
      CA(4)=(U1*CA(3)+U2*CA(2)+U3*CA(1))
     1  /((S+3)*(S+2)-U0)
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      Q1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      P2P1=S*(S-1)*CA(1)+(S+1)*S*CA(2)+(S+2)*(S+1)*CA(3)
     1  +(S+3)*(S+2)*CA(4)
C
      DO I=5,60
        K=I-1
        CA(I)=(U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)+U4*CA(I-4))
     1    /((S+K)*(S+K-1)-U0)
        P1=P1+CA(I)
        DQ1=(S+K)*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+(S+K-1)*DQ1
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*(P2P1-UT*P1))
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 1
      ENDDO
C  ****  Renormalization. (ME-4.18)
 1    CONTINUE
      NSUM=K+1
      Q1=Q1/(ABS(P1)*R1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (ME-4.10 to ME-4.13)
C
 2    CONTINUE
      S=0.0D0
      H=R1-R0
      H2=H*H
C
      RHO=H/R0
      U0=AL*(AL+1)+2*R0*(RV0+R0*(RVE+R0*(RV2+R0*RV3)))
      U1=2*(RV0+R0*(2*RVE+R0*(3*RV2+R0*4*RV3)))*H
      U2=2*(RVE+R0*(3*RV2+R0*6*RV3))*H2
      U3=2*(RV2+R0*4*RV3)*H2*H
      U4=2*RV3*H2*H2
      UT=U0+U1+U2+U3+U4
C
      CA(1)=P0
      CA(2)=Q0*H
      CA(3)=RHO*RHO*U0*CA(1)/2
      CA(4)=RHO*(RHO*(U0*CA(2)+U1*CA(1))-4*CA(3))/6
      CAK=(U0-2)*CA(3)+U1*CA(2)+U2*CA(1)
      CA(5)=RHO*(RHO*CAK-12*CA(4))/12
      CAK=(U0-6)*CA(4)+U1*CA(3)+U2*CA(2)+U3*CA(1)
      CA(6)=RHO*(RHO*CAK-24*CA(5))/20
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)+CA(5)+CA(6)
      Q1=CA(2)+2*CA(3)+3*CA(4)+4*CA(5)+5*CA(6)
      P2P1=2*CA(3)+6*CA(4)+12*CA(5)+20*CA(6)
C
      DO I=7,60
        K=I-1
        CAK=(U0-(K-2)*(K-3))*CA(I-2)+U1*CA(I-3)+U2*CA(I-4)
     1    +U3*CA(I-5)+U4*CA(I-6)
        CA(I)=RHO*(RHO*CAK-2*(K-1)*(K-2)*CA(K))/(K*(K-1))
        P1=P1+CA(I)
        DQ1=K*CA(I)
        Q1=Q1+DQ1
        P2P1=P2P1+K*(K-1)*CA(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(P2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(CA(I))
        T2=ABS(R1*R1*P2P1-H2*UT*P1)
        TST1=EPS*MAX(ABS(P1),ABS(Q1)/I)
        TST2=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 3
      ENDDO
C
 3    NSUM=K+1
      Q1=Q1/H
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR
C  *********************************************************************
      SUBROUTINE DIR(E,AK,EPS,ISUM)
C
C     This subroutine solves the Dirac radial equation for a central
C  potential V(R) such that
C             R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AK .................... relativistic angular momentum quantum
C                             number,
C     ISUM .................. normalization flag.
C
C  Input (common POTEN):
C     RV0, RV1, RV2, RV3 .... potential parameters.
C
C  Input-output (common DINOUT):
C     RA, RB ................ interval end points (input),
C     PPI, QQI .............. values of the large and small radial
C                             functions at RA (input),
C     PF, QF ................ values of the large and small radial
C                             functions at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     EPS ................... estimate of the global error in PF and QF,
C     RSUM .................. normalization integral,
C     NSTEP ................. number of steps,
C     NCHS .................. number of zeros of P(R) in (RA,RB).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PPI,QQI,PF,QF,RA,RB,RLN,RSUM,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
      NCHS=0
      RLN=0.0D0
      RSUM=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      P1=PPI
      Q1=QQI
 1    CONTINUE
      R0=R1
      P0=P1
      Q0=Q1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL DIR0(E,AK,EPS)
      HP=H
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C  ****  Normalization integral (ME-5.18)
      IF(ISUM.EQ.1) THEN
        RSUMI=0.0D0
        CDEMS=2.0D0*S-1.0D0
        CDEMT=2.0D0*(S+T)-1.0D0
        DO I1=1,NSUM
          DO I2=1,NSUM
            RSUMI=RSUMI+CA(I1)*CA(I2)/(CDEMS+DBLE(I1+I2))
     1        +CB(I1)*CB(I2)/(CDEMT+DBLE(I1+I2))
          ENDDO
        ENDDO
        RSUM=RSUM+DIRECT*HP*RSUMI
      ENDIF
C
      NSTEP=NSTEP+1
      TST=ABS(P1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        P1=P1/TST
        Q1=Q1/TST
        RSUM=RSUM/TST**2
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      PF=P1
      QF=Q1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DIR0
C  *********************************************************************
      SUBROUTINE DIR0(E,AK,EPS)
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,S,T,NSUM
C
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      RVE=RV1-E
C
      IF(R0.GT.1.0D-10) GO TO 4
C
C  ****  First interval. (ME-4.28 to ME-4.39)
C
      U0=RV0/SL
      U1=RVE*R1/SL
      U2=RV2*R1**2/SL
      U3=RV3*R1**3/SL
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
      UH=U1-2*SL*R1
      IF(ABS(U0).LT.1.0D-10) GO TO 1
C
C  ****  U0.NE.0. (ME-4.30 to ME-4.34)
      S=SQRT(AK*AK-U0*U0)
      T=0.0D0
      DS=S+S
      CA(1)=1.0D0
      CB(1)=-(S+AK)/U0
      CAI=U1*CA(1)
      CBI=UH*CB(1)
      CA(2)=(-U0*CAI-(S+1-AK)*CBI)/(DS+1)
      CB(2)=((S+1+AK)*CAI-U0*CBI)/(DS+1)
      CAI=U1*CA(2)+U2*CA(1)
      CBI=UH*CB(2)+U2*CB(1)
      CA(3)=(-U0*CAI-(S+2-AK)*CBI)/(2*(DS+2))
      CB(3)=((S+2+AK)*CAI-U0*CBI)/(2*(DS+2))
      P1=CA(1)+CA(2)+CA(3)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)
      Q1=CB(1)+CB(2)+CB(3)
      QP1=S*CB(1)+(S+1)*CB(2)+(S+2)*CB(3)
C
      DO I=4,60
        K=I-1
        CAI=U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)
        CBI=UH*CB(K)+U2*CB(I-2)+U3*CB(I-3)
        CA(I)=(-U0*CAI-(S+K-AK)*CBI)/(K*(DS+K))
        CB(I)=((S+K+AK)*CAI-U0*CBI)/(K*(DS+K))
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=1. (ME-4.35, ME-4.36)
 1    CONTINUE
      IF(ISIG.LT.0) GO TO 2
      S=ABS(AK)
      T=1.0D0
      DS1=S+S+1
      CA(1)=1.0D0
      CB(1)=U1*CA(1)/DS1
      CA(2)=0.0D0
      CB(2)=U2*CA(1)/(DS1+1)
      CA(3)=-UH*CB(1)/2
      CB(3)=(U1*CA(3)+U3*CA(1))/(DS1+2)
      CA(4)=-(UH*CB(2)+U2*CB(1))/3
      CB(4)=(U1*CA(4)+U2*CA(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S+1)*CB(1)+(S+2)*CB(2)+(S+3)*CB(3)
C
      DO I=5,60
        K=I-1
        CA(I)=-(UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K
        CB(I)=(U1*CA(I)+U2*CA(K)+U3*CA(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+I)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=-1. (ME-4.37, ME-4.38)
 2    CONTINUE
      S=ABS(AK)+1
      T=-1.0D0
      DS1=S+ABS(AK)
      IF(UH.GT.0.0D0) THEN
        CB(1)=-1.0D0
      ELSE
        CB(1)=1.0D0
      ENDIF
      CA(1)=-UH*CB(1)/DS1
      CB(2)=0.0D0
      CA(2)=-U2*CB(1)/(DS1+1)
      CB(3)=U1*CA(1)/2
      CA(3)=-(UH*CB(3)+U3*CB(1))/(DS1+2)
      CB(4)=(U1*CA(2)+U2*CA(1))/3
      CA(4)=-(UH*CB(4)+U2*CB(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S-1)*CB(1)+S*CB(2)+(S+1)*CB(3)
C
      DO I=5,60
        K=I-1
        CB(I)=(U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K
        CA(I)=-(UH*CB(I)+U2*CB(K)+U3*CB(I-2))/(DS1+K)
        P1=P1+CA(I)
        PP1=PP1+(S+K)*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+(S+K-1)*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
C  ****  Renormalization. (ME-4.39)
 3    CONTINUE
      NSUM=K+1
      Q1=Q1/ABS(P1)
      P1=P1/ABS(P1)
      RETURN
C
C  ****  Middle region. (ME-4.23 to ME-4.27)
C
 4    CONTINUE
      S=0.0D0
      T=0.0D0
      RHO=H/R0
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/SL
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/SL
      U2=(RV2+R0*3*RV3)*H2/SL
      U3=RV3*H*H2/SL
      UB=U0-2*SL*R0
      UH=U1-2*SL*H
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
C
      CA(1)=P0
      CB(1)=Q0
      CA(2)=-RHO*(AK*CA(1)+UB*CB(1))
      CB(2)=RHO*(AK*CB(1)+U0*CA(1))
      CA(3)=-RHO*((AK+1)*CA(2)+UB*CB(2)+UH*CB(1))/2
      CB(3)=RHO*((AK-1)*CB(2)+U0*CA(2)+U1*CA(1))/2
      CA(4)=-RHO*((AK+2)*CA(3)+UB*CB(3)+UH*CB(2)+U2*CB(1))/3
      CB(4)=RHO*((AK-2)*CB(3)+U0*CA(3)+U1*CA(2)+U2*CA(1))/3
C
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=CA(2)+2*CA(3)+3*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=CB(2)+2*CB(3)+3*CB(4)
C
      DO I=5,60
        K=I-1
        CA(I)=-RHO*((AK+K-1)*CA(K)+UB*CB(K)+UH*CB(I-2)+U2*CB(I-3)
     1       +U3*CB(I-4))/K
        CB(I)=RHO*((AK-K+1)*CB(K)+U0*CA(K)+U1*CA(I-2)+U2*CA(I-3)
     1       +U3*CA(I-4))/K
        P1=P1+CA(I)
        PP1=PP1+K*CA(I)
        Q1=Q1+CB(I)
        QP1=QP1+K*CB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(P1),ABS(Q1),ABS(PP1),ABS(QP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*PP1+H*(AK*P1+UQ*Q1))
        T1B=ABS(R1*QP1-H*(AK*Q1+UT*P1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(CA(I)),ABS(CB(I)))
        TST=EPS*MAX(ABS(P1),ABS(Q1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 5
      ENDDO
C
 5    CONTINUE
      NSUM=K+1
      RETURN
      END


C  *********************************************************************
C                       SUBROUTINE ZVINT
C  *********************************************************************
      SUBROUTINE ZVINT(R,RV,RW,NV)
C
C     Natural cubic spline interpolation for R*V(R) from the input radii
C  and potential values. (ME-4.3)
C
C  ****  Complex potential; R*V(R)=RV+SQRT(-1)*RW
C        It is assumed that the imaginary part RW has a finite range.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/STORE/X(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      COMMON/ZSTORE/RT(NPTG),VT(NPTG),WT(NPTG),Y1(NPTG),Y2(NPTG),
     1  AUX1(NPTG),AUX2(NPTG),AUX3(NPTG)
      DIMENSION R(NV),RV(NV),RW(NV)
C
      IF(NV.GT.NDIM) THEN
        WRITE(6,2101) NV,NDIM
 2101   FORMAT(1X,'*** Error in ZVINT: input potential grid with NV = ',
     1    I5,' data points.',/5X,'NV must be less than NDIM = ',I5,'.')
        STOP
      ENDIF
      IF(NV.LT.4) THEN
        WRITE(6,2102) NV
 2102   FORMAT(1X,'*** Error in ZVINT: the input potential grid must ',
     1    /5X,'have more than 4 data points. NV =',I5,'.')
        STOP
      ENDIF
      IF(R(1).LT.0.0D0) THEN
        WRITE(6,2103)
 2103   FORMAT(1X,'*** Error in ZVINT: R(1).LT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
      IF(R(1).GT.1.0D-15) THEN
        WRITE(6,2104)
 2104   FORMAT(1X,'*** Error in ZVINT: R(1).GT.0.')
        WRITE(6,'(5X,''R(1) = '',1PE14.7)') R(1)
        STOP
      ENDIF
C
      RT(1)=0.0D0
      VT(1)=RV(1)
      WT(1)=RW(1)
      DO I=2,NV
        RT(I)=R(I)
        VT(I)=RV(I)
        WT(I)=RW(I)
        IF(RT(I-1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I-1)))) THEN
          WRITE(6,2105)
 2105     FORMAT(1X,'*** Error in ZVINT: X values in',
     1      'decreasing order.',
     2      /5X,'Details in file ''ZVINT-error.dat''.')
          OPEN(33,FILE='ZVINT-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,NV
              WRITE(33,'(I5,1P,3E18.10)') J,R(J),RV(J),RW(J)
            ENDDO
          CLOSE(33)
          STOP 'ZVINT: X values in decreasing order.'
        ENDIF
      ENDDO
C
C  ****  Coulomb tail of the real part of the potential.
C
      ZINF=VT(NV)
      TOL=MAX(ABS(ZINF)*1.0D-10,1.0D-10)
      NVI=NV
      DO I=NV,4,-1
        IF(ABS(VT(I-1)-ZINF).GT.TOL) THEN
          NVE=I
          IF(RT(I)-RT(I-1).LT.EPS*MAX(ABS(RT(I-1)),ABS(RT(I)))) THEN
            RT(I)=RT(I-1)
            WT(I)=WT(I-1)
            VT(I)=ZINF
            NVE=I-1
          ELSE
            NVI=NVI+1  ! Add a discontinuity.
            DO J=NVI,NVE+1,-1
              RT(J)=RT(J-1)
              WT(J)=WT(J-1)
            ENDDO
            VT(NVE+1)=ZINF
          ENDIF
          GO TO 10
        ELSE
          VT(I)=ZINF
        ENDIF
      ENDDO
      NVE=4
 10   CONTINUE
C
C  ****  Tail of the imaginary part of the potential.
C
      TOL=1.0D-8
      DO I=NVI,4,-1
        IF(ABS(WT(I-1)).GT.TOL) THEN
          WT(I)=0.0D0
          NWE=I-1
          IF(ABS(NWE-NVE).LT.3) NWE=NVE
          GO TO 20
        ELSE
          WT(I)=0.0D0
        ENDIF
      ENDDO
      NWE=4
 20   CONTINUE
      NVE=MAX(NVE,NWE)
C
C  ****  Natural cubic spline interpolation, piecewise.
C
      IO=0
      I=0
      K=0
    1 I=I+1
      K=K+1
      X(K)=RT(I)
      Y1(K)=VT(I)
      Y2(K)=WT(I)
      IF(I.EQ.NVE) GO TO 2
C  ****  Duplicated points are considered as discontinuities.
      IF(RT(I+1)-RT(I).GT.EPS*MAX(ABS(RT(I)),ABS(RT(I+1)))) GO TO 1
    2 CONTINUE
C
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y1,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y1,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      IOO=IO
      DO J=1,K-1
        IO=IO+1
        RG(IO)=X(J)
        RVG(IO)=Y1(J)
        VA(IO)=A(J)
        VB(IO)=B(J)
        VC(IO)=C(J)
        VD(IO)=D(J)
      ENDDO
      IF(K.GT.3) THEN
        CALL SPLIN0(X,Y2,A,B,C,D,0.0D0,0.0D0,K)
      ELSE
        CALL SPLINE(X,Y2,A,B,C,D,0.0D0,0.0D0,K)
      ENDIF
      IO=IOO
      DO J=1,K-1
        IO=IO+1
        RWG(IO)=Y2(J)
        WA(IO)=A(J)
        WB(IO)=B(J)
        WC(IO)=C(J)
        WD(IO)=D(J)
      ENDDO
      IF(I.LT.NVE) THEN
        K=0
        GO TO 1
      ENDIF
C  ****  The last sets of coefficients of the splines are replaced by
C        those of the potential tails.
      IO=IO+1
      NVT=IO
      RG(IO)=X(K)
      RVG(IO)=ZINF
      RWG(IO)=0.0D0
      VA(IO)=ZINF
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
      WA(IO)=0.0D0
      WB(IO)=0.0D0
      WC(IO)=0.0D0
      WD(IO)=0.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZRVSPL
C  *********************************************************************
      SUBROUTINE ZRVSPL(R,RVS,RWS)
C
C     This function gives the (natural cubic spline) interpolated values
C  of R*V(R) and R*W(R) at R, i.e., the potential functions that are
C  effectively used in the numerical solution.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
C
      IF(R.LT.0.0D0) THEN
        RVS=0.0D0
        RWS=0.0D0
      ELSE
        RVS=SPLVAL(R,RG,VA,VB,VC,VD,NVT)
        RWS=SPLVAL(R,RG,WA,WB,WC,WD,NVT)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSFREE
C  *********************************************************************
      SUBROUTINE ZSFREE(E,EPS,PHASER,PHASEI,L,IRWF)
C
C     This subroutine solves the Schrodinger radial equation for free
C  states of a complex optical potential.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/ZSTORE/ZPA(NPTG),ZQA(NPTG),ZPB(NPTG),ZQB(NPTG)
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
      ZI=DCMPLX(0.0D0,1.0D0)
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in ZSFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in ZSFREE: L.LT.0.')
        STOP
      ENDIF
      FL1=0.5D0*L*(L+1)
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in ZSFREE: E.LT.0.0001')
        RETURN
      ENDIF
      RK=SQRT(E+E)
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in ZSFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      RZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
 2    CONTINUE
      ILAST=NRT+1
      IF(ABS(RZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          T=EPS*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(ZRVN).GT.T) GO TO 3
          BNL1=SBESJN(2,L+1,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL=SBESJN(2,L,X)
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+1,X)
          ILAST=IL
          ZPA(ILAST)=X*BJL
          ZPB(ILAST)=-X*BNL
          ZQA(ILAST)=RK*((L+1.0D0)*BJL-X*BJL1)
          ZQB(ILAST)=-RK*((L+1.0D0)*BNL-X*BNL1)
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(RZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          IF(ABS(ZRVN-RZINF).GT.TAS) GO TO 3
          CALL SCOULF(RZINF,E,L,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          ZPA(ILAST)=P0
          ZPB(ILAST)=P1
          ZQA(ILAST)=Q0
          ZQB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (ZSFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in ZSFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL ZSOUTW(E,EPS,L,ILAST)
C
C  ****  Phase shift. (ME-6.62)
C
      ZPO=ZP(ILAST)
      ZPOP=ZQ(ILAST)
      ZPIA=ZPA(ILAST)
      ZPIAP=ZQA(ILAST)
      ZPIB=ZPB(ILAST)
      ZPIBP=ZQB(ILAST)
C
      ZPHASE=(ZPO*(ZPIAP+ZI*ZPIBP)-ZPOP*(ZPIA+ZI*ZPIB))
     1      /(ZPOP*(ZPIA-ZI*ZPIB)-ZPO*(ZPIAP-ZI*ZPIBP))
C  ****  The real phase shift is reduced to the interval (-PI/2,PI/2).
      ZPH=-ZI*CDLOG(ZPHASE)*0.5D0
      PHASER=ZPH
      PHASEI=-ZI*ZPH
      TT=ABS(PHASER)
      IF(TT.GT.PIH) PHASER=PHASER*(1.0D0-PI/TT)
C
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.63 and ME-6.64, ME-6.65)
C
      ZPHASE=CDEXP(ZI*PHASER-PHASEI)
      ZCD=(ZPHASE+1.0D0/ZPHASE)/2.0D0
      ZSD=-ZI*(ZPHASE-1.0D0/ZPHASE)/2.0D0
      IF(ABS(ZPO).GT.EPS) THEN
        ZNORM=(ZCD*ZPIA+ZSD*ZPIB)/ZPO
      ELSE
        ZNORM=(ZCD*ZPIAP+ZSD*ZPIBP)/ZPOP
      ENDIF
C
      DO I=1,ILAST
        ZP(I)=ZNORM*ZP(I)
        ZQ(I)=ZNORM*ZQ(I)
        IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
        IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          ZP(I)=ZCD*ZPA(I)+ZSD*ZPB(I)
          ZQ(I)=ZCD*ZQA(I)+ZSD*ZQB(I)
          IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
          IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=ZP(J)
          Q(I)=ZQ(J)
          PIM(I)=-ZI*ZP(J)
          QIM(I)=-ZI*ZQ(J)
        ELSE
          P(I)=ZP(J+1)
          Q(I)=ZQ(J+1)
          PIM(I)=-ZI*ZP(J+1)
          QIM(I)=-ZI*ZQ(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSOUTW
C  *********************************************************************
      SUBROUTINE ZSOUTW(E,EPS,L,IOTP)
C
C     Outward solution of the Schrodinger radial equation for a complex
C  piecewise cubic potential. Power series method, free states only.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      AL=L
      N1=IOTP-1
C
      ZP(1)=0.0D0
      ZQ(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        RW0=WA(IN)
        RW1=WB(IN)
        RW2=WC(IN)
        RW3=WD(IN)
        ZPI=ZP(I)
        ZQI=ZQ(I)
        CALL ZSCH(E,AL,EPS)
        ZP(I+1)=ZPF
        ZQ(I+1)=ZQF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO K=1,I
          ZP(K)=ZP(K)*FACT
          ZQ(K)=ZQ(K)*FACT
          ENDDO
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSCH
C  *********************************************************************
      SUBROUTINE ZSCH(E,AL,EPS)
C
C     This subroutine solves the Schrodinger radial equation for a
C  central COMPLEX potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C                     +ZI*(RW0+RW1*R+RW2*R**2+RW3*R**3)
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AL .................... orbital angular momentum quantum number.
C  Output argument:
C     EPS ................... estimate of the global error in ZPF
C                             and ZQF.
C
C  Input (common ZPOTEN):
C     RV0, RV1, RV2, RV3 .... real potential parameters,
C     RW0, RW1, RW2, RW3 .... imaginary potential parameters.
C
C  Input-output (common ZSINOU):
C     RA, RB ................ interval end points (input),
C     ZPI, ZQI .............. values of the radial function and its
C                             derivative at RA (input),
C     ZPF, ZQF .............. values of the radial function and its
C                             derivative at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     NSTEP ................. number of steps.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      COMMON/ZSSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),R0,R1,NSUM
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
        DIRECT=-1.0D0
      ELSE
        DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      ZP1=ZPI
      ZQ1=ZQI
 1    CONTINUE
      R0=R1
      ZP0=ZP1
      ZQ0=ZQ1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL ZSCH0(E,AL,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C
      NSTEP=NSTEP+1
      TST=ABS(ZP1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        ZP1=ZP1/TST
        ZQ1=ZQ1/TST
      ENDIF
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      ZPF=ZP1
      ZQF=ZQ1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZSCH0
C  *********************************************************************
      SUBROUTINE ZSCH0(E,AL,EPS)
C
C  Power series solution of the Schrodinger eq. for a central potential
C  with an absorptive imaginary component.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZSSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),R0,R1,NSUM
C
      ZI=DCMPLX(0.0D0,1.0D0)
C
      RVE=RV1-E
      IF(R0.GT.1.0D-10) GO TO 2
C
C  ****  First interval. (ME-4.15 to ME-4.17)
C
      S=AL+1
      ZU0=AL*S
      ZU1=2*(RV0+ZI*RW0)*R1
      ZU2=2*(RVE+ZI*RW1)*R1**2
      ZU3=2*(RV2+ZI*RW2)*R1**3
      ZU4=2*(RV3+ZI*RW3)*R1**4
      ZUT=ZU0+ZU1+ZU2+ZU3+ZU4
C
      ZCA(1)=1.0D0
      ZCA(2)=ZU1*ZCA(1)/((S+1)*S-ZU0)
      ZCA(3)=(ZU1*ZCA(2)+ZU2*ZCA(1))/((S+2)*(S+1)-ZU0)
      ZCA(4)=(ZU1*ZCA(3)+ZU2*ZCA(2)+ZU3*ZCA(1))
     1  /((S+3)*(S+2)-ZU0)
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZQ1=S*ZCA(1)+(S+1)*ZCA(2)+(S+2)*ZCA(3)+(S+3)*ZCA(4)
      ZP2P1=S*(S-1)*ZCA(1)+(S+1)*S*ZCA(2)+(S+2)*(S+1)*ZCA(3)
     1  +(S+3)*(S+2)*ZCA(4)
C
      DO I=5,60
        K=I-1
        ZCA(I)=(ZU1*ZCA(K)+ZU2*ZCA(I-2)+ZU3*ZCA(I-3)+ZU4*ZCA(I-4))
     1    /((S+K)*(S+K-1)-ZU0)
        ZP1=ZP1+ZCA(I)
        ZDQ1=(S+K)*ZCA(I)
        ZQ1=ZQ1+ZDQ1
        ZP2P1=ZP2P1+(S+K-1)*ZDQ1
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZP2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(ZCA(I))
        T2=ABS(R1*R1*(ZP2P1-ZUT*ZP1))
        TST1=EPS*MAX(ABS(ZP1),ABS(ZQ1)/I)
        TST2=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 1
      ENDDO
C  ****  Renormalization. (ME-4.18)
 1    CONTINUE
      NSUM=K+1
      ZQ1=ZQ1/(ABS(ZP1)*R1)
      ZP1=ZP1/ABS(ZP1)
      RETURN
C
C  ****  Middle region. (ME-4.10 to ME-4.13)
C
 2    CONTINUE
      S=0.0D0
      H=R1-R0
      H2=H*H
C
      ZV0=RV0+ZI*RW0
      ZV1=RVE+ZI*RW1
      ZV2=RV2+ZI*RW2
      ZV3=RV3+ZI*RW3
C
      RHO=H/R0
      ZU0=AL*(AL+1)+2*R0*(ZV0+R0*(ZV1+R0*(ZV2+R0*ZV3)))
      ZU1=2*(ZV0+R0*(2*ZV1+R0*(3*ZV2+R0*4*ZV3)))*H
      ZU2=2*(ZV1+R0*(3*ZV2+R0*6*ZV3))*H2
      ZU3=2*(ZV2+R0*4*ZV3)*H2*H
      ZU4=2*ZV3*H2*H2
      ZUT=ZU0+ZU1+ZU2+ZU3+ZU4
C
      ZCA(1)=ZP0
      ZCA(2)=ZQ0*H
      ZCA(3)=RHO*RHO*ZU0*ZCA(1)/2
      ZCA(4)=RHO*(RHO*(ZU0*ZCA(2)+ZU1*ZCA(1))-4*ZCA(3))/6
      ZCAK=(ZU0-2)*ZCA(3)+ZU1*ZCA(2)+ZU2*ZCA(1)
      ZCA(5)=RHO*(RHO*ZCAK-12*ZCA(4))/12
      ZCAK=(ZU0-6)*ZCA(4)+ZU1*ZCA(3)+ZU2*ZCA(2)+ZU3*ZCA(1)
      ZCA(6)=RHO*(RHO*ZCAK-24*ZCA(5))/20
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)+ZCA(5)+ZCA(6)
      ZQ1=ZCA(2)+2*ZCA(3)+3*ZCA(4)+4*ZCA(5)+5*ZCA(6)
      ZP2P1=2*ZCA(3)+6*ZCA(4)+12*ZCA(5)+20*ZCA(6)
C
      DO I=7,60
        K=I-1
        ZCAK=(ZU0-(K-2)*(K-3))*ZCA(I-2)+ZU1*ZCA(I-3)+ZU2*ZCA(I-4)
     1    +ZU3*ZCA(I-5)+ZU4*ZCA(I-6)
        ZCA(I)=RHO*(RHO*ZCAK-2*(K-1)*(K-2)*ZCA(K))/(K*(K-1))
        ZP1=ZP1+ZCA(I)
        ZDQ1=K*ZCA(I)
        ZQ1=ZQ1+ZDQ1
        ZP2P1=ZP2P1+K*(K-1)*ZCA(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZP2P1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1=ABS(ZCA(I))
        T2=ABS(R1*R1*ZP2P1-H2*ZUT*ZP1)
        TST1=EPS*MAX(ABS(ZP1),ABS(ZQ1)/I)
        TST2=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST1.AND.T2.LT.TST2) GO TO 3
      ENDDO
C
 3    NSUM=K+1
      ZQ1=ZQ1/H
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDFREE
C  *********************************************************************
      SUBROUTINE ZDFREE(E,EPS,PHASER,PHASEI,K,IRWF)
C
C     This subroutine solves the Dirac radial equation for free states
C  of a complex optical potential.
C     When IRWF=0, the radial function is not returned.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI,PIH=0.5D0*PI)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
      COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/ZSTORE/ZPA(NPTG),ZQA(NPTG),ZPB(NPTG),ZQB(NPTG)
      COMMON/OCOUL/RK,ETA,DELTA
      ETA=0.0D0
      DELTA=0.0D0
      IER=0
      ZI=DCMPLX(0.0D0,1.0D0)
C
      IF(EPS.LT.1.0D-15) THEN
        WRITE(6,2100) EPS
 2100   FORMAT(1X,'*** Error in ZDFREE: EPS =',1P,E13.6,
     1    ' is too small.')
        STOP
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in ZDFREE: K.EQ.0.')
        STOP
      ENDIF
C
      IF(E.LT.0.0001D0) THEN
        IER=7
        WRITE(6,1007)
 1007   FORMAT(1X,'*** Error 7 in ZDFREE: E.LT.0.0001')
        RETURN
      ENDIF
C  ****  Orbital angular momentum quantum number. (ME-2.19d)
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
C  ****  Merge the 'RG' and 'RAD' grids.
C
      IF(NGP.GT.NDIM) THEN
        WRITE(6,2102) NDIM
 2102   FORMAT(1X,'*** Error in ZDFREE: User radial grid with',
     1    ' more than ',I5,' data points.')
        STOP
      ENDIF
      T=MAX(0.5D0*EPS,1.0D-10)
      DO I=1,NVT
        R(I)=RG(I)
        IND(I)=I
      ENDDO
      NRT=NVT
      RZINF=RV(NVT)
C
      DO 1 I=1,NGP
        RLOC=RAD(I)
        DO J=1,NRT
          IF(ABS(RLOC-R(J)).LT.T) GO TO 1
        ENDDO
        NRT=NRT+1
        CALL FINDI(RLOC,RG,NVT,J)
        R(NRT)=RLOC
        IND(NRT)=J
 1    CONTINUE
C  ****  ... and sort the resulting R-grid in increasing order.
      DO I=1,NRT-1
        RMIN=1.0D35
        IMIN=I
        DO J=I,NRT
          IF(R(J).LT.RMIN) THEN
            RMIN=R(J)
            IMIN=J
          ENDIF
        ENDDO
        IF(IMIN.NE.I) THEN
          RMIN=R(I)
          R(I)=R(IMIN)
          R(IMIN)=RMIN
          INDMIN=IND(I)
          IND(I)=IND(IMIN)
          IND(IMIN)=INDMIN
        ENDIF
      ENDDO
C
C  ****  Asymptotic solution.
C
      IWARN=0
  2   CONTINUE
      ILAST=NRT+1
      IF(ABS(RZINF).LT.EPS) THEN
        ETA=0.0D0
        DELTA=0.0D0
C  ****  Finite range potentials.
        FACTOR=SQRT(E/(E+2.0D0*SL*SL))
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1      +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          T=EPS*RN*ABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(ZRVN).GT.T) GO TO 3
          BNL=SBESJN(2,L,X)
          IF(ABS(BNL).GT.1.0D6) GO TO 3  ! Test cutoff.
          BNL1=SBESJN(2,L+KSIGN,X)
          IF(ABS(BNL1).GT.1.0D6) GO TO 3  ! Test cutoff.
          BJL=SBESJN(1,L,X)
          BJL1=SBESJN(1,L+KSIGN,X)
          ILAST=IL
          ZPA(ILAST)=X*BJL
          ZPB(ILAST)=-X*BNL
          ZQA(ILAST)=-FACTOR*KSIGN*X*BJL1
          ZQB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
      ELSE
C  ****  Coulomb potentials.
        TAS=MAX(1.0D-11,EPS)*ABS(RZINF)
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
     1        +(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          IF(ABS(ZRVN-RZINF).GT.TAS) GO TO 3
          CALL DCOULF(RZINF,E,K,RN,P0,Q0,P1,Q1,ERRF,ERRG)
          ERR=MAX(ERRF,ERRG)
          IF(ERR.GT.EPS.OR.ABS(P1).GT.1.0D6) GO TO 3  ! Test cutoff.
          ILAST=IL
          ZPA(ILAST)=P0
          ZPB(ILAST)=P1
          ZQA(ILAST)=Q0
          ZQB(ILAST)=Q1
        ENDDO
      ENDIF
 3    CONTINUE
      IF(ILAST.EQ.NRT+1) THEN
C  ****  Move R(NRT) outwards, seeking a possible matching point.
        R(NRT)=1.2D0*R(NRT)
        CALL FINDI(R(NRT),RG,NVT,J)
        IND(NRT)=J
        IF(IWARN.EQ.0) THEN
          WRITE(6,1008)
 1008     FORMAT(1X,'*** Warning (ZDFREE): RAD(NGP) is too small.'
     1      /5X,'Tentatively, it is moved outwards to')
          IWARN=1
        ENDIF
        WRITE(6,'(7X,''R(NRT) ='',1P,E13.6)') R(NRT)
        IF(R(NRT).LT.1.0D4) GO TO 2
      ENDIF
C
      IF(IWARN.EQ.1.AND.IRWF.NE.0) THEN
        IER=8
        WRITE(6,1009) R(NRT)
 1009   FORMAT(1X,'*** Error 8 in ZDFREE: RAD(NGP) is too small.'
     1    /5X,'Extend the grid to radii larger than',1P,E13.6)
        RETURN
      ENDIF
C
C  ****  Outward solution.
C
      CALL ZDOUTW(E,EPS,K,ILAST)
C
C  ****  Phase shift.  (ME-6.62)
C
      RM=R(ILAST)
      IL=IND(ILAST-1)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
     1  +(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPO=ZP(ILAST)
      ZPOP=-K*ZPO/RM+ZFG*ZQ(ILAST)
      IL=IND(ILAST)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
     1  +(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPIA=ZPA(ILAST)
      ZPIAP=-K*ZPIA/RM+ZFG*ZQA(ILAST)
      ZPIB=ZPB(ILAST)
      ZPIBP=-K*ZPIB/RM+ZFG*ZQB(ILAST)
C
      ZPHASE=(ZPO*(ZPIAP+ZI*ZPIBP)-ZPOP*(ZPIA+ZI*ZPIB))
     1      /(ZPOP*(ZPIA-ZI*ZPIB)-ZPO*(ZPIAP-ZI*ZPIBP))
C  ****  The real phase shift is reduced to the interval (-PI/2,PI/2).
      ZPH=-ZI*CDLOG(ZPHASE)*0.5D0
      PHASER=ZPH
      PHASEI=-ZI*ZPH
      TT=ABS(PHASER)
      IF(TT.GT.PIH) PHASER=PHASER*(1.0D0-PI/TT)
      IF(IRWF.EQ.0) RETURN
C
C  ****  Normalized wave function. (ME-6.63 and ME-6.64, ME-6.65)
C
      ZPHASE=CDEXP(ZI*PHASER-PHASEI)
      ZCD=(ZPHASE+1.0D0/ZPHASE)/2.0D0
      ZSD=-ZI*(ZPHASE-1.0D0/ZPHASE)/2.0D0
      IF(ABS(ZPO).GT.EPS) THEN
        ZNORM=(ZCD*ZPIA+ZSD*ZPIB)/ZPO
      ELSE
        ZNORM=(ZCD*ZPIAP+ZSD*ZPIBP)/ZPOP
      ENDIF
C
      DO I=1,ILAST
        ZP(I)=ZNORM*ZP(I)
        ZQ(I)=ZNORM*ZQ(I)
        IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
        IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
      ENDDO
      IF(ILAST.LT.NRT) THEN
        DO I=ILAST+1,NRT
          ZP(I)=ZCD*ZPA(I)+ZSD*ZPB(I)
          ZQ(I)=ZCD*ZQA(I)+ZSD*ZQB(I)
          IF(ABS(ZP(I)).LT.1.0D-99) ZP(I)=0.0D0
          IF(ABS(ZQ(I)).LT.1.0D-99) ZQ(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  Extract the 'RAD' grid...
C
      DO I=1,NGP
        RLOC=RAD(I)
        CALL FINDI(RLOC,R,NRT,J)
        IF(J.EQ.NRT) J=NRT-1
        IF(RLOC-R(J).LT.R(J+1)-RLOC) THEN
          P(I)=ZP(J)
          Q(I)=ZQ(J)
          PIM(I)=-ZI*ZP(J)
          QIM(I)=-ZI*ZQ(J)
        ELSE
          P(I)=ZP(J+1)
          Q(I)=ZQ(J+1)
          PIM(I)=-ZI*ZP(J+1)
          QIM(I)=-ZI*ZQ(J+1)
        ENDIF
      ENDDO
C
      RLOC=R(ILAST)
      CALL FINDI(RLOC,RAD,NGP,ILAST)
      ILAST=MIN(ILAST+1,NGP)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDOUTW
C  *********************************************************************
      SUBROUTINE ZDOUTW(E,EPS,K,IOTP)
C
C     Outward solution of the Dirac radial equation for a complex
c  piecewise cubic potential. Power series method, free states only.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/VGRIDI/RW(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      AK=K
      N1=IOTP-1
C
      ZP(1)=0.0D0
      ZQ(1)=0.0D0
      DO 1 I=1,N1
        RA=R(I)
        RB=R(I+1)
        IN=IND(I)
        RV0=VA(IN)
        RV1=VB(IN)
        RV2=VC(IN)
        RV3=VD(IN)
        RW0=WA(IN)
        RW1=WB(IN)
        RW2=WC(IN)
        RW3=WD(IN)
        ZPI=ZP(I)
        ZQI=ZQ(I)
        CALL ZDIR(E,AK,EPS)
        ZP(I+1)=ZPF
        ZQ(I+1)=ZQF
        IF(I.EQ.1) GO TO 1
C  ****  Renormalization.
        IF(RLN.GT.0.0D0) THEN
          FACT=EXP(-RLN)
          DO J=1,I
            ZP(J)=ZP(J)*FACT
            ZQ(J)=ZQ(J)*FACT
          ENDDO
        ENDIF
 1    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDIR
C  *********************************************************************
      SUBROUTINE ZDIR(E,AK,EPS)
C
C     This subroutine solves the Dirac radial equation for a central
C  COMPLEX potential V(R) such that
C              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
C                     +ZI*(RW0+RW1*R+RW2*R**2+RW3*R**3)
C  Given the boundary conditions (i.e. the value of the large and small
C  radial functions) at RA, the solution in the interval between RA and
C  RB is generated by using a piecewise power series expansion for a
C  partition of the interval, suitably chosen to allow fast convergence
C  of the series.
C
C  Input arguments:
C     E ..................... particle kinetic energy,
C     AK .................... relativistic angular momentum quantum
C                             number.
C  Output argument:
C     EPS ................... estimate of the global error in ZPF
C                             and ZQF.
C  Input (common ZPOTEN):
C     RV0, RV1, RV2, RV3 .... real potential parameters,
C     RW0, RW1, RW2, RW3 .... imaginary potential parameters.
C
C  Input-output (common ZDINOU):
C     RA, RB ................ interval end points (input),
C     ZPI, ZQI .............. values of the large and small radial
C                             functions at RA (input),
C     ZPF, ZQF .............. values of the large and small radial
C                             functions at RB (output),
C     RLN ................... LOG of the re-normalizing factor,
C     NSTEP ................. number of steps.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
      RLN=0.0D0
C
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
C
      R1=RA
      ZP1=ZPI
      ZQ1=ZQI
 1    CONTINUE
      R0=R1
      ZP0=ZP1
      ZQ0=ZQ1
 2    CONTINUE
      IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
        R1=RB
        H=RB-R0
        IOUT=1
      ENDIF
      CALL ZDIR0(E,AK,EPS)
C
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
 3    CONTINUE
      IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
 4    CONTINUE
C
      NSTEP=NSTEP+1
      TST=ABS(ZP1)
      IF(TST.GT.1.0D2) THEN
C  ****  Renormalization.
        RLN=RLN+LOG(TST)
        ZP1=ZP1/TST
        ZQ1=ZQ1/TST
      ENDIF
      IF(IOUT.EQ.0) GO TO 1
C  ****  Output.
      ZPF=ZP1
      ZQF=ZQ1
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE ZDIR0
C  *********************************************************************
      SUBROUTINE ZDIR0(E,AK,EPS)
C
C  Power series solution of the Dirac eq. for a central potential
C  with an imaginary component (negative for absorptive interac-
C  tions).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER*4 (I-N),
     1  COMPLEX*16 (Z)
      PARAMETER (OVER=1.0D15)  ! Overflow level.
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
C
      ZI=DCMPLX(0.0D0,1.0D0)
C
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      ZRVE=RV1+ZI*RW1-E
      ZRV0=RV0+ZI*RW0
      ZRV2=RV2+ZI*RW2
      ZRV3=RV3+ZI*RW3
C
      IF(R0.GT.1.0D-10) GO TO 4
C
C  ****  First interval. (ME-4.28 to ME-4.39)
C
      ZU0=ZRV0/SL
      ZU1=ZRVE*R1/SL
      ZU2=ZRV2*R1**2/SL
      ZU3=ZRV3*R1**3/SL
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
      ZUH=ZU1-2*SL*R1
      IF(ABS(ZU0).LT.1.0D-10) GO TO 1
C
C  ****  U0.NE.0. (ME-4.30 to ME-4.34)
      ZS=(AK*AK-ZU0*ZU0)**0.5D0
      ZDS=ZS+ZS
      ZCA(1)=1.0D0
      ZCB(1)=-(ZS+AK)/ZU0
      ZCAI=ZU1*ZCA(1)
      ZCBI=ZUH*ZCB(1)
      ZCA(2)=(-ZU0*ZCAI-(ZS+1-AK)*ZCBI)/(ZDS+1)
      ZCB(2)=((ZS+1+AK)*ZCAI-ZU0*ZCBI)/(ZDS+1)
      ZCAI=ZU1*ZCA(2)+ZU2*ZCA(1)
      ZCBI=ZUH*ZCB(2)+ZU2*ZCB(1)
      ZCA(3)=(-ZU0*ZCAI-(ZS+2-AK)*ZCBI)/(2*(ZDS+2))
      ZCB(3)=((ZS+2+AK)*ZCAI-ZU0*ZCBI)/(2*(ZDS+2))
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)
      ZQP1=ZS*ZCB(1)+(ZS+1)*ZCB(2)+(ZS+2)*ZCB(3)
C
      DO I=4,60
        K=I-1
        ZCAI=ZU1*ZCA(K)+ZU2*ZCA(I-2)+ZU3*ZCA(I-3)
        ZCBI=ZUH*ZCB(K)+ZU2*ZCB(I-2)+ZU3*ZCB(I-3)
        ZCA(I)=(-ZU0*ZCAI-(ZS+K-AK)*ZCBI)/(K*(ZDS+K))
        ZCB(I)=((ZS+K+AK)*ZCAI-ZU0*ZCBI)/(K*(ZDS+K))
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(ZS+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(ZS+K)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=1. (ME-4.35, ME-4.36)
 1    CONTINUE
      IF(ISIG.LT.0) GO TO 2
      ZS=ABS(AK)
      ZDS1=ZS+ZS+1
      ZCA(1)=1.0D0
      ZCB(1)=ZU1*ZCA(1)/ZDS1
      ZCA(2)=0.0D0
      ZCB(2)=ZU2*ZCA(1)/(ZDS1+1)
      ZCA(3)=-ZUH*ZCB(1)/2
      ZCB(3)=(ZU1*ZCA(3)+ZU3*ZCA(1))/(ZDS1+2)
      ZCA(4)=-(ZUH*ZCB(2)+ZU2*ZCB(1))/3
      ZCB(4)=(ZU1*ZCA(4)+ZU2*ZCA(3))/(ZDS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)+(ZS+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(ZS+1)*ZCB(1)+(ZS+2)*ZCB(2)+(ZS+3)*ZCB(3)
C
      DO I=5,60
        K=I-1
        ZCA(I)=-(ZUH*ZCB(I-2)+ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
        ZCB(I)=(ZU1*ZCA(I)+ZU2*ZCA(K)+ZU3*ZCA(I-2))/(ZDS1+K)
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(ZS+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(ZS+I)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
      GO TO 3
C
C  ****  U0.EQ.0 and SIG=-1. (ME-4.37, ME-4.38)
 2    CONTINUE
      S=ABS(AK)+1
      DS1=S+ABS(AK)
      RZUH=ZUH
      IF(RZUH.GT.0.0D0) THEN
        ZCB(1)=-1.0D0
      ELSE
        ZCB(1)=1.0D0
      ENDIF
      ZCA(1)=-ZUH*ZCB(1)/DS1
      ZCB(2)=0.0D0
      ZCA(2)=-ZU2*ZCB(1)/(DS1+1)
      ZCB(3)=ZU1*ZCA(1)/2
      ZCA(3)=-(ZUH*ZCB(3)+ZU3*ZCB(1))/(DS1+2)
      ZCB(4)=(ZU1*ZCA(2)+ZU2*ZCA(1))/3
      ZCA(4)=-(ZUH*ZCB(4)+ZU2*ZCB(3))/(DS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=S*ZCA(1)+(S+1)*ZCA(2)+(S+2)*ZCA(3)+(S+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(S-1)*ZCB(1)+S*ZCB(2)+(S+1)*ZCB(3)
C
      DO I=5,60
        K=I-1
        ZCB(I)=(ZU1*ZCA(I-2)+ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
        ZCA(I)=-(ZUH*ZCB(I)+ZU2*ZCB(K)+ZU3*ZCB(I-2))/(DS1+K)
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+(S+K)*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+(S+K-1)*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 3
      ENDDO
C  ****  Renormalization. (ME-4.39)
 3    CONTINUE
      NSUM=K+1
      ZQ1=ZQ1/ABS(ZP1)
      ZP1=ZP1/ABS(ZP1)
      RETURN
C
C  ****  Middle region. (ME-4.23 to ME-4.27)
C
 4    CONTINUE
      RHO=H/R0
      ZU0=(ZRV0+R0*(ZRVE+R0*(ZRV2+R0*ZRV3)))/SL
      ZU1=(ZRVE+R0*(2*ZRV2+R0*3*ZRV3))*H/SL
      ZU2=(ZRV2+R0*3*ZRV3)*H2/SL
      ZU3=ZRV3*H*H2/SL
      ZUB=ZU0-2*SL*R0
      ZUH=ZU1-2*SL*H
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
C
      ZCA(1)=ZP0
      ZCB(1)=ZQ0
      ZCA(2)=-RHO*(AK*ZCA(1)+ZUB*ZCB(1))
      ZCB(2)=RHO*(AK*ZCB(1)+ZU0*ZCA(1))
      ZCA(3)=-RHO*((AK+1)*ZCA(2)+ZUB*ZCB(2)+ZUH*ZCB(1))/2
      ZCB(3)=RHO*((AK-1)*ZCB(2)+ZU0*ZCA(2)+ZU1*ZCA(1))/2
      ZCA(4)=-RHO*((AK+2)*ZCA(3)+ZUB*ZCB(3)+ZUH*ZCB(2)
     1      +ZU2*ZCB(1))/3
      ZCB(4)=RHO*((AK-2)*ZCB(3)+ZU0*ZCA(3)+ZU1*ZCA(2)
     1      +ZU2*ZCA(1))/3
C
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZCA(2)+2*ZCA(3)+3*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=ZCB(2)+2*ZCB(3)+3*ZCB(4)
C
      DO I=5,60
        K=I-1
        ZCA(I)=-RHO*((AK+K-1)*ZCA(K)+ZUB*ZCB(K)+ZUH*ZCB(I-2)
     1     +ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
        ZCB(I)=RHO*((AK-K+1)*ZCB(K)+ZU0*ZCA(K)+ZU1*ZCA(I-2)
     1     +ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
        ZP1=ZP1+ZCA(I)
        ZPP1=ZPP1+K*ZCA(I)
        ZQ1=ZQ1+ZCB(I)
        ZQP1=ZQP1+K*ZCB(I)
C  ****  Check overflow limit.
        TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
        IF(TST.GT.OVER) THEN
          NSUM=100
          RETURN
        ENDIF
        T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
        T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
        T1=MAX(T1A,T1B)
        T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
        TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
        IF(T1.LT.TST.AND.T2.LT.TST) GO TO 5
      ENDDO
C
 5    CONTINUE
      NSUM=K+1
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC      Quantum defect and high-energy Dirac phase shift     CCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE QNTDEF
C  *********************************************************************
      SUBROUTINE QNTDEF(K,QD0,A,B,EPS,ERRM)
C
C     This subroutine determines the quantum defect MU for states with
C  relativistic angular-momentum quantum number K (kappa) as a function
C  of the energy eigenvalue E. For free states, the quantum defect is
C  defined as the inner phase shift divided by PI. MU is calculated
C  explicitly for bound states with principal quantum number N from 20
C  up to about 35, and for free states with energies between 1.0D-4 and
C  1.0D-3. The quantum defect is expressed as MU(E)=QD0+A*E+B*E*E, with
C  the parameters QD0, A and B determined from a least-squares fit to
C  the calculated data; ERRM is the largest relative error (%) of the
C  fit. Note that quantum defects are defined only for partially-
C  screened Coulomb potentials with an attractive tail.
C  (ME-6.12 to ME-6.15)
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (SL2=SL**2,TSL2=2.0D0*SL2)
      COMMON/RADWF/RAD(NDIM),P(NDIM),Q(NDIM),NGP,ILAST,IER
C ****  Potential table.
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RP(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      DIMENSION DIFR(NDIM)
C  ****  Auxiliary arrays.
      PARAMETER (NP=3)
      DIMENSION X(100),Y(100),AF(NP,NP),BF(NP),PAR(NP),NN(100)
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I6)') K
        STOP 'QNTDEF: The K value is not allowed'
      ENDIF
C
      ZINF=RVG(NVT)
      IF(ZINF.GT.-0.5D0) THEN
        WRITE(6,'('' ZINF ='',1P,E16.8)') ZINF
        STOP 'QNTDEF: The Coulomb tail is too weak.'
      ENDIF
C
C ****  Bound states.
C
      NPR=MIN(5000,NDIM)
      RMAX=3000.0D0
      CALL SGRID(RAD,DIFR,RMAX,1.0D-6,1.0D0,NPR,NDIM,IERS)
      IF(IERS.NE.0) STOP 'QNTDEF: Radial grid error (1).'
      NGP=NPR
      GAMMA=SQRT(K**2-(ZINF/SL)**2)
C
      IF(K.LT.0) THEN
        N=ABS(K)
      ELSE
        N=K+1
      ENDIF
      N=MAX(N-1,19)
C
      E=0.0D0
      RMU=0.0D0
      NDAT=0
 100  CONTINUE
      N=N+1
      CALL DBOUND(E,EPS,N,K)
      IF(IER.NE.0) THEN
        WRITE(6,*) '    This error is harmless.'
        WRITE(6,*) ' '
        WRITE(6,*) '    The following calculation may be quite slow.'
        WRITE(6,*) '    Please, be patient...'
        GO TO 101
      ENDIF
      RNEF=ABS(ZINF/SL)*(E+SL2)/(SQRT(-E*(E+TSL2)))
      RMU=N+GAMMA-ABS(K)-RNEF
      NDAT=NDAT+1
      X(NDAT)=E
      Y(NDAT)=RMU
      NN(NDAT)=N
      WRITE(6,'(1P,2E16.8,I4)') E,RMU,N
      IF(N.LT.45) GO TO 100
 101  CONTINUE
      IF(ABS(RMU).LT.1.0D-9) THEN
        QD0=0.0D0
        A=0.0D0
        B=0.0D0
        RETURN
      ENDIF
C
C ****  Free states.
C
      DE=(X(NDAT)-X(1))/4
      E=MAX(1.0D-4,-X(NDAT))-DE
      NF=0
 200  CONTINUE
      NF=NF+1
      E=E+DE
      WAVEL=2.0D0*PI/SQRT(E*(2.0D0+E/SL**2))
      DRN=WAVEL/20.0D0
      NPR=MIN(5000,NDIM)
      RMAX=DRN*DBLE(NPR-300)
      CALL SGRID(RAD,DIFR,RMAX,1.0D-6,DRN,NPR,NDIM,IERS)
      IF(IERS.NE.0) STOP 'QNTDEF: Radial grid error (2).'
      NGP=NPR
      CALL DFREE(E,EPS,PHASE,K,0)
      IF(IER.NE.0) STOP 'QNTDEF: Fatal error in DFREE.'
      DEL=PHASE/PI
      IF(DEL.LT.RMU-0.5D0) THEN
        ISH=INT(RMU-DEL+0.5D0)
        DEL=DEL+ISH
      ELSE IF(DEL.GT.RMU+0.5D0) THEN
        ISH=INT(DEL-RMU+0.5D0)
        DEL=DEL-ISH
      ENDIF
      NDAT=NDAT+1
      X(NDAT)=E
      Y(NDAT)=DEL
      NN(NDAT)=0
      WRITE(6,'(1P,2E16.8)') E,DEL
      IF(NF.LT.5) GO TO 200
C
C  ****  Quantum-defect function. Least-squares fit.
C
      NPAR=3
 300  CONTINUE
      DO I=1,NP
        DO J=1,NP
          AF(I,J)=0.0D0
        ENDDO
        BF(I)=0.0D0
        PAR(I)=0.0D0
      ENDDO
      DO II=1,NDAT
        DO I=1,NPAR
          DO J=1,NPAR
            AF(I,J)=AF(I,J)+X(II)**(I+J-2)
          ENDDO
          BF(I)=BF(I)+Y(II)*X(II)**(I-1)
        ENDDO
      ENDDO
      CALL SLQS(AF,BF,PAR,DET,NPAR,NP,IER)
      IF(IER.NE.0) THEN
        NPAR=NPAR-1
        GO TO 300
      ENDIF
      QD0=PAR(1)
      A=PAR(2)
      B=PAR(3)
C
      OPEN(33,FILE='qntdef.dat')
      WRITE(33,'('' # Quantum defect function.  K = '',I4)') K
      WRITE(33,'('' #'',4X,''QD0 = '',1P,E16.8)') QD0
      WRITE(33,'('' #'',4X,''  A = '',1P,E16.8)') A
      WRITE(33,'('' #'',4X,''  B = '',1P,E16.8)') B
      WRITE(33,1001)
 1001 FORMAT(/' #',4X,'N',8X,'E',15X,'MU',13X,'FIT',12X,'ERR (%)',
     1   /' #',2X,69('-'))
      ERRM=0.0D0
      DO II=1,NDAT
        FIT=0.0D0
        DO I=1,NPAR
          FIT=FIT+PAR(I)*X(II)**(I-1)
        ENDDO
        ERR=100.0D0*(FIT-Y(II))/MAX(ABS(Y(II)),1.0D-35)
        ERRM=MAX(ERRM,ABS(ERR))
        WRITE(33,'(3X,I4,1P,5E16.8)') NN(II),X(II),Y(II),FIT,ERR
      ENDDO
      WRITE(33,'(/'' #   Largest error (%) ='',1P,E9.2)') ERRM
      CLOSE(UNIT=33)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SLQS
C  *********************************************************************
      SUBROUTINE SLQS(A,B,X,DET,N,NP,IER)
C
C     Solution of the system of linear equations A*X=B, with N equations
C  and N unknowns, by Gauss-Jordan elimination with partial pivoting. On
C  input, the matrix A(1:N,1:N) is stored in an array of physical dimen-
C  sions NP by NP; the vector B(1:N) is stored in an array of physical
C  dimension NP. The solution X(1:N) is returned in an array of physical
C  dimension NP. The output value of DET is the determinant of the
C  matrix A. IER is an error flag, IER=0 means that the calculation has
C  been completed successfully, a value IER=1 is returned when
C  DET=0.0D0. The input matrix A and the vector B are destroyed.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION A(NP,NP),B(NP),X(NP)
      IER=1
      DET=1.0D0
      NM1=N-1
      IF(NM1.LT.0) RETURN
      IF(NM1.EQ.0) GO TO 2
C  ****  Gauss ordering.
      DO I=1,NM1
        I1=I+1
        IM=I
C  ****  Maximum pivot.
        TST=ABS(A(I,I))
        DO J=I1,N
          PST=ABS(A(J,I))
          IF(PST.GT.TST) THEN
            TST=PST
            IM=J
          ENDIF
        ENDDO
        IF(TST.LT.1.D-15) RETURN
        IF(IM.EQ.I) GO TO 1
C  ****  Re-ordering of rows.
        DO K=I,N
          SAVE=A(IM,K)
          A(IM,K)=A(I,K)
          A(I,K)=SAVE
        ENDDO
        SAVE=B(IM)
        B(IM)=B(I)
        B(I)=SAVE
        DET=-DET
C  ****  Renormalization.
 1      AUX=1.0D0/A(I,I)
        DO J=I1,N
          SAVE=A(J,I)*AUX
          DO K=I,N
            A(J,K)=A(J,K)-A(I,K)*SAVE
          ENDDO
          B(J)=B(J)-B(I)*SAVE
        ENDDO
        DET=DET*A(I,I)
        DO J=I,N
          A(I,J)=A(I,J)*AUX
        ENDDO
        B(I)=B(I)*AUX
      ENDDO
 2    DET=DET*A(N,N)
C  ****  Solution.
      IF(ABS(A(N,N)).LT.1.0D-15) RETURN
      X(N)=B(N)/A(N,N)
      IER=0
      IF(NM1.EQ.0) RETURN
      DO I1=1,NM1
        I=N-I1
        SUM=0.0D0
        ID=I+1
        DO J=ID,N
          SUM=SUM+A(I,J)*X(J)
        ENDDO
        X(I)=B(I)-SUM
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DELINF
C  *********************************************************************
      SUBROUTINE DELINF(HEDEL)
C
C     This subroutine calculates the high-energy limit HEDEL of the
C  Dirac inner phase shifts. HEDEL = -SUMV/SL, where SUMV is the
C  integral over R, from zero to infinity, of the short-range part of
C  the potential V(R), that is, excluding the Coulomb tail. (ME-6.11)
C
C     The output value HEDEL = +- 1.0D35 indicates that the short-range
C  potential has a pole at R=0 (the integral diverges).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
C
      ZINF=RVG(NVT)
      DO I=1,NVT
        A(I)=VA(I)-ZINF
      ENDDO
      IF(A(1).LT.-1.0D-12) THEN
        HEDEL=1.0D35
      ELSE IF(A(1).GT.1.0D-12) THEN
        HEDEL=-1.0D35
      ELSE
C  ****  NB: The lower limit of the integral is set equal to 1.0D-34 to
C  pass a consistency check, which protects against taking the logarithm
C  of zero. This does not affect the result.
        SUMV=SPLINT(RG,A,VB,VC,VD,1.0D-34,RG(NVT),NVT,-1)
        HEDEL=-SUMV/SL
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC         Coulomb and Bessel functions         CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SCOULF
C  *********************************************************************
      SUBROUTINE SCOULF(Z,E,L,R,F,FP,G,GP,ERRF,ERRG)
C
C     This subroutine computes radial Schrodinger-Coulomb wave functions
C  for free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     L ........ orbital angular momentum quantum number (.GE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     F, FP .... regular Schrodinger-Coulomb function and its
C                derivative,
C     G, GP .... irregular Schrodinger-Coulomb function and its
C                derivative,
C     ERRF, ERRG ... accuracy of the computed functions (relative
C                    uncertainty).
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalized so that, for large R, they
C  oscillate with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PTOL=1.0D-10)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C
C  ****  Parameters.
C
      WAVNUM=SQRT(E+E)
      IF(ABS(Z).GT.0.00001D0) THEN
        ETA=Z/WAVNUM
        ICAL=0
      ELSE
        ETA=0.0D0
        DELTA=0.0D0
        ICAL=1
      ENDIF
      RLAMB=L
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.L.LT.0) THEN
        F=0.0D0
        FP=0.0D0
        G=1.0D35
        GP=-1.0D35
        DELTA=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in SCOULF: E is too small.')
        IF(L.LT.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in SCOULF: L.LT.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ************  Coulomb functions and phase shift.
C
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
      FP=FP*WAVNUM
      GP=GP*WAVNUM
      DELTA=DELTAC(ETA,RLAMB)
      ERRF=ERR
      ERRG=ERR
      IF(ERR.GE.PTOL.OR.ABS(G).GT.9.999999999D34) THEN
C  ****  Very small radii.
        RLAMB1=RLAMB+1.0D0
        CALL FCRS(ETA,RLAMB,X,F,ERR1)
        CALL FCRS(ETA,RLAMB1,X,FP1,ERR2)
        ERR=MAX(ERR1,ERR2)
        IF(ERR.LT.PTOL.AND.X.GT.1.0D-20) THEN
          SX=RLAMB1+X*ETA/RLAMB1
          RR=SQRT(RLAMB1**2+ETA**2)/RLAMB1
          FP=SX*(F/X)-RR*FP1
          ERRF=ERR
        ELSE
          F=0.0D0
          FP=0.0D0
          ERRF=1.0D0
        ENDIF
        G=1.0D35
        GP=-1.0D35
        ERRG=1.0D0
      ENDIF
      RETURN
C
C  ************  Z=0. Spherical Bessel functions.
C
 1    CONTINUE
      F=X*SBESJN(1,L,X)
      FP=((L+1)*SBESJN(1,L,X)-X*SBESJN(1,L+1,X))*WAVNUM
      G=-X*SBESJN(2,L,X)
      GP=-((L+1)*SBESJN(2,L,X)-X*SBESJN(2,L+1,X))*WAVNUM
      DELTA=0.0D0
      ERRF=0.0D0
      IF(ABS(G).GT.1.0D30) THEN
        ERRG=1.0D0
      ELSE
        ERRG=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCOULF
C  *********************************************************************
      SUBROUTINE DCOULF(Z,E,K,R,FU,FL,GU,GL,ERRF,ERRG)
C
C     This subroutine computes radial Dirac-Coulomb wave functions for
C  free states.
C
C  **** All quantities in atomic units.
C
C  Input arguments:
C     Z ........ potential strength, i.e. value of R*V(R) (assumed
C                constant),
C     E ........ particle kinetic energy (positive),
C     K ........ angular momentum quantum number kappa (.NE.0),
C     R ........ radial distance (positive).
C
C  Output arguments:
C     FU, FL ... upper and lower components of the regular Dirac-
C                Coulomb function,
C     GU, GL ... upper and lower components of the irregular Dirac-
C                Coulomb function,
C     ERRF, ERRG ... accuracy of the computed functions (relative
C                    uncertainty),
C
C  Output through common /OCOUL/:
C     WAVNUM ... wave number,
C     ETA ...... Sommerfeld's parameter,
C     DELTA .... Coulomb phase shift (modulus 2*PI).
C
C     Radial functions are normalized so that, for large r, the upper
C  component function oscillates with unit amplitude.
C
C     Other subprograms required: subroutines FCOUL and SUM2F0,
C                                 and function CLGAM.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI)
      PARAMETER (PTOL=1.0D-10)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
C
      IF(ABS(Z).GT.0.00001D0) THEN
        ZETA=Z*ALPHA
        ICAL=0
      ELSE
        ZETA=0.0D0
        ICAL=1
      ENDIF
      RLAMBS=K*K-ZETA*ZETA
      RLAMB=SQRT(RLAMBS)
      PC=SQRT(E*(E+TSL2))
      WAVNUM=PC/SL
      X=WAVNUM*R
C
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERRF=1.0D0
        ERRG=1.0D0
        DELTA=0.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** Error in DCOULF: E is too small.')
        IF(K.EQ.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DCOULF: K.EQ.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
C
C  ****  Parameters.
C
      RLAMB1=RLAMB-1.0D0
      W=E+SL2
      ETA=ZETA*W/PC
      RLA=SQRT(RLAMBS+ETA*ETA)
      P1=K+RLAMB
      P2=RLAMB*SL2-K*W
      RNUR=ZETA*(W+SL2)
      RNUI=-P1*PC
      RNU=ATAN2(RNUI,RNUR)
      RNORM=1.0D0/(SQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)
C
C  ****  Coulomb phase shift.
C
      IF(K.GT.0) THEN
        L=K
      ELSE
        L=-K-1
      ENDIF
      DELTA0=DELTAC(ETA,RLAMB1)
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN
        RNORM=-RNORM
        DELTA=DELTA-PI
      ENDIF
      IF(DELTA.GE.0.0D0) THEN
        DELTA=MOD(DELTA,TPI)
      ELSE
        DELTA=-MOD(-DELTA,TPI)
      ENDIF
C
C  ****  Coulomb functions.
C
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR0)
C
      Q2=P1*P2*RNORM
      Q1=RLA*PC*RNORM
      P1=P1*Q1
      Q1=ZETA*Q1
      P2=ZETA*P2*RNORM
C  ****  Very small radii.
      IF(ERR0.GE.PTOL) THEN
        CALL FCRS(ETA,RLAMB,X,F,ERR1)
        CALL FCRS(ETA,RLAMB1,X,FM1,ERR2)
        ERR=MAX(ERR1,ERR2)
        ERRF=ERR
        ERRG=ERR
        IF(ERR.LT.PTOL) THEN
          FU=P1*F+P2*FM1
          FL=-Q1*F-Q2*FM1
        ELSE
          FU=0.0D0
          FL=0.0D0
          ERRF=1.0D0
        ENDIF
        GU=1.0D35
        GL=-1.0D35
        ERRG=1.0D0
        RETURN
      ENDIF
      SLA=(RLAMB/X)+(ETA/RLAMB)
      F=RLAMB*(SLA*FM1-FPM1)/RLA
      G=RLAMB*(SLA*GM1-GPM1)/RLA
C
      FU=P1*F+P2*FM1
      GU=P1*G+P2*GM1
      FL=-Q1*F-Q2*FM1
      GL=-Q1*G-Q2*GM1
      ERRF=ERR0
      ERRG=ERR0
      RETURN
C
C  ****  Z=0. Spherical Bessel functions.
C
 1    CONTINUE
      RLAMB=ABS(K)
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)
      DELTA=0.0D0
      IF(ERR.GE.PTOL) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERRF=1.0D0
        ERRG=1.0D0
        RETURN
      ENDIF
      FM1=(RLAMB*F/X)+FP
      GM1=(RLAMB*G/X)+GP
      FACT=SQRT(E/(E+TSL2))
      IF(K.LT.0) THEN
        FU=FM1
        FL=-FACT*F
        GU=GM1
        GL=-FACT*G
      ELSE
        FU=F
        FL=FACT*FM1
        GU=G
        GL=FACT*GM1
      ENDIF
      ERRF=ERR
      IF(ABS(GL).GT.1.0D30) THEN
        ERRG=ERR
      ELSE
        ERRG=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SCOULB
C  *********************************************************************
      SUBROUTINE SCOULB(Z,N,L,E,R,P)
C
C     This subroutine computes the Schrodinger radial functions of bound
C  states in an attractive Coulomb field.
C
C  Input arguments:
C     Z ..... field strength (it must be negative),
C     N ..... principal quantum number,
C     L ..... orbital angular momentum quantum number
C             (Note: 0. LE. L .LE. N-1),
C     R ..... radial distance (positive).
C
C  Output arguments:
C     E ..... binding energy,
C     P ..... radial function at R.
C
C     Radial functions are normalized to unity.
C
C     Other subprograms required: function RLGAMA.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(Z.GE.0.0D0) THEN
        WRITE(6,2100) Z
 2100   FORMAT(1X,'*** Error in SCOULB: Z =',1P,E13.6,'  .GE.0.0')
        STOP
      ENDIF
      IF(N.LE.0) THEN
        WRITE(6,2101) DBLE(N)
 2101   FORMAT(1X,'*** Error in SCOULB: N =',1P,E13.6,'  .LE.0')
        STOP
      ENDIF
      IF(L.LT.0) THEN
        WRITE(6,2102) DBLE(L)
 2102   FORMAT(1X,'*** Error in SCOULB: L =',1P,E13.6,'  .LT.0')
        STOP
      ENDIF
      IF(L.GT.N-1) THEN
        WRITE(6,2103) DBLE(L),DBLE(N)
 2103   FORMAT(1X,'*** Error in SCOULB: L.GT.N-1',/5X,
     1    1P,'L =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      ZZ=ABS(Z)
C
      E=-ZZ*ZZ/(2.0D0*N*N)
      NR=N-L-1
      LP1=L+1
      A=ZZ/DBLE(N)
      B=L+L+2.0D0
      CN1=EXP(0.5D0*(RLGAMA(N+L+1.0D0)-RLGAMA(N-L*1.0D0))
     1   -RLGAMA(B))*SQRT(ZZ)/N
C
C  ****  Evaluation of the Kummer function and other variable
C        factors.
C
      FKUM=1.0D0
      TERM=1.0D0
      X=2.0D0*A*R
      IF(NR.GT.0) THEN
        DO I=0,NR-1
        TERM=TERM*(I-NR)*X/((I+B)*(I+1))
        FKUM=FKUM+TERM
        ENDDO
      ENDIF
      P=CN1*(X**LP1)*EXP(-A*R)*FKUM
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DCOULB
C  *********************************************************************
      SUBROUTINE DCOULB(Z,N,K,E,R,P,Q)
C
C     This subroutine computes the Dirac radial functions of bound
C  states in an attractive Coulomb field.
C
C  Input arguments:
C     Z ..... field strength (it must be negative),
C     N ..... principal quantum number,
C     K ..... relativistic angular momentum quantum number
C             (Note: -N .LE. K .LE. N-1, K .NE. 0),
C     R ..... radial distance (positive).
C
C  Output arguments:
C     E ..... binding energy,
C     P ..... large radial function at R,
C     Q ..... small radial function at R.
C
C     Radial functions are normalized to unity.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(Z.GE.0.0D0) THEN
        WRITE(6,2100) Z
 2100   FORMAT(1X,'*** Error in DCOULB: Z =',1P,E13.6,'  .GE.0.0')
        STOP
      ENDIF
      IF(N.LE.0) THEN
        WRITE(6,2101) DBLE(N)
 2101   FORMAT(1X,'*** Error in DCOULB: N =',1P,E13.6,'  .LE.0')
        STOP
      ENDIF
      IF(K.EQ.0) THEN
        WRITE(6,*) ' K =',K
        WRITE(6,2102)
 2102   FORMAT(1X,'*** Error in DCOULB: K = 0.')
        STOP
      ENDIF
      IF(K.LT.-N) THEN
        WRITE(6,2103) DBLE(K),DBLE(N)
 2103   FORMAT(1X,'*** Error in DCOULB: K.LT.-N',/5X,
     1    1P,'K =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      IF(K.GT.N-1) THEN
        WRITE(6,2104) DBLE(K),DBLE(N)
 2104   FORMAT(1X,'*** Error in DCOULB: K.GT.N-1',/5X,
     1    1P,'K =',E13.6,',  N =',E13.6)
        STOP
      ENDIF
      ZZ=ABS(Z)
      NR=N-ABS(K)
      ZETA=ZZ/SL
      RLAMB=SQRT(K*K-ZETA*ZETA)
      XE=(ZETA/(NR+RLAMB))**2
      IF(XE.GT.5.0D-4) THEN
        E=SL**2*((1.0D0/SQRT(1.0D0+XE))-1.0D0)
      ELSE
        E=SL**2*(-XE/2.0D0+3.0D0*XE**2/8.0D0-15.0D0*XE**3/48.0D0
     1   +105.0D0*XE**4/384.0D0)
      ENDIF
      TAU=ZETA/(RLAMB+NR)
      A=SL*SQRT(TAU*TAU/(1.0D0+TAU*TAU))
      X=2.0D0*A*R
C
      IF(NR.EQ.0) THEN
        FACT=SQRT(2.0D0*A/(EXP(RLGAMA(2.0D0*RLAMB+1.0D0))
     1      *(ZETA**2+(K+RLAMB)**2)))
        AUX1=X**RLAMB*EXP(-A*R)
        P=FACT*ZETA*AUX1
        Q=FACT*(K+RLAMB)*AUX1
        RETURN
      ENDIF
C
      NR2=NR-1
      AUX1=2.0D0*RLAMB*(2.0D0*RLAMB+1.0D0)
      AUX2=K+RLAMB*SQRT(1.0D0+TAU*TAU)
      P1=ZETA*AUX1
      P2=(K+RLAMB)*AUX2/TAU
      Q1=(K+RLAMB)*AUX1
      Q2=ZETA*AUX2/TAU
C
      B1=RLAMB+RLAMB
      B2=B1+2.0D0
      CN0=(K+RLAMB)*(RLAMB+NR)*AUX2*(1.0D0+1.0D0/(TAU*TAU))
      CN1=EXP(0.5D0*(RLGAMA(B1+NR+1)-RLGAMA(NR*1.0D0))
     1   -RLGAMA(B2))*SQRT(A/CN0)/B1
C
C  ****  Evaluation of the Kummer functions and other variable factors.
C
      FKUM1=1.0D0
      TERM=1.0D0
      IF(NR.GT.0) THEN
        DO I=0,NR-1
        TERM=TERM*(I-NR)*X/((I+B1)*(I+1))
        FKUM1=FKUM1+TERM
        ENDDO
      ENDIF
      FKUM2=1.0D0
      TERM=1.0D0
      IF(NR2.GT.0) THEN
        DO I=0,NR2-1
        TERM=TERM*(I-NR2)*X/((I+B2)*(I+1))
        FKUM2=FKUM2+TERM
        ENDDO
      ENDIF
      FACT=CN1*X**RLAMB*EXP(-A*R)
      P=FACT*(P1*FKUM1+P2*X*FKUM2)
      Q=FACT*(Q1*FKUM1+Q2*X*FKUM2)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RLGAMA
C  *********************************************************************
      FUNCTION RLGAMA(R)
C
C     This function gives LOG(GAMMA(R)) for real positive arguments.
C
C   Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C         Functions'. Dover, New York (1974). pp 255-257.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      ZA=R
      RLGAMA=80.5D0
      IF(ZA.LT.1.0D-16) RETURN
C
      ZFAC=1.0D0
      ZFL=0.0D0
 1    ZFAC=ZFAC/ZA
      IF(ZFAC.GT.1.0D8) THEN
        ZFL=ZFL+LOG(ZFAC)
        ZFAC=1.0D0
      ENDIF
      ZA=ZA+1.0D0
      IF(ZA.GT.15.0D0) GO TO 2
      GO TO 1
C  ****  Stirling's expansion of LOG(GAMMA(ZA)).
 2    ZI2=1.0D0/(ZA*ZA)
      ZS=(43867.0D0/244188.0D0)*ZI2
      ZS=(ZS-3617.0D0/122400.0D0)*ZI2
      ZS=(ZS+1.0D0/156.0D0)*ZI2
      ZS=(ZS-691.0D0/360360.0D0)*ZI2
      ZS=(ZS+1.0D0/1188.0D0)*ZI2
      ZS=(ZS-1.0D0/1680.0D0)*ZI2
      ZS=(ZS+1.0D0/1260.0D0)*ZI2
      ZS=(ZS-1.0D0/360.0D0)*ZI2
      ZS=(ZS+1.0D0/12.0D0)/ZA
      RLGAMA=(ZA-0.5D0)*LOG(ZA)-ZA+9.1893853320467274D-1+ZS
     1     +ZFL+LOG(ZFAC)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FCRS
C  *********************************************************************
      SUBROUTINE FCRS(ETA,RLAMB,X,F,ERR)
C
C     Regular Coulomb functions for real ETA, RLAMB.GT.-1 and small X.
C  Evaluated in terms of Kummer's hypergeometric series.
C
C  Input arguments:
C     ETA ...... Sommerfels's parameter,
C     RLAMB .... angular momentum,
C     X ........ variable (=wave number times radial distance).
C
C  Output arguments:
C     F ........ regular function,
C     ERR ...... relative numerical uncertainty. A value of the order of
C                10**(-N) means that the calculated functions are
C                accurate to N decimal figures. The maximum accuracy
C                attainable with double precision arithmetic is about
C                1.0D-15.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0)
C
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,2100) RLAMB
 2100   FORMAT(1X,'*** Error in FCRS: RLAMB =',1P,E13.6,'  .LT.-0.999')
        STOP
      ENDIF
C
C  ****  Numerical constants.
C
      CI=DCMPLX(0.0D0,1.0D0)
      B=2.0D0*(RLAMB+1.0D0)
      CA=RLAMB+1.0D0+CI*ETA
      CX=-2.0D0*CI*X
C
C  ****  Normalization constant.
C
      RCL=2**RLAMB*EXP(-0.5D0*ETA*PI-RLGAMA(B))
     1   *CDABS(CDEXP(CLGAM(CA)))
C
C  ****  Evaluation of Kummer's function. (ME-3.18)
C
      ERR=1.0D-15
      CFK=1.0D0
      CTERM=1.0D0
      ICONV=0
      DO I=0,200
        CTERM=CTERM*(I+CA)*CX/((I+B)*(I+1))
        CFK=CFK+CTERM
        TERM=CDABS(CTERM)
        IF(TERM.GT.1.0D35) THEN
          F=1.01D35
          ERR=1.0D0
          RETURN
        ENDIF
        IF(TERM.LT.CDABS(CFK)*1.0D-15) THEN
          ICONV=ICONV+1
          IF(ICONV.GT.4) GO TO 1
        ELSE
          ICONV=0
        ENDIF
      ENDDO
      ERR=TERM/MAX(CDABS(CFK),1.0D-16)
 1    CONTINUE
      CF=RCL*X**(RLAMB+1.0D0)*CDEXP(CI*X)*CFK
      F=CF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FCOUL
C  *********************************************************************
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
C
C     Calculation of (Schrodinger) Coulomb functions for real ETA,
C  RLAMB.GT.-1 and X larger than, or of the order of XTP0 (the turning
C  point for RLAMB=0). Steed's continued fraction method is combined
C  with several recursion relations and an asymptotic expansion. The
C  output value ERR=1.0D0 indicates that the evaluation algorithm is not
C  applicable (X is too small).
C
C  Input arguments:
C     ETA ...... Sommerfeld's parameter,
C     RLAMB .... angular momentum,
C     X ........ variable (=wave number times radial distance).
C
C  Output arguments:
C     F, FP .... regular function and its derivative,
C     G, GP .... irregular function and its derivative,
C     ERR ...... relative numerical uncertainty. A value of the
C                order of 10**(-N) means that the calculated
C                functions are accurate to N decimal figures.
C                The maximum accuracy attainable with double
C                precision arithmetic is about 1.0D-15.
C
C     Other subprograms required: subroutine SUM2F0 and
C                                 functions DELTAC and CLGAM.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)
      PARAMETER (PTOL=1.0D-10)
C
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,2100) RLAMB
 2100   FORMAT(1X,'*** Error in FCOUL: RLAMB =',1P,E13.6,'  .LT.-0.999')
        STOP
      ENDIF
      IF(X.LT.EPS) GO TO 5
C
C  ****  Numerical constants.
C
      CI=DCMPLX(0.0D0,1.0D0)
      CI2=2.0D0*CI
      CIETA=CI*ETA
      X2=X*X
      ETA2=ETA*ETA
C
C  ****  Turning point (XTP). (ME-3.30)
C
      IF(RLAMB.GE.0.0D0) THEN
        XTP=ETA+SQRT(ETA2+RLAMB*(RLAMB+1.0D0))
      ELSE
        XTP=EPS
      ENDIF
      ERRS=10.0D0
      IF(X.LT.XTP) GO TO 1
C
C  ************  Asymptotic expansion. (ME-3.50 to ME-3.54)
C
C  ****  Coulomb phase-shift.
      DELTA=DELTAC(ETA,RLAMB)
C
      CPA=CIETA-RLAMB
      CPB=CIETA+RLAMB+1.0D0
      CPZ=CI2*X
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)
      CQA=CPA+1.0D0
      CQB=CPB+1.0D0
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)
C  ****  Functions.
      THETA=X-ETA*LOG(2.0D0*X)-RLAMB*PIH+DELTA
      IF(THETA.GT.1.0D4) THETA=MOD(THETA,TPI)
      CEITH=CDEXP(CI*THETA)
      CGIF=C2F0*CEITH
      G=CGIF
      F=-CI*CGIF
C  ****  Derivatives.
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH
      GP=CGIFP
      FP=-CI*CGIFP
C  ****  Global uncertainty. The Wronskian may differ from 1 due
C        to truncation and round off errors.
      ERR=MAX(ERR1,ERR2,ABS(G*FP-F*GP-1.0D0))
      IF(ERR.LE.EPS) RETURN
      ERRS=ERR
C
C  ************  Steed's continued fraction method.
C
 1    CONTINUE
      CIETA2=CIETA+CIETA
      ETAX=ETA*X
C
C  ****  Continued fraction for F. (ME-3.40 to ME-3.49)
C
      INULL=0
      RLAMBN=RLAMB+1.0D0
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN
      B0=(RLAMBN/X)+(ETA/RLAMBN)
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
      FA3=B0
      FA2=B0*B1+A1
      FB3=1.0D0
      FB2=B1
C  ************  Error identified and corrected by Prof. A. Stauffer.
C  The next line originally was RF=FA3. When BN=0 with N=2, this caused
C  premature convergence to an incorrect value.
      RF=FA2/FB2
C
      ICONV=0
      DO N=2,NTERM
        RFO=RF
        DAF=ABS(RF)
        RLAMBN=RLAMB+N
        AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2
        BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
        FA1=FA2*BN+FA3*AN
        FB1=FB2*BN+FB3*AN
C
        TST=ABS(FB1)
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) THEN
            WRITE(6,2200)
 2200       FORMAT(1X,'*** Warning (FCOUL): multiple null factors (1).')
          ENDIF
          INULL=1
          FA3=FA2
          FA2=FA1
          FB3=FB2
          FB2=FB1
          RF=RFO
        ELSE
          FA3=FA2/TST
          FA2=FA1/TST
          FB3=FB2/TST
          FB2=FB1/TST
          RF=FA2/FB2
          IF(ABS(RF-RFO).LT.EPS*DAF) THEN
            ICONV=ICONV+1
            IF(ICONV.GT.4) GO TO 2
          ELSE
            ICONV=0
          ENDIF
        ENDIF
      ENDDO
 2    CONTINUE
      IF(DAF.GT.1.0D-25) THEN
        ERRF=ABS(RF-RFO)/DAF
      ELSE
        ERRF=EPS
      ENDIF
      IF(ERRF.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.PTOL) GO TO 5
        RETURN
      ENDIF
C
C  ****  Downward recursion for F and FP. Only if RLAMB.GT.1 and
C        X.LT.XTP. (ME-3.31a,b)
C
      RLAMB0=RLAMB
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN
        ISHIFT=0
        XTPC=XTP
        RFM=0.0D0
      ELSE
        FT=1.0D0
        FTP=RF
        IS0=RLAMB0+PTOL
        TST=X*(X-2.0D0*ETA)
        RL1T=0.0D0
        DO I=1,IS0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SXL=(RLAMB0/X)+ETARL0
          RLAMB0=RLAMB0-1.0D0
          FTO=FT
          FT=(SXL*FT+FTP)/RL
          FTP=SXL*FT-RL*FTO
          IF(FT.GT.1.0D10) THEN
            FTP=FTP/FT
            FT=1.0D0
          ENDIF
          RL1T=RLAMB0*(RLAMB0+1.0D0)
          IF(TST.GT.RL1T) THEN
            ISHIFT=I
            GO TO 3
          ENDIF
        ENDDO
        ISHIFT=IS0
 3      CONTINUE
        XTPC=ETA+SQRT(ETA2+RL1T)
        RFM=FTP/FT
      ENDIF
C
C  ****  Continued fraction for P+CI*Q with RLAMB0. (ME-3.55 to ME-3.58)
C
      INULL=0
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)
      CB0=X-ETA
      CBN=2.0D0*(X-ETA+CI)
      CFA3=CB0
      CFA2=CB0*CBN+CAN
      CFB3=1.0D0
      CFB2=CBN
C  ************  Error identified and corrected by Prof. A. Stauffer.
C  The next line originally was CPIQ=CFA3. When BN=0 with N=2, this
C  caused premature convergence to an incorrect value.
      CPIQ=CFA2/CFB2
C
      DO N=2,NTERM
        CPIQO=CPIQ
        DAPIQ=CDABS(CPIQ)
        CAN=CAN+CIETA2+(N+N-2)
        CBN=CBN+CI2
        CFA1=CFA2*CBN+CFA3*CAN
        CFB1=CFB2*CBN+CFB3*CAN
        TST=CDABS(CFB1)
C
        IF(TST.LT.1.0D-25) THEN
          IF(INULL.GT.0) THEN
            WRITE(6,2300)
 2300       FORMAT(1X,'*** Warning (FCOUL): multiple null factors (2).')
          ENDIF
          INULL=1
          CFA3=CFA2
          CFA2=CFA1
          CFB3=CFB2
          CFB2=CFB1
          CPIQ=CPIQO
        ELSE
          CFA3=CFA2/TST
          CFA2=CFA1/TST
          CFB3=CFB2/TST
          CFB2=CFB1/TST
          CPIQ=CFA2/CFB2
          IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 4
        ENDIF
      ENDDO
 4    CONTINUE
      IF(DAPIQ.GT.1.0D-25) THEN
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ
      ELSE
        ERRPIQ=EPS
      ENDIF
      IF(ERRPIQ.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.PTOL) GO TO 5
        RETURN
      ENDIF
      CPIQ=CI*CPIQ/X
C
      RP=CPIQ
      RQ=-CI*CPIQ
      IF(RQ.LE.1.0D-25) GO TO 5
      ERR=MAX(ERRF,ERRPIQ)
C
C  ****  Inverting Steed's transformation. (ME-3.36, ME-3.37)
C
      IF(ISHIFT.LT.1) THEN
        RFP=RF-RP
        F=SQRT(RQ/(RFP**2+RQ**2))
        IF(FB2.LT.0.0D0) F=-F
        FP=RF*F
        G=RFP*F/RQ
        GP=(RP*RFP-RQ**2)*F/RQ
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 5
      ELSE
        RFP=RFM-RP
        FM=SQRT(RQ/(RFP**2+RQ**2))
        G=RFP*FM/RQ
        GP=(RP*RFP-RQ**2)*FM/RQ
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 5
C  ****  Upward recursion for G and GP (if ISHIFT.GT.0). (ME-3.32a,b)
        DO I=1,ISHIFT
          RLAMB0=RLAMB0+1.0D0
          ETARL0=ETA/RLAMB0
          RL=SQRT(1.0D0+ETARL0**2)
          SXL=(RLAMB0/X)+ETARL0
          GO=G
          G=(SXL*GO-GP)/RL
          GP=RL*GO-SXL*G
          IF(G.GT.1.0D35) GO TO 5
        ENDDO
        W=RF*G-GP
        F=1.0D0/W
        FP=RF/W
      ENDIF
C  ****  The Wronskian may differ from 1 due to round off errors.
      ERR=MAX(ERR,ABS(FP*G-F*GP-1.0D0))
      IF(ERR.LT.PTOL) RETURN
C
 5    CONTINUE
      F=0.0D0
      FP=0.0D0
      G=1.0D35
      GP=-1.0D35
      ERR=1.0D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SUM2F0
C  *********************************************************************
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)
C
C     Summation of the 2F0(CA,CB;1/CZ) hypergeometric asymptotic series.
C  The positive and negative contributions to the real and imaginary
C  parts are added separately to obtain an estimate of rounding errors.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)
      RRP=1.0D0
      RRN=0.0D0
      RIP=0.0D0
      RIN=0.0D0
      CDF=1.0D0
      ERR2=0.0D0
      ERR3=1.0D0
      AR=0.0D0
      AF=0.0D0
C  ****  Asymptotic series. (ME-3.54)
      DO I=1,NTERM
        J=I-1
        CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
        ERR1=ERR2
        ERR2=ERR3
        ERR3=CDABS(CDF)
        IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 1
        AR=CDF
        IF(AR.GT.0.0D0) THEN
          RRP=RRP+AR
        ELSE
          RRN=RRN+AR
        ENDIF
        AI=DCMPLX(0.0D0,-1.0D0)*CDF
        IF(AI.GT.0.0D0) THEN
          RIP=RIP+AI
        ELSE
          RIN=RIN+AI
        ENDIF
        CF=DCMPLX(RRP+RRN,RIP+RIN)
        AF=CDABS(CF)
        IF(AF.GT.1.0D25) THEN
          CF=0.0D0
          ERR=1.0D0
          RETURN
        ENDIF
        IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN
           ERR=EPS
           RETURN
        ENDIF
      ENDDO
C  ****  Round off error.
 1    CONTINUE
      TR=ABS(RRP+RRN)
      IF(TR.GT.1.0D-25) THEN
        ERRR=(RRP-RRN)*ACCUR/TR
      ELSE
        ERRR=1.0D0
      ENDIF
      TI=ABS(RIP+RIN)
      IF(TI.GT.1.0D-25) THEN
        ERRI=(RIP-RIN)*ACCUR/TI
      ELSE
        ERRI=1.0D0
      ENDIF
C  ****  ... and truncation error.
      IF(AF.GT.1.0D-25) THEN
        ERR=MAX(ERRR,ERRI)+ERR2/AF
      ELSE
        ERR=MAX(ERRR,ERRI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DELTAC
C  *********************************************************************
      FUNCTION DELTAC(ETA,RLAMB)
C
C     Calculation of Coulomb phase shift (modulus 2*PI). (ME-3.21)
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      CI=DCMPLX(0.0D0,1.0D0)
C  ****  Coulomb phase-shift.
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)
      IF(DELTAC.GE.0.0D0) THEN
        DELTAC=MOD(DELTAC,TPI)
      ELSE
        DELTAC=-MOD(-DELTAC,TPI)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION CLGAM
C  *********************************************************************
      FUNCTION CLGAM(CZ)
C
C     This function gives LOG(GAMMA(CZ)) for complex arguments.
C
C   Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C         Functions'. Dover, New York (1974). PP 255-257.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1  INTEGER*4 (I-N)
      CZA=CZ
      ICONJ=0
      AR=CZA
      CLGAM=36.84136149D0
      IF(CDABS(CZA).LT.1.0D-16) RETURN
C
      AI=CZA*DCMPLX(0.0D0,-1.0D0)
      IF(AI.GT.0.0D0) THEN
        ICONJ=0
      ELSE
        ICONJ=1
        CZA=DCONJG(CZA)
      ENDIF
C
      CZFAC=1.0D0
      CZFL=0.0D0
 1    CONTINUE
      CZFAC=CZFAC/CZA
      IF(CDABS(CZFAC).GT.1.0D8) THEN
        CZFL=CZFL+CDLOG(CZFAC)
        CZFAC=1.0D0
      ENDIF
      CZA=CZA+1.0D0
      AR=CZA
      IF(CDABS(CZA).LT.1.0D-16) RETURN
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2
      GO TO 1
C  ****  Stirling's expansion of CDLOG(GAMMA(CZA)).
 2    CONTINUE
      CZI2=1.0D0/(CZA*CZA)
      CZS=(43867.0D0/244188.0D0)*CZI2
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2
      CZS=(CZS+1.0D0/156.0D0)*CZI2
      CZS=(CZS-691.0D0/360360.0D0)*CZI2
      CZS=(CZS+1.0D0/1188.0D0)*CZI2
      CZS=(CZS-1.0D0/1680.0D0)*CZI2
      CZS=(CZS+1.0D0/1260.0D0)*CZI2
      CZS=(CZS-1.0D0/360.0D0)*CZI2
      CZS=(CZS+1.0D0/12.0D0)/CZA
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS
     1     +CZFL+CDLOG(CZFAC)
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SBESJN
C  *********************************************************************
      FUNCTION SBESJN(JN,N,X)
C
C     This function computes the spherical Bessel functions of the first
C  kind and spherical Bessel functions of the second kind (also known as
C  spherical Neumann functions) for real positive arguments.
C
C  Input arguments:
C        JN ...... kind: 1(Bessel) or 2(Neumann).
C        N ....... order (integer).
C        X ....... argument (real and positive).
C
C  Ref.: M. Abramowitz and I.A. Stegun, 'Handbook of Mathematical
C        Functions'. Dover, New York (1974), pp 435-478.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      IF(X.LT.0) THEN
        WRITE(6,1000)
 1000   FORMAT(1X,'*** Negative argument in function SBESJN.')
        STOP
      ENDIF
C  ****  Order and phase correction for Neumann functions.
C        Abramowitz and Stegun, Eq. 10.1.15.
      IF(JN.EQ.2) THEN
        NL=-N-1
        IPH=2*MOD(ABS(N),2)-1
      ELSE
        NL=N
        IPH=1
      ENDIF
C  ****  Selection of calculation mode.
      IF(NL.LT.0) GO TO 5
      IF(X.GT.1.0D0*NL) GO TO 3
      XI=X*X
      IF(XI.GT.NL+NL+3.0D0) GO TO 2
C  ****  Power series for small arguments and positive orders.
C        Abramowitz and Stegun, Eq. 10.1.2.
      F1=1.0D0
      IP=1
      IF(NL.NE.0) THEN
        DO I=1,NL
          IP=IP+2
          F1=F1*X/IP
        ENDDO
      ENDIF
      XI=0.5D0*XI
      SBESJN=1.0D0
      PS=1.0D0
      DO I=1,1000
        IP=IP+2
        PS=-PS*XI/(I*IP)
        SBESJN=SBESJN+PS
        IF(ABS(PS).LT.1.0D-18*ABS(SBESJN)) GO TO 1
      ENDDO
 1    SBESJN=IPH*F1*SBESJN
      RETURN
C  ****  Miller's method for positive orders and intermediate arguments.
C        Abramowitz and Stegun, Eq. 10.1.19.
 2    XI=1.0D0/X
      F2=0.0D0
      F3=1.0D-35
      IP=2*(NL+31)+3
      DO I=1,31
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
      SBESJN=F3
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D30) THEN
          SBESJN=SBESJN/F3
          F2=F2/F3
          F3=1.0D0
        ENDIF
      ENDDO
C
      IF(MAX(ABS(F2),ABS(F3)).LT.1.0D-99) THEN
        STOP '*** Error in SBESJN: Inconsistent low-order values.'
      ENDIF
      FJ1=XI*SIN(X)
      IF(MIN(ABS(FJ1),ABS(F3)).GT.1.0D-15) THEN
        FACT=FJ1/F3
      ELSE
        FJ2=XI*(FJ1-COS(X))
        FACT=FJ2/F2
      ENDIF
      SBESJN=IPH*SBESJN*FACT
      RETURN
C  ****  Recurrence relation for arguments greater than order.
C        Abramowitz and Stegun, Eq. 10.1.19.
 3    XI=1.0D0/X
      F3=XI*SIN(X)
      IF(NL.EQ.0) GO TO 4
      F2=F3
      F3=XI*(F2-COS(X))
      IF(NL.EQ.1) GO TO 4
      IP=1
      DO I=2,NL
        F1=F2
        F2=F3
        IP=IP+2
        F3=IP*XI*F2-F1
      ENDDO
 4    SBESJN=IPH*F3
      RETURN
C  ****  Recurrence relation for negative orders.
C        Abramowitz and Stegun, Eq. 10.1.19.
 5    NL=ABS(NL)
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN
        SBESJN=-1.0D35
        RETURN
      ENDIF
      XI=1.0D0/X
      F3=XI*SIN(X)
      F2=XI*(F3-COS(X))
      IP=3
      DO I=1,NL
        F1=F2
        F2=F3
        IP=IP-2
        F3=IP*XI*F2-F1
        IF(ABS(F3).GT.1.0D35) THEN
          SBESJN=-1.0D35
          RETURN
        ENDIF
      ENDDO
      SBESJN=IPH*F3
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC         Asymptotic expansions          CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SBAS0
C  *********************************************************************
      SUBROUTINE SBAS0(Z,E,L,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Schrodinger radial functions for bound states of
C  modified Coulomb potentials. The expansion is valid only for radii
C  beyond the start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential, must be negative.
C    E ........ energy of the state, must be negative.
C    L ........ orbital angular momentum quantum number.
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CSBAS/):
C    AS ....... 'a' parameter (imaginary part of the wavenumber),
C    ETAS ..... Sommerfeld parameter (imaginary part).
C    RACN ..... radius of convergence of the asymptotic series.
C    SP(1:NTERM) ... coefficients of the asymptotic expansion of P(R)
C               in powers of (RACN/R).
C    SQ(1:NTERM+1) ... coefficients of the asymptotic expansion of
C               Q(R)=P'(R) in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSBAS/SP(MNT),SQ(MNT),AS,ETAS,RACN,NTERM
      DIMENSION BPP(MNT)
C
      IF(E.GT.-1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'SBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(Z.GT.0.0D0) THEN
        WRITE(6,'('' Z ='',1P,E16.8)') Z
        STOP 'SBAS0: The Coulomb charge, Z, must be negative.'
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,'('' L ='',I6)') L
        STOP 'SBAS0: L.LT.0'
      ENDIF
C
      TOL=0.1D0*EPS
      IF(TOL.LT.1.0D-13) THEN
        TOL=1.0D-13
      ELSE IF(TOL.GT.1.0D-7) THEN
        TOL=1.0D-7
      ENDIF
C
      RLAMB=L
      EE=E
 1    CONTINUE
      AS=SQRT(-2.0D0*EE)
      ETAS=-Z/AS
      PA=RLAMB+1.0D0-ETAS
      PB=-RLAMB-ETAS
      R=MAX(RCI,1.0D0)
C
 11   CONTINUE
      IRL=0
      IRU=0
      RL=0.0D0
      RU=10.0D0
      IEIGEN=0
      IEND=0
C
C  ****  Asymptotic expansion. Coefficients and convergence.
C        (ME-7.34 to ME-7.36)
C
 2    CONTINUE
      DO I=1,MNT
        SP(I)=0.0D0
        SQ(I)=0.0D0
        BPP(I)=0.0D0
      ENDDO
C
      RA=R
      TWOARA=-2.0D0*AS*RA
      SP(1)=1.0D0
      SMIN=SP(1)
      IM=0
      P=SP(1)
      DO I=2,MNT
        SP(I)=SP(I-1)*(PA+DBLE(I-2))*(PB+DBLE(I-2))/(DBLE(I-1)*TWOARA)
C
        IF(MAX(ABS(SP(I)),ABS(SP(I-1))).LT.1.0D-90) THEN
          IF(I.LT.MNT-5) THEN  ! E seems to be an exact eigenvalue.
            NTERM=I-2
            IEIGEN=1
            GO TO 3
          ENDIF
        ENDIF
C
        IF(ABS(SP(I)).LT.SMIN) THEN
          IF(SMIN.GT.TOL) THEN
            SMIN=ABS(SP(I))
            IM=I
          ENDIF
        ELSE
          IF(IM.GT.2) THEN
            IF(ABS(SP(I)).GT.ABS(SP(I-1)).AND.
     1         ABS(SP(I-1)).LT.ABS(SP(I-2))) THEN
              SP(I)=0.0D0
              NTERM=I-1
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
C
        P=P+SP(I)
        IF(ABS(P).GT.1.0D35.OR.ABS(P).LT.1.0D-35) THEN
          R=1.1D0*R
          GO TO 11
        ENDIF
        IF(ABS(SP(I)).LT.TOL*ABS(P).AND.ABS(SP(I)).GT.1.0D-90) THEN
          NTERM=I
          IF(R.LT.1.0001D0*RCI) IEND=1
          GO TO 3
         ENDIF
      ENDDO
      NTERM=MNT
C
 3    CONTINUE
      SQ(1)=-AS*SP(1)
      BPP(1)=-AS*SQ(1)
      P=SP(1)
      Q=SQ(1)
      PPP=BPP(1)
      TST=ABS(-0.5D0*PPP+(0.5D0*RLAMB*(RLAMB+1.0D0)/R**2+Z/R-E)*P)
      DO I=2,MIN(NTERM+2,MNT)
        SQ(I)=(ETAS-(I-2))*SP(I-1)/RA-AS*SP(I)
        BPP(I)=(ETAS-(I-2))*SQ(I-1)/RA-AS*SQ(I)
        P=P+SP(I)
        Q=Q+SQ(I)
        PPP=PPP+BPP(I)
        TST=ABS(-0.5D0*PPP+(0.5D0*RLAMB*(RLAMB+1.0D0)/R**2+Z/R-E)*P)
        IF(NTERM.GT.MNT-1) THEN
          IF(ABS(SP(I)).LT.TOL*ABS(P).AND.TST.LT.EPS*ABS(P)) THEN
            NTERM=I
            GO TO 4
          ENDIF
        ENDIF
      ENDDO
 4    CONTINUE
      TST=TST/ABS(P)
C
      IF(IEIGEN.EQ.1) THEN
        IF(TST.GT.EPS) THEN
          WRITE(6,1000)
 1000     FORMAT(/1X,'*** Warning (SBAS0): Eigenstate radial function',
     1      ' did not converge.')
          EE=0.99D0*EE
          IEIGEN=0
          GO TO 1
        ELSE
          IEND=1
          GO TO 5
        ENDIF
      ELSE
        IF(IEND.EQ.1) GO TO 5
      ENDIF
C
      IF(TST.LT.EPS) THEN
        IF(R.LT.1.0001D0*RCI) GO TO 5
        RU=R
        IRU=1
        IF(IRL.EQ.0) THEN
          R=0.5D0*R
          GO TO 2
        ENDIF
      ELSE
        RL=R
        IRL=1
        IF(IRU.EQ.0) THEN
          R=2.0D0*R
          GO TO 2
        ENDIF
      ENDIF
C
      IF(ABS(RU-RL).GT.1.0D-4*RL) THEN
        R=0.5D0*(RL+RU)
        GO TO 2
      ENDIF
      R=RU
C
C  ****  Output convergence radius.
C
 5    CONTINUE
      IF(IEND.EQ.0) THEN
        IEND=1
        GO TO 2
      ENDIF
      RACN=R
      RACONV=RACN
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SBAS
C  *********************************************************************
      SUBROUTINE SBAS(R,P,Q,IER)
C
C     This subroutine calculates Schrodinger radial functions of bound
C  states of modified Coulomb potentials using asymptotic expansions,
C  whose coefficients are precalculated by subroutine SBAS0. The
C  expansions are valid only for radii beyond the start of the Coulomb
C  tail. (ME-7.34, ME-7.36)
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... _unnormalized_ radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.ge.RACN, the asymptotic series converge.
C               =1, for R.lt.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSBAS/SP(MNT),SQ(MNT),AS,ETAS,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        IER=0
        RAOR=RACN/R
        P=0.0D0
        Q=0.0D0
        RPOW=1.0D0
        DO I=1,MIN(NTERM+1,MNT)
          P=P+SP(I)*RPOW
          Q=Q+SQ(I)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        FACT=(2.0D0*AS*R)**ETAS*EXP(-AS*R)
        P=FACT*P
        Q=FACT*Q
        IF(ABS(P).LT.1.0D-75) P=0.0D0
        IF(ABS(Q).LT.1.0D-75) Q=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFAS0
C  *********************************************************************
      SUBROUTINE SFAS0(Z,E,L,PHASE,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of radial Schrodinger functions for free states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential.
C    E ........ kinetic energy.
C    L ........ orbital angular momentum quantum number.
C    PHASE .... inner phase shift (zero for pure Coulomb fields).
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CSFAS/):
C    WAVNUM ... wave number.
C    ETA ...... Sommerfeld parameter.
C    PHASE0 ... asymptotic R-independent phase (in ATOMIC units).
C    RACN ..... radius of convergence of the asymptotic series.
C    A(1:2,1:NTERM) ... matrix of coefficients of the asymptotic
C               expansion of P(R)  in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), INTEGER*4 (I-N),
     1  COMPLEX*16 (C)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/CSFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
      DIMENSION F(2,MNT),G(2,MNT)
C
      IF(E.LT.1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'SFAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(L.LT.0) THEN
        WRITE(6,'('' L ='',I6)') L
        STOP 'SFAS0: L.LT.0'
      ENDIF
C
      IF(EPS.LT.1.0D-15) THEN
        TOL=1.0D-15
      ELSE IF(EPS.GT.1.0D-7) THEN
        TOL=1.0D-7
      ELSE
        TOL=EPS
      ENDIF
C
      NTERM=MNT
      DO I=1,MNT
        A(1,I)=0.0D0
        A(2,I)=0.0D0
        A(3,I)=0.0D0
        A(4,I)=0.0D0
      ENDDO
C
      RLAMB=L
      WAVNUM=SQRT(E+E)
      ETA=Z/WAVNUM
C  ****  Classical turning point. (ME-7.27)
      RTURN=(ETA+SQRT(ETA**2+RLAMB*(RLAMB+1.0D0)))/WAVNUM
      RACN=MAX(0.75D0*RTURN,RCI)
      IF(RACN.LT.1.0D-6) THEN
        WRITE(6,'('' RTURN ='',1P,E16.8)') RTURN
        WRITE(6,'(''   RCI ='',1P,E16.8)') RCI
        STOP 'SFAS0: RCI is less than 1.0D-6 (?).'
      ENDIF
C
 1    CONTINUE
      DO I=1,MNT
        F(1,I)=0.0D0
        G(1,I)=0.0D0
      ENDDO
      RACNI=1.0D0/RACN
C
C  ****  Asymptotic expansion. (ME-7.9, ME-7.10)
C
      CI=DCMPLX(0.0D0,1.0D0)
      CA=CI*ETA-RLAMB
      CB=CI*ETA+RLAMB+1.0D0
      F(1,1)=1.0D0
      G(1,1)=0.0D0
      CTERM=1.0D0
      CSUM1=CTERM
      TMIN1=1.0D0
      NTERM=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM1=CSUM1+CTERM
        TABS1=CDABS(CTERM)
        IF(TABS1.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(1,I+1)=CTERM
        G(1,I+1)=-CI*CTERM
        IF(TABS1.LT.TMIN1) THEN
          TMIN1=TABS1
          NTERM=I+1
          IF(TMIN1.LT.TOL*CDABS(CSUM1)) GO TO 2
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 2    CONTINUE
C
C  ****  R-independent phase shift. (ME-7.8)
C
      DELTC=DELTAC(ETA,RLAMB)
      PHASE0=PHASE+DELTC-ETA*LOG(2.0D0*WAVNUM)-RLAMB*PIH
      RCK=WAVNUM*RACN
C
C  ****  Coefficients of the asymptotic expansion. (ME-7.9, ME-7.10)
C
      DO I=1,NTERM
        A(1,I)=G(1,I)
        A(2,I)=F(1,I)
      ENDDO
      A(3,1)=RCK*A(2,1)/RACN
      A(4,1)=-RCK*A(1,1)/RACN
      A(3,2)=(RCK*A(2,2)-ETA*A(2,1))/RACN
      A(4,2)=(-RCK*A(1,2)+ETA*A(1,1))/RACN
      DO I=3,NTERM
        A(3,I)=(-(I-2)*A(1,I-1)+RCK*A(2,I)-ETA*A(2,I-1))/RACN
        A(4,I)=(-(I-2)*A(2,I-1)-RCK*A(1,I)+ETA*A(1,I-1))/RACN
      ENDDO
      RACONV=RACN
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SFAS
C  *********************************************************************
      SUBROUTINE SFAS(R,P,Q,IER)
C
C     This subroutine calculates radial Schrodinger functions of free
C  states of modified Coulomb potentials using asymptotic expansions,
C  whose coefficients are precalculated by subroutine SFAS0. The
C  expansions are valid only for radii beyond the start of the Coulomb
C  tail.
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.GE.RACN, the asymptotic series converge.
C               =1, for R.LT.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CSFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        RAOR=RACN/R
        A1=0.0D0
        A2=0.0D0
        A3=0.0D0
        A4=0.0D0
        RPOW=1.0D0
        DO N=1,NTERM
          A1=A1+A(1,N)*RPOW
          A2=A2+A(2,N)*RPOW
          A3=A3+A(3,N)*RPOW
          A4=A4+A(4,N)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        PHI=WAVNUM*R-ETA*LOG(R)+PHASE0
        P=A1*COS(PHI)+A2*SIN(PHI)
        Q=A3*COS(PHI)+A4*SIN(PHI)
        IER=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS0
C  *********************************************************************
      SUBROUTINE DBAS0(Z,E,K,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Dirac radial functions for bound states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail. (ME-7.43, ME-7.44)
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential, must be negative.
C    E ........ energy of the state, must be negative.
C    K ........ angular momentum quantum number, kappa.
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CDBAS/):
C    AD ....... 'a' parameter (imaginary part of the wavenumber),
C    ETAD ..... Sommerfeld parameter (imaginary part).
C    RACN ..... radius of convergence of the asymptotic series.
C    DP(1:NTERM) ... coefficients of the asymptotic expansion of P(R)
C               in powers of (RACN/R).
C    DQ(1:NTERM) ... coefficients of the asymptotic expansion of Q(R)
C               in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2)
      COMMON/CDBAS/DP(MNT),DQ(MNT),AD,ETAD,RACN,NTERM
      COMMON/CDBAS1/SP(MNT),NTERM1
      DIMENSION DPP(MNT),DQP(MNT)
C
      IF(E.GT.-1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'DBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(Z.GT.0.0D0) THEN
        WRITE(6,'('' Z ='',1P,E16.8)') Z
        STOP 'DBAS0: The Coulomb charge, Z, must be negative.'
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I4)') K
        STOP 'DBAS0: The K (kappa) value cannot be zero.'
      ENDIF
C
      AK=DBLE(K)
      ZETA=Z/SL
      RLAMB2=AK*AK-ZETA*ZETA
      RLAMB=SQRT(RLAMB2)
      AD=SQRT(-E*(E+TSL2))/SL
      ETAD=SQRT((ZETA*(E+SL2))**2/(-E*(E+TSL2)))
C
      C11=(AK+RLAMB)*(AD*(AK+RLAMB)+E*ZETA/SL)
      C21=-ZETA*(AD*(AK+RLAMB)+E*ZETA/SL)
      C12=-ZETA*(AD*ZETA+E*(AK+RLAMB)/SL)
      C22=(AK+RLAMB)*(AD*ZETA+E*(AK+RLAMB)/SL)
C
      RCID=RCI
 1    CONTINUE
      DO I=1,MNT
        DP(I)=0.0D0
        DQ(I)=0.0D0
        DPP(I)=0.0D0
        DQP(I)=0.0D0
      ENDDO
C  ****  S-Coulomb function with RLAMB.
      CALL DBAS01(AD,ETAD,RLAMB,RCID,RAC1,EPS)
      NTERM=NTERM1
      DO I=1,NTERM
        DP(I)=C11*SP(I)
        DQ(I)=C21*SP(I)
      ENDDO
C  ****  S-Coulomb function with RLAMB-1.0. Usually converges faster.
      CALL DBAS01(AD,ETAD,RLAMB-1.0D0,RAC1,RAC2,EPS)
      IF(RAC2.GT.RAC1) THEN
        RCID=RAC2
        GO TO 1
      ENDIF
      RACN=RAC1
      RACONV=RAC1
      NTERM=MAX(NTERM,NTERM1)
      DO I=1,NTERM
        DP(I)=DP(I)+C12*SP(I)
        DQ(I)=DQ(I)+C22*SP(I)
      ENDDO
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS01
C  *********************************************************************
      SUBROUTINE DBAS01(AS,ETAS,RLAMB,RCI,RACONV,EPS)
C
C     Asymptotic expansion of Schrodinger-Coulomb radial function.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDBAS1/SP(MNT),NTERM1
      DIMENSION SQ(MNT),BPP(MNT)
C
      TOL=0.1D0*EPS
      IF(TOL.LT.1.0D-13) THEN
        TOL=1.0D-13
      ELSE IF(TOL.GT.1.0D-7) THEN
        TOL=1.0D-7
      ENDIF
C
      PA=RLAMB+1.0D0-ETAS
      PB=-RLAMB-ETAS
      R=MAX(RCI,1.0D0)
C
 1    CONTINUE
      IRL=0
      IRU=0
      RL=0.0D0
      RU=10.0D0
      IEND=0
C
C  ****  Asymptotic expansion. Coefficients and convergence.
C
 2    CONTINUE
      DO I=1,MNT
        SP(I)=0.0D0
        SQ(I)=0.0D0
        BPP(I)=0.0D0
      ENDDO
C
      RA=R
      TWOARA=-2.0D0*AS*RA
      SP(1)=1.0D0
      SMIN=SP(1)
      IM=0
      P=SP(1)
      DO I=2,MNT
        SP(I)=SP(I-1)*(PA+DBLE(I-2))*(PB+DBLE(I-2))/(DBLE(I-1)*TWOARA)
C
        IF(MAX(ABS(SP(I)),ABS(SP(I-1))).LT.1.0D-90) THEN
          IF(I.LT.MNT-5) THEN  ! E seems to be an exact eigenvalue.
            NTERM1=I-2
            SP(I)=0.0D0
            SP(I-1)=0.0D0
            GO TO 3
          ENDIF
        ENDIF
C
        IF(ABS(SP(I)).LT.SMIN) THEN
          IF(SMIN.GT.TOL) THEN
            SMIN=ABS(SP(I))
            IM=I
          ENDIF
        ELSE
          IF(IM.GT.2) THEN
            IF(ABS(SP(I)).GT.ABS(SP(I-1)).AND.
     1         ABS(SP(I-1)).LT.ABS(SP(I-2))) THEN
              SP(I)=0.0D0
              NTERM1=I-1
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
C
        P=P+SP(I)
        IF(ABS(P).GT.1.0D35.OR.ABS(P).LT.1.0D-35) THEN
          R=1.1D0*R
          GO TO 1
        ENDIF
        IF(ABS(SP(I)).LT.TOL*ABS(P).AND.ABS(SP(I)).GT.1.0D-90) THEN
          NTERM1=I
          IF(R.LT.1.0001D0*RCI) IEND=1
          GO TO 3
         ENDIF
      ENDDO
      NTERM1=MNT
C
 3    CONTINUE
      SQ(1)=-AS*SP(1)
      BPP(1)=-AS*SQ(1)
      P=SP(1)
      Q=SQ(1)
      PPP=BPP(1)
      TST=ABS(PPP/P-AS**2+2.0D0*AS*ETAS/R-RLAMB*(RLAMB+1.0D0)/R**2)
      DO I=2,MIN(NTERM1+2,MNT)
        SQ(I)=(ETAS-(I-2))*SP(I-1)/RA-AS*SP(I)
        BPP(I)=(ETAS-(I-2))*SQ(I-1)/RA-AS*SQ(I)
        P=P+SP(I)
        Q=Q+SQ(I)
        PPP=PPP+BPP(I)
        TST=ABS(PPP/P-AS**2+2.0D0*AS*ETAS/R-RLAMB*(RLAMB+1.0D0)/R**2)
        IF(NTERM1.GT.MNT-1) THEN
          IF(ABS(SP(I)).LT.TOL*ABS(P).AND.TST.LT.EPS) THEN
            NTERM1=I
            GO TO 4
          ENDIF
        ENDIF
      ENDDO
 4    CONTINUE
      IF(IEND.EQ.1) GO TO 5
C
      IF(TST.LT.EPS) THEN
        IF(R.LT.1.0001D0*RCI) GO TO 5
        RU=R
        IRU=1
        IF(IRL.EQ.0) THEN
          R=0.5D0*R
          GO TO 2
        ENDIF
      ELSE
        RL=R
        IRL=1
        IF(IRU.EQ.0) THEN
          R=2.0D0*R
          GO TO 2
        ENDIF
      ENDIF
C
      IF(ABS(RU-RL).GT.1.0D-4*RL) THEN
        R=0.5D0*(RL+RU)
        GO TO 2
      ENDIF
      R=RU
C
C  ****  Output convergence radius.
C
 5    CONTINUE
      IF(IEND.EQ.0) THEN
        IEND=1
        GO TO 2
      ENDIF
      RACONV=R
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DBAS
C  *********************************************************************
      SUBROUTINE DBAS(R,P,Q,IER)
C
C     This subroutine calculates Dirac radial functions of bound states
C  of modified Coulomb potentials using asymptotic expansions, whose
C  coefficients are precalculated by subroutine DBAS0. The expansions
C  are valid only for radii beyond the start of the Coulomb tail.
C  (ME-7.43, ME-7.44)
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... _unnormalized_ radial functions at R; Q=P'.
C    IER ...... error flag:
C               =0, for R.GE.RACN, the asymptotic series converge.
C               =1, for R.LT.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDBAS/DP(MNT),DQ(MNT),AD,ETAD,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        IER=0
        RAOR=RACN/R
        P=0.0D0
        Q=0.0D0
        RPOW=1.0D0
        DO I=1,NTERM
          P=P+DP(I)*RPOW
          Q=Q+DQ(I)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        FACT=(2.0D0*AD*R)**ETAD*EXP(-AD*R)
        P=FACT*P
        Q=FACT*Q
        IF(ABS(P).LT.1.0D-75) P=0.0D0
        IF(ABS(Q).LT.1.0D-75) Q=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFAS0
C  *********************************************************************
      SUBROUTINE DFAS0(Z,E,K,PHASE,RCI,RACONV,EPS)
C
C     This subroutine determines the coefficients of the asymptotic
C  expansion of Dirac radial functions for free states of modified
C  Coulomb potentials. The expansion is valid only for radii beyond the
C  start of the Coulomb tail.
C
C  Input parameters:
C    Z ........ strength of the Coulomb potential.
C    E ........ kinetic energy.
C    K ........ angular momentum quantum number kappa.
C    PHASE .... inner phase shift (zero for pure Coulomb fields).
C    RCI ...... user cutoff radius. The asymptotic expansion is needed
C               only for radii larger than RCI.
C    EPS ...... global tolerance, i.e. allowed relative error in the
C               summation of the asymptotic series.
C
C  Output argument:
C    RACONV ... effective convergence radius.
C               A value of RACONV larger than RCI indicates that the
C               asymptotic expansion converges only for radii larger
C               than RACONV.
C
C  Output parameters (through common block /CDFAS/):
C    WAVNUM ... wave number.
C    ETA ...... Sommerfeld parameter.
C    PHASE0 ... asymptotic R-independent phase (in ATOMIC units).
C    RACN ..... radius of convergence of the asymptotic series.
C    A(1:4,1:NTERM) ... Matrix of coefficients of the asymptotic
C               expansions in powers of (RACN/R).
C    NTERM .... number of terms in the asymptotic series (.LE.MNT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), INTEGER*4 (I-N),
     1  COMPLEX*16 (C)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/CDFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
      DIMENSION F(2,MNT),G(2,MNT)
C
      IF(E.LT.1.0D-6) THEN
        WRITE(6,'('' E ='',1P,E16.8)') E
        STOP 'DBAS0: The energy value, E, must be less than -1.0D-4.'
      ENDIF
C
      IF(K.EQ.0) THEN
        WRITE(6,'('' K ='',I4)') K
        STOP 'DBAS0: The K (kappa) value cannot be zero.'
      ENDIF
C
      IF(EPS.LT.1.0D-15) THEN
        TOL=1.0D-15
      ELSE IF(EPS.GT.1.0D-7) THEN
        TOL=1.0D-7
      ELSE
        TOL=EPS
      ENDIF
C
      NTERM=MNT
      DO I=1,MNT
        A(1,I)=0.0D0
        A(2,I)=0.0D0
        A(3,I)=0.0D0
        A(4,I)=0.0D0
      ENDDO
C
      ZETA=Z/SL
      RLAMB2=K*K-ZETA*ZETA
      RLAMB=SQRT(RLAMB2)
      PC2=E*(E+TSL2)
      PC=SQRT(PC2)
      WAVNUM=PC/SL
      W=E+SL2
      ETA=ZETA*W/PC
C  ****  Classical turning point. (ME-7.27)
      RTURN=(ETA+SQRT(ETA**2+RLAMB*(RLAMB+1.0D0)))/WAVNUM
      RACN=MAX(0.75D0*RTURN,RCI)
      IF(RACN.LT.1.0D-6) THEN
        WRITE(6,'('' RTURN ='',1P,E16.8)') RTURN
        WRITE(6,'(''   RCI ='',1P,E16.8)') RCI
        STOP 'DFAS0: RCI is less than 1.0D-6 (?).'
      ENDIF
C
 1    CONTINUE
      DO I=1,MNT
        F(1,I)=0.0D0
        G(1,I)=0.0D0
        F(2,I)=0.0D0
        G(2,I)=0.0D0
      ENDDO
      RACNI=1.0D0/RACN
C
C  ****  Asymptotic expansion of the first Coulomb function.
C        (ME-7.4, ME-7.5)
C
      CI=DCMPLX(0.0D0,1.0D0)
      CA=CI*ETA-RLAMB
      CB=CI*ETA+RLAMB+1.0D0
      F(1,1)=1.0D0
      G(1,1)=0.0D0
      CTERM=1.0D0
      CSUM1=CTERM
      TMIN1=1.0D0
      NT1=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM1=CSUM1+CTERM
        TABS1=CDABS(CTERM)
        IF(TABS1.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(1,I+1)=CTERM
        G(1,I+1)=-CI*CTERM
        IF(TABS1.LT.TMIN1) THEN
          TMIN1=TABS1
          NT1=I+1
          IF(TMIN1.LT.TOL*CDABS(CSUM1)) GO TO 2
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 2    CONTINUE
C
C  ****  Asymptotic expansion of the second Coulomb function.
C        (ME-7.4, ME-7.5)
C
      CA=CA+1.0D0
      CB=CB-1.0D0
      F(2,1)=1.0D0
      G(2,1)=0.0D0
      CTERM=1.0D0
      CSUM2=CTERM
      TMIN2=1.0D0
      NT2=1
      DO I=1,MNT-1
        N=I-1
        CTERM=CTERM*(CA+N)*(CB+N)*RACNI/(DBLE(I)*2*CI*WAVNUM)
        CSUM2=CSUM2+CTERM
        TABS2=CDABS(CTERM)
        IF(TABS2.GT.1.0D35) THEN
          RACN=1.05D0*RACN
          GO TO 1
        ENDIF
        F(2,I+1)=CTERM
        G(2,I+1)=-CI*CTERM
        IF(TABS2.LT.TMIN2) THEN
          TMIN2=TABS2
          NT2=I+1
          IF(TMIN2.LT.TOL*CDABS(CSUM2)) GO TO 3
        ENDIF
      ENDDO
      RACN=1.05D0*RACN
      GO TO 1
 3    CONTINUE
C
C  ****  R-independent phase shift. (ME-7.21)
C
      DELTC=-CI*CLGAM(RLAMB+CI*ETA)
      PHASE0=PHASE+DELTC-ETA*LOG(2.0D0*WAVNUM)-(RLAMB-1.0D0)*PIH
C
C  ****  Coefficients of the asymptotic expansion. (ME-7.17 to ME-7.25)
C
      DTHETA=-PIH+ATAN2(ETA,RLAMB)
      RNORM=SQRT((ZETA*(W+SL2))**2+(K+RLAMB)**2*PC2)
      IF(ABS(RNORM).LT.1.0D-14) THEN
C  ****  Spherical Bessel functions if Z=0 and K.lt.0.
        TCOS=0.0D0
        TSIN=1.0D0
        A11=0.0D0
        A12=1.0D0
        A31=SQRT(E/(W+SL2))
        A32=0.0D0
      ELSE
        TCOS=COS(DTHETA)
        TSIN=SIN(DTHETA)
        RNORM=1.0D0/(RLAMB*RNORM)
        IF(ZETA.LT.0.0D0.AND.K.LT.0) RNORM=-RNORM
        RLA=SQRT(RLAMB2+ETA*ETA)
        A11=RNORM*(K+RLAMB)*RLA*PC
        A12=RNORM*ZETA*(RLAMB*SL2-K*W)
        A31=-RNORM*ZETA*RLA*PC
        A32=-RNORM*(K+RLAMB)*(RLAMB*SL2-K*W)
      ENDIF
C
      NTERM=MAX(NT1,NT2)
      DO I=1,NTERM
        IF(I.LE.NT1) THEN
          A(1,I)=A11*(F(1,I)*TSIN+G(1,I)*TCOS)
          A(2,I)=A11*(F(1,I)*TCOS-G(1,I)*TSIN)
          A(3,I)=A31*(F(1,I)*TSIN+G(1,I)*TCOS)
          A(4,I)=A31*(F(1,I)*TCOS-G(1,I)*TSIN)
        ENDIF
        IF(I.LE.NT2) THEN
          A(1,I)=A(1,I)+A12*G(2,I)
          A(2,I)=A(2,I)+A12*F(2,I)
          A(3,I)=A(3,I)+A32*G(2,I)
          A(4,I)=A(4,I)+A32*F(2,I)
        ENDIF
      ENDDO
      RACONV=RACN
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DFAS
C  *********************************************************************
      SUBROUTINE DFAS(R,P,Q,IER)
C
C     This subroutine calculates Dirac radial functions of free states
C  of modified Coulomb potentials using asymptotic expansions, whose
C  coefficients are precalculated by subroutine DFAS0. The expansions
C  are valid only for radii beyond the start of the Coulomb tail.
C
C  Input parameters:
C    R ........ radius.
C
C  Output parameters:
C    P,Q ...... radial functions at R.
C    IER ...... error flag:
C               =0, for R.ge.RACN, the asymptotic series converge.
C               =1, for R.lt.RACN, the series do not converge.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CDFAS/A(4,MNT),WAVNUM,ETA,PHASE0,RACN,NTERM
C
      IF(R.LT.RACN) THEN
        P=0.0D0
        Q=0.0D0
        IER=1
      ELSE
        RAOR=RACN/R
        A1=0.0D0
        A2=0.0D0
        A3=0.0D0
        A4=0.0D0
        RPOW=1.0D0
        DO N=1,NTERM
          A1=A1+A(1,N)*RPOW
          A2=A2+A(2,N)*RPOW
          A3=A3+A(3,N)*RPOW
          A4=A4+A(4,N)*RPOW
          RPOW=RPOW*RAOR
        ENDDO
        PHI=WAVNUM*R-ETA*LOG(R)+PHASE0
        P=A1*COS(PHI)+A2*SIN(PHI)
        Q=A3*COS(PHI)+A4*SIN(PHI)
        IER=0
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC          Cubic spline interpolation          CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C     SUBROUTINE SPLIN0(X,Y,A,B,C,D,S1,SN,N)
C     SUBROUTINE FINDI(XC,X,N,I)
C       FUNCTION SPLVAL(XC,X,A,B,C,D,N)
C     SUBROUTINE SPLERR(X,Y,S1,SN,ERR,N,IWR)
C     SUBROUTINE SPLSET(FUNC,XL,XU,X,Y,TOL,ERR,NPM,NFIX,NU,N)
C       FUNCTION SPLINT(X,A,B,C,D,XL,XU,N,NPOW)
C
C  *********************************************************************
C                       SUBROUTINE SPLINE
C  *********************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
C
C     This subroutine determines the coefficients of a piecewise cubic
C  spline that interpolates the input table (X,Y) of function values.
C  Duplicated abscissas are considered as discontinuities; a separate
C  spline is used for each interval between consecutive discontinuities,
C  with 'natural' end-point shape (null second derivative) at the
C  discontinuities.
C
C  Input:
C     X(1:N) ..... grid points (must be in non-decreasing order).
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N).
C     N .......... number of grid points.
C
C  Output:
C     A,B,C,D(1:N) ... spline coefficients.
C
C  Other subprograms used: subroutine SPLIN0.
C
C  The interpolating cubic polynomial in the I-th interval, from X(I) to
C  X(I+1), is
C               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      PARAMETER (NM=25000)  ! Auxiliary storage.
      DIMENSION XP(NM),YP(NM),AP(NM),BP(NM),CP(NM),DP(NM)
C
      SS1=S1
      SSN=SN
C
      IO=0
      I=0
      NP=0
 1    I=I+1
      NP=NP+1
      XP(NP)=X(I)
      YP(NP)=Y(I)
      IF(I.EQ.N) GO TO 2
C
      IF(ABS(X(I+1)-X(I)).GT.EPS*MAX(ABS(X(I)),ABS(X(I+1)))) THEN
        GO TO 1
      ELSE
        X(I)=X(I+1)
      ENDIF
 2    CONTINUE
C
      IF(NP.LT.2) THEN
        WRITE(6,10)
 10     FORMAT(1X,'*** Error in SPLINE: More than 2 coinciding ',
     1    'abscissas.',/5X,'Interpolation is not possible.')
        OPEN(33,FILE='SPLINE-error.dat')
          WRITE(33,*) '# Error at I =',I
          DO J=1,N
            WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
          ENDDO
        CLOSE(33)
        STOP 'SPLINE: More than 2 coinciding abscissas.'
      ELSE IF(NP.EQ.2) THEN  ! Linear interpolation.
        AP(1)=(XP(2)*YP(1)-XP(1)*YP(2))/(XP(2)-XP(1))
        BP(1)=(YP(2)-YP(1))/(XP(2)-XP(1))
        CP(1)=0.0D0
        DP(1)=0.0D0
        AP(2)=AP(1)
        BP(2)=BP(1)
        CP(2)=CP(1)
        DP(2)=DP(1)
      ELSE IF(NP.EQ.3) THEN  ! Quadratic interpolation.
        BB=((XP(1)-XP(2))**2*(YP(3)-YP(2))
     1    -(XP(3)-XP(2))**2*(YP(1)-YP(2)))
     2    /((XP(3)-XP(1))*(XP(3)-XP(2))*(XP(2)-XP(1)))
        CC=((XP(3)-XP(2))*(YP(1)-YP(2))
     1    -(XP(1)-XP(2))*(YP(3)-YP(2)))
     2    /((XP(3)-XP(1))*(XP(3)-XP(2))*(XP(2)-XP(1)))
        AP(1)=YP(2)-BB*XP(2)+CC*XP(2)**2
        BP(1)=BB-2.0D0*CC*XP(2)
        CP(1)=CC
        DP(1)=0.0D0
        AP(2)=AP(1)
        BP(2)=BP(1)
        CP(2)=CP(1)
        DP(2)=DP(1)
        AP(3)=AP(1)
        BP(3)=BP(1)
        CP(3)=CP(1)
        DP(3)=DP(1)
      ELSE
        IF(NP.GT.NM) THEN
          WRITE(6,11)
 11       FORMAT(1X,'*** Error in SPLINE: too many grid points ',
     1      'between two',/5X,'discontinuities of Y(X). ',
     2      /5X,'Details in file ''SPLIN0-error.dat''.')
          WRITE(6,*) '     NP =',NP
          STOP 'SPLINE: NP is larger than NPM.'
        ENDIF
        IF(IO.EQ.1) THEN
          SP1=SS1
        ELSE
          SP1=0.0D0
        ENDIF
        IF(I.EQ.N) THEN
          SPN=SSN
        ELSE
          SPN=0.0D0
        ENDIF
        CALL SPLIN0(XP,YP,AP,BP,CP,DP,SP1,SPN,NP)
      ENDIF
C
      DO J=1,NP
        IO=IO+1
        A(IO)=AP(J)
        B(IO)=BP(J)
        C(IO)=CP(J)
        D(IO)=DP(J)
      ENDDO
      IF(I.LT.N) THEN
        NP=0
        GO TO 1
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPLIN0
C  *********************************************************************
      SUBROUTINE SPLIN0(X,Y,A,B,C,D,S1,SN,N)
C
C     Initialization of cubic spline interpolation of tabulated data.
C  It is assumed that the function and its first two derivatives are
C  continuous.
C
C  Input:
C     X(1:N) ..... grid points (the X values must be in strictly
C                      increasing order).
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N). The natural
C                  spline corresponds to taking S1=SN=0.
C     N .......... number of grid points.
C  Output:
C     A,B,C,D(1:N) ... spline coefficients.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF(N.LT.4) THEN
        WRITE(6,10) N
 10     FORMAT(1X,'*** Error in SPLIN0: interpolation cannot be ',
     1    'performed with',I4,' points.',
     2    /5X,'Details in file ''SPLIN0-error.dat''.')
        OPEN(33,FILE='SPLIN0-error.dat')
          DO J=1,N
            WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
          ENDDO
        CLOSE(33)
        STOP 'SPLIN0: N is less than 4.'
      ENDIF
      N1=N-1
      N2=N-2
C  ****  Auxiliary arrays H(=A) and DELTA(=D).
      DO I=1,N1
        A(I)=X(I+1)-X(I)  ! h_i
        IF(A(I).LT.EPS*MAX(ABS(X(I)),ABS(X(I+1)))) THEN
          WRITE(6,11)
 11       FORMAT(1X,'*** Error in SPLIN0: X values not in',
     1      'increasing order.',
     2      /5X,'Details in file ''SPLIN0-error.dat''.')
          OPEN(33,FILE='SPLIN0-error.dat')
            WRITE(33,'(A,I5)') 'Order error at I =',I
            DO J=1,N
              WRITE(33,'(I5,1P,2E18.10)') J,X(J),Y(J)
            ENDDO
          CLOSE(33)
          STOP 'SPLIN0: X values not in increasing order.'
        ENDIF
        D(I)=(Y(I+1)-Y(I))/A(I)  ! delta_i
      ENDDO
C  ****  Symmetric coefficient matrix (augmented).
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))  ! C_i
      ENDDO
      DO K=N1,2,-1
        D(K)=6.0D0*(D(K)-D(K-1))  ! D_i
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
C  ****  Gauss solution of the tridiagonal system.
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
C  ****  The SIGMA coefficients are stored in array D.
      D(N)=SN
      D(N1)=D(N1)/B(N2)
      DO K=N2,2,-1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
C  ****  Spline coefficients.
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)
     1      +HI*(Y(I)*X(I+1)-Y(I+1)*X(I))
     2      +(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
        B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)
     1      +HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
        C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
        D(I)=(HI/6.0D0)*(SI1-SI)
      ENDDO
C  ****  Quadratic extrapolation for X.GT.X(N). Natural spline if SN=0.
      A(N)=A(N1)+D(N1)*X(N)**3
      B(N)=B(N1)-3.0D0*D(N1)*X(N)**2
      IF(ABS(SN).LT.1.0D-16) THEN
        C(N)=0.0D0
      ELSE
        C(N)=C(N1)+3.0D0*D(N1)*X(N)
      ENDIF
      D(N)=0.0D0
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE FINDI
C  *********************************************************************
      SUBROUTINE FINDI(XC,X,N,I)
C
C     This subroutine finds the interval (X(I),X(I+1)) that contains the
C  value XC by using the binary search algorithm.
C
C  Input:
C     XC ............. point to be located.
C     X(1:N) ......... grid points.
C     N  ............. number of grid points.
C  Output:
C     I .............. interval index.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N)
C
      IF(XC.GT.X(N)) THEN
        I=N
      ELSE IF(XC.LT.X(1)) THEN
        I=1
      ELSE
        I=1
        I1=N
 1      IT=(I+I1)/2
        IF(XC.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SPLVAL
C  *********************************************************************
      FUNCTION SPLVAL(XC,X,A,B,C,D,N)
C
C     This function gives the value of a cubic spline at the point XC;
C  quadratic extrapolation is used for points outside the interval
C  (X(1),X(N)).
C
C  Input:
C     XC ............. spline argument.
C     N .............. number of grid points.
C     X(1:N) ......... grid points.
C     A,B,C,D(1:N) ... spline coefficients.
C
C  Output:
C     SPLVAL ......... value of the spline function at XC.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C
      IF(XC.LT.X(1)) THEN
C  ****  Quadratic extrapolation for X.LT.X(1). Natural spline if S1=0.
        A0=A(1)+D(1)*X(1)**3
        B0=B(1)-3.0D0*D(1)*X(1)**2
        C0=C(1)+3.0D0*D(1)*X(1)
        SPLVAL=A0+XC*(B0+XC*C0)
      ELSE IF(XC.GT.X(N)) THEN
        SPLVAL=A(N)+XC*(B(N)+XC*C(N))
      ELSE
        I=1
        I1=N
 1      IT=(I+I1)/2
        IF(XC.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
        SPLVAL=A(I)+XC*(B(I)+XC*(C(I)+XC*D(I)))
      ENDIF
      RETURN
      END
C  *********************************************************************
C                        SUBROUTINE SPLERR
C  *********************************************************************
      SUBROUTINE SPLERR(X,Y,S1,SN,ERR,N,IWR)
C
C     This subroutine estimates the error introduced by the cubic spline
C  interpolation of a table X(I),Y(I) (I=1:N). The interpolation error
c  in the vicinity of X(K) is approximated by the difference between
C  Y(K) and the value obtained from the spline that interpolates the
C  table with the K-th point removed. ERR is the largest relative error
C  in the table.
C
C     When IWR is greater than zero, a table of relative differences is
C  written in a file named 'SPLERR.dat', which is open as UNIT=IWR.
C
C  Input:
C     X(1:N) ..... grid points.
C     Y(1:N) ..... corresponding function values.
C     S1,SN ...... second derivatives at X(1) and X(N). The natural
C                  spline corresponds to taking S1=SN=0.
C     N .......... number of grid points.
C     IWR ........ printing flag.
C  Output:
C     ERR ........ estimated largest relative error of the spline.
C
C  Other subprograms used: subroutines SPLINE and SPLIN0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10,RES=2.5D-5)
      PARAMETER (NM=25000)
      DIMENSION X(N),Y(N),XT(NM),YT(NM)
      DIMENSION XS(NM),YS(NM),A(NM),B(NM),C(NM),D(NM)
      COMMON/SPLCOM/IME  ! Grid point with the largest error.
C
      ERR=0.0D0
      IME=0
C
      IF(N.LT.5) THEN
        WRITE(6,10) N
 10     FORMAT(1X,'*** SPLERR: input table with N = ',I5,
     1    ' data points.',/5X,'N must be greater than 4.')
        RETURN
      ENDIF
C
      IF(N.GT.NM) THEN
        WRITE(6,11) N,NM
 11     FORMAT(1X,'*** Error in SPLERR: input table with N = ',I5,
     1    ' data points.',/5X,'N must be less than NM = ',I5,'.')
        STOP 'SPLERR: Too many data points.'
      ENDIF
C
      IF(IWR.GT.0) THEN
        OPEN(33,FILE='SPLERR.dat')
        WRITE(33,12)
      ENDIF
 12   FORMAT(1X,'# Test of the spline interpolation (output fro',
     1  'm subroutine SPLERR)',/1X,'#',8X,'x',16X,'y(x)',10X,
     1  'y_spline(x)',6X,'rel.dif.')
C
      YMAX=0.0D0
      DO I=1,N
        YMAX=MAX(YMAX,ABS(Y(I)))
      ENDDO
      YCUT=YMAX*1.0D-35
C
C  ****  Locating the discontinuities.
C
      IT=0
      NP=0
      SS1=S1
      SSN=0.0D0
 1    IT=IT+1
      IF(IT.EQ.N) SSN=SN
      NP=NP+1
      XT(NP)=X(IT)
      YT(NP)=Y(IT)
      IF(IT.EQ.N) GO TO 2
      IF(X(IT+1)-X(IT).GT.EPS*MAX(ABS(X(IT)),ABS(X(IT+1)))) GO TO 1
 2    CONTINUE
C
      IF(NP.LT.3) THEN
        NP=0
        IF(IT.EQ.N) GO TO 3
        SS1=0.0D0
        SSN=0.0D0
        GO TO 1
      ENDIF
C  ****  Removing individual grid points to estimate the error.
      NP1=NP-1
      DO K=2,NP1
        IF(XT(K+1)-XT(K-1).GT.MAX(RES*ABS(XT(K)),1.0D-10)) THEN
          DO J=1,NP1
            IF(J.LT.K) THEN
              XS(J)=XT(J)
              YS(J)=YT(J)
            ELSE
              XS(J)=XT(J+1)
              YS(J)=YT(J+1)
            ENDIF
          ENDDO
C  ****  Points near the zeros of the spline are not analyzed.
          IF(ABS(YT(K)).GT.YCUT) THEN
            IF(NP1.GT.3) THEN
              CALL SPLIN0(XS,YS,A,B,C,D,SS1,SSN,NP1)
            ELSE
              CALL SPLINE(XS,YS,A,B,C,D,SS1,SSN,NP1)
            ENDIF
            YK=A(K-1)+XT(K)*(B(K-1)+XT(K)*(C(K-1)+XT(K)*D(K-1)))
            ERRP=ABS(YK/YT(K)-1.0D0)
            IF(ERRP.GT.ERR) THEN
              ERR=ERRP
              IME=IT-NP+K
            ENDIF
          ELSE
            YK=0.0D0
            ERRP=Y(K)
          ENDIF
          IF(IWR.GT.0) THEN
            IF(ERRP.GT.1.0D-35) THEN
               WRITE(33,13) XT(K),YT(K),YK,ERRP
            ELSE
               WRITE(33,13) XT(K),YT(K),YK
            ENDIF
 13         FORMAT(1X,1P,3E18.10,E12.4)
          ENDIF
        ENDIF
      ENDDO
C
      IF(IT.LT.N) THEN
        NP=0
        SS1=0.0D0
        SSN=0.0D0
        GO TO 1
      ENDIF
C
 3    CONTINUE
      IF(IWR.GT.0) CLOSE(33)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPLSET
C  *********************************************************************
      SUBROUTINE SPLSET(FUNC,XL,XU,X,Y,TOL,ERR,NPM,NFIX,NU,N)
C
C     This subroutine determines a table (X,Y) of the external function
C  FUNC(X) suited for natural cubic spline approximation in the interval
C  (XL,XU).
C
C  The X grid is built starting from an initial subgrid consisting of NU
C  uniformly spaced points in (XL,XU). If NU is negative, the initial
C  grid points are logarithmically spaced within that interval. The grid
C  is refined by adding new points in regions where the interpolation
C  seems to have the larger errors. Optionally a number NFIX of X-values
C  can be fixed; they are to be entered in the first NFIX positions of
C  the input array X(). TOL is the tolerance; the subroutine ends when
C  the largest relative error is estimated to be less than TOL.
C
C  Input arguments:
C    FUNC ..... name of the external function.
C    XL,XU .... end points of the considered interval.
C    TOL ...... tolerance, desired relative error of the interpolation.
C    NPM ...... physical dimension of arrays X and Y.
C    NFIX ..... number of fixed points. Their abscissas must be entered
C               as the first NFIX elements of the array X. A doubled
C               value is considered as a discontinuity.
C    NU ....... number of points in the initial 'uniform' subgrid.
C    N ........ desired number of points in the table. Must be larger
C               than NFIX+ABS(NU) and less than NPM.
C
C  Output arguments:
C    X(1:N),Y(1:N) ... generated arrays of abscisas and function values.
C    ERR ...... estimate of the interpolation error.
C    N ........ number of points in the generated grid. It may differ
C               from the input value because the subroutine stops
C               adding grid points as soon as the required tolerance
C               or the allowed minimum spacing (RES) is attained.
C
C  Other subprograms used: subroutines SPLINE and SPLERR.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10, RES=2.5D-5)
      DIMENSION X(NPM),Y(NPM)
      PARAMETER (LH=8, L=LH+LH)
      DIMENSION XS(L),YS(L),AS(L),BS(L),CS(L),DS(L)
      COMMON/SPLCOM/IME  ! Grid point with the largest error.
C
      TTOL=MAX(1.0D-8,TOL)
      NFFIX=MAX(NFIX,0)
      IF(NPM.LE.NFFIX+ABS(NU).OR.N.GT.NPM) THEN
        WRITE(6,10)
 10     FORMAT(1X,'*** Error in SPLSET: The physical dimension ',
     1    'NPM is too small.')
        WRITE(6,*) ' NPM =',NPM
        WRITE(6,*) 'NFIX =',NFFIX
        WRITE(6,*) '  NU =',ABS(NU)
        WRITE(6,*) ' NPR =',MAX(NFFIX+ABS(NU),N)
        STOP 'SPLSET: NPM must be larger than NPR.'
      ENDIF
      IF(XU.LT.XL) THEN
        DX=XU
        XU=XL
        XL=DX
      ENDIF
      XUC=XL+1.0D5*EPS*MAX(ABS(XL),ABS(XU))
      IF(XU.LT.XUC) THEN
        WRITE(6,11)
 11     FORMAT(1X,'*** Error in SPLSET: the interval endpoints ',
     1    ' are not sufficiently spaced.')
        WRITE(6,*) ' XL =',XL
        WRITE(6,*) ' XU =',XU
        WRITE(6,*) 'XUC =',XUC
        STOP 'SPLSET: XU must be large than XUC.'
      ENDIF
C
      NP=N
      N=0
      IF(NFFIX.GT.0) THEN
        DO I=1,NFFIX
          Y(I)=X(I)
        ENDDO
        DO I=1,NFFIX
          IF(Y(I).GT.XL.AND.Y(I).LT.XU) THEN
            N=N+1
            X(N)=Y(I)
          ENDIF
        ENDDO
      ENDIF
C
      X(N+1)=XL
      X(N+2)=XU
      N=N+2
C
      IF(NU.NE.0) THEN
        IF(NU.LT.0) THEN
          XXL=MAX(XL,MIN(1.0D-7,XU*0.1D0))
          FACT=(XU/XXL)**(1.0D0/DBLE(ABS(NU)+1))
          DO I=1,ABS(NU)
            XX=XXL*FACT**I
            DMIN=1.0D99
            DO J=1,N
              DMIN=MIN(DMIN,ABS(XX-X(J)))
            ENDDO
            IF(DMIN.GT.RES) THEN
              N=N+1
              X(N)=XX
            ENDIF
          ENDDO
        ELSE IF(NU.GT.0) THEN
          DX=(XU-XL)/DBLE(NU+1)
          DO I=1,NU
            XX=XL+DX*I
            DMIN=1.0D99
            DO J=1,N
              DMIN=MIN(DMIN,ABS(XX-X(J)))
            ENDDO
            IF(DMIN.GT.RES) THEN
              N=N+1
              X(N)=XX
            ENDIF
          ENDDO
        ENDIF
      ENDIF
C
C  ****  Clean the initial grid.
C
      IFILL=1
 1    CONTINUE
C  ****  Sort grid points in increasing order.
      N1=N-1
      DO I=1,N1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            SAVE=X(I)
            X(I)=X(J)
            X(J)=SAVE
          ENDIF
        ENDDO
      ENDDO
C  ****  Remove doubled grid points that are not discontinuities.
      IF(X(2)-X(1).LT.MAX(EPS*ABS(X(1)),EPS*ABS(X(2)),1.0D-16)) THEN
        X(2)=X(N)
        N=N-1
        GO TO 1
      ENDIF
C
      IF(X(N)-X(N-1).LT.MAX(EPS*ABS(X(N-1)),EPS*ABS(X(N)),1.0D-16)) THEN
        X(N-1)=X(N)
        N=N-1
        GO TO 1
      ENDIF
C
      DO I=1,N-2
        IF((X(I+1)-X(I).LT.EPS*ABS(X(I))).AND.
     1    (X(I+2)-X(I+1).LT.EPS*ABS(X(I)))) THEN
          X(I+2)=X(N)
          N=N-1
          GO TO 1
        ENDIF
      ENDDO
C  ****  Ensure that duplicated abscissas are equal.
      DO I=2,N-2
        IF(X(I).GT.X(I+1)-EPS*ABS(X(I+1))) X(I+1)=X(I)
      ENDDO
C  ****  Add 4 points between each pair of consecutive discontinuities.
      IF(IFILL.EQ.1) THEN
        XN=X(N)
        XN1=X(N1)
        DO I=2,N1-1
          IF((X(I)-X(I-1).LT.EPS*MAX(ABS(X(I)),ABS(X(I-1)))).AND.
     1      (X(I+2)-X(I+1).LT.EPS*MAX(ABS(X(I+1)),ABS(X(I+2))))) THEN
            DX=(X(I+1)-X(I))/5.0D0
            DO J=1,4
              N=N+1
              X(N)=X(I)+J*DX
            ENDDO
          ENDIF
        ENDDO
C
        IF(XN-XN1.GT.5.0D0*RES) THEN
          DX=(XN-XN1)/5.0D0
          DO J=1,4
            N=N+1
            X(N)=XN1+J*DX
          ENDDO
        ENDIF
C
        IF(N1.GT.1) THEN
          IF(X(2)-X(1).GT.5.0D0*RES) THEN
            DX=(X(2)-X(1))/5.0D0
            DO J=1,4
              N=N+1
              X(N)=X(1)+J*DX
            ENDDO
          ENDIF
        ENDIF
        IFILL=0
        GO TO 1
      ENDIF
C
C  ****  Function values at the grid points.
C
      DO I=1,N
        IF(I.GT.1.AND.I.LT.N) THEN
          IF(X(I)-X(I-1).LT.EPS*ABS(X(I))) THEN
            IF(X(I).GT.0.0D0) THEN
              XX=X(I)*(1.0D0+1.0D-14)
            ELSE
              XX=X(I)*(1.0D0-1.0D-14)
            ENDIF
          ELSE IF(X(I+1)-X(I).LT.EPS*ABS(X(I))) THEN
            IF(X(I).GT.0.0D0) THEN
              XX=X(I)*(1.0D0-1.0D-14)
            ELSE
              XX=X(I)*(1.0D0+1.0D-14)
            ENDIF
          ELSE
            XX=X(I)
          ENDIF
        ELSE
          XX=X(I)
        ENDIF
        Y(I)=FUNC(XX)
      ENDDO
C
C  ****  Adding new grid points adaptively.
C
 2    CONTINUE
      IF(N.GT.4.AND.N.LT.L+1) THEN
        CALL SPLERR(X,Y,0.0D0,0.0D0,ERR,N,0)
        IF(ERR.LT.TTOL.OR.IME.EQ.0) GO TO 6
      ELSE
        IME=0
        ERR=0.0D0
        DO 5 I=2,N-1
          IF(MIN(X(I)-X(I-1),X(I+1)-X(I)).LT.EPS*ABS(X(I))
     1      .OR.X(I+1)-X(I-1).LT.RES*ABS(X(I))) GO TO 5
          IL=MAX(1,I-L)
          DO J=I,MAX(2,I-L),-1
            IF(X(J)-X(J-1).LT.EPS*ABS(X(J))) THEN
              IL=J
              GO TO 3
            ENDIF
          ENDDO
 3        CONTINUE
C
          IU=MIN(I+L,N)
          DO J=I,MIN(I+L,N-1)
            IF(X(J+1)-X(J).LT.EPS*ABS(X(J))) THEN
              IU=J
              GO TO 4
            ENDIF
          ENDDO
 4        CONTINUE
C
          IF(I-IL.LT.LH) THEN
            IU=MIN(IU,IL+L-1)
          ELSE IF(IU-I.LT.LH) THEN
            IL=MAX(IL,IU-L+1)
          ELSE
            IL=I-LH
            IU=I+LH
          ENDIF
          IF(IU.LT.IL+2) GO TO 5
C
          NS=0
          DO J=IL,I-1
            NS=NS+1
            XS(NS)=X(J)
            YS(NS)=Y(J)
          ENDDO
          JT=NS
          DO J=I+1,IU
            NS=NS+1
            XS(NS)=X(J)
            YS(NS)=Y(J)
          ENDDO
C
          CALL SPLINE(XS,YS,AS,BS,CS,DS,0.0D0,0.0D0,NS)
          YI=AS(JT)+X(I)*(BS(JT)+X(I)*(CS(JT)+X(I)*DS(JT)))
          IF(ABS(Y(I)).GT.1.0D-60) THEN
            ERRP=ABS(YI/Y(I)-1.0D0)
          ELSE
            ERRP=ABS(YI-Y(I))
          ENDIF
          IF(ERRP.GT.ERR) THEN
            ERR=ERRP
            IME=I
          ENDIF
 5      CONTINUE
        IF(ERR.LT.TTOL.OR.IME.EQ.0) GO TO 6
      ENDIF
C
      TST=0.5D0*RES*ABS(X(IME))
      IF(X(IME+1)-X(IME-1).LT.2.01D0*TST) GO TO 6
C
      IF(X(IME+1)-X(IME).GT.TST) THEN
        XX=(X(IME)+X(IME+1))*0.5D0
        N=N+1
        DO I=N,IME+1,-1
          X(I)=X(I-1)
          Y(I)=Y(I-1)
        ENDDO
        X(IME+1)=XX
        Y(IME+1)=FUNC(XX)
      ENDIF
      IF(N.EQ.NP) GO TO 6
C
      IF(X(IME)-X(IME-1).GT.TST) THEN
        XX=(X(IME)+X(IME-1))*0.5D0
        N=N+1
        DO I=N,IME+1,-1
          X(I)=X(I-1)
          Y(I)=Y(I-1)
        ENDDO
        X(IME)=XX
        Y(IME)=FUNC(XX)
      ENDIF
      IF(N.LT.NP) GO TO 2
C
 6    CONTINUE
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SPLINT
C  *********************************************************************
      FUNCTION SPLINT(X,A,B,C,D,XL,XU,N,NPOW)
C
C     This function calculates the integral of a cubic spline function
C  SPL(X) multiplied by a power of X, i.e.,
C         SPLINT = INTEGRAL (from XL to XU) of SPL(X)*X**NPOW
C
C  Input:
C     X(1:N) ........ grid points.
C     A,B,C,D(1:N) ... spline coefficients.
C     XL, XU ........ lower and upper limits of the integral.
C     N ............. number of grid points.
C     NPOW .......... power of X in the integrand.
C
C  Output:
C     SPLINT ......... value of the integral.
C
C  Other subprograms used: subroutines FINDI and SPLIN1.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION X(N),A(N),B(N),C(N),D(N)
      DIMENSION S(4)
C
C  ****  Set integration limits in increasing order.
C
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
      SPLINT=0.0D0
C
C  ****  MIN(XL,XU) is less than X(1).
C
      IF(XLL.LT.X(1)) THEN
        A0=A(1)+D(1)*X(1)**3
        B0=B(1)-3.0D0*D(1)*X(1)**2
        C0=C(1)+3.0D0*D(1)*X(1)
C
        X1=XLL
        IF(XUU.LT.X(1)) THEN
          X2=XUU
        ELSE
          X2=X(1)
        ENDIF
C
        CALL SPLIN1(X1,X2-X1,NPOW,S)
        SPLINT=SPLINT+A0*S(1)+B0*S(2)+C0*S(3)
C
        IF(XUU.LT.X(1)) THEN
          SPLINT=SIGN*SPLINT
          RETURN
        ENDIF
        IL=1
        XLL=X(1)
      ELSE
        CALL FINDI(XLL,X,N,IL)
      ENDIF
      CALL FINDI(XUU,X,N,IU)
C  ****  Contributions from different intervals.
      DO I=IL,IU
        IF(I.EQ.IL) THEN
          X1=XLL
        ELSE
          X1=X(I)
        ENDIF
        IF(I.EQ.IU) THEN
          X2=XUU
        ELSE
          X2=X(I+1)
        ENDIF
C
        CALL SPLIN1(X1,X2-X1,NPOW,S)
        SPLINT=SPLINT+A(I)*S(1)+B(I)*S(2)+C(I)*S(3)+D(I)*S(4)
      ENDDO
C
      SPLINT=SIGN*SPLINT
      RETURN
      END
C  *********************************************************************
      SUBROUTINE SPLIN1(X,DX,N,S)
C
C  Integrals from X to X+DX of X**(N+I) with I=0, 1, 2, and 3.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (F1O3=1.0D0/3.0D0)
      DIMENSION S(4)
C
      IF(N.GT.-4.AND.N.LT.0) THEN
        IF(MIN(ABS(X),ABS(X+DX)).LT.1.0D-60) THEN
          WRITE(6,10)
 10       FORMAT(1X,'*** Error in SPLINT: negative or zero argumen',
     1      't to LOG.')
          WRITE(6,'(A,1P,E14.6)') '  X =',X
          WRITE(6,'(A,1P,E14.6)') ' DX =',DX
          STOP 'SPLINT: LOG of a small or negative number.'
        ENDIF
      ENDIF
C
      IF(ABS(X).LT.1.0D-16) THEN
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=((X+DX)**NI-X**NI)/DBLE(NI)
          ELSE
            S(I)=LOG(1.0D0+DX/X)
          ENDIF
        ENDDO
        RETURN
      ENDIF
C
      DEL=DX/X
      IF(ABS(DEL).GT.1.0D-5) THEN
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=X**NI*((1.0D0+DEL)**NI-1.0D0)/DBLE(NI)
          ELSE
            S(I)=LOG(1.0D0+DEL)
          ENDIF
        ENDDO
      ELSE
        DO I=1,4
          NI=N+I
          IF(NI.NE.0) THEN
            S(I)=X**NI*DEL*(1.0D0+DEL*0.5D0*(NI-1)
     1        *(1.0D0+DEL*F1O3*(NI-2)*(1.0D0+DEL*0.25D0*(NI-3))))
          ELSE
            S(I)=DEL*(1.0D0-DEL*(0.5D0-DEL*(F1O3-DEL*0.25D0)))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC                 Radial grid                  CCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  *********************************************************************
C                       SUBROUTINE SGRID
C  *********************************************************************
      SUBROUTINE SGRID(R,DR,RN,R2,DRN,N,NMAX,IER)
C
C     This subroutine sets up a radial grid R(I) (I=1:N) such that
C  1) A*R(I)+B*LOG(R(I)+C)+D=I.
C  2) R(1)=0, R(N)=RN, R(2)=R2 and R(N-1)=RN-DRN (approximately).
C  3) The spacing between consecutive grid points, R(I+1)-R(I),
C     increases with I and is always less than DRN.
C
C  Input arguments:
C    RN ..... outer grid point (the grid extends from 0 up to RN).
C             RN must be greater than 1.0D-5
C    R2  .... approximately =R(2) (controls the grid spacing at small
C             radii). R2 must be less than 1.0D-2 and less than RN.
C    DRN .... R(N)-R(N-1) (controls the grid spacing at large radial
C             distances).
C    N ...... tentative number of grid points (it may be increased to
C             meet conditions 2 and 3).
C    NMAX ... physical dimensions of the arrays R(.) and DR(.); N cannot
C             exceed NMAX.
C
C  Output arguments:
C    N ...... number of grid points. It may be greater than the input
C             value.
C    R(1:N) ... radial grid points.
C    DR(1:N) ... values of the derivative of R with respect to I, which
C             are required to evaluate integrals using, e.g., Simpson's
C             or Lagrange's quadrature formulas.
C    IER .... error flag:
C             IER=0, the grid has been successfully defined.
C             IER>0, the grid could not be defined.
C
C     To describe the radial wave functions of bound electrons of an
C  atom or positive ion in its ground state configuration, the following
C  values of the input parameters should be adequate: RN of the order of
C  50, R2 about 1.0D-5 or smaller, DRN about 0.5, N=750 or larger.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(NMAX),DR(NMAX)
      IER=0
C
 1000 FORMAT(' RN=',1P,E13.6,', R2=',E13.6,', DRN=',E13.6,', N=',I6)
      RRN=RN
      IF(RRN.LT.1.0D-5) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: RN is .LT. 1.0E-5.'')')
        WRITE(6,1000) RRN,R2,DRN,N
        IER=1
        RETURN
      ENDIF
C
      RR2=R2
      IF(RR2.LT.1.0D-8) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): R2 is .LT. 1.0E-8.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        RR2=1.0D-8
        WRITE(6,1000) RRN,RR2,DRN,N
      ENDIF
C
      IF(RR2.GE.1.0D-2*RRN) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): R2 is .GE. 1.0E-2*RN.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        RR2=1.0D-2
        WRITE(6,1000) RRN,RR2,DRN,N
      ENDIF
C
      IF(RR2.GE.RRN) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: R2 is .GE. RN.'')')
        WRITE(6,1000) RRN,RR2,DRN,N
        IER=2
        RETURN
      ENDIF
C
      RDRN=DRN
      IF(RDRN.LE.RR2) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): DRN is .LE. R2.'')')
        WRITE(6,1000) RRN,RR2,RDRN,N
        RDRN=MAX(2.0D0*(RRN-RR2)/MAX(N,10),RR2)
        WRITE(6,1000) RRN,RR2,RDRN,N
      ENDIF
C
      NR=MAX(N,10)
      TST=5.0D0*(RRN-RR2)/NR
      IF(RDRN.GT.TST) RDRN=TST
      RLOW=(RRN/RDRN)+10.0D0
      IF(RLOW.GT.DBLE(NMAX)) THEN
        NLOW=NMAX
      ELSE
        NLOW=RLOW
      ENDIF
      IF(NR.LT.NLOW) THEN
        WRITE(6,'(/1X,''*** Warning (SGRID): NR is .LT. NLOW.'')')
        WRITE(6,'('' NLOW ='',I8)') NLOW
        WRITE(6,1000) RRN,RR2,RDRN,NR
        NR=NLOW
        WRITE(6,1000) RRN,RR2,RDRN,NR
      ENDIF
C
      IF(NR.GT.NMAX) THEN
        WRITE(6,'(/1X,''*** Error in SGRID: NR is .GT. NMAX.'')')
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' NMAX ='',I8)') NMAX
        IER=3
        RETURN
      ENDIF
C
      HIGH=0.5D0*((RRN/RR2)+(RRN/RDRN))-20.0D0
      IF(HIGH.GT.DBLE(NMAX)) THEN
        NHIGH=NMAX
      ELSE
        NHIGH=HIGH
      ENDIF
      IF(NR.GT.NHIGH) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,'(/1X,''*** Error in SGRID: NR is .GT. NHIGH.'')')
        WRITE(6,'('' NHIGH ='',I8)') NHIGH
        WRITE(6,1000) RRN,RR2,RDRN,NR
        IER=4
        RETURN
      ENDIF
C
C  ****  Grid parameters. (ME-8.6 to ME-8.10)
C
      AG=(RRN-(NR-1)*RR2)*RDRN/(RRN*(RDRN-RR2))
      IF(AG.GT.1.0D0.OR.AG.LT.0.5002D0) THEN
        A=(NR-1)/RRN
        AA=1.0D0/A
        IF(AA.LT.RDRN.AND.AA.LT.RR2) THEN
          B=0.0D0  ! Linear grid.
          C=1.0D0
          D=0.0D0
          GO TO 3
        ENDIF
        WRITE(6,'(/1X,''*** Error in SGRID: AG is out of range.'')')
        WRITE(6,1000) RRN,RR2,RDRN,NR
        WRITE(6,'('' AG ='',1P,E13.6)') AG
        IER=5
        RETURN
      ENDIF
C
      XL=0.0D0
      XU=1000.0D0
 1    X=0.5D0*(XL+XU)
      F=(1.0D0+X)*(1.0D0-X*LOG((1.0D0+X)/X))
      IF(F.GT.AG) THEN
        XL=X
      ELSE
        XU=X
      ENDIF
      IF(XU-XL.GT.1.0D-15) GO TO 1
C
      C=X*RRN
      B=X*(C+RRN)*(RDRN-RR2)/(RDRN*RR2)
      A=(C-B*RR2)/(C*RR2)
      D=1.0D0-B*LOG(C)
C
      R(1)=0.0D0
      DR(1)=(R(1)+C)/(A*(R(1)+C)+B)
      RR=1.0D-35
      DO I=2,NR
        RL=RR
        RU=RRN
 2      RR=0.5D0*(RU+RL)
        FR=A*RR+B*LOG(RR+C)+D-DBLE(I)
        IF(FR.GT.0.0D0) THEN
          RU=RR
        ELSE
          RL=RR
        ENDIF
        IF(RU-RL.GT.1.0D-15*RR) GO TO 2
        R(I)=RR
        DR(I)=(RR+C)/(A*(RR+C)+B)
        IF(DR(I).LT.DR(I-1)) THEN
C  **** The grid spacing does not increase with I.
          WRITE(6,'(/1X,''*** Error in SGRID: non-increasing grid '',
     1      ''spacing.'')')
          WRITE(6,1000) RRN,RR2,RDRN,NR
          IER=6
          RETURN
        ENDIF
      ENDDO
      N=NR
      R(N)=RRN
      RETURN
C
 3    CONTINUE
      R(1)=0.0D0
      DR(1)=AA
      DO I=2,NR
        R(I)=(I-1)*AA
        DR(I)=AA
      ENDDO
      N=NR
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SLAG6
C  *********************************************************************
      SUBROUTINE SLAG6(H,Y,S,N)
C
C     Piecewise six-point Lagrange integration of a uniformly tabulated
C  function. (ME-8.13)
C
C  Input arguments:
C     H ............ grid spacing.
C     Y(I) (1:N) ... array of function values (ordered abscissas).
C     N ............ number of data points.
C
C  Output argument:
C     S(I) (1:N) ... array of integral values defined as
C                S(I)=INTEGRAL(Y) from X(1) to X(I)=X(1)+(I-1)*H.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION Y(N),S(N)
      IF(N.LT.6) STOP 'SLAG6: too few data points.'
      HR=H/1440.0D0
      Y1=0.0D0
      Y2=Y(1)
      Y3=Y(2)
      Y4=Y(3)
      S(1)=0.0D0
      S(2)=HR*(475.0D0*Y2+1427.0D0*Y3-798.0D0*Y4+482.0D0*Y(4)
     1   -173.0D0*Y(5)+27.0D0*Y(6))
      S(3)=S(2)+HR*(-27.0D0*Y2+637.0D0*Y3+1022.0D0*Y4-258.0D0*Y(4)
     1   +77.0D0*Y(5)-11.0D0*Y(6))
      DO I=4,N-2
        Y1=Y2
        Y2=Y3
        Y3=Y4
        Y4=Y(I)
        S(I)=S(I-1)+HR*(11.0D0*(Y1+Y(I+2))-93.0D0*(Y2+Y(I+1))
     1    +802.0D0*(Y3+Y4))
      ENDDO
      Y5=Y(N-1)
      Y6=Y(N)
      S(N-1)=S(N-2)+HR*(-27.0D0*Y6+637.0D0*Y5+1022.0D0*Y4-258.0D0*Y3
     1  +77.0D0*Y2-11.0D0*Y1)
      S(N)=S(N-1)+HR*(475.0D0*Y6+1427.0D0*Y5-798.0D0*Y4+482.0D0*Y3
     1  -173.0D0*Y2+27.0D0*Y1)
      RETURN
      END
C  NOTE: The present subroutine package uses I/O units 33, 98 and 99.
C        Do not use these unit numbers in your main program.
C
C  *********************************************************************
C                       SUBROUTINE ELSEPA
C  *********************************************************************
      SUBROUTINE ELSEPA(IELEC,EV,IZ,NELEC,MNUCL,MELEC,MUFIN,RMUF,VMOL,
     1  MEXCH,MCPOL,VPOLA,VPOLB,MABS,VABSA,VABSD,IHEF,IW)
C
C
C                       F. Salvat, D. Bote, A. Jablonski and C.J. Powell
C                              September 25, 2008
C
C  Updated in September 2020 by F. Salvat, to run with the subroutine
C  package RADIAL.
C
C  Ref.: F. Salvat and J. M. Fernandez-Varea,
C        'RADIAL: a Fortran subroutine package for the solution of
C        radial Schrodinger and Dirac wave equations',
C        Comput. Phys. Commun. 240 (2019) 165-177.
C
C     This subroutine computes scattering amplitudes, differential cross
C  sections and total (integrated) cross sections for ELastic Scattering
C  of Electrons and Positrons by neutral Atoms and positive ions.
C
C     The interaction is described through a static (central) field,
C  which consists of the electrostatic potential and, for projectile
C  electrons, an approximate local exchange potential. For slow
C  projectiles, a correlation-polarization potential and an absorptive
C  imaginary potential can optionally be included. The differential
C  cross section is evaluated by means of relativistic (Dirac) partial-
C  wave analysis, or from approximate high-energy factorizations.
C
C  Input arguments:
C    IELEC ..... electron-positron flag;
C                =-1 for electrons,
C                =+1 for positrons.
C    EV ........ projectile's kinetic energy (in eV).
C    IZ ........ atomic number of the target atom or ion.
C    NELEC ..... number of bound atomic electrons.
C    MNUCL ..... nuclear charge density model.
C                  1 --> point nucleus (P),
C                  2 --> uniform distribution (U),
C                  3 --> Fermi distribution (F),
C                  4 --> Helm's uniform-uniform distribution (Uu).
C    MELEC ..... electron density model.
C                  1 --> TFM analytical density,
C                  2 --> TFD analytical density,
C                  3 --> DHFS analytical density,
C                  4 --> DF numerical density, read from 'Z_zzz.DEN',
C                  5 --> density read from file 'density.usr'.
C    MUFIN ..... Aggregation effects...
C                  0 --> free atom,
C                  1 --> muffin-tin model.
C      RMUF .... Muffin-tin radius (in cm).
C      VMOL .... Number of atoms per unit volume (in 1/cm**3).
C                    Reciprocal of the Wigner-Seitz cell volume.
C    MEXCH ..... exchange correction for electrons.
C                  0 --> no exchange correction,
C                  1 --> Furness-McCarthy (FM),
C                  2 --> Thomas-Fermi (TF),
C                  3 --> Riley-Truhlar (RT).
C    MCPOL ..... correlation-polarization correction.
C                  0 --> no correlation-polarization correction,
C                  1 --> Buckingham potential (B),
C                  2 --> Local density approximation (LDA).
C      VPOLA ... atomic polarizability (in cm**3).
C      VPOLB ... cutoff radius parameter b_pol
C                    (used only when MCPOL>0).
C    MABS ...... absorption correction (imaginary potential).
C                  0 --> no absorption correction,
C                  1 --> LDA-I (electron-hole excitations only).
C                  2 --> LDA-II (full Lindhard dielectric function).
C      VABSA ... strength of the absorption potential.
C      VABSD ... energy gap, DELTA (eV).
C                    (used only when MABS is different from 0).
C    IHEF ...... =0: phase shifts are computed for the electrostatic
C                    field of the whole atom (nucleus+electron cloud)
C                    with optional exchange, polarization and absorption
C                    corrections.
C                =1: the differential cross section is obtained from a
C                    high-energy factorization. The phase shifts are
C                    evaluated for the bare nucleus. The screening of
C                    the nuclear charge by the atomic electrons is
C                    accounted for by means of a pre-evaluated high-
C                    energy correction factor, which is read from file
C                    'Z_zzz.DFS'.
C                =2: when the energy is larger than 100 MeV, the DCS is
C                    obtained as the product of the Mott DCS for a point
C                    nucleus, the Helm Uu nuclear form factor (with an
C                    empirical Coulomb correction) and the electron
C                    screening factor.
C    IW ........ output unit (to be defined in the main program).
C
C  The electrostatic potential and electron density of the target atom
C  or ion are calculated by subroutine EFIELD and delivered through the
C  the named common block /CFIELD/.
C
C  Output (through the common block /DCSTAB/):
C     ECS ........ total cross section (cm**2).
C     TCS1 ....... 1st transport cross section (cm**2).
C     TCS2 ....... 2nd transport cross section (cm**2).
C     TH(I) ...... scattering angles (in deg)
C     XT(I) ...... values of (1-COS(TH(I)))/2.
C     DCST(I) .... differential cross section per unit solid
C                  angle at TH(I) (in cm**2/sr).
C     SPOL(I) .... Sherman spin-polarization function at TH(I).
C     ERROR(I) ... relative uncertainty of the computed DCS
C                  values. Estimated from the convergence of the
C                  series.
C     NTAB ....... number of angles in the table.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
C
C  ****  The parameter IWR defines the amount of information printed on
C  output files:
C  IWR>0 => the scattering potential is printed on file 'scfield.dat'.
C  IWR>1 => the scattering amplitudes are printed on file 'scatamp.dat'.
      PARAMETER (IWR=2)
C
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (TREV=REV+REV)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C
      CHARACTER*120 SCFILE,FILE1,CS120,NULL
      CHARACTER*1 LIT10(10),LIT1,LIT2,LIT3
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
C
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
      DIMENSION Q2T(NGT),FQ(NGT)
C  ****  Link with the RADIAL package.
      COMMON/RADWF/RRR(NDIM),P(NDIM),Q(NDIM),NRT,ILAST,IER
C
      COMMON/CFIELD/R(NDIM),RVN(NDIM),DEN(NDIM),RVST(NDIM),NPOT
      COMMON/FIELD/RAD(NDIM),RV(NDIM),NP
      COMMON/FIELDI/RADI(NDIM),RVI(NDIM),RW(NDIM),IAB,NPI
      DIMENSION RVEX(NDIM),RVPOL(NDIM)
C
      PARAMETER (NPC=1500,NDM=25000)
      DIMENSION SA(NDM),SB(NDM),SC(NDM),SD(NDM)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CFX,CGX,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
      DIMENSION DENA(NPC),DENB(NPC),DENC(NPC),DEND(NPC)
C
C  ****  Path to the ELSEPA database
      CHARACTER*100 PATHE
      PATHE='./database/'
C
C  ----  Mott DCS and spin polarization (point unscreened nucleus)
      IF(NELEC.EQ.0.AND.MNUCL.EQ.1) THEN
        CALL MOTTSC(IELEC,IZ,EV,IW)
        RETURN
      ENDIF
C  ----  High-energy Mott-Born approximation for neutral atoms.
      IF(EV.GT.100.0D6.AND.IHEF.EQ.2.AND.IZ.EQ.NELEC) THEN
        CALL HEBORN(IELEC,IZ,MNUCL,EV,IW)
        RETURN
      ENDIF
C
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine ELSEPA. Elastic scattering of ',
     1  'electrons and positrons',/1X,'#',20X,
     2  'by neutral atoms and positive ions')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',
     1       E12.5,' a.u.')
      IF(MUFIN.EQ.1.AND.NELEC.NE.IZ) THEN
        WRITE(IW,*) '   IZ =',IZ
        WRITE(IW,*) 'NELEC =',NELEC
        WRITE(IW,'(''  STOP. For muffin-tin atoms, NELEC '',
     1    ''must equal IZ.'')')
        WRITE(6,*) '   IZ =',IZ
        WRITE(6,*) 'NELEC =',NELEC
        STOP 'ELSEPA: For muffin-tin atoms, NELEC must equal IZ.'
      ENDIF
C
C  ****  You may wish to comment off the next condition to run the
C  program for kinetic energies less that 5 eV. However, the results
C  for these energies may be highly inaccurate.
C
      IF(EV.LT.4.9990D0) THEN
        WRITE(IW,'(''  STOP. The kinetic energy is too small.'')')
        STOP 'ELSEPA: The kinetic energy is too small.'
      ENDIF
      IF(EV.LT.100.0D0) THEN
        WRITE(IW,'(1X,''#'',/1X,''#  ***  WARNING: Energy is '',
     1   ''too low.'',/1X,''#'',16X,''The reliability of the '',
     2   ''results is questionable.'')')
      ENDIF
      ERE0=0.0D0
      EIM=0.0D0
      VMOL1=-1.0D-23  ! Serves only to prevent compiler warnings.
C
C  ****  Electrostatic field.
C
      IF(MNUCL.LT.1.OR.MNUCL.GT.4) THEN
        WRITE(IW,*) 'ELSEPA: incorrect MNUCL value.'
        STOP 'ELSEPA: incorrect MNUCL value.'
      ENDIF
      IF(MELEC.LT.1.OR.MELEC.GT.5) THEN
        WRITE(IW,*) 'ELSEPA: incorrect MELEC value.'
        STOP 'ELSEPA: incorrect MELEC value.'
      ENDIF
      CALL EFIELD(IZ,NELEC,MNUCL,MELEC,IW,0)
C
      IHEF0=0
      IAB=0
      IF(NELEC.EQ.0) THEN  ! Bare nuclei treated separately.
        IHEF0=2
        GO TO 100
      ENDIF
      IF(NELEC.EQ.IZ.AND.MABS.EQ.0.AND.EV.GT.20.1D3*IZ) THEN
        IF(IHEF.GT.0) THEN  ! Neutral atoms, high-energy factorization.
          IHEF0=1
          GO TO 100
        ENDIF
      ENDIF
C  ****  Caution: Partially stripped atoms are difficult to compute...
      ZINF=IELEC*DBLE(IZ-NELEC)
      DO I=1,NPOT
        RV(I)=DBLE(IELEC)*RVST(I)
      ENDDO
      RV(NPOT)=ZINF
C
C  ************  Muffin-tin model for scattering in solids.
C
      NMT=0
      IF(MUFIN.EQ.1) THEN
        WRITE(IW,1600) RMUF
 1600   FORMAT(1X,'#',/1X,'# Muffin-tin model:    Rmt =',1P,E12.5,' cm')
        IF(ABS(ZINF).GT.1.0D-6) THEN
          WRITE(IW,*) 'ELSEPA: muffin-tin model does not work for ions.'
          STOP 'ELSEPA: muffin-tin model does not work for ions.'
        ENDIF
        CALL SPLINE(R,RV,SA,SB,SC,SD,0.0D0,0.0D0,NPOT)
        CALL SPLINE(R,DEN,DENA,DENB,DENC,DEND,0.0D0,0.0D0,NPOT)
        RMT=RMUF/A0B
        IF(RMT.LT.R(NPOT)) THEN
          CALL FINDI(RMT,R,NPOT,J)
          IF(J.LT.5) STOP 'ELSEPA: The muffin-tin radius is too small.'
          DENRMT=DENA(J)+RMT*(DENB(J)+RMT*(DENC(J)+RMT*DEND(J)))
        ELSE
          RMT=R(NPOT)
          DENRMT=DEN(NPOT)
        ENDIF
C
        IF(VMOL.GT.0.0D0) THEN
          VMOL1=VMOL*A0B**3
          WRITE(IW,1601) VMOL
        ELSE
          VMOL1=1.0D0/(2.0D0*RMT)**3
          WRITE(IW,1601) VMOL1/A0B**3
        ENDIF
 1601   FORMAT(1X,'#   Atomic density:   Vmol =',1P,E12.5,' 1/cm**3')
        DO I=1,NPOT
          IF(R(I).GT.RMT) THEN
            IF(RAD(I-1).LT.RMT*0.9999999D0) THEN
              NP=I
              RAD(NP)=RMT
            ELSE
              NP=I-1
              RAD(NP)=RMT
            ENDIF
            RC1=RMT
            CALL FINDI(RC1,R,NPOT,J)
            V1=SA(J)+RC1*(SB(J)+RC1*(SC(J)+RC1*SD(J)))
            DEN1=DENA(J)+RC1*(DENB(J)+RC1*(DENC(J)+RC1*DEND(J)))
            RV(NP)=2.0D0*V1
            DEN(NP)=2.0D0*DEN1
            RVST(NP)=RV(NP)/DBLE(IELEC)
            GO TO 1
          ELSE
            RAD(I)=R(I)
C
            RC1=R(I)
            FD1=FOURPI*RC1**2
            V1=RV(I)
            DEN1=DEN(I)
C
            RC2=2.0D0*RMT-R(I)
            FD2=FOURPI*RC2**2
            CALL FINDI(RC2,R,NPOT,J)
            V2=SA(J)+RC2*(SB(J)+RC2*(SC(J)+RC2*SD(J)))
            DEN2=DENA(J)+RC2*(DENB(J)+RC2*(DENC(J)+RC2*DEND(J)))
C
            IF(I.GT.1) THEN
              RV(I)=V1+RC1*(V2/RC2)
              DEN(I)=DEN1+FD1*(DEN2/FD2)
            ELSE
              RV(I)=V1
              DEN(I)=DEN1
            ENDIF
          ENDIF
          RVST(I)=RV(I)/DBLE(IELEC)
        ENDDO
        NP=NPOT
 1      CONTINUE
C
C  ****  Ensure proper normalization of the muffin-tin electron density.
C
        SUM=SMOMLL(RAD,DEN,RAD(1),RAD(NP),NP,0,0)
        RHOU=(NELEC-SUM)/(FOURPI*RMT**3/3.0D0)
        DO I=1,NP
          DEN(I)=DEN(I)+RHOU*FOURPI*RAD(I)**2
        ENDDO
        SUM=SMOMLL(RAD,DEN,RAD(1),RAD(NP),NP,0,0)
        WRITE(6,*) 'Electron density normalization =',SUM
C
        NMT=NP
        ERE0=-RV(NMT)/RAD(NMT)
        E=E-RV(NMT)/RAD(NMT)
        DO I=2,NMT
          RV(I)=RV(I)-RV(NMT)*RAD(I)/RAD(NMT)
          RVST(I)=RVST(I)-RVST(NMT)*RAD(I)/RAD(NMT)
        ENDDO
        NPP=NP+10
        IF(NP.LT.NPP) THEN
          NP=NP+1
          DO I=NP,NPP
            IF(I.EQ.NP) THEN
              RAD(I)=RAD(I-1)
            ELSE
              RAD(I)=RAD(I-1)+0.01D0*RMT
            ENDIF
            RV(I)=0.0D0
            RVST(I)=0.0D0
            DEN(I)=0.0D0
          ENDDO
          NP=NPP
        ENDIF
      ELSE
        NP=NPOT
        DO I=1,NPOT
          RAD(I)=R(I)
        ENDDO
      ENDIF
C
C  ************  Exchange correction for electrons.
C
      IF(IELEC.EQ.-1.AND.MEXCH.NE.0) THEN
        IF(MEXCH.EQ.1) THEN
C  ****  Furness-McCarthy exchange potential.
          WRITE(IW,1500)
 1500     FORMAT(1X,'#',/1X,'# Furness-McCarthy exchange',
     1      ' potential')
          DO I=2,NP
            AUX=RAD(I)*E*(1.0D0+EV/TREV)+RVST(I)
            AUX2=AUX*AUX
            IF(DEN(I).GT.1.0D-5*AUX2) THEN
              RVEX(I)=0.5D0*(AUX-SQRT(AUX2+DEN(I)))
            ELSE
              T=DEN(I)/AUX2
              RVEX(I)=-0.5D0*AUX*T*(0.5D0-T*(0.125D0-T*0.065D0))
            ENDIF
            RV(I)=RV(I)+RVEX(I)
          ENDDO
        ELSE IF(MEXCH.EQ.2) THEN
C  ****  Thomas-Fermi exchange potential.
          WRITE(IW,1400)
 1400     FORMAT(1X,'#',/1X,'# Thomas-Fermi exchange potential')
          DO I=1,NP
            RHO=DEN(MAX(2,I))/(FOURPI*RAD(MAX(2,I))**2)
            SKF=(3.0D0*PI*PI*RHO)**3.333333333333333D-1
            EF=0.5D0*SKF*SKF
            SKL=SQRT(2.0D0*(E*(1.0D0+EV/TREV)+EF))
            X=SKF/SKL
            IF(X.LT.0.001D0) THEN
              FX=(2.0D0/3.0D0)*X**3
            ELSE
              FX=X-0.5D0*(1.0D0-X*X)*LOG(ABS((1.0D0+X)/(1.0D0-X)))
            ENDIF
            RVEX(I)=-(SKL/PI)*FX*RAD(I)
            RV(I)=RV(I)+RVEX(I)
          ENDDO
        ELSE IF(MEXCH.EQ.3) THEN
C  ****  Riley-Truhlar exchange potential.
          WRITE(IW,1501)
 1501     FORMAT(1X,'#',/1X,'# Riley-Truhlar exchange potential')
          DO I=1,NP
            AUX=4.0D0*(RAD(I)*E*(1.0D0+EV/TREV)+RVST(I))
            IF(AUX.GT.1.0D-16*DEN(I)) THEN
              RVEX(I)=-DEN(I)/AUX
              RV(I)=RV(I)+RVEX(I)
            ENDIF
          ENDDO
        ELSE
          WRITE(IW,*) 'ELSEPA: incorrect MEXCH value.'
          STOP 'ELSEPA: incorrect MEXCH value.'
        ENDIF
        IF(NMT.GT.1) THEN
          ERE0=ERE0-RV(NMT)/RAD(NMT)
          DO I=2,NMT
            RV(I)=RV(I)-RAD(I)*(RV(NMT)/RAD(NMT))
            RVEX(I)=RVEX(I)-RAD(I)*(RVEX(NMT)/RAD(NMT))
          ENDDO
        ENDIF
      ELSE
        IF(IELEC.EQ.-1) WRITE(IW,1511)
 1511   FORMAT(1X,'#',/1X,'# No exchange potential')
        DO I=1,NP
          RVEX(I)=0.0D0
        ENDDO
      ENDIF
C
C  ********  Absorption potential.
C
      IF(MABS.EQ.1.AND.VABSA.GT.1.0D-12.AND.EV.LT.1.001D6) THEN
C
C  ****  LDA-I model
C
        IAB=1
        WRITE(IW,1502) VABSD,VABSA
 1502   FORMAT(1X,'#',/1X,'# LDA-I absorption potential (only electr',
     1    'on-hole excitations):',/1X,'#',
     2    27X,'Delta =',1P,E12.5,' eV',/1X,'#',28X,'Aabs =',E12.5)
        DELTA=VABSD/HREV
        AABS=VABSA
C
        RW(1)=0.0D0
        RVPOL(1)=0.0D0
        DO I=2,NP
          RHO=DEN(I)/(FOURPI*RAD(I)**2)
C  ****  Local kinetic energy.
          EKIN=E-RV(I)/RAD(I)
          IF(RHO.GT.1.0D-16.AND.EKIN.GT.DELTA) THEN
            VEL=SQRT(2.0D0*EKIN)
C  ****  Relativistic correction.
            EKEV=EKIN*HREV
            FREL=SQRT(2.0D0*(EKEV+REV)**2/(REV*(EKEV+2.0D0*REV)))
C  ****  Only electron-hole excitations.
            CALL XSFEG(RHO,DELTA,IELEC,EKIN,0,XSEC,2)
            RW(I)=-0.5D0*VEL*RHO*XSEC*RAD(I)*AABS*FREL
          ELSE
            RW(I)=0.0D0
          ENDIF
          RVPOL(I)=0.0D0
          WRITE(6,1503) I,RAD(I),RW(I)
 1503     FORMAT(1X,'i=',I4,',   r=',1P,E12.5,',   r*Wabs= ',E12.5)
        ENDDO
      ELSE IF(MABS.EQ.2.AND.VABSA.GT.1.0D-12.AND.EV.LT.1.001D6) THEN
C
C  ****  LDA-II model.
C
        IAB=1
        WRITE(IW,1504) VABSD,VABSA
 1504   FORMAT(1X,'#',/1X,'# LDA-II absorption potential (plasmon and',
     1    ' e-h excitations):',/1X,'#',
     2    27X,'Delta =',1P,E12.5,' eV',/1X,'#',28X,'Aabs =',E12.5)
        DELTA=VABSD/HREV
        AABS=VABSA
C
        RW(1)=0.0D0
        RVPOL(1)=0.0D0
        DO I=2,NP
          RHO=DEN(I)/(FOURPI*RAD(I)**2)
          IF(RHO.GT.1.0D-16) THEN
            EFERMI=0.5D0*(3.0D0*PI**2*RHO)**0.666666666666666D0
C  ****  Local kinetic energy.
            IF(IELEC.EQ.-1) THEN
              EKIN=E+EFERMI
            ELSE
              EKIN=E-EFERMI
            ENDIF
            IF(EKIN.GT.DELTA) THEN
              VEL=SQRT(2.0D0*EKIN)
C  ****  Relativistic correction.
              EKEV=EKIN*HREV
              FREL=SQRT(2.0D0*(EKEV+REV)**2/(REV*(EKEV+2.0D0*REV)))
C  ****  Plasmon and electron-hole excitations.
              CALL XSFEG(RHO,DELTA,IELEC,EKIN,0,XSEC,1)
              RW(I)=-0.5D0*VEL*RHO*XSEC*RAD(I)*AABS*FREL
            ELSE
              RW(I)=0.0D0
            ENDIF
          ELSE
            RW(I)=0.0D0
          ENDIF
          RVPOL(I)=0.0D0
          WRITE(6,1503) I,RAD(I),RW(I)
        ENDDO
      ELSE
        WRITE(IW,1514)
 1514   FORMAT(1X,'#',/1X,'# No absorption potential')
        DO I=1,NP
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
        ENDDO
      ENDIF
C
C  ****  We add a 'constant' tail to the potential to extend the grid
C  up to a point where irregular Coulomb functions can be be calculated.
C
      IF((MCPOL.NE.1.AND.MCPOL.NE.2).OR.EV.GT.1.0D4) THEN
        WRITE(IW,1710)
 1710   FORMAT(1X,'#',/1X,'# No correlation-polarization potential')
        IF(NP.LT.NDIM-10) THEN
          IF(RAD(NP)-RAD(NP-1).LT.1.0D-16) THEN
            I=NP
            RAD(I)=RAD(I-1)
          ELSE
            I=NP+1
            RAD(I)=RAD(I-1)
          ENDIF
          RV(I)=ZINF
          RVST(I)=ZINF/DBLE(IELEC)
          DEN(I)=0.0D0
          RVEX(I)=0.0D0
          RVPOL(I)=0.0D0
          RW(I)=0.0D0
          IST=I+1
          NADD=1
          DO I=IST,NDIM
            NADD=NADD+1
            RAD(I)=2.0D0*RAD(I-1)
            RV(I)=ZINF
            RVST(I)=ZINF/DBLE(IELEC)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RVPOL(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            IF(RAD(I).GT.1.0D4.AND.NADD.GT.4) GO TO 2
          ENDDO
        ELSE
          STOP 'ELSEPA: Not enough memory space 1.'
        ENDIF
 2      CONTINUE
      ELSE IF(MCPOL.EQ.1) THEN
C
C  ************  Atomic polarizability correction.
C
C  ****  Buckingham empirical potential.
        WRITE(IW,1700) VPOLA,VPOLB
 1700   FORMAT(1X,'#',/1X,'# Correlation-polarization potential (Buc',
     1    'kingham):',/1X,'#',27X,'Alpha =',1P,E12.5,' cm**3',
     2     /1X,'#',28X,'Bpol =',E12.5)
        IF(VPOLB.LT.0.01D0) THEN
          WRITE(IW,*) 'ELSEPA: VPOLB cannot be less than 0.01.'
          STOP 'ELSEPA: VPOLB cannot be less than 0.01.'
        ENDIF
        ALPHA=VPOLA/A0B**3
        D2=SQRT(0.5D0*ALPHA*VPOLB**2/DBLE(IZ)**3.333333333333333D-1)
        NPOL=NP
        DO I=1,NPOL
          VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
          RVPOL(I)=VPOL*RAD(I)
          RV(I)=RV(I)+RVPOL(I)
        ENDDO
        IF(NPOL.LT.NDIM-10) THEN
          DO I=NPOL+1,NDIM-10
            RAD(I)=1.25D0*RAD(I-1)
            VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
            RVPOL(I)=VPOL*RAD(I)
            RVST(I)=ZINF
            RV(I)=ZINF+RVPOL(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            IF(ABS(VPOL).LT.1.0D-8*MAX(E,1.0D1*ABS(ZINF)/RAD(I))
     1        .AND.RAD(I).GT.50.0D0) GO TO 3
          ENDDO
        ENDIF
        STOP 'ELSEPA: Not enough memory space 2.'
 3      CONTINUE
        IF(NP.LT.NDIM-10) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RVST(I)=ZINF
          RV(I)=ZINF
          DEN(I)=0.0D0
          RVEX(I)=0.0D0
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
          NDIN=NP+10
          DO I=NP+2,NDIN
            RAD(I)=2.0D0*RAD(I-1)
            RVST(I)=ZINF
            RV(I)=ZINF
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            RVPOL(I)=0.0D0
            NP=I
            IF(RAD(I).GT.1.0D4) GO TO 33
          ENDDO
        ELSE
          STOP 'ELSEPA: Not enough memory space 3.'
        ENDIF
 33     CONTINUE
C
      ELSE IF(MCPOL.EQ.2) THEN
C  ****  LDA correlation-polarization potential.
        IF(MUFIN.NE.1) THEN
          WRITE(IW,1701) VPOLA,VPOLB
 1701     FORMAT(1X,'#',/1X,'# Correlation-polarization potential (LDA',
     1      '):',/1X,'#',27X,'Alpha =',1P,E12.5,' cm**3',
     2      /1X,'#',28X,'Bpol =',E12.5)
          IF(VPOLB.LT.0.01D0) THEN
            WRITE(IW,*) 'ELSEPA: VPOLB cannot be less than 0.01.'
            STOP 'ELSEPA: VPOLB cannot be less than 0.01.'
          ENDIF
        ELSE
          WRITE(IW,1711)
 1711     FORMAT(1X,'#',/1X,'# Correlation potential (muffin-tin',
     1      ' model): LDA or Lindhard')
        ENDIF
        NPOL=NP
        IMODE=0
        IF(MUFIN.NE.1) THEN
          ALPHA=VPOLA/A0B**3
          D2=SQRT(0.5D0*ALPHA*VPOLB**2/DBLE(IZ)**3.333333333333333D-1)
          DO I=NPOL,1,-1
            RIP=RAD(MAX(2,I))
            RHO=DEN(MAX(2,I))/(FOURPI*RIP**2)
            VCO=VCPOL(IELEC,RHO)
C
            VPAS=-0.5D0*ALPHA/(RIP**2+D2)**2
            IF(IMODE.EQ.0) THEN
              VPOL=VPAS
              IF(VCO.LT.VPAS) IMODE=1
            ELSE
              VPOL=MAX(VCO,VPAS)
            ENDIF
            RVPOL(I)=VPOL*RAD(I)
            RV(I)=RV(I)+RVPOL(I)
          ENDDO
          IF(IMODE.EQ.0) THEN
            WRITE(IW,1702)
            WRITE(6,1702)
 1702       FORMAT(1X,'#',/1X,'# ERROR: The correlation and pol',
     1        'arization potentials do not cross.')
            STOP 'ELSEPA: V_corr and V_pol do not cross.'
          ENDIF
          VCOUT=-1.0D16
        ELSE
          VCOUT=0.0D0  ! Serves only to prevent compiler warnings.
          DO I=1,NPOL
            RIP=RAD(MAX(2,I))
            RHO=DEN(MAX(2,I))/(FOURPI*RIP**2)
            IF(I.LE.NMT) THEN
              VCO=VCPOL(IELEC,RHO)
C  ****  Lindhard high-E correlation potential.
              EKIN=E-DBLE(IELEC)*RVST(I)/RIP
              IF(EKIN.GT.1.0D-12) THEN
                VCOL=-SQRT((PI/2.0D0)**3*RHO/EKIN)
              ELSE
                VCOL=-1.0D35
              ENDIF
              VCO=MAX(VCO,VCOL)
              IF(I.EQ.NMT) VCOUT=VCO
            ELSE
              VCO=VCOUT
            ENDIF
            RVPOL(I)=VCO*RAD(I)
            RV(I)=RV(I)+RVPOL(I)
          ENDDO
          ERE0=ERE0-VCOUT
          DO I=2,NPOL
            RV(I)=RV(I)-RAD(I)*VCOUT
            RVPOL(I)=RVPOL(I)-RAD(I)*VCOUT
          ENDDO
          GO TO 34
        ENDIF
C
        IF(NPOL.LT.NDIM-10) THEN
          DO I=NPOL+1,NDIM-10
            IF(IMODE.EQ.0) THEN
              RAD(I)=RAD(I-1)+0.05D0
            ELSE
              RAD(I)=1.25D0*RAD(I-1)
            ENDIF
            VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
            IF(IMODE.EQ.0) THEN
              IF(VPOL.GT.VCOUT) IMODE=1
              VPOL=MAX(VCOUT,VPOL)
            ENDIF
            RVPOL(I)=VPOL*RAD(I)
            RVST(I)=ZINF
            RV(I)=ZINF+RVPOL(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            IF(ABS(VPOL).LT.1.0D-8*MAX(E,1.0D1*ABS(ZINF)/RAD(I))
     1        .AND.RAD(I).GT.50.0D0) GO TO 34
          ENDDO
        ENDIF
        STOP 'ELSEPA: Not enough memory space 4.'
 34     CONTINUE
        IF(NP.LT.NDIM-10) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RVST(I)=ZINF
          RV(I)=ZINF
          DEN(I)=0.0D0
          RVEX(I)=0.0D0
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
          NDIN=NP+10
          DO I=NP+2,NDIN
            RAD(I)=2.0D0*RAD(I-1)
            RVST(I)=ZINF
            RV(I)=ZINF
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            RVPOL(I)=0.0D0
            NP=I
            IF(RAD(I).GT.1.0D4) GO TO 35
          ENDDO
        ELSE
          STOP 'ELSEPA: Not enough memory space 5.'
        ENDIF
 35     CONTINUE
      ENDIF
C
C  ****  Muffin-tin absorption tail.
C
      IF(NMT.GT.0.AND.MABS.NE.0) THEN
        WAR=RW(NMT)/RAD(NMT)
        EIM=-WAR
        DO I=1,NP
          IF(RAD(I).LT.RAD(NMT)) THEN
            RW(I)=RW(I)-WAR*RAD(I)
          ELSE
            RW(I)=0.0D0
          ENDIF
        ENDDO
      ENDIF
C
C  ****  At high energies, we compute the DCS for scattering by the bare
C  nucleus and multiply it by a pre-evaluated screening factor.
C
 100  CONTINUE
      IF(IHEF0.EQ.1) THEN
        WRITE(IW,1800)
 1800   FORMAT(1X,'#',/1X,'# WARNING: High-energy factorization',
     1    ' with free-atom DF screening.',/1X,'#',
     2    10X,'Absorption, polarization and exchange corrections are',
     3    /1X,'#',10X,'switched off.',/1X,'#',10X,
     4    'Phase shifts are calculated for the bare nucleus.'/1X,'#',
     5    10X,'Scattering amplitudes are not evaluated.')
C  ****  Read screening function from data files.
        JT=IZ
        J1=JT-10*(JT/10)
        JT=(JT-J1)/10
        J2=JT-10*(JT/10)
        JT=(JT-J2)/10
        J3=JT-10*(JT/10)
        LIT1=LIT10(J1+1)
        LIT2=LIT10(J2+1)
        LIT3=LIT10(J3+1)
        FILE1=PATHE//'z_'//LIT3//LIT2//LIT1//'.dfs'
        NC=120
        CS120=' '
        I=0
        DO J=1,NC
          IF(FILE1(J:J).NE.' ') THEN
            I=I+1
            CS120(I:I)=FILE1(J:J)
          ENDIF
        ENDDO
        SCFILE=CS120
        WRITE(6,'(A)') SCFILE
C
        OPEN(99,FILE=SCFILE,STATUS='OLD',ERR=4)
        READ(99,'(1X,A1)') NULL
        READ(99,'(1X,A1)') NULL
        READ(99,'(1X,A1)') NULL
        NQS=0
        DO I=1,NGT
          READ(99,*,END=4) Q2T(I),FQ(I)
          NQS=I
        ENDDO
 4      CONTINUE
        CLOSE(UNIT=99)
        IF(NQS.EQ.0) STOP 'ELSEPA: I/O error. SCFILE does not exist.'
      ELSE IF(IHEF0.EQ.2) THEN
        WRITE(IW,1801)
 1801   FORMAT(1X,'#',/1X,'# WARNING: Scattering by the bare ',
     1    'nucleus.')
      ENDIF
C
      IF(IHEF0.GT.0) THEN
C  ----  The Coulomb tail of the potential is removed.
        NP=NPOT
        ZTAIL=RVN(NP)
        DO I=NPOT,1,-1
          IF(ABS(RVN(I)-ZTAIL).GT.1.0D-10) THEN
            NP=I+4
            GO TO 41
          ENDIF
        ENDDO
 41     CONTINUE
        DO I=1,NP
          RAD(I)=R(I)
          RV(I)=DBLE(IELEC)*RVN(I)
          RVST(I)=RVN(I)
          RVEX(I)=0.0D0
          RVPOL(I)=0.0D0
          RW(I)=0.0D0
        ENDDO
        RV(NP)=DBLE(IELEC)*DBLE(IZ)
C  ----  A Coulomb tail is appended to the potential table...
        IF(NP.LT.NDIM-4) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RVST(I)=RVST(I-1)
          RV(I)=RV(I-1)
          DO I=NP+2,NDIM
            RAD(I)=2.0D0*RAD(I-1)
            RV(I)=RV(I-1)
            RVST(I)=RVST(I-1)
            RVEX(I)=0.0D0
            RVPOL(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            IF(RAD(I).GT.1.0D4) GO TO 5
          ENDDO
 5        CONTINUE
        ELSE
          STOP 'ELSEPA: Not enough memory space 6'
        ENDIF
      ENDIF
C
      EEV=EV+ERE0*HREV
C
      OPEN(99,FILE='scfield.dat')
      WRITE(99,3001)
 3001 FORMAT(1X,'#  Scattering field.',/1X,'#  All quantities in',
     1  ' atomic units (a.u.), unless otherwise indicated.')
      IF(IELEC.EQ.-1) THEN
        WRITE(99,3002) IZ,NELEC
 3002   FORMAT(1X,'#  Z =',I4,', NELEC =',I4,
     1    ',   projectile: electron')
      ELSE
        WRITE(99,3003) IZ,NELEC
 3003   FORMAT(1X,'#  Z =',I4,', NELEC =',I4,
     1    ',   projectile: positron')
      ENDIF
      WRITE(99,3004) EV/HREV,EV
 3004 FORMAT(1X,'#  Kinetic energy =',1P,E12.5,' a.u. =',
     1    E12.5,' eV',/1X,'#')
      IF(NMT.GT.0) THEN
        WRITE(99,3005) NMT,RAD(NMT)
 3005   FORMAT(1X,'#  Muffin-tin radius = RAD(',I3,') =',
     1     1P,E12.5,' a.u.')
        WRITE(99,3105) VMOL
 3105   FORMAT(1X,'#  Atomic density    =',1P,E12.5,' 1/cm**3')
        WRITE(99,3006) ERE0,ERE0*HREV
 3006   FORMAT(1X,'#  Zero-energy shift = ',1P,E12.5,
     1     ' a.u. = ',E12.5,' eV')
        WRITE(99,3106) EIM,EIM*HREV
 3106   FORMAT(1X,'#  Background absorption potential = ',1P,E12.5,
     1     ' a.u. = ',E12.5,' eV')
        WRITE(99,3007) EEV/HREV,EEV
 3007   FORMAT(1X,'#  Effective kinetic energy =',1P,E12.5,
     1     ' a.u. =',E12.5,' eV',/1X,'#')
      ENDIF
      WRITE(99,3008)
 3008 FORMAT(1X,'#',3X,'i',7X,'r',11X,'r*V',9X,'r*Vst',8X,'r*Vex',7X,
     1      'r*Vpol',7X,'r*Wabs',8X,'rho_e',/1X,'#',96('-'))
      DO I=1,NP
        IF(IHEF0.GT.0) THEN
          RHO=0.0D0
        ELSE
          RHO=DEN(MAX(2,I))/(FOURPI*RAD(MAX(2,I))**2)
        ENDIF
        WRITE(99,'(2X,I4,1P,7E13.5)') I,RAD(I),RV(I),
     1    IELEC*RVST(I),RVEX(I),RVPOL(I),RW(I),RHO
      ENDDO
      CLOSE(99)
C
C  ************  Partial-wave analysis.
C
      NRT=NP
      NPI=NP
      DO I=1,NRT
        RRR(I)=RAD(I)
        RADI(I)=RAD(I)
        RVI(I)=RV(I)
      ENDDO
C
      NDELTA=NDM
      IF(IAB.EQ.0) THEN
        IF(IHEF0.GT.0) THEN
          CALL DPWA0(EEV,NDELTA,2)
        ELSE
          IF(EEV.LT.1.001D3) THEN
            ISCH=1
          ELSE
            ISCH=2
          ENDIF
          IF(MUFIN.EQ.1) ISCH=1
          CALL DPWA0(EEV,NDELTA,ISCH)
        ENDIF
        TOTCS=ECS
        ABCS=0.0D0
      ELSE
        IF(EEV.LT.1.001D3) THEN
          ISCH=1
        ELSE
          ISCH=2
        ENDIF
        IF(MUFIN.EQ.1) ISCH=1
        CALL DPWAI0(EEV,TOTCS,ABCS,NDELTA,ISCH)
      ENDIF
C
C  ************  DCS table.
C
      IF(IHEF0.EQ.1) THEN
        NTABT=NTAB
        DO I=1,NTAB
C  ****  Screening correction and check for numerical artifacts.
          Q2=4.0D0*RK2*XT(I)
          IF(Q2.LT.Q2T(NQS)) THEN
            CALL FINDI(Q2,Q2T,NQS,J)
            F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
          ELSE
            F=1.0D0
          ENDIF
          IF(TH(I).GT.1.0D0.AND.ERROR(I).LT.1.0D-2) THEN
            DCST(I)=DCST(I)*(Q2*F/(1.0D0+Q2))**2
          ELSE IF(TH(I).LT.10.0001D0) THEN
            THRAD=TH(I)*PI/180.0D0
            RMR=DPWAC(THRAD)
            DCST(I)=RUTHC*F*F*RMR/(1.0D0+Q2)**2
            ERROR(I)=1.0D-5
          ELSE
            IF(ERROR(I).GT.1.0D-2) THEN
              NTABT=I-1
              GO TO 6
            ENDIF
          ENDIF
          SPOL(I)=0.0D0
        ENDDO
 6      CONTINUE
        IF(NTABT.LT.NTAB) THEN
          DO I=NTABT,NTAB
            DCST(I)=1.0D-45
            ERROR(I)=1.0D0
          ENDDO
        ENDIF
      ENDIF
C
C  ****  Small-angle DCS for ions.
C
      XTL=0.0D0
      IMATCH=1
      IF(IZ.NE.NELEC) THEN
        DO I=1,NTAB-2
          IF(MAX(ERROR(I),ERROR(I+1),ERROR(I+2)).LT.5.0D-4) THEN
            IMATCH=I
            GO TO 7
          ENDIF
        ENDDO
 7      CONTINUE
        THL=MAX(0.5D0,TH(IMATCH))
        WRITE(IW,2012) THL
 2012   FORMAT(1X,'#',/1X,'# WARNING: DCSs are calculated and inte',
     1    'grated only for angles',/1X,'#',10X,'THETA .gt.',1P,E10.3,
     2    ' deg')
        XTL=SIN(THL*PI/360.0D0)**2
        DO I=1,IMATCH-1
          DCST(I)=0.0D0
          SPOL(I)=0.0D0
          ERROR(I)=1.0D0
        ENDDO
      ENDIF
C
C  ****  Integrated cross sections.
C
      ECS0=FOURPI*SMOMLL(XT,DCST,XTL,XT(NTAB),NTAB,0,0)
      ECS1=FOURPI*SMOMLL(XT,DCST,XTL,XT(NTAB),NTAB,1,0)
      ECS2=FOURPI*SMOMLL(XT,DCST,XTL,XT(NTAB),NTAB,2,0)
      ECS=ECS0
      TCS1=2.0D0*ECS1
      TCS2=6.0D0*(ECS1-ECS2)
      WRITE(IW,2212)
 2212 FORMAT(1X,'#',/1X,'# Integrated cross sections:')
      WRITE(IW,2013) ECS,ECS/A0B2
 2013 FORMAT(1X,'# Total elastic cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2')
      WRITE(IW,2014) TCS1,TCS1/A0B2
 2014 FORMAT(1X,'# 1st transport cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2')
      WRITE(IW,2015) TCS2,TCS2/A0B2
 2015 FORMAT(1X,'# 2nd transport cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2',/1X,'#')
      IF(IAB.EQ.1.AND.IZ.EQ.NELEC) THEN
        IF(NMT.GT.0) THEN
C  ****  Muffin-tin, atoms in solids. Contribution of the constant
C  imaginary potential to the absorption cross section.
          BETA=SQRT(EV*(EV+2.0D0*REV)/(EV+REV)**2)
          ABCSO=A0B2*2.0D0*EIM/(BETA*SL*VMOL1)
          ABCS=ABCS+ABCSO
          TOTCS=TOTCS+ABCSO
        ENDIF
C
        WRITE(IW,2016) ABCS,ABCS/A0B2
 2016   FORMAT(1X,'#    Absorption cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2')
        WRITE(IW,2017) TOTCS,TOTCS/A0B2
 2017   FORMAT(1X,'#   Grand total cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2',/1X,'#')
      ENDIF
C
      ECUT=MAX(20.0D3*IZ,2.0D6)
      WRITE(IW,'(1X,''#'')')
      WRITE(IW,'(1X,''# Differential cross section:'',23X,
     1  ''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',
     1  10X,''DCS'',8X,''Sherman'',6X,''error''/1X,''#  (deg)'',
     2  17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',4X,''function'',
     3  /1X,''#'',70(''-''))')
C
      DO I=1,NTAB
        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I),ERROR(I)
      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,3E13.5,E9.1)
C
C  ****  Scattering amplitudes.
C
      OPEN(99,FILE='scatamp.dat')
      WRITE(99,'(1X,''#  Scattering amplitudes (in cm)'',
     1  6X,''MU=(1-cos(TH))/2'')')
      IF(IELEC.EQ.-1) THEN
        WRITE(99,3010) IZ
 3010   FORMAT(1X,'#  Z =',I4,',   projectile: electron')
      ELSE
        WRITE(99,3011) IZ
 3011   FORMAT(1X,'#  Z =',I4,',   projectile: positron')
      ENDIF
      WRITE(99,3012) EV
 3012 FORMAT(1X,'#  Kinetic energy =',1P,E12.5,' eV',/1X,'#')
      WRITE(99,'(1X,''# TH (deg)'',6X,''MU'',10X,''Re(F)'',8X,
     1  ''Im(F)'',8X,''Re(G)'',8X,''Im(G)'',/1X,''#'',
     2  74(''-''))')
      IF(IHEF0.NE.0) RETURN
      DO I=1,NTAB
        THRAD=TH(I)*PI/180.0D0
        IF(EV.LT.ECUT) THEN
          CALL DPWA(THRAD,CF,CG,DCS,SPL,ERRF,ERRG)
        ELSE
          CF=0.0D0
          CG=0.0D0
        ENDIF
        WRITE(99,3013) TH(I),XT(I),CF,CG
      ENDDO
 3013 FORMAT(1X,1P,E10.3,E13.5,4E13.5)
      CLOSE(99)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE EFIELD
C  *********************************************************************
      SUBROUTINE EFIELD(IZ,NELEC,MNUCL,MELEC,IW,IWR)
C
C     Electrostatic potential of atoms and ions.
C
C  Input parameters:
C     IZ ....... atomic number (INTEGER).
C     NELEC .... number of electrons (INTEGER, .GE.0 and .LE.IZ).
C     MNUCL .... nuclear density model (INTEGER).
C                 1 --> point nucleus,
C                 2 --> uniform distribution,
C                 3 --> Fermi distribution,
C                 4 --> Helm's uniform-uniform distribution.
C     MELEC .... electron density model (INTEGER).
C                 1 --> TFM analytical density,
C                 2 --> TFD analytical density,
C                 3 --> DHFS analytical density,
C                 4 --> DF density from pre-evaluated files,
C                 5 --> density read from file 'density.usr'.
C     IW ....... output unit (to be defined in the main program).
C     IWR ...... if >0, the potential is written in a file named
C                'esfield.dat'.
C
C  Output (through common block /CFIELD/):
C     R(I) ..... radial grid points. R(1)=0.0D0.
C     RVN(I) ... nuclear potential times R.
C     DEN(I) ... radial electron density, i.e. the electron density
C                multiplied by 4*PI*R**2.
C     RVST(I) ... atomic electrostatic potential (nuclear+electronic)
C                times R.
C     NPOT ..... number of grid points where the potential function is
C                tabulated. For I.GT.NPOT, RVST(I) is set equal to
C                RVST(NPOT).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*1 LIT10(10),LIT1,LIT2,LIT3
      CHARACTER*120 ELFILE,FILE1,CS120,NULL
      CHARACTER*2 LSYMBL(103)
C
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
C
      DIMENSION DIFR(NDIM),DENN(NDIM),RVE(NDIM),ELAW(103)
      COMMON/CFIELD/R(NDIM),RVN(NDIM),DEN(NDIM),RVST(NDIM),NPOT
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      DIMENSION AUX(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
C
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
C
      DATA LSYMBL       /' H','He','Li','Be',' B',' C',' N',' O',
     1    ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K',
     2    'Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3    'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb',
     4    'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
     5    ' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     6    'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W',
     7    'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',
     8    'Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm',
     9    'Bk','Cf','Es','Fm','Md','No','Lr'/
      DATA ELAW     /1.007900D0,4.002600D0,6.941000D0,9.012200D0,
     1    1.081100D1,1.201070D1,1.400670D1,1.599940D1,1.899840D1,
     2    2.017970D1,2.298980D1,2.430500D1,2.698150D1,2.808550D1,
     3    3.097380D1,3.206600D1,3.545270D1,3.994800D1,3.909830D1,
     4    4.007800D1,4.495590D1,4.786700D1,5.094150D1,5.199610D1,
     5    5.493800D1,5.584500D1,5.893320D1,5.869340D1,6.354600D1,
     6    6.539000D1,6.972300D1,7.261000D1,7.492160D1,7.896000D1,
     7    7.990400D1,8.380000D1,8.546780D1,8.762000D1,8.890590D1,
     8    9.122400D1,9.290640D1,9.594000D1,9.890630D1,1.010700D2,
     9    1.029055D2,1.064200D2,1.078682D2,1.124110D2,1.148180D2,
     1    1.187100D2,1.217600D2,1.276000D2,1.269045D2,1.312900D2,
     1    1.329054D2,1.373270D2,1.389055D2,1.401160D2,1.409076D2,
     2    1.442400D2,1.449127D2,1.503600D2,1.519640D2,1.572500D2,
     3    1.589253D2,1.625000D2,1.649303D2,1.672600D2,1.689342D2,
     4    1.730400D2,1.749670D2,1.784900D2,1.809479D2,1.838400D2,
     5    1.862070D2,1.902300D2,1.922170D2,1.950780D2,1.969666D2,
     6    2.005900D2,2.043833D2,2.072000D2,2.089804D2,2.089824D2,
     7    2.099871D2,2.220176D2,2.230197D2,2.260254D2,2.270277D2,
     8    2.320381D2,2.310359D2,2.380289D2,2.370482D2,2.440642D2,
     9    2.430614D2,2.470000D2,2.470000D2,2.510000D2,2.520000D2,
     1    2.570000D2,2.580000D2,2.590000D2,2.620000D2/
C
C  ****  Path to the ELSEPA database
      CHARACTER*100 PATHE
      PATHE='./database/'
C
      IF(IZ.LE.0) STOP 'EFIELD: Negative atomic number.'
      IF(IZ.GT.103) STOP 'EFIELD: Atomic number larger than 103.'
      IF(NELEC.GT.IZ) STOP 'EFIELD: Negative ion.'
      IF(NELEC.LT.0) STOP 'EFIELD: Negative number of electrons.'
C
      Z=DBLE(IZ)
      AW=ELAW(IZ)
      IF(IW.GT.0) WRITE(IW,1001) LSYMBL(IZ),IZ,AW
 1001 FORMAT(1X,'#',/1X,'# Element: ',A2,',  Z = ',I3,
     1  ',  atomic weight =',1P,E12.5,' g/mol')
      NDIN=1000
C
C  ************  Nuclear electrostatic potential (times R).
C
      IF(MNUCL.EQ.1) THEN
        IF(IW.GT.0) WRITE(IW,1002)
 1002   FORMAT(1X,'#',/1X,'# Nuclear model: point charge')
C
        RTN=100.0D0
        RT2=2.0D-8
        DRN=RTN/DBLE(NDIN/3)
        CALL SGRID(R,DIFR,RTN,RT2,DRN,NDIN,NDIM,IER)
C
        DO I=1,NDIN
          RVN(I)=Z
          DENN(I)=0.0D0
        ENDDO
        GO TO 10
      ELSE IF (MNUCL.EQ.2) THEN
C  ****  The factor F2BOHR=1.889726D-5 transforms from fm to Bohr.
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=2.0D0*F2BOHR
        R1=R1*SQRT((1.0D0+2.5D0*(R2/R1)**2)
     1             /(1.0D0+0.75D0*(R2/R1)**2))
        R0=R1/10.0D0
C
        RTN=100.0D0
        RT2=MAX(0.1D0*R0,1.01D-8)
        DRN=RTN/DBLE(NDIN/3)
        CALL SGRID(R,DIFR,RTN,RT2,DRN,NDIN,NDIM,IER)
C
        IF(IW.GT.0) WRITE(IW,1004) R1*A0B
 1004   FORMAT(1X,'#',/1X,'# Nuclear model: uniform spherical d',
     1    'istribution',/1X,'#',16X,'Nuclear radius =',1P,E12.5,
     2    ' cm')
        DO I=1,NDIN
          X=R(I)
          IF(X.LT.R1) THEN
            RVN(I)=Z*(1.5D0-0.5D0*(X/R1)**2)*(X/R1)
            DENN(I)=Z/(FOURPI*R1**3/3.0D0)
          ELSE
            RVN(I)=Z
            DENN(I)=0.0D0
          ENDIF
        ENDDO
        GO TO 10
      ELSE IF (MNUCL.EQ.3) THEN
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=0.546D0*F2BOHR
        R0=R1/10.0D0
C
        RTN=100.0D0
        RT2=MAX(0.1D0*R0,1.01D-8)
        DRN=RTN/DBLE(NDIN/3)
        CALL SGRID(R,DIFR,RTN,RT2,DRN,NDIN,NDIM,IER)
C
        IF(IW.GT.0) WRITE(IW,1005) R1*A0B,R2*A0B
 1005   FORMAT(1X,'#',/1X,'# Nuclear model: Fermi distribution',
     1    /1X,'#',16X,'Average radius =',1P,E12.5, ' cm',
     2    /1X,'#',16X,'Skin thickness =',E12.5,' cm')
C  ****  The array DENN contains the nuclear charge density.
C        (unnormalized).
        DO I=1,NDIN
          X=R(I)
          XX=EXP((R1-X)/R2)
          DENN(I)=XX/(1.0D0+XX)
          RVN(I)=DENN(I)*X**2*DIFR(I)
        ENDDO
      ELSE
        RNUC=1.070D0*F2BOHR*AW**0.3333333333333333D0
        R1=0.96219D0*RNUC+0.435D0*F2BOHR
        R2=2.0D0*F2BOHR
        IF(IW.GT.0) WRITE(IW,1105) R1*A0B,R2*A0B
 1105   FORMAT(1X,'#',/1X,'# Nuclear model: Helm''s Uu distribu',
     1    'tion',/1X,'#',16X,'  Inner radius =',1P,E12.5, ' cm',
     2    /1X,'#',16X,'Skin thickness =',E12.5,' cm')
        IF(R2.GT.R1) THEN
          STORED=R1
          R1=R2
          R2=STORED
        ENDIF
        R0=R1/10.0D0
C
        RTN=100.0D0
        RT2=MAX(0.1D0*R0,1.01D-8)
        DRN=RTN/DBLE(NDIN/3)
        CALL SGRID(R,DIFR,RTN,RT2,DRN,NDIN,NDIM,IER)
C
        DO I=1,NDIN
          RR=R(I)
          IF(RR.LT.R1-R2) THEN
            V=1.0D0
          ELSE IF(RR.GT.R1+R2) THEN
            V=0.0D0
          ELSE
            T=RR*RR+R1*R1-R2*R2
            V1=(T+4.0D0*RR*R1)*(T-2.0D0*RR*R1)**2
            T=RR*RR+R2*R2-R1*R1
            V2=(T+4.0D0*RR*R2)*(T-2.0D0*RR*R2)**2
            V=(V1+V2)/(32.0D0*(R2*RR)**3)
          ENDIF
          DENN(I)=V
          RVN(I)=DENN(I)*RR**2*DIFR(I)
        ENDDO
      ENDIF
C
      CALL SLAG6(1.0D0,RVN,RVN,NDIN)
      NDIN1=NDIN+1
      DO I=1,NDIN
        K=NDIN1-I
        AUX(I)=DENN(K)*R(K)*DIFR(K)
      ENDDO
      CALL SLAG6(1.0D0,AUX,AUX,NDIN)
      FNORM=Z/RVN(NDIN)
      DO I=1,NDIN
        RVN(I)=FNORM*(RVN(I)+AUX(NDIN1-I)*R(I))
        DENN(I)=FNORM*DENN(I)/FOURPI
        IF(DENN(I).LT.1.0D-35) DENN(I)=0.0D0
      ENDDO
 10   CONTINUE
C
C  ************  Electronic electrostatic potential (times R).
C
      IF(IW.GT.0) WRITE(IW,1006) NELEC
 1006 FORMAT(1X,'#',/1X,'# Number of electrons =',I3)
      IF(NELEC.EQ.0) THEN
        DO I=1,NDIN
          DEN(I)=0.0D0
          RVE(I)=0.0D0
        ENDDO
        GO TO 2
      ENDIF
C
      IF(MELEC.LT.4.OR.MELEC.GT.5) THEN
C  ****  Analytical electron density models.
        IF(MELEC.EQ.1) THEN
          CALL TFM(IZ,A1,A2,A3,AL1,AL2,AL3)
          IF(IW.GT.0) WRITE(IW,1007)
 1007     FORMAT(1X,'#',/1X,
     1      '# Electron density: analytical TFM model')
        ELSE IF(MELEC.EQ.2) THEN
          CALL TFD(IZ,A1,A2,A3,AL1,AL2,AL3)
          IF(IW.GT.0) WRITE(IW,1008)
 1008     FORMAT(1X,'#',/1X,
     1      '# Electron density: analytical TFD model')
        ELSE
          CALL DHFS(IZ,A1,A2,A3,AL1,AL2,AL3)
          IF(IW.GT.0) WRITE(IW,1009)
 1009     FORMAT(1X,'#',/1X,
     1      '# Electron density: analytical DHFS model')
        ENDIF
        IF(IW.GT.0) WRITE(IW,1010) A1,AL1,A2,AL2,A3,AL3
 1010   FORMAT(1X,'#',19X,'A1 = ',1P,D12.5,' ,   ALPHA1 =',D12.5,
     1      /1X,'#',19X,'A2 = ',D12.5,' ,   ALPHA2 =',D12.5,
     2      /1X,'#',19X,'A3 = ',D12.5,' ,   ALPHA3 =',D12.5)
        XN=DBLE(NELEC)
        DO I=1,NDIN
          DEN(I)=(A1*AL1*AL1*EXP(-AL1*R(I))
     1           +A2*AL2*AL2*EXP(-AL2*R(I))
     2           +A3*AL3*AL3*EXP(-AL3*R(I)))*XN
        ENDDO
      ELSE
C  ****  Electron density read from a file.
        NE=0
        IF(MELEC.EQ.4) THEN
          JT=IZ
          J1=JT-10*(JT/10)
          JT=(JT-J1)/10
          J2=JT-10*(JT/10)
          JT=(JT-J2)/10
          J3=JT-10*(JT/10)
          LIT1=LIT10(J1+1)
          LIT2=LIT10(J2+1)
          LIT3=LIT10(J3+1)
          FILE1=PATHE//'z_'//LIT3//LIT2//LIT1//'.den'
          NC=120
          CS120=' '
          I=0
          DO J=1,NC
            IF(FILE1(J:J).NE.' ') THEN
              I=I+1
              CS120(I:I)=FILE1(J:J)
            ENDIF
          ENDDO
          ELFILE=CS120
          WRITE(6,'(A)') ELFILE
        ELSE
          ELFILE='density.usr'
        ENDIF
C
        IF(IW.GT.0) WRITE(IW,1011) ELFILE
 1011   FORMAT(1X,'#',/1X,'# Electron density: Read from file ',A120)
        OPEN(99,FILE=ELFILE,STATUS='OLD',ERR=1)
        READ(99,'(A12)') NULL
        READ(99,'(A12)') NULL
        READ(99,'(A12)') NULL
        DO I=1,NDIN
          READ(99,*,END=1) AUX(I),RVE(I)
          NE=I
          RVE(I)=LOG(RVE(I))
        ENDDO
        STOP 'EFIELD: File is too large.'
 1      CONTINUE
        IF(NE.EQ.0) STOP 'EFIELD: I/O error in EFIELD.'
        CLOSE(99)
C
        IF(IW.GT.0) WRITE(IW,1012) NE
 1012   FORMAT(1X,'#',19X,'Number of data points = ',I4)
        IF(NE.LT.4) STOP 'EFIELD: SPLINE needs more than 4 points.'
C  ****  ... and interpolated (lin-log cubic spline).
        CALL SPLINE(AUX,RVE,A,B,C,D,0.0D0,0.0D0,NE)
        B(NE)=(RVE(NE)-RVE(NE-1))/(AUX(NE)-AUX(NE-1))
        A(NE)=RVE(NE-1)-B(NE)*AUX(NE-1)
        C(NE)=0.0D0
        D(NE)=0.0D0
        DO I=1,NDIN
          X=R(I)
          IF(X.GT.AUX(NE)) THEN
            DEN(I)=0.0D0
          ELSE
            CALL FINDI(X,AUX,NE,J)
            DEN(I)=EXP(A(J)+X*(B(J)+X*(C(J)+X*D(J))))*X*FOURPI
          ENDIF
        ENDDO
      ENDIF
C  ****  Calculation of the electrostatic potential.
      DO I=1,NDIN
        RVE(I)=DEN(I)*R(I)*DIFR(I)
      ENDDO
      CALL SLAG6(1.0D0,RVE,RVE,NDIN)
      NDIN1=NDIN+1
      DO I=1,NDIN
        K=NDIN1-I
        AUX(I)=DEN(K)*DIFR(K)
      ENDDO
      CALL SLAG6(1.0D0,AUX,AUX,NDIN)
      IF(IW.GT.0) WRITE(IW,1013) RVE(NDIN)
 1013 FORMAT(1X,'#',19X,'Volume integral =',1P,E12.5)
      FNORM=DBLE(NELEC)/RVE(NDIN)
      DO I=1,NDIN
        RVE(I)=FNORM*(RVE(I)+AUX(NDIN1-I)*R(I))
        DEN(I)=FNORM*DEN(I)*R(I)
      ENDDO
C
 2    CONTINUE
      ZINF=DBLE(IZ-NELEC)
      DO I=1,NDIN
        RVST(I)=RVN(I)-RVE(I)
      ENDDO
      NPOT=NDIN
      DO I=NDIN,6,-1
        IF((ABS(RVST(I)-ZINF).LT.5.0D-12).
     1       AND.(ABS(DEN(I)).LT.5.0D-12)) THEN
          RVST(I)=ZINF
          NPOT=I
          IF(R(I).LT.1.0D0) GO TO 3
        ELSE
          GO TO 3
        ENDIF
      ENDDO
 3    CONTINUE
C
      IF(IWR.GT.0) THEN
        OPEN(99,FILE='esfield.dat')
        WRITE(99,'(1X,''#  Electrostatic field (a.u.)'')')
        WRITE(99,2001) LSYMBL(IZ),IZ,AW
 2001   FORMAT(1X,'#  Element: ',A2,',  Z = ',I3,
     1      ',  atomic weight =',1P,E14.7,' g/mol')
        WRITE(99,2002) NELEC
 2002 FORMAT(1X,'#  Number of electrons =',I3)
        WRITE(99,2003) MNUCL,MELEC
 2003 FORMAT(1X,'#  MNUCL =',I3,',   MELEC =',I3,/1X,'#')
        WRITE(99,2004)
 2004 FORMAT(1X,'#',3X,'I',7X,'R(I)',11X,'RHON(I)',9X,'RVN(I)',10X,
     1         'DEN(I)',10X,'RVST(I)',/1X,'#',84('-'))
        DO I=1,NPOT
         WRITE(99,'(2X,I4,1P,5E16.8)')
     1     I,R(I),DENN(I),RVN(I),DEN(I),RVST(I)
        ENDDO
        CLOSE(99)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE TFM
C  *********************************************************************
      SUBROUTINE TFM(IZ,A1,A2,A3,AL1,AL2,AL3)
C
C     Parameters in Moliere's analytical approximation (three Yukawa
C  terms) to the Thomas-Fermi atomic screening function.
C     Ref.: G. Moliere, Z. Naturforsch. 2a (1947) 133.
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      Z=DBLE(IZ)
      RTF=0.88534D0/Z**0.33333333333D0
      AL1=6.0D0/RTF
      AL2=1.2D0/RTF
      AL3=0.3D0/RTF
      A1=0.10D0
      A2=0.55D0
      A3=0.35D0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE TFD
C  *********************************************************************
      SUBROUTINE TFD(IZ,A1,A2,A3,AL1,AL2,AL3)
C
C     Parameters in the analytical approximation (three Yukawa terms)
C  for the Thomas-Fermi-Dirac atomic screening function.
C     Ref.: R.A. Bonham and T.G. Strand, J. Chem. Phys. 39 (1963) 2200.
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION AA1(5),AA2(5),AA3(5),AAL1(5),AAL2(5),AAL3(5)
      DATA AA1/1.26671D-2,-2.61047D-2,2.14184D-2,-2.35686D-3,
     12.10672D-5/
      DATA AA2/5.80612D-2,2.93077D-2,8.57135D-2,-2.23342D-2,
     11.64675D-3/
      DATA AA3/9.27968D-1,-1.64643D-3,-1.07685D-1,2.47998D-2,
     1-1.67822D-3/
      DATA AAL1/1.64564D2,-1.52192D2,6.23879D1,-1.15005D1,
     18.08424D-1/
      DATA AAL2/1.13060D1,-6.31902D0,2.26025D0,-3.70738D-1,
     12.61151D-2/
      DATA AAL3/1.48219D0,-5.57601D-2,1.64387D-2,-4.39703D-3,
     19.97225D-4/
C
      IF(IZ.LE.0) THEN
        WRITE(6,100)
 100    FORMAT(5X,'*** TFD: Negative atomic number. STOP.')
        STOP 'TFD: Negative atomic number.'
      ENDIF
C
      X=LOG(DBLE(IZ))
      A1=AA1(1)+X*(AA1(2)+X*(AA1(3)+X*(AA1(4)+X*AA1(5))))
      A2=AA2(1)+X*(AA2(2)+X*(AA2(3)+X*(AA2(4)+X*AA2(5))))
      A3=AA3(1)+X*(AA3(2)+X*(AA3(3)+X*(AA3(4)+X*AA3(5))))
      AL1=AAL1(1)+X*(AAL1(2)+X*(AAL1(3)+X*(AAL1(4)+X*AAL1(5))))
      AL2=AAL2(1)+X*(AAL2(2)+X*(AAL2(3)+X*(AAL2(4)+X*AAL2(5))))
      AL3=AAL3(1)+X*(AAL3(2)+X*(AAL3(3)+X*(AAL3(4)+X*AAL3(5))))
      A3=1.0D0-A1-A2
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DHFS
C  *********************************************************************
      SUBROUTINE DHFS(IZ,A1,A2,A3,AL1,AL2,AL3)
C
C     DHFS analytical screening function parameters for free neutral
C  atoms. The input argument is the atomic number.
C
C     Ref.: F. Salvat et al., Phys. Rev. A36 (1987) 467-474.
C     Elements from Z=93 to 103 added in march 1992.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION B1(103),B2(103),BL1(103),BL2(103),BL3(103)
      DATA B1/-7.05665D-6,-2.25920D-1,6.04537D-1,3.27766D-1,
     1   2.32684D-1,1.53676D-1,9.95750D-2,6.25130D-2,3.68040D-2,
     2   1.88410D-2,7.44440D-1,6.42349D-1,6.00152D-1,5.15971D-1,
     3   4.38675D-1,5.45871D-1,7.24889D-1,2.19124D+0,4.85607D-2,
     4   5.80017D-1,5.54340D-1,1.11950D-2,3.18350D-2,1.07503D-1,
     5   4.97556D-2,5.11841D-2,5.00039D-2,4.73509D-2,7.70967D-2,
     6   4.00041D-2,1.08344D-1,6.09767D-2,2.11561D-2,4.83575D-1,
     7   4.50364D-1,4.19036D-1,1.73438D-1,3.35694D-2,6.88939D-2,
     8   1.17552D-1,2.55689D-1,2.69313D-1,2.20138D-1,2.75057D-1,
     9   2.71053D-1,2.78363D-1,2.56210D-1,2.27100D-1,2.49215D-1,
     A   2.15313D-1,1.80560D-1,1.30772D-1,5.88293D-2,4.45145D-1,
     B   2.70796D-1,1.72814D-1,1.94726D-1,1.91338D-1,1.86776D-1,
     C   1.66461D-1,1.62350D-1,1.58016D-1,1.53759D-1,1.58729D-1,
     D   1.45327D-1,1.41260D-1,1.37360D-1,1.33614D-1,1.29853D-1,
     E   1.26659D-1,1.28806D-1,1.30256D-1,1.38420D-1,1.50030D-1,
     F   1.60803D-1,1.72164D-1,1.83411D-1,2.23043D-1,2.28909D-1,
     G   2.09753D-1,2.70821D-1,2.37958D-1,2.28771D-1,1.94059D-1,
     H   1.49995D-1,9.55262D-2,3.19155D-1,2.40406D-1,2.26579D-1,
     I   2.17619D-1,2.41294D-1,2.44758D-1,2.46231D-1,2.55572D-1,
     J   2.53567D-1,2.43832D-1,2.41898D-1,2.44050D-1,2.40237D-1,
     K   2.34997D-1,2.32114D-1,2.27937D-1,2.29571D-1/
      DATA B2/-1.84386D+2,1.22592D+0,3.95463D-1,6.72234D-1,
     1   7.67316D-1,8.46324D-1,9.00425D-1,9.37487D-1,9.63196D-1,
     2   9.81159D-1,2.55560D-1,3.57651D-1,3.99848D-1,4.84029D-1,
     3  5.61325D-1,-5.33329D-1,-7.54809D-1,-2.2852D0,7.75935D-1,
     4   4.19983D-1,4.45660D-1,6.83176D-1,6.75303D-1,7.16172D-1,
     5   6.86632D-1,6.99533D-1,7.14201D-1,7.29404D-1,7.95083D-1,
     6   7.59034D-1,7.48941D-1,7.15671D-1,6.70932D-1,5.16425D-1,
     7   5.49636D-1,5.80964D-1,7.25336D-1,7.81581D-1,7.20203D-1,
     8   6.58088D-1,5.82051D-1,5.75262D-1,5.61797D-1,5.94338D-1,
     9   6.11921D-1,6.06653D-1,6.50520D-1,6.15496D-1,6.43990D-1,
     A   6.11497D-1,5.76688D-1,5.50366D-1,5.48174D-1,5.54855D-1,
     B   6.52415D-1,6.84485D-1,6.38429D-1,6.46684D-1,6.55810D-1,
     C   7.05677D-1,7.13311D-1,7.20978D-1,7.28385D-1,7.02414D-1,
     D   7.42619D-1,7.49352D-1,7.55797D-1,7.61947D-1,7.68005D-1,
     E   7.73365D-1,7.52781D-1,7.32428D-1,7.09596D-1,6.87141D-1,
     F   6.65932D-1,6.46849D-1,6.30598D-1,6.17575D-1,6.11402D-1,
     G   6.00426D-1,6.42829D-1,6.30789D-1,6.21959D-1,6.10455D-1,
     H   6.03147D-1,6.05994D-1,6.23324D-1,6.56665D-1,6.42246D-1,
     I   6.24013D-1,6.30394D-1,6.29816D-1,6.31596D-1,6.49005D-1,
     J   6.53604D-1,6.43738D-1,6.48850D-1,6.70318D-1,6.76319D-1,
     K   6.65571D-1,6.88406D-1,6.94394D-1,6.82014D-1/
      DATA BL1/ 4.92969D+0,5.52725D+0,2.81741D+0,4.54302D+0,
     1   5.99006D+0,8.04043D+0,1.08122D+1,1.48233D+1,2.14001D+1,
     2   3.49994D+1,4.12050D+0,4.72663D+0,5.14051D+0,5.84918D+0,
     3   6.67070D+0,6.37029D+0,6.21183D+0,5.54701D+0,3.02597D+1,
     4   6.32184D+0,6.63280D+0,9.97569D+1,4.25330D+1,1.89587D+1,
     5   3.18642D+1,3.18251D+1,3.29153D+1,3.47580D+1,2.53264D+1,
     6   4.03429D+1,2.01922D+1,2.91996D+1,6.24873D+1,8.78242D+0,
     7   9.33480D+0,9.91420D+0,1.71659D+1,5.52077D+1,3.13659D+1,
     8   2.20537D+1,1.42403D+1,1.40442D+1,1.59176D+1,1.43137D+1,
     9   1.46537D+1,1.46455D+1,1.55878D+1,1.69141D+1,1.61552D+1,
     A   1.77931D+1,1.98751D+1,2.41540D+1,3.99955D+1,1.18053D+1,
     B   1.65915D+1,2.23966D+1,2.07637D+1,2.12350D+1,2.18033D+1,
     C   2.39492D+1,2.45984D+1,2.52966D+1,2.60169D+1,2.54973D+1,
     D   2.75466D+1,2.83460D+1,2.91604D+1,2.99904D+1,3.08345D+1,
     E   3.16806D+1,3.13526D+1,3.12166D+1,3.00767D+1,2.86302D+1,
     F   2.75684D+1,2.65861D+1,2.57339D+1,2.29939D+1,2.28644D+1,
     G   2.44080D+1,2.09409D+1,2.29872D+1,2.37917D+1,2.66951D+1,
     H   3.18397D+1,4.34890D+1,2.00150D+1,2.45012D+1,2.56843D+1,
     I   2.65542D+1,2.51930D+1,2.52522D+1,2.54271D+1,2.51526D+1,
     J   2.55959D+1,2.65567D+1,2.70360D+1,2.72673D+1,2.79152D+1,
     K   2.86446D+1,2.93353D+1,3.01040D+1,3.02650D+1/
      DATA BL2/ 2.00272D+0,2.39924D+0,6.62463D-1,9.85154D-1,
     1   1.21347D+0,1.49129D+0,1.76868D+0,2.04035D+0,2.30601D+0,
     2   2.56621D+0,8.71798D-1,1.00247D+0,1.01529D+0,1.17314D+0,
     3   1.34102D+0,2.55169D+0,3.38827D+0,4.56873D+0,3.12426D+0,
     4   1.00935D+0,1.10227D+0,4.12865D+0,3.94043D+0,3.06375D+0,
     5   3.78110D+0,3.77161D+0,3.79085D+0,3.82989D+0,3.39276D+0,
     6   3.94645D+0,3.47325D+0,4.12525D+0,4.95015D+0,1.69671D+0,
     7   1.79002D+0,1.88354D+0,3.11025D+0,4.28418D+0,4.24121D+0,
     8   4.03254D+0,2.97020D+0,2.86107D+0,3.36719D+0,2.73701D+0,
     9   2.71828D+0,2.61549D+0,2.74124D+0,3.08408D+0,2.88189D+0,
     A   3.29372D+0,3.80921D+0,4.61191D+0,5.91318D+0,1.79673D+0,
     B   2.69645D+0,3.45951D+0,3.46574D+0,3.48193D+0,3.50982D+0,
     C   3.51987D+0,3.55603D+0,3.59628D+0,3.63834D+0,3.73639D+0,
     D   3.72882D+0,3.77625D+0,3.82444D+0,3.87344D+0,3.92327D+0,
     E   3.97271D+0,4.09040D+0,4.20492D+0,4.24918D+0,4.24261D+0,
     F   4.23412D+0,4.19992D+0,4.14615D+0,3.73461D+0,3.69138D+0,
     G   3.96429D+0,3.24563D+0,3.62172D+0,3.77959D+0,4.25824D+0,
     H   4.92848D+0,5.85205D+0,2.90906D+0,3.55241D+0,3.79223D+0,
     I   4.00437D+0,3.67795D+0,3.63966D+0,3.61328D+0,3.43021D+0,
     J   3.43474D+0,3.59089D+0,3.59411D+0,3.48061D+0,3.50331D+0,
     K   3.61870D+0,3.55697D+0,3.58685D+0,3.64085D+0/
      DATA BL3/ 1.99732D+0,1.00000D+0,1.00000D+0,1.00000D+0,
     1   1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,
     2   1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,
     3   1.00000D+0,1.67534D+0,1.85964D+0,2.04455D+0,7.32637D-1,
     4   1.00000D+0,1.00000D+0,1.00896D+0,1.05333D+0,1.00137D+0,
     5   1.12787D+0,1.16064D+0,1.19152D+0,1.22089D+0,1.14261D+0,
     6   1.27594D+0,1.00643D+0,1.18447D+0,1.35819D+0,1.00000D+0,
     7   1.00000D+0,1.00000D+0,7.17673D-1,8.57842D-1,9.47152D-1,
     8   1.01806D+0,1.01699D+0,1.05906D+0,1.15477D+0,1.10923D+0,
     9   1.12336D+0,1.43183D+0,1.14079D+0,1.26189D+0,9.94156D-1,
     A   1.14781D+0,1.28288D+0,1.41954D+0,1.54707D+0,1.00000D+0,
     B   6.81361D-1,8.07311D-1,8.91057D-1,9.01112D-1,9.10636D-1,
     C   8.48620D-1,8.56929D-1,8.65025D-1,8.73083D-1,9.54998D-1,
     D   8.88981D-1,8.96917D-1,9.04803D-1,9.12768D-1,9.20306D-1,
     E   9.28838D-1,1.00717D+0,1.09456D+0,1.16966D+0,1.23403D+0,
     F   1.29699D+0,1.35350D+0,1.40374D+0,1.44284D+0,1.48856D+0,
     G   1.53432D+0,1.11214D+0,1.23735D+0,1.25338D+0,1.35772D+0,
     H   1.46828D+0,1.57359D+0,7.20714D-1,8.37599D-1,9.33468D-1,
     I   1.02385D+0,9.69895D-1,9.82474D-1,9.92527D-1,9.32751D-1,
     J   9.41671D-1,1.01827D+0,1.02554D+0,9.66447D-1,9.74347D-1,
     K   1.04137D+0,9.90568D-1,9.98878D-1,1.04473D+0/
C
      IIZ=IABS(IZ)
      IF(IIZ.GT.103) IIZ=103
      IF(IIZ.EQ.0) IIZ=1
      A1=B1(IIZ)
      A2=B2(IIZ)
      A3=1.0D0-(A1+A2)
      IF(ABS(A3).LT.1.0D-15) A3=0.0D0
      AL1=BL1(IIZ)
      AL2=BL2(IIZ)
      AL3=BL3(IIZ)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE MOTTCS
C  *********************************************************************
      SUBROUTINE MOTTSC(IELEC,IZ,EV,IW)
C
C     Mott cross section for elastic scattering of high-energy electrons
C  and positrons by unscreened point nuclei.
C
C  Input parameters:
C    IELEC ..... electron-positron flag;
C                =-1 for electrons,
C                =+1 for positrons.
C    IZ ........ atomic number of the target atom.
C    EV ........ projectile's kinetic energy (in eV).
C    IW ........ output unit (to be defined in the main program).
C
C  The Mott DCS and spin polarization function are printed on unit IW.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
C
      PARAMETER (SL=137.035999074D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B=5.2917721092D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.21138505D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
C
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
C
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine MOTTSC. Elastic scattering of ',
     1  'electrons and positrons',/1X,'#',20X,
     2  'by unscreened Coulomb fields')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',
     1       E12.5,' a.u.')
C
      IF(IZ.LE.0) STOP 'MOTTCS: Negative atomic number.'
      WRITE(IW,1001) IZ
 1001 FORMAT(1X,'#',/1X,'# Z = ',I3)
C
      Z=DBLE(IZ*IELEC)
      CALL DPWAC0(Z,EV)
C
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
 10   CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
C
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
        IF(TH(I).GT.1.0D-5) THEN
          Q2=4.0D0*RK2*XT(I)
          RMR=DPWAC(THR)
          DCST(I)=RUTHC*RMR/Q2**2
C  ****  Spin polarization (Sherman) function.
          CF=CF*A0B
          CG=CG*A0B
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          DCS=ACF+ACG
          IF(DCS.GT.1.0D-45) THEN
            ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/DCS
          ELSE
            ERR=1.0D0
          ENDIF
          ERROR(I)=ERR
        ELSE
          DCST(I)=1.0D-45
          ERROR(I)=1.0D0
        ENDIF
      ENDDO
C
      WRITE(IW,'(1X,''#'')')
      WRITE(IW,'(1X,''# Differential cross section'',6X,
     1  ''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',
     1  10X,''DCS'',7X,''McKinl-Fesh'',4X,''error''/1X,''#  (deg)'',
     2  17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',3X,''(a0**2/sr)'',
     3  /1X,''#'',71(''-''))')
C
      BETA2=E*(E+2.0D0*SL**2)/(E+SL**2)**2
      BETA=SQRT(BETA2)
      GAMMA=1.0D0+E/SL**2
      CONS=(0.5D0*Z/(BETA2*GAMMA*SL**2))**2
C
      DO I=1,NTAB
        STH2=SIN(0.5D0*MAX(TH(I),1.0D-5)*PI/180.0D0)
        XSMOT=(CONS/STH2**4)*(1.0D0-BETA2*STH2**2
     1    -PI*Z*(BETA/SL)*STH2*(1.0D0-STH2))
        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,XSMOT,ERROR(I)
      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,3E13.5,2X,E8.1)
C
      ECS=1.0D35
      TCS1=1.0D35
      TCS2=1.0D35
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE HEBORN
C  *********************************************************************
      SUBROUTINE HEBORN(IELEC,IZ,MNUCL,EV,IW)
C
C     Mott-Born cross section for elastic scattering of high-energy
C  electrons and positrons by neutral atoms.
C
C    The DCS is obtained as the product of the Mott DCS for a point
C  nucleus, the Helm uniform-uniform nuclear form factor (with an
C  empirical Coulomb correction) and the high-energy DF screening
C  factor.
C
C  Input parameters:
C    IELEC ..... electron-positron flag;
C                =-1 for electrons,
C                =+1 for positrons.
C    IZ ........ atomic number of the target atom.
C    MNUCL ..... nuclear charge density model.
C                  1 --> point nucleus (P),
C                  2 --> uniform distribution (U),
C                3,4 --> Helm's uniform-uniform distribution (Uu).
C    EV ........ projectile's kinetic energy (in eV).
C    IW ........ output unit (to be defined in the main program).
C
C  Output (through the common block /DCSTAB/):
C     ECS ........ total cross section (cm**2).
C     TCS1 ....... 1st transport cross section (cm**2).
C     TCS2 ....... 2nd transport cross section (cm**2).
C     TH(I) ...... scattering angles (in deg)
C     XT(I) ...... values of (1-COS(TH(I)))/2.
C     DCST(I) .... differential cross section per unit solid angle at
C                  TH(I) (in cm**2/sr).
C     ERROR(I) ... relative uncertainty of the computed DCS values.
C                  Estimated from the convergence of the series.
C     NTAB ....... number of angles in the table.
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
C
      PARAMETER (SL=137.035999074D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B=5.2917721092D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.21138505D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
      COMMON/CDCSHE/Q2T(NGT),FQ(NGT),U1,U2,NQS,MOM
C
      DIMENSION ELAW(103)
C
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
C
      CHARACTER*120 SCFILE,FILE1,CS120,NULL
      CHARACTER*1 LIT10(10),LIT1,LIT2,LIT3
      CHARACTER*2 LSYMBL(103)
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
C
      DATA LSYMBL       /' H','He','Li','Be',' B',' C',' N',' O',
     1    ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K',
     2    'Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3    'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb',
     4    'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',
     5    ' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     6    'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W',
     7    'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At',
     8    'Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm',
     9    'Bk','Cf','Es','Fm','Md','No','Lr'/
C
      DATA ELAW     /1.007900D0,4.002600D0,6.941000D0,9.012200D0,
     1    1.081100D1,1.201070D1,1.400670D1,1.599940D1,1.899840D1,
     2    2.017970D1,2.298980D1,2.430500D1,2.698150D1,2.808550D1,
     3    3.097380D1,3.206600D1,3.545270D1,3.994800D1,3.909830D1,
     4    4.007800D1,4.495590D1,4.786700D1,5.094150D1,5.199610D1,
     5    5.493800D1,5.584500D1,5.893320D1,5.869340D1,6.354600D1,
     6    6.539000D1,6.972300D1,7.261000D1,7.492160D1,7.896000D1,
     7    7.990400D1,8.380000D1,8.546780D1,8.762000D1,8.890590D1,
     8    9.122400D1,9.290640D1,9.594000D1,9.890630D1,1.010700D2,
     9    1.029055D2,1.064200D2,1.078682D2,1.124110D2,1.148180D2,
     1    1.187100D2,1.217600D2,1.276000D2,1.269045D2,1.312900D2,
     1    1.329054D2,1.373270D2,1.389055D2,1.401160D2,1.409076D2,
     2    1.442400D2,1.449127D2,1.503600D2,1.519640D2,1.572500D2,
     3    1.589253D2,1.625000D2,1.649303D2,1.672600D2,1.689342D2,
     4    1.730400D2,1.749670D2,1.784900D2,1.809479D2,1.838400D2,
     5    1.862070D2,1.902300D2,1.922170D2,1.950780D2,1.969666D2,
     6    2.005900D2,2.043833D2,2.072000D2,2.089804D2,2.089824D2,
     7    2.099871D2,2.220176D2,2.230197D2,2.260254D2,2.270277D2,
     8    2.320381D2,2.310359D2,2.380289D2,2.370482D2,2.440642D2,
     9    2.430614D2,2.470000D2,2.470000D2,2.510000D2,2.520000D2,
     1    2.570000D2,2.580000D2,2.590000D2,2.620000D2/
C
      EXTERNAL DCSHB
C
C  ****  Path to the ELSEPA database
      CHARACTER*100 PATHE
      PATHE='./database/'
C
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine HEBORN. Elastic scattering of ',
     1  'electrons and positrons',/1X,'#',20X,
     2  'by neutral atoms')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',
     1       E12.5,' a.u.')
C
      WRITE(IW,'(1X,''#'',/1X,''#  ***  WARNING: High-energy '',
     1  ''Mott-Born approximation. Neutral atom.'')')
C
      IF(IZ.LE.0) STOP 'HEBORN: Negative atomic number.'
      IF(IZ.GT.103) STOP 'HEBORN: Atomic number larger than 103.'
      AW=ELAW(IZ)
      WRITE(IW,1001) LSYMBL(IZ),IZ,AW
 1001 FORMAT(1X,'#',/1X,'# Element: ',A2,',  Z = ',I3,
     1  ',  atomic weight =',1P,E12.5,' g/mol')
C
      Z=DBLE(IZ*IELEC)
      CALL DPWAC0(Z,EV)
C
C  ****  Read screening function from data files.
      JT=IZ
      J1=JT-10*(JT/10)
      JT=(JT-J1)/10
      J2=JT-10*(JT/10)
      JT=(JT-J2)/10
      J3=JT-10*(JT/10)
      LIT1=LIT10(J1+1)
      LIT2=LIT10(J2+1)
      LIT3=LIT10(J3+1)
      FILE1=PATHE//'z_'//LIT3//LIT2//LIT1//'.dfs'
      NC=120
      CS120=' '
      I=0
      DO J=1,NC
        IF(FILE1(J:J).NE.' ') THEN
          I=I+1
          CS120(I:I)=FILE1(J:J)
        ENDIF
      ENDDO
      SCFILE=CS120
      WRITE(6,'(A)') SCFILE
      OPEN(99,FILE=SCFILE,STATUS='OLD',ERR=4)
      READ(99,'(1X,A1)') NULL
      READ(99,'(1X,A1)') NULL
      READ(99,'(1X,A1)') NULL
      NQS=0
      DO I=1,NGT
        READ(99,*,END=4) Q2T(I),FQ(I)
        NQS=I
      ENDDO
 4    CONTINUE
      CLOSE(UNIT=99)
      IF(NQS.EQ.0) THEN
        WRITE(IW,*) 'HEBORN: I/O error. SCFILE does not exist.'
        STOP 'HEBORN: I/O error. SCFILE does not exist.'
      ENDIF
C
C  ****  Nuclear charge density parameters.
C
      IF(MNUCL.EQ.1) THEN
C  ****  Point nucleus.
        WRITE(IW,1002)
 1002   FORMAT(1X,'#',/1X,'# Nuclear model: point charge')
        U1=0.0D0
        U2=0.0D0
      ELSE IF (MNUCL.EQ.2) THEN
C  ****  Uniform distribution..
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=2.00D0*F2BOHR
        R1=R1*SQRT((1.0D0+2.5D0*(R2/R1)**2)
     1             /(1.0D0+0.75D0*(R2/R1)**2))
        WRITE(IW,1004) R1*A0B
 1004   FORMAT(1X,'#',/1X,'# Nuclear model: uniform spherical d',
     1    'istribution',/1X,'#',16X,'Nuclear radius =',1P,E12.5,
     2    ' cm')
        U1=R1**2
        U2=0.0D0
      ELSE IF(MNUCL.EQ.4.OR.MNUCL.EQ.3) THEN
C  ****  Helm's Uu distribution.
        RNUC=1.070D0*F2BOHR*AW**0.3333333333333333D0
        R1=0.962D0*RNUC+0.435D0*F2BOHR
        R2=2.0D0*F2BOHR
        IF(R2.GT.R1) THEN
          STORE=R1
          R1=R2
          R2=STORE
        ENDIF
        WRITE(IW,1105) R1*A0B,R2*A0B
 1105   FORMAT(1X,'#',/1X,'# Nuclear model: Helm''s Uu distribu',
     1    'tion',/1X,'#',16X,'  Inner radius =',1P,E12.5, ' cm',
     2    /1X,'#',16X,'Skin thickness =',E12.5,' cm')
        U1=R1**2
        U2=R2**2
      ELSE
        WRITE(IW,1003)
 1003   FORMAT(1X,'#',/1X,'# Undefined nuclear charge density model.',
     1    /1X,'# The calculation was aborted by subroutine HEBORN.')
        STOP 'HEBORN: Undefined nuclear charge density model.'
      ENDIF
C
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
 10   CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
C
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
C  ****  Screening correction.
        Q2=4.0D0*RK2*XT(I)
        IF(Q2.LT.Q2T(NQS)) THEN
          CALL FINDI(Q2,Q2T,NQS,J)
          F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
        ELSE
          F=1.0D0
        ENDIF
C  ****  Nuclear form factor.
        QR2=Q2*U1
        QR=SQRT(QR2)
        IF(QR2.LT.1.0D-8) THEN
          FR=1.0D0+QR2*(-0.1D0+QR2*3.5714285714285714D-3)
        ELSE
          FR=3.0D0*(SIN(QR)-QR*COS(QR))/(QR*QR2)
        ENDIF
        QU2=Q2*U2
        QU=SQRT(QU2)
        IF(QU2.LT.1.0D-8) THEN
          FU=1.0D0+QU2*(-0.1D0+QU2*3.5714285714285714D-3)
        ELSE
          FU=3.0D0*(SIN(QU)-QU*COS(QU))/(QU*QU2)
        ENDIF
        FN=FR*FU
C
        RMR=DPWAC(THR)
        DCST(I)=RUTHC*(F*FN)**2*RMR/(1.0D0+Q2)**2
      ENDDO
C
C  ****  Integrated cross sections.
C
      SUM0=0.0D0
      SUM1=0.0D0
      SUM2=0.0D0
      RMUL=0.0D0
      RMUU=1.0D-16
 20   CONTINUE
      TOL=1.0D-9
      MOM=0
      SUMP0=SUMGA(DCSHB,RMUL,RMUU,TOL)
      MOM=1
      SUMP1=SUMGA(DCSHB,RMUL,RMUU,TOL)
      MOM=2
      SUMP2=SUMGA(DCSHB,RMUL,RMUU,TOL)
      SUM0=SUM0+SUMP0
      SUM1=SUM1+SUMP1
      SUM2=SUM2+SUMP2
      RMUL=RMUU
      RMUU=MIN(2.0D0*RMUL,1.0D0)
      IF(RMUL.LT.0.9999999D0) GO TO 20
      ECS0=FOURPI*SUM0
      ECS1=FOURPI*SUM1
      ECS2=FOURPI*SUM2
C
      ECS=ECS0
      TCS1=2.0D0*ECS1
      TCS2=6.0D0*(ECS1-ECS2)
      WRITE(IW,2013) ECS,ECS/A0B2
 2013 FORMAT(1X,'#',/1X,'# Total elastic cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2')
      WRITE(IW,2014) TCS1,TCS1/A0B2
 2014 FORMAT(1X,'# 1st transport cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2')
      WRITE(IW,2015) TCS2,TCS2/A0B2
 2015 FORMAT(1X,'# 2nd transport cross section =',1P,
     1     E12.5,' cm**2 =',E12.5,' a0**2',/1X,'#')
C
      WRITE(IW,'(1X,''#'')')
      WRITE(IW,'(1X,''# Differential cross section'',6X,
     1  ''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',
     1  10X,''DCS'',8X,''Sherman'',7X,''error''/1X,''#  (deg)'',
     2  17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',4X,''function'',
     3  /1X,''#'',71(''-''))')
      DO I=1,NTAB
        SPOL(I)=0.0D0
        ERROR(I)=1.0D-5
        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I),ERROR(I)
      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,3E13.5,2X,E8.1)
C
      RETURN
      END
C  *********************************************************************
      FUNCTION DCSHB(RMU)
C     Mott-Born DCS for elastic scattering of high-energy electrons and
C  positrons by neutral atoms.
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
      PARAMETER (NGT=650)
      COMMON/CDCSHE/Q2T(NGT),FQ(NGT),U1,U2,NQS,MOM
C  ****  Screening correction.
      Q2=4.0D0*RK2*RMU
      IF(Q2.LT.Q2T(NQS)) THEN
        CALL FINDI(Q2,Q2T,NQS,J)
        F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
      ELSE
        F=1.0D0
      ENDIF
C  ****  (nuclear form factor)**2.
      QR2=Q2*U1
      QR=SQRT(QR2)
      IF(QR2.LT.1.0D-8) THEN
        FR=1.0D0+QR2*(-0.1D0+QR2*3.5714285714285714D-3)
      ELSE
        FR=3.0D0*(SIN(QR)-QR*COS(QR))/(QR*QR2)
      ENDIF
      QU2=Q2*U2
      QU=SQRT(QU2)
      IF(QU2.LT.1.0D-8) THEN
        FU=1.0D0+QU2*(-0.1D0+QU2*3.5714285714285714D-3)
      ELSE
        FU=3.0D0*(SIN(QU)-QU*COS(QU))/(QU*QU2)
      ENDIF
      FN=FR*FU
C
      RMR=DPWAC(ACOS(1.0D0-2.0D0*RMU))
      DCSHB=(RUTHC*(F*FN)**2*RMR/(1.0D0+Q2)**2)*RMU**MOM
      RETURN
      END
C  *********************************************************************
C                       FUNCTION VCPOL
C  *********************************************************************
      FUNCTION VCPOL(IELEC,DEN)
C
C     This function gives the correlation potential of an electron
C  (IELEC=-1) or positron (IELEC=+1) in an homogeneous electron gas of
C  density DEN (electrons per unit volume).
C
C  ****  All quantities are in atomic units.
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C
      IF(DEN.LT.1.0D-12) THEN
        VCPOL=0.0D0
        RETURN
      ENDIF
      RS=(3.0D0/(FOURPI*DEN))**3.333333333333D-1
      RSL=LOG(RS)
C
      IF(IELEC.EQ.-1) THEN
C  ****  Electron exchange-correlation potential.
C        Ref:  Padial and Norcross, Phys. Rev. A 29(1984)1742.
C              Perdew and Zunger, Phys. Rev. B 23(1981)5048.
        IF(RS.LT.1.0D0) THEN
          VCPOL=0.0311D0*RSL-0.0584D0+0.00133D0*RS*RSL-0.0084D0*RS
        ELSE
          GAM=-0.1423D0
          BET1=1.0529D0
          BET2=0.3334D0
          RSS=SQRT(RS)
          VCPOL=GAM*(1.0D0+(7.0D0/6.0D0)*BET1*RSS
     1         +(4.0D0/3.0D0)*BET2*RS)/(1.0D0+BET1*RSS+BET2*RS)**2
        ENDIF
      ELSE
C  ****  Positron correlation potential.
C        Ref:  Jain, Phys. Rev. A 41(1990)2437.
        IF(RS.LT.0.302D0) THEN
          VCPOL=(-1.82D0/SQRT(RS))+(0.051D0*RSL-0.115D0)*RSL+1.167D0
        ELSE IF(RS.LT.0.56D0) THEN
          VCPOL=-0.92305D0-0.09098D0/RS**2
        ELSE IF(RS.LT.8.0D0) THEN
          RSD=1.0D0/(RS+2.5D0)
          VCPOL=(-8.7674D0*RS*RSD**3)+(-13.151D0+0.9552D0*RS)*RSD**2
     1         +2.8655D0*RSD-0.6298D0
        ELSE
          VCPOL=-179856.2768D0*3.0D0*DEN**2+186.4207D0*2.0D0*DEN
     1         -0.524D0
        ENDIF
        VCPOL=0.5D0*VCPOL
      ENDIF
      RETURN
      END
C  *********************************************************************
C                        SUBROUTINE XSFEG
C  *********************************************************************
      SUBROUTINE XSFEG(DEN,DELTA,IELEC,EK,MORD,XSEC,IMODE)
C
C     This subroutine computes restricted (W>DELTA) total cross sections
C  for interactions of electrons (IELEC=-1) or positrons (IELEC=+1) with
C  a degenerate free electron gas (per electron in the gas). The DCS
C  is obtained from Lindhard's dielectric function (i.e. within the
C  first Born approximation), with the Ochkur exchange correction for
C  electrons.
C
C  Ref.: F. Salvat, Phys. Rev. A 68 (2003) 012708.
C
C
C  Input arguments:
C     DEN ...... density of the electron gas (electrons per unit
C                volume).
C     DELTA .... energy gap (or minimum energy loss).
C     IELEC .... kind of projectile.
C                =-1, electron; =+1, positron.
C     EK ....... kinetic energy of the projectile.
C     MORD ..... order of the calculated cross section;
C                =0, total cross section,
C                =1, stopping cross section,
C                =2, energy straggling cross section.
C     XSEC ..... total integrated cross section.
C     IMODE .... =1, the complete DCS (for electron-hole and plasmon
C                    excitations) is calculated.
C                =2, the output value XSEC corresponds to electron-hole
C                    excitations (binary collisions) only.
C
C                                  (All quantities in atomic units).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(F2O3=2.0D0/3.0D0,F3O16=3.0D0/16.0D0)
      PARAMETER(PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      PARAMETER(NPM=50)
      COMMON/CXSPL0/Z
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,
     1              AX(NPM),BX(NPM),CX(NPM),DX(NPM),
     2              AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
C  ****  Contributions from electron-hole and plasmon excitations.
      COMMON/XSFEGO/XSEH,XSPL
C
      COMMON/CSUMGA/ERR,IERGA,NCALL  ! Error, code, function calls.
      EXTERNAL PLSTR
C  ****  Plasmon excitation functions are printed if IWR=1.
      IWR=0
      TOL=1.0D-6
C
      IF(MORD.LT.0.OR.MORD.GT.2) THEN
        STOP 'XSFEG: Wrong MORD value.'
      ENDIF
C
C  ****  Constants and energy-independent parameters.
C
      DENM=MAX(DEN,1.0D-5)
      IELPO=IELEC
      EP2=FOURPI*DENM
      EP=SQRT(EP2)
      EF=0.5D0*(3.0D0*PI**2*DENM)**F2O3
      XE=EK/EF
      IF(XE.LT.1.0D-3) THEN
        XSEC=0.0D0
        RETURN
      ENDIF
      SXE=SQRT(XE)
      XP=EP/EF
      CHI2=F3O16*XP*XP
C  ****  Plasmon cutoff momentum.
        ZL=XP-1.0D-3
 1      CONTINUE
        FL=ZL*ZL+CHI2*F1(ZL,4.0D0*ZL*(ZL+1.0D0))
        IF(FL.GT.0.0D0) THEN
          ZL=0.5D0*ZL
          GO TO 1
        ENDIF
        ZU=XP+1.0D-2
 2      CONTINUE
        FU=ZU*ZU+CHI2*F1(ZU,4.0D0*ZU*(ZU+1.0D0))
        IF(FU.LT.0.0D0) THEN
          ZU=ZU+ZU
          GO TO 2
        ENDIF
 3      ZC=0.5D0*(ZL+ZU)
        FT=ZC*ZC+CHI2*F1(ZC,4.0D0*ZC*(ZC+1.0D0))
        IF(FT.LT.0.0D0) THEN
          ZL=ZC
        ELSE
          ZU=ZC
        ENDIF
        IF(ABS(ZL-ZU).GT.1.0D-15*ZC) GO TO 3
        XC=4.0D0*ZC*(ZC+1.0D0)
C
C  ************  Electron-hole contribution.
C
      CALL SEH0(XSEH,DELTA,MORD,IWR)
C
      IF(IMODE.EQ.2) THEN
        XSEC=XSEH
        RETURN
      ENDIF
C
C  ************  Plasmon contribution.
C
      IF(XE.LT.XP) THEN
        XSPL=0.0D0
      ELSE
C  ****  Plasmon line.
        ZPT(1)=0.0D0
        XPT(1)=XP
        FPL(1)=1.0D0
        ZANAL=0.0D0
        XANAL=0.0D0
C  ****  Varying step: 2*DZ for I<NPH.
        NPH=2*NPM/3
        DFZ=0.999999999D0*ZC/DBLE(NPM+NPH-3)
        DO I=2,NPM
          IF(I.LT.NPH) THEN
            Z=ZPT(I-1)+DFZ*2.0D0
          ELSE
            Z=ZPT(I-1)+DFZ
          ENDIF
          IF(Z.GT.0.02D0*ZC) THEN
C  The starting endpoints must be outside the Lindhard continuum.
            XL=MAX(4.0D0*Z*(Z+1.0D0)+1.0D-9,0.9D0*XP)
            XU=1.1D0*XC
 4          X=0.5D0*(XL+XU)
            FT=Z*Z+CHI2*F1(Z,X)
            IF(FT.GT.0.0D0) THEN
              XU=X
            ELSE
              XL=X
            ENDIF
C           WRITE(6,'('' X,FT ='',1P,3E18.11)') X,FT
            IF(FT.GT.1.0D-6) GO TO 4
            IF(ABS(XL-XU).GT.1.0D-13*X) GO TO 4
          ELSE
            X=SQRT(XP**2+(48.0D0/5.0D0)*Z**2+16.0D0*Z**4)
          ENDIF
          XPT(I)=X
          ZPT(I)=Z
        ENDDO
        DO I=2,NPM-1
          Z=ZPT(I)
          XUP=4.0D0*Z*(Z+1.0D0)-1.0D-9
          SUM=SUMGA(PLSTR,1.0D-10,XUP,TOL)
          IF(IERGA.EQ.1) THEN
            WRITE(6,*) 'SUMGA error in XSFEG.'
            STOP 'XSFEG: SUMGA error (1).'
          ENDIF
          FPL(I)=1.0D0-SUM*(6.0D0/(16.0D0*PI))
          XAP=XP+(24.0D0/5.0D0)*Z*Z/XP
          IF(ABS(XAP-XPT(I)).LT.1.0D-3*XPT(I).AND.
     1      FPL(I).GT.0.999D0) THEN
            ZANAL=ZPT(I)
            XANAL=XPT(I)
          ENDIF
        ENDDO
        FPL(NPM)=FPL(NPM-1)+(FPL(NPM-1)-FPL(NPM-2))
     1      *(XPT(NPM)-XPT(NPM-1))/(XPT(NPM-1)-XPT(NPM-2))
C
        IF(IWR.EQ.1) THEN
          OPEN(99,FILE='plasma.dat')
          Z=1.1D0*ZC
          XLOW=MAX(0.0D0,4.0D0*Z*(Z-1.0D0))+1.0D-9
          XUP=4.0D0*Z*(Z+1.0D0)-1.0D-9
          SUM=SUMGA(PLSTR,XLOW,XUP,TOL)
          IF(IERGA.EQ.1) THEN
            WRITE(6,*) 'SUMGA error in XSFEG.'
            STOP 'XSFEG: SUMGA error (2).'
          ENDIF
          BETHE=SUM*(6.0D0/(16.0D0*PI))
          WRITE(99,*) '#  BETHE SUM =',BETHE
          WRITE(99,*) '#  AN. APPROX. VALID FOR Z <',ZANAL
          WRITE(99,*) '#  AN. APPROX. VALID FOR X <',XANAL
          DO I=1,NPM
            Z=ZPT(I)
            XAP=XP+(24.0D0/5.0D0)*Z*Z/XP
            WRITE(99,'(I4,1P,5E14.6)') I,ZPT(I),XPT(I),FPL(I),XAP
          ENDDO
          CLOSE(99)
        ENDIF
C
        CALL SPLINE(ZPT,XPT,AX,BX,CX,DX,0.0D0,0.0D0,NPM)
        CALL SPLINE(ZPT,FPL,AF,BF,CF,DF,0.0D0,0.0D0,NPM)
        CALL SPL0(XSPL,DELTA,MORD)
      ENDIF
C
      XSEC=XSEH+XSPL
      RETURN
      END
C  *********************************************************************
C                       FUNCTION PLSTR
C  *********************************************************************
      FUNCTION PLSTR(X)
C
C     Integrand of the DDCS for a point (Z,X) within the Lindhard
C  continuum.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSPL0/Z
C
      PLSTR=0.0D0
      IF(Z.LT.1.0D-8) RETURN
C  ****  F2 function.
      IF(X.LT.1.0D0) THEN
        IF(X.LT.4.0D0*Z*(1.0D0-Z)) THEN
          F2=PI*X*0.125D0/Z  ! Region a.
        ELSE
          ZIN=1.0D0/Z
          ZM=Z-X*ZIN*0.25D0
          F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
        ENDIF
      ELSE
        ZIN=1.0D0/Z
        ZM=Z-X*ZIN*0.25D0
        F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
      ENDIF
C
      PLSTR=X*Z*Z*F2/((Z*Z+CHI2*F1(Z,X))**2+(CHI2*F2)**2)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SPL0
C  *********************************************************************
      SUBROUTINE SPL0(XSPL,DELTA,MORD)
C
C     Restricted total cross sections for plasmon excitation.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      PARAMETER(NPM=50)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,
     1              AX(NPM),BX(NPM),CX(NPM),DX(NPM),
     2              AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
      COMMON/CSUMGA/ERR,IERGA,NCALL  ! Error, code, function calls.
      EXTERNAL SPL1
C
      TOL=1.0D-6
      XSPL=0.0D0
      IF(XE.LT.XP+1.0D-8) THEN
        WRITE(6,*) 'WARNING: X is less than XP (SPL0).'
        RETURN
      ENDIF
C  ****  Minimum and maximum allowed Z-values.
      I1=0
      IN=NPM
      DO I=2,NPM
        IF(IELPO.EQ.-1) THEN
          XUP=MIN(4.0D0*ZPT(I)*(SXE-ZPT(I)),XE-1.0D0)
        ELSE
          XUP=4.0D0*ZPT(I)*(SXE-ZPT(I))
        ENDIF
        IF(XUP.GT.XPT(I)) THEN
          IF(I1.EQ.0) I1=I
        ENDIF
        IF(XUP.LT.XPT(I).AND.I1.GT.0) THEN
          IN=I
          GO TO 1
        ENDIF
      ENDDO
 1    CONTINUE
      IF(I1.EQ.0) RETURN
C
      I=I1-1
      ZL=ZPT(I)
      ZU=ZPT(I+1)
 2    Z=0.5D0*(ZL+ZU)
      X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
      IF(IELPO.EQ.-1) THEN
        XMIN=MIN(4.0D0*Z*(SXE-Z),XE-1.0D0)
      ELSE
        XMIN=4.0D0*Z*(SXE-Z)
      ENDIF
      IF(XMIN.GT.X) THEN
        ZU=Z
      ELSE
        ZL=Z
      ENDIF
C       WRITE(6,'('' Z1,X-XCON ='',1P,3E18.11)') Z,X-XMIN
      IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 2
      ZMIN=Z
C
      IF(IN.LT.NPM) THEN
        I=IN-1
        ZL=ZPT(I)
        ZU=ZPT(I+1)
 3      Z=0.5D0*(ZL+ZU)
        X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
        IF(IELPO.EQ.-1) THEN
          XMAX=MIN(4.0D0*Z*(SXE-Z),XE-1.0D0)
        ELSE
          XMAX=4.0D0*Z*(SXE-Z)
        ENDIF
        IF(XMAX.LT.X) THEN
          ZU=Z
        ELSE
          ZL=Z
        ENDIF
C         WRITE(6,'('' Z2,X-XCON ='',1P,3E18.11)') Z,X-XMAX
        IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 3
        ZMAX=Z
      ELSE
        XMAX=XC
        ZMAX=ZC
      ENDIF
C
      XDEL=DELTA/EF
      IF(XDEL.GT.XMAX) RETURN
      IF(XDEL.GT.XMIN) THEN
        CALL FINDI(XDEL,XPT,NPM,I)
        ZL=ZPT(I)
        ZU=ZPT(I+1)
 4      Z=0.5D0*(ZL+ZU)
        X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
        IF(XDEL.LT.X) THEN
          ZU=Z
        ELSE
          ZL=Z
        ENDIF
C         WRITE(6,'('' Z1,X-XCON ='',1P,3E18.11)') Z,X-XDEL
        IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 4
        ZMIN=Z
        XMIN=XDEL
      ENDIF
C
      IF(XMIN.GT.XMAX) RETURN
C
C  ****  Soft plasmon excitation.
C
      FACT= 3.0D0/(16.0D0*CHI2)
      SUMP=0.0D0
      IF(XMIN.LT.XANAL.AND.XMAX.GT.XANAL) THEN
        IF(MORD.EQ.0) THEN
          X=XANAL
          S0U=X+(XP/2.0D0)*LOG((X-XP)/(X+XP))
          X=XMIN
          S0L=X+(XP/2.0D0)*LOG((X-XP)/(X+XP))
          SUMP=FACT*(S0U-S0L)
          ZMIN=ZANAL
        ELSE IF(MORD.EQ.1) THEN
          X=XANAL
          S1U=(X**2/2.0D0)+(XP**2/2.0D0)*LOG(X*X-XP*XP)
          X=XMIN
          S1L=(X**2/2.0D0)+(XP**2/2.0D0)*LOG(X*X-XP*XP)
          SUMP=FACT*(S1U-S1L)
          ZMIN=ZANAL
        ELSE IF(MORD.EQ.2) THEN
          X=XANAL
          S2U=(X**3/3.0D0)+XP**2*X+(XP**3/2.0D0)
     1       *LOG((X-XP)/(X+XP))
          X=XMIN
          S2L=(X**3/3.0D0)+XP**2*X+(XP**3/2.0D0)
     1       *LOG((X-XP)/(X+XP))
          SUMP=FACT*(S2U-S2L)
          ZMIN=ZANAL
        ELSE
          STOP 'SPL0: Wrong MORD value.'
        ENDIF
      ENDIF
C
      IF(ZMIN.LT.ZMAX) THEN
        MOM=MORD
        SUM=SUMGA(SPL1,ZMIN,ZMAX,TOL)
        IF(IERGA.NE.0) THEN
          OPEN(99,FILE='plasma.dat')
          DO I=1,NPM
            WRITE(99,'(I4,1P,5E14.6)') I,XPT(I),ZPT(I),FPL(I)
          ENDDO
          CLOSE(99)
          WRITE(6,*) 'Accumulated numerical errors...'
          WRITE(6,*) 'SUMGA error in SPL0.'
          STOP 'SPL0: SUMGA error.'
        ENDIF
      ELSE
        SUM=0.0D0
      ENDIF
      XSPL=(2.0D0*PI/(XE*EF*EF))*(SUM+SUMP)*EF**MORD
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SPL1
C  *********************************************************************
      FUNCTION SPL1(Z)
C
C     DCS for plasmon excitations.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(R96O5=96.0D0/5.0D0,R3O4=3.0D0/4.0D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      PARAMETER(NPM=50)
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,
     1              AX(NPM),BX(NPM),CX(NPM),DX(NPM),
     2              AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
C
      CALL FINDI(Z,ZPT,NPM,I)
      X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
      FP=AF(I)+Z*(BF(I)+Z*(CF(I)+Z*DF(I)))
      SPL1=FP*X**MOM/(Z*X)
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE SEH0
C  *********************************************************************
      SUBROUTINE SEH0(XSEH,DELTA,MORD,IWR)
C
C  Restricted total cross sections for electron-hole excitations.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      PARAMETER(NHM=150)
      DIMENSION XT(NHM),DW(NHM)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
C
      IF(IELPO.EQ.-1) THEN
        XMAX=0.5D0*(XE-1.0D0)
C       XMAX=XE-1.0D-8  !!!! NO EXCHANGE !!!!
      ELSE
        XMAX=XE-1.0D-8
      ENDIF
      XMIN=MAX(DELTA/EF,1.0D-10)
      IF(XMIN.GT.XMAX) THEN
        XSEH=0.0D0
        RETURN
      ENDIF
C
      FACTL=6.0D0/(16.0D0*PI)
      FACTR=0.5D0
      NP=1
      XT(1)=XMIN
      DW(1)=SEH1(XT(1))
      IF(XMIN.LT.1.2D0*XC) THEN
        NS1=2*NHM/3
        DX=(MIN(1.2D0*XC,XMAX)-XMIN)/DBLE(NS1-1)
        DO I=2,NS1
          NP=I
          XT(I)=XT(I-1)+DX
          DW(I)=FACTL*SEH1(XT(I))
        ENDDO
      ENDIF
      IF(XT(NP).LT.XMAX-1.0D-10) THEN
        DFX=EXP(LOG((XMAX)/XT(NP))/DBLE(NHM-NP))
        NP1=NP+1
        ICALC=0
        DO I=NP1,NHM
          NP=I
          XT(I)=XT(I-1)*DFX
          IF(ICALC.EQ.0) THEN
            DW(I)=FACTL*SEH1(XT(I))
            DWA=FACTR/XT(I)**2
            IF(IELPO.EQ.-1) THEN  ! Exchange correction.
              FEXP=XT(I)/(XE-XT(I))
C             FEXP=1.0D0  !!!! NO EXCHANGE !!!!
              DWA=DWA*(1.0D0-FEXP*(1.0D0-FEXP))
            ENDIF
            IF(ABS(DW(I)-DWA).LT.1.0D-4*DWA) ICALC=1
          ELSE
C  ****  High-Z electron-hole excitations. Moller or Rutherford
C        differential cross section.
            DW(I)=FACTR/XT(I)**2
            IF(IELPO.EQ.-1) THEN  ! Exchange correction.
              FEXP=XT(I)/(XE-XT(I))
C             FEXP=1.0D0  !!!! NO EXCHANGE !!!!
              DW(I)=DW(I)*(1.0D0-FEXP*(1.0D0-FEXP))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF(NP.LT.3) THEN
        XSEH=0.0D0
        WRITE(6,*) 'WARNING: NP is too small (SEH0).'
        RETURN
      ENDIF
      DW(NP)=EXP(LOG(DW(NP-1))+LOG(DW(NP-1)/DW(NP-2))
     1      *(XT(NP)-XT(NP-1))/(XT(NP-1)-XT(NP-2)))
C
      IF(IWR.EQ.1) THEN
        OPEN(99,FILE='ehdcs.dat')
        DO I=1,NP
          WRITE(99,'(1X,1P,5E14.6)') XT(I)/XC,DW(I)
        ENDDO
        CLOSE(99)
      ENDIF
C
      FACT=2.0D0*PI/(XE*EF*EF)
      XSEH=FACT*SMOMLL(XT,DW,XT(1),XT(NP),NP,MORD,0)*EF**MORD
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SEH1
C  *********************************************************************
      FUNCTION SEH1(X)
C
C     Integral of the DDCS over Z within the Lindhard continuum.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSEH1/XX
      COMMON/CSUMGA/ERR,IERGA,NCALL  ! Error, code, function calls.
      EXTERNAL SEH2
C
      TOL=1.0D-6
      SEH1=0.0D0
      SXP1=SQRT(X+1.0D0)
      SXEX=SQRT(XE-X)
      ZMIN=MAX(0.5D0*(SXE-SXEX),0.5D0*(SXP1-1.0D0))+1.0D-10
      ZMAX=MIN(0.5D0*(SXE+SXEX),0.5D0*(SXP1+1.0D0))-1.0D-10
      IF(ZMIN.GT.ZMAX) RETURN
C
      XX=X
      IF(ABS(X-XC).LT.2.0D-2*XC) THEN
        ZMINM=ZMIN+1.0D-7*(ZMAX-ZMIN)
        SUM=SUMGA(SEH2,ZMINM,ZMAX,TOL)
        IF(IERGA.EQ.1) THEN
          WRITE(6,*) 'SUMGA error in SEH1 (1).'
          STOP 'SUMGA error in SEH1 (1).'
        ENDIF
        SEH1=SUM
      ELSE
        SUM=SUMGA(SEH2,ZMIN,ZMAX,TOL)
        IF(IERGA.EQ.1) THEN
          WRITE(6,*) 'SUMGA error in SEH1 (2).'
          STOP 'SUMGA error in SEH1 (2).'
        ENDIF
        SEH1=SUM
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SEH2
C  *********************************************************************
      FUNCTION SEH2(Z)
C
C     Integrand of the DDCS for a point (Z,X) within the Lindhard
C  continuum.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSEH1/X
C
      SEH2=0.0D0
      IF(Z.LT.1.0D-8) RETURN
C  ****  F2 function.
      IF(X.LT.1.0D0) THEN
        IF(X.LT.4.0D0*Z*(1.0D0-Z)) THEN
          F2=PI*X*0.125D0/Z  ! Region a.
        ELSE
          ZIN=1.0D0/Z
          ZM=Z-X*ZIN*0.25D0
          F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
        ENDIF
      ELSE
        ZIN=1.0D0/Z
        ZM=Z-X*ZIN*0.25D0
        F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
      ENDIF
C
      SEH2=Z*F2/((Z*Z+CHI2*F1(Z,X))**2+(CHI2*F2)**2)
      IF(IELPO.EQ.-1) THEN  ! Exchange correction for electrons.
        FEXP=4.0D0*Z*Z/(XE-X)
C       FEXP=1.0D0  !!!! NO EXCHANGE !!!!
        SEH2=SEH2*(1.0D0-FEXP*(1.0D0-FEXP))
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       FUNCTION F1
C  *********************************************************************
      FUNCTION F1(Z,X)
C
C     Lindhard's f_1(z,x) function.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      IF(Z.LT.1.0D-5*X) THEN
        R=(Z/X)**2
        F1=-((16.0D0/3.0D0)+(256.0D0/5.0D0)*R)*R
        RETURN
      ENDIF
C
      ZIN=1.0D0/Z
      ZM=Z-X*ZIN*0.25D0
      IF(ABS(ZM).LT.1.0D-8) THEN
        AUX1=2.0D0*ZM-(4.0D0/3.0D0)*ZM**3-(4.0D0/15.0D0)*ZM**5
      ELSE
        ARGL=ABS((1.0D0+ZM)/(1.0D0-ZM))
        IF(ARGL.LT.1.0D-25.OR.ARGL.GT.1.0D25) THEN
          AUX1=0.0D0
        ELSE
          AUX1=(1.0D0-ZM**2)*LOG(ARGL)
        ENDIF
      ENDIF
C
      ZP=Z+X*ZIN*0.25D0
      IF(ABS(ZP).LT.1.0D-8) THEN
        AUX2=2.0D0*ZP-(4.0D0/3.0D0)*ZP**3-(4.0D0/15.0D0)*ZP**5
      ELSE
        ARGL=ABS((1.0D0+ZP)/(1.0D0-ZP))
        IF(ARGL.LT.1.0D-25.OR.ARGL.GT.1.0D25) THEN
          AUX2=0.0D0
        ELSE
          AUX2=(1.0D0-ZP**2)*LOG(ARGL)
        ENDIF
      ENDIF
C
      F1=0.5D0+0.125D0*(AUX1+AUX2)*ZIN
      RETURN
      END


C  *********************************************************************
C                       FUNCTION SUMGA
C  *********************************************************************
      FUNCTION SUMGA(FCT,XL,XU,TOL)
C
C     This function calculates the value SUMGA of the integral of the
C  (external) function FCT over the interval (XL,XU) using the 20-point
C  Gauss-Legendre quadrature method with an adaptive-bisection scheme.
C
C  TOL is the tolerance, i.e. maximum allowed relative error; it should
C  not be less than 1.0D-13. A warning message is written in unit 6 when
C  the required accuracy is not attained. The common block CSUMGA can be
C  used to transfer to the calling program the error ERR, the error flag
C  IERGA, and the number NCALL of calculated function values.
C
C                                    Francesc Salvat. 17 February, 2020.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NP=10, NP2=2*NP, NP4=4*NP, NOIT=256, NCALLT=100000)
      DIMENSION X(NP),W(NP),XM(NP),XP(NP)
      DIMENSION S(NOIT),SN(NOIT),XR(NOIT),XRN(NOIT)
C  Output error codes:
C     IERGA = 0, no problem, the calculation has converged.
C           = 1, too many open subintervals.
C           = 2, too many function calls.
      COMMON/CSUMGA/ERR,IERGA,NCALL  ! Error, code, function calls.
      DATA IWR/6/
C
C  ****  Gauss 20-point quadrature formula.
C  Abscissas.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,
     1       3.7370608871541956D-01,5.1086700195082710D-01,
     2       6.3605368072651503D-01,7.4633190646015079D-01,
     3       8.3911697182221882D-01,9.1223442825132591D-01,
     4       9.6397192727791379D-01,9.9312859918509492D-01/
C  Weights.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,
     1       1.4209610931838205D-01,1.3168863844917663D-01,
     2       1.1819453196151842D-01,1.0193011981724044D-01,
     3       8.3276741576704749D-02,6.2672048334109064D-02,
     4       4.0601429800386941D-02,1.7614007139152118D-02/
C
      DO I=1,NP
        XM(I)=1.0D0-X(I)
        XP(I)=1.0D0+X(I)
      ENDDO
C  ****  Global and partial tolerances.
      TOLG=MIN(MAX(TOL,1.0D-13),1.0D-5)  ! Global tolerance.
      SUMGA=0.0D0
      IERGA=0
      ERRP=0.0D0
C  ****  Straight integration from XL to XU.
      H=XU-XL
      HH=0.5D0*H
      X1=XL
      SP=W(1)*(FCT(X1+XM(1)*HH)+FCT(X1+XP(1)*HH))
      DO J=2,NP
        SP=SP+W(J)*(FCT(X1+XM(J)*HH)+FCT(X1+XP(J)*HH))
      ENDDO
      S(1)=SP*HH
      XR(1)=X1
      NCALL=NP2
      NOI=1
      IDONE=1  ! To prevent a compilation warning.
C
C  ****  Adaptive-bisection scheme.
C
 1    CONTINUE
      H=HH  ! Subinterval length.
      HH=0.5D0*H
      SUMR=0.0D0
      NOIP=NOI
      NOI=0
      ERRPA=ERRP
      ERRP=0.0D0
      DO I=1,NOIP
        SI=S(I)  ! Bisect the I-th open interval.
C
        X1=XR(I)
        SP=W(1)*(FCT(X1+XM(1)*HH)+FCT(X1+XP(1)*HH))
        DO J=2,NP
          SP=SP+W(J)*(FCT(X1+XM(J)*HH)+FCT(X1+XP(J)*HH))
        ENDDO
        S1=SP*HH
C
        X2=X1+H
        SP=W(1)*(FCT(X2+XM(1)*HH)+FCT(X2+XP(1)*HH))
        DO J=2,NP
          SP=SP+W(J)*(FCT(X2+XM(J)*HH)+FCT(X2+XP(J)*HH))
        ENDDO
        S2=SP*HH
C
        IDONE=I
        NCALL=NCALL+NP4
        S12=S1+S2  ! Sum of integrals on the two subintervals.
        IF(ABS(S12-SI).LT.MAX(TOLG*ABS(S12),1.0D-35)) THEN
C  ****  The integral over the parent interval has converged.
          SUMGA=SUMGA+S12
        ELSE
          ERRP=ERRP+ABS(S12-SI)
          SUMR=SUMR+S12
          NOI=NOI+2
          IF(NOI.LT.NOIT) THEN
C  ****  Store open intervals.
            SN(NOI-1)=S1
            XRN(NOI-1)=X1
            SN(NOI)=S2
            XRN(NOI)=X2
          ELSE
C  ****  Too many open intervals.
            IERGA=1
            GO TO 2
          ENDIF
        ENDIF
        IF(NCALL.GT.NCALLT) THEN
C  ****  Too many calls to FCT.
          IERGA=2
          GO TO 2
        ENDIF
      ENDDO
C
C  ****  Analysis of partial results and error control.
C
      IF(IERGA.EQ.0) THEN
        IF(ABS(SUMR).LT.MAX(TOLG*ABS(SUMGA+SUMR),1.0D-35).
     1    OR.NOI.EQ.0) THEN
          ERR=TOLG
          SUMGA=SUMGA+SUMR
          RETURN
        ELSE
          DO I=1,NOI
            S(I)=SN(I)
            XR(I)=XRN(I)
          ENDDO
          GO TO 1
        ENDIF
      ENDIF
C
C  ****  Warning (low accuracy) message.
C
 2    CONTINUE
      IF(IDONE.LT.NOIP) THEN
        DO I=IDONE+1,NOIP
          SUMR=SUMR+S(I)
        ENDDO
        NOI=NOI+(NOIP-IDONE)
      ENDIF
      ERR=ERRPA+TOLG*ABS(SUMGA)
      SUMGA=SUMGA+SUMR
      IF(ERR.LT.10.0D0*TOLG*ABS(SUMGA)) THEN
        IF(ABS(SUMGA).GT.1.0D-16) ERR=ERR/ABS(SUMGA)
        IERGA=0
        RETURN
      ENDIF
      IF(IWR.GT.0) WRITE(IWR,11)
 11   FORMAT(/2X,'>>> SUMGA. Gauss adaptive-bisection quadrature.')
      IF(IWR.GT.0) WRITE(IWR,12) XL,XU,TOL
 12   FORMAT(2X,'XL =',1P,E15.8,', XU =',E15.8,', TOL =',E8.1)
      IF(ABS(SUMGA).GT.1.0D-16) THEN
        ERR=ERR/ABS(SUMGA)
        IF(IWR.GT.0) WRITE(IWR,13) SUMGA,ERR
 13     FORMAT(2X,'SUMGA =',1P,E22.15,', relative error =',E8.1)
      ELSE
        IF(IWR.GT.0) WRITE(IWR,14) SUMGA,ERR
 14     FORMAT(2X,'SUMGA =',1P,E22.15,', absolute error =',E8.1)
      ENDIF
      IF(IWR.GT.0) WRITE(IWR,15) NCALL,NOI,HH
 15   FORMAT(2X,'NCALL =',I6,', open subintervals =',I4,', H =',
     1  1P,E10.3)
      IF(IERGA.EQ.1) THEN
        IF(IWR.GT.0) WRITE(IWR,16)
 16     FORMAT(2X,'IERGA = 1, too many open subintervals.')
      ELSE IF(IERGA.EQ.2) THEN
        IF(IWR.GT.0) WRITE(IWR,17)
 17     FORMAT(2X,'IERGA = 2, too many function calls.')
      ELSE IF(IERGA.EQ.3) THEN
        IF(IWR.GT.0) WRITE(IWR,18)
 18     FORMAT(2X,'IERGA = 3, subintervals are too narrow.')
      ENDIF
      IF(IWR.GT.0) WRITE(IWR,19)
 19   FORMAT(2X,'WARNING: the required accuracy has not been ',
     1  'attained.'/)
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                  *****************************
C                  *  SUBROUTINE PACKAGE DPWA  *
C                  *****************************
C
C
C                                          Francesc Salvat.
C                                          Universitat de Barcelona.
C                                          January 2, 2017
C
C
C     Dirac Partial Wave Analysis for elastic scattering of electrons
C  and positrons by Coulomb fields with short-range central
C  modifications. The radial Dirac equation is solved by using the
C  Fortran subroutine package RADIAL described in
C
C      F. Salvat and J.M. Fernandez-Varea,
C      'RADIAL: a FORTRAN subroutine package for the solution of the
C      radial Schrodinger and Dirac wave equations'.
C      Internal report, University of Barcelona, 2017.
C
C     The calling sequence from the main program is:
C
C****   CALL DPWA0(EV,NDELTA,ISCH)
C
C     This subroutine determines the phase shifts. It acts as the
C  initialization routine for the evaluation of scattering amplitudes
C  and differential cross sections.
C
C****   CALL DPWA(TH,CF,CG,DCS,SPL,ERRF,ERRG)
C
C     Subroutine DPWA gives elastic scattering functions at the
C  scattering angle TH obtained from the phase shifts calculated
C  previously by subroutine DPWA0.
C
C
C            ****  All I/O energies and lengths in eV and cm, resp.
C
C  *********************************************************************
C                      SUBROUTINE DPWA0
C  *********************************************************************
      SUBROUTINE DPWA0(EV,NDELTA,ISCH)
C
C     This subroutine computes Dirac phase shifts, differential cross
C  sections and scattering amplitudes for elastic scattering of
C  electrons in central fields.
C
C  Input arguments:
C     EV ....... effective kinetic energy of the projectile (eV).
C     NDELTA ... number of required phase shifts (LT.25000).
C     ISCH ..... =1: all phase shifts are computed by solving the radial
C                    equation.
C                =2: only phase shifts of selected orders are computed
C                    from the solution of the radial equation, the
C                    others are obtained by lin-log natural cubic spline
C                    interpolation. For high energies, ISCH=2 leads to a
C                    considerable reduction of the calculation time.
C
C  Input (through the common block /FIELD/):
C     R(I) .... radial grid points (radii in increasing order). The
C               first point in the grid must be the origin, i.e. R(1)=0.
C               Repeated values are interpreted as discontinuities.
C     RV(I).... R(I) times the potential energy at R=R(I). The last
C               component, RV(NP), is assumed to be equal to the
C               asymptotic value.
C     NP ...... number of input grid points.
C
C *** NOTE: The radii and potential values, R(I) and RV(I), are in
C           atomic units.
C
C  Output (through the common block /DCSTAB/):
C     ECS ........ total cross section (cm**2)
C                    (only for finite range fields).
C     TCS1 ....... 1st transport cross section (cm**2)
C                    (only for finite range fields).
C     TCS2 ....... 2nd transport cross section (cm**2)
C                    (only for finite range fields).
C     TH(I) ...... scattering angles (in deg)
C     XT(I) ...... values of (1-COS(TH(I)))/2.0D0.
C     DCST(I) .... differential cross section per unit solid angle at
C                    TH(I) (cm**2/sr).
C     ERROR(I) ... estimated relative uncertainty of the computed DCS
C                    value.
C     NTAB ....... number of angles in the table.
C
C  NOTE: The values of ECS, TCS1 and TCS2 are computed from the DCS
C  table. This introduces a certain error (of the order of 0.01 per
C  cent) but ensures consistency of multiple scattering simulations
C  using the DCS table.
C
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C  ****  Input-output.
      COMMON/FIELD/R(NDIM),RV(NDIM),NP
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
C  ****  Link with the RADIAL package.
      COMMON/RADWF/RRR(NDIM),P(NDIM),Q(NDIM),NRT,ILAST,IER
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C  ****  Phase shifts and partial wave series coefficients.
      PARAMETER (NPC=1500,NDM=25000)
      DIMENSION XL(NDM),DPI(NDM),DMI(NDM)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      DIMENSION X(NDM),Y(NDM),SA(NDM),SB(NDM),SC(NDM),SD(NDM)
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),
     1              CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
C
      CI=DCMPLX(0.0D0,1.0D0)
      ISUMP=0
C
      OPEN(98,FILE='dpwa.dat')
      WRITE(98,2000)
 2000 FORMAT(//2X,'**** Partial wave analysis (DPWA0) ',
     1  42('*')/)
C
      NDELT=NDELTA
      IF(NDELT.GT.NDM) THEN
        WRITE(98,2001)
 2001   FORMAT(/2X,'Warning: NDELTA is too large.')
        NDELT=NDM
      ENDIF
      IF(NDELT.LT.6) NDELT=6
      EPS=1.0D-15
      EPSCUT=1.0D-9
C
      E=EV/HREV
C
C  ****  Initialization of the RADIAL package.
C
      CALL VINT(R,RV,NP)
      ZINF=RVG(NVT)
      IF(ABS(ZINF).GT.1.0D-10) THEN
        IZINF=1
        CALL DPWAC0(ZINF,EV)
        IF(NDELT.GT.NPC) NDELT=NPC
      ELSE
        IZINF=0
      ENDIF
C
      WRITE(98,2002) EV
 2002 FORMAT(/2X,'Kinetic energy =',1P,E12.5,' eV')
      IF(E.LT.0.0D0) THEN
        WRITE(6,2002) EV
        WRITE(98,2003)
        WRITE(6,2003)
 2003   FORMAT(//2X,'Negative energy. Stop.')
        STOP 'DPWA0: Negative energy.'
      ENDIF
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
      IF(IZINF.EQ.1) WRITE(98,2004)
 2004 FORMAT(/2X,'Only inner phase shifts are tabulated')
      WRITE(98,2005)
 2005 FORMAT(/6X,'L',7X,'Phase(spin UP)',5X,'Phase(spin DOWN)',
     1       /2X,47('-'))
      ISCH0=ISCH
      IF(ISCH0.EQ.2.AND.EV.GT.1000.0D0) GO TO 1
C
C  ****  ISCH0=1, all phase shifts are computed by solving the radial
C        equation.
C
      L=0
      CALL DFREE(E,EPS,PHP,-1,0)
      IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (1).'
      PHM=0.0D0
      WRITE(98,2006) L,PHP,PHM
      WRITE(6,2006) L,PHP,PHM
 2006 FORMAT(3X,I5,4X,1P,E16.8,4X,E16.8)
      DP(1)=PHP
      DM(1)=0.0D0
      NPH=1
C
      IFIRST=2
 33   CONTINUE
      ISUMP=1
      TST=0.0D0
      DO I=IFIRST,NDELT
        L=I-1
        CALL DFREE(E,EPS,PHP,-L-1,0)
        IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (2).'
        CALL DFREE(E,EPS,PHM,L,0)
        IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (3).'
        DP(I)=PHP
        DM(I)=PHM
        TST=MAX(ABS(PHP),ABS(PHM),ABS(DP(I-1)))
        NPH=I
        WRITE(98,2006) L,PHP,PHM
        WRITE(6,2006) L,PHP,PHM
        IF(TST.LT.EPSCUT.AND.L.GT.10) GO TO 6
C  ****  When the last phase shift (spin up) differs in more than 20 per
C  cent from the quadratic extrapolation, accumulated roundoff errors
C  may be important and the calculation of phase shifts is discontinued.
        IF(I.GT.500) THEN
          DPEXT=DP(I-3)+3.0D0*(DP(I-1)-DP(I-2))
          DPMAX=MAX(ABS(DP(I-3)),ABS(DP(I-2)),ABS(DP(I-1)),ABS(DP(I)))
          IF(ABS(DP(I)-DPEXT).GT.0.20D0*DPMAX) THEN
            NPH=I-1
            WRITE(98,2107)
            WRITE(6,2107)
 2107 FORMAT(/2X,'WARNING: Possible accumulation of round-off errors.')
            GO TO 6
          ENDIF
        ENDIF
      ENDDO
      WRITE(98,2007) TST
      WRITE(6,2007) TST
 2007 FORMAT(/2X,'WARNING: TST =',1P,E11.4,'. Check convergence.')
      GO TO 6
C
C  ****  ISCH0=2, only inner phase shifts of orders L in a given grid
C        are computed from the solution of the radial equation. Phase
C        shifts of orders not included in this grid are obtained by
C        lin-log cubic spline interpolation.
C          The adopted grid is: 0(1)100(5)300(10) ...
C
C        This is a somewhat risky procedure, which is based on the
C        observed variation of the calculated phase shifts with L for
C        atomic scattering fields. When a change of sign is found, all
C        the phases are recalculated.
C
 1    CONTINUE
      L=0
      CALL DFREE(E,EPS,PHP,-1,0)
      IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (1).'
      PHM=0.0D0
      WRITE(98,2006) L,PHP,PHM
      WRITE(6,2006) L,PHP,PHM
      DP(1)=PHP
      DM(1)=0.0D0
C
      LMAX=NDELT-1
      IND=0
      IADD=1
      NADD=0
      LPP=1
 2    CONTINUE
      L=LPP
      CALL DFREE(E,EPS,PHP,-L-1,0)
      IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (2).'
      CALL DFREE(E,EPS,PHM,L,0)
      IF(IER.NE.0) STOP 'DPWA0: Error in DFREE (3).'
      WRITE(6,2006) L,PHP,PHM
C
      DP(L+1)=PHP
      DM(L+1)=PHM
C
      IF(L.LT.95) THEN
        WRITE(98,2006) L,PHP,PHM
      ELSE
        IND=IND+1
        NADD=NADD+1
        XL(IND)=L
        DPI(IND)=PHP
        DMI(IND)=PHM
        IF(IND.GT.1) THEN
          TST1=DPI(IND)*DPI(IND-1)
          TST2=DMI(IND)*DMI(IND-1)
          IF(TST1.LT.0.0D0.OR.TST2.LT.0.0D0) THEN
            IF(L.LT.600) THEN
              ISCH0=1
              IFIRST=MIN(L,94)
              GO TO 33
            ELSE
              IND=IND-1
              L=XL(IND)+0.5D0
              GO TO 3
            ENDIF
          ENDIF
C
          IF(L.GT.600.AND.NADD.GT.3) THEN
            I=IND
            DPEXT=DPI(I-3)+3.0D0*(DPI(I-1)-DPI(I-2))
            DPMAX=MAX(ABS(DPI(I-3)),ABS(DPI(I-2)),ABS(DPI(I-1)),
     1            ABS(DPI(I)))
            IF(ABS(DPI(I)-DPEXT).GT.0.20D0*DPMAX) THEN
              IND=I-1
              L=XL(IND)+0.5D0
              WRITE(98,2107)
              WRITE(6,2107)
              GO TO 3
            ENDIF
            DMEXT=DMI(I-3)+3.0D0*(DMI(I-1)-DMI(I-2))
            DMMAX=MAX(ABS(DMI(I-3)),ABS(DMI(I-2)),ABS(DMI(I-1)),
     1            ABS(DMI(I)))
            IF(ABS(DMI(I)-DMEXT).GT.0.20D0*DMMAX) THEN
              IND=I-1
              L=XL(IND)+0.5D0
              WRITE(98,2107)
              WRITE(6,2107)
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      TST=MAX(ABS(PHP),ABS(PHM))
      IF(TST.LT.EPSCUT) THEN
        IF(L.LT.100) THEN
          NPH=L+1
          GO TO 6
        ELSE
          GO TO 3
        ENDIF
      ENDIF
      IF(L.GE.LMAX) GO TO 3
C
      IADDO=IADD
      IF(L.GT.99) IADD=5
      IF(L.GT.299) IADD=10
      IF(L.GT.599) IADD=20
      IF(L.GT.1199) IADD=50
      IF(L.GT.2999) IADD=100
      IF(L.GT.9999) IADD=250
      IF(IADD.NE.IADDO) NADD=0
      LPP=L+IADD
      IF(LPP.GT.LMAX) LPP=LMAX
      GO TO 2
C
C  ****  Check consistency of sparsely tabulated phase shifts.
C        A discontinuity larger than 0.25*PI is considered as
C        a symptom of numerical inconsistencies.
C
 3    CONTINUE
      IF(IND.LT.5.OR.IADD.EQ.1) GO TO 6
      NPH=XL(IND)+1.5D0
      TST=0.0D0
      DO I=1,IND
        WRITE(98,2008) INT(XL(I)+0.5D0),DPI(I),DMI(I)
 2008   FORMAT(3X,I5,4X,1P,E16.8,4X,E16.8,'  i')
        IF(I.GT.1) THEN
          TST=MAX(TST,ABS(DPI(I)-DPI(I-1)),ABS(DMI(I)-DMI(I-1)))
        ENDIF
      ENDDO
      IF(TST.GT.0.25D0*PI) THEN
        WRITE(98,2009)
        WRITE(6,2009)
 2009   FORMAT(/2X,'ERROR: Directly computed phase shifts show',
     1    ' large discontinuities.')
        STOP 'DPWA0: Phase shifts do not vary continuously with L.'
      ENDIF
C
C  ****  Interpolated phase shifts (lin-log cubic spline).
C
      IF(DPI(4).GT.0.0D0) THEN
        ITRAN=+1
      ELSE
        ITRAN=-1
      ENDIF
      JT=0
      DO I=1,IND
        IF(DPI(I)*ITRAN.LT.0.0D0.OR.ABS(DPI(I)).LT.1.0D-12) THEN
          GO TO 4
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=LOG(ABS(DPI(I)))
        ENDIF
      ENDDO
 4    CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DP(I)=ITRAN*EXP(SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J))))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DP(I)=0.0D0
        ENDDO
      ENDIF
C
      IF(DMI(4).GT.0.0D0) THEN
        ITRAN=+1
      ELSE
        ITRAN=-1
      ENDIF
      JT=0
      DO I=1,IND
        IF(DMI(I)*ITRAN.LT.0.0D0.OR.ABS(DMI(I)).LT.1.0D-12) THEN
          GO TO 5
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=LOG(ABS(DMI(I)))
        ENDIF
      ENDDO
 5    CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DM(I)=ITRAN*EXP(SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J))))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DM(I)=0.0D0
        ENDDO
      ENDIF
C
      TST=MAX(ABS(DP(NPH)),ABS(DM(NPH)))
      IF(TST.GT.10.0D0*EPSCUT) THEN
        WRITE(98,2007) TST
        WRITE(6,2007) TST
      ENDIF
C
C  ************  Coefficients in the partial-wave expansion.
C
 6    CONTINUE
      CFACT=1.0D0/(2.0D0*CI*RK)
      IF(IZINF.EQ.1) THEN
        CXP=CDEXP(2.0D0*CI*DP(1))
        CXPC=CDEXP(2.0D0*CI*DPC(1))
        CFL(1)=CXPC*(CXP-1)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXP=CDEXP(2.0D0*CI*DP(I))
          CXM=CDEXP(2.0D0*CI*DM(I))
          CXPC=CDEXP(2.0D0*CI*DPC(I))
          CXMC=CDEXP(2.0D0*CI*DMC(I))
          CFL(I)=((L+1)*CXPC*(CXP-1)+L*CXMC*(CXM-1))*CFACT
          CGL(I)=(CXMC*(CXM-1)-CXPC*(CXP-1))*CFACT
        ENDDO
      ELSE
        CXP=CDEXP(2*CI*DP(1))
        CFL(1)=(CXP-1.0D0)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXP=CDEXP(2.0D0*CI*DP(I))
          CXM=CDEXP(2.0D0*CI*DM(I))
          CFL(I)=((L+1)*(CXP-1)+L*(CXM-1))*CFACT
          CGL(I)=CXM*(1.0D0-CDEXP(2.0D0*CI*(DP(I)-DM(I))))*CFACT
        ENDDO
      ENDIF
C
C  ****  Reduced series (two iterations).
C
      IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
        DO I=1,NPH
          CFM(I)=CFL(I)
          CGM(I)=CGL(I)
        ENDDO
C
        NPHM=NPH
        DO 7 NTR=1,2
          NPHM=NPHM-1
          CFC=0.0D0
          CFP=CFM(1)
          CGC=0.0D0
          CGP=CGM(1)
          DO I=1,NPHM
            RL=I-1
            CFA=CFC
            CFC=CFP
            CFP=CFM(I+1)
            CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
            CGA=CGC
            CGC=CGP
            CGP=CGM(I+1)
            CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
          ENDDO
 7      CONTINUE
      ENDIF
C
*     OPEN(99, file='pwa-coefs0.dat')
*     WRITE(99,'(A)') '# Coefficients in the partial-wave expansions'
*     WRITE(99,'(A,I6)') '# NPH =',NPH
*     WRITE(99,'(A)') '# L, CFL(L),CGL(L) all in a.u.'
*     DO I=1,NPH
*       WRITE(99,'(I5,1P,2(2X,E14.6,E14.6))') I-1,CFL(I),CGL(I)
*     ENDDO
*     CLOSE(99)
*     IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
*       OPEN(99, file='pwa-coefs2.dat')
*       WRITE(99,'(A)') '# Coefficients in the partial-wave expansions'
*       WRITE(99,'(A)') '#   after the reduced-series transformation'
*       WRITE(99,'(A,I6)') '# NPHM =',NPHM
*       WRITE(99,'(A)') '# L, CFM(L),CGM(L) all in a.u.'
*       DO I=1,NPHM
*         WRITE(99,'(I5,1P,2(2X,E14.6,E14.6))') I-1,CFM(I),CGM(I)
*       ENDDO
*       CLOSE(99)
*     ENDIF
C
C  ****  Scattering amplitudes and DCS.
C
      WRITE(98,2010)
 2010 FORMAT(//2X,'*** Scattering amplitudes and different',
     1  'ial cross section ***')
      WRITE(98,2011)
 2011 FORMAT(/4X,'Angle',6X,'DCS',7X,'Asymmetry',4X,'Direct amplitu',
     1  'de',7X,'Spin-flip amplitude',5X,'error',/4X,'(deg)',3X,
     2  '(cm**2/sr)',22X,'(cm)',20X,'(cm)',/2X,91('-'))
C
C  ****  Angular grid (TH in deg).
C
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
 10   CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(I.GT.NGT) STOP 'DPWA0. The NGT parameter is too small.'
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
C
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
        CALL DPWA(THR,CF,CG,DCS,SPL,ERRF,ERRG)
        IF(MAX(ERRF,ERRG).GT.0.95D0) THEN
          ERR=1.0D0
        ELSE
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/MAX(DCS,1.0D-45)
        ENDIF
        DCST(I)=DCS
        ERROR(I)=MAX(ERR,1.0D-7)
        SPOL(I)=SPL
        WRITE(98,2012) TH(I),DCST(I),SPOL(I),CF,CG,ERROR(I)
 2012   FORMAT(1X,1P,E10.3,E12.5,1X,E10.3,2(1X,'(',E10.3,',',
     1    E10.3,')'),E10.2)
      ENDDO
C
C  ************  Total and momentum transfer cross sections.
C                Convergence test (only for finite range fields).
C
      IF(IZINF.EQ.0) THEN
        INC=5
        IF(ISUMP.EQ.1) INC=1
        TST1=0.0D0
        TST2=0.0D0
        ECS=4.0D0*PI*CFL(1)*DCONJG(CFL(1))
        TCS=0.0D0
        ECSO=ECS
        TCSO=TCS
        DO I=2,NPH
          L=I-1
          RL=L
          DECS=CFL(I)*DCONJG(CFL(I))+RL*(L+1)*CGL(I)*DCONJG(CGL(I))
          DECS=4.0D0*PI*DECS/(L+L+1)
          DTCS=CFL(L)*DCONJG(CFL(I))+DCONJG(CFL(L))*CFL(I)
     1        +(L-1)*(RL+1)*(CGL(L)*DCONJG(CGL(I))
     2        +DCONJG(CGL(L))*CGL(I))
          DTCS=4.0D0*PI*DTCS*L/((RL+L-1)*(L+L+1))
          ECS=ECS+DECS
          TCS=TCS+DTCS
C  ****  Convergence test.
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST1=ABS(ECS-ECSO)/(ABS(ECS)+1.0D-35)
            TST2=ABS(TCS-TCSO)/(ABS(TCS)+1.0D-35)
            ECSO=ECS
            TCSO=TCS
          ENDIF
        ENDDO
        TST=MAX(TST1,TST2)
        TCS=ECS-TCS
        IF(TST.GT.1.0D-5.AND.NPH.GT.40) THEN
          WRITE(98,2007) TST
          WRITE(6,2007) TST
        ENDIF
        ECS=ECS*A0B2
        TCS=TCS*A0B2
C
C  ****  ECS and TCSs are evaluated from the DCS table.
C
        ECS0=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,0,0)
        ECS1=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,1,0)
        ECS2=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,2,0)
        TST1=ABS(ECS-ECS0)/(ABS(ECS)+1.0D-35)
        WRITE(98,2013) ECS,ECS0,TST1
        WRITE(6,2013) ECS,ECS0,TST1
 2013   FORMAT(/2X,'Total elastic cross section =',1P,E13.6,' cm**2',
     1         /2X,'             from DCS table =',E13.6,
     2         '  (rel. dif. =',E9.2,')')
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TST2=ABS(TCS-TCS1)/(ABS(TCS)+1.0D-35)
        WRITE(98,2014) TCS,TCS1,TST2
        WRITE(6,2014) TCS,TCS1,TST2
 2014   FORMAT(/2X,'1st transport cross section =',1P,E13.6,' cm**2',
     1         /2X,'             from DCS table =',E13.6,
     2         '  (rel. dif. =',E9.2,')')
        WRITE(98,2015) TCS2
        WRITE(6,2015) TCS2
 2015   FORMAT(/2X,'2nd transport cross section =',1P,E13.6,' cm**2')
        TST=MAX(TST1,TST2)
        IF(TST.GT.2.0D-3) THEN
          WRITE(98,2016)
          WRITE(6,2016)
        ENDIF
      ENDIF
 2016 FORMAT(/2X,'WARNING: relative differences are too large.',
     1       /11X,'The dcs table is not consistent.')
C
      WRITE(98,2017)
 2017 FORMAT(/2X,'**** DPWA0 ended ',60('*')/)
      CLOSE(UNIT=98)
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DPWA
C  *********************************************************************
      SUBROUTINE DPWA(TH,CF,CG,DCS,SPL,ERRF,ERRG)
C
C    This subroutine gives various elastic scattering functions at the
C  scattering angle TH (in radians) computed from Dirac phase shifts.
C  It should be previously initialized by calling subroutine DPWA0.
C
C  Input argument:
C     TH ....... scattering angle (in rad)
C
C  Output arguments:
C     CF ....... F scattering amplitude (cm).
C     CG ....... G scattering amplitude (cm).
C     DCS ...... differential cross section per unit solid angle for
C                unpolarized beams.
C     SPL ...... asymmetry function.
C     ERRF ..... relative uncertainty of CF.
C     ERRG ..... relative uncertainty of CG.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (NPC=1500,NDM=25000,TOL=5.0D-8)
C  ****  Phase shifts and partial wave series coefficients.
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),
     1              CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
C
      X=COS(TH)
      Y=SIN(TH)
      TST=1.0D35
      CTERM=0.0D0
C
C  ************  Reduced series method. Only when TH is greater
C                than 0.5 deg and NPH.ge.250.
C
      IF(TH.LT.0.008726D0.OR.NPH.LT.250.OR.ISUMP.EQ.1) THEN
        CFO=0.0D0
        CGO=0.0D0
        ERRFO=1.0D10
        ERRGO=1.0D10
        GO TO 10
      ENDIF
      FACT=1.0D0/(1.0D0-X)**2
C
C  ****  F scattering amplitude.
C
      P2=1.0D0
      P3=X
      CFS=CFM(1)
      CFSO=CFS
      CFS=CFS+CFM(2)*P3
      DO I=3,NPHM
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CTERM=CFM(I)*P3
        CFS=CFS+CTERM
C  ****  Convergence test.
        IF(L.LT.149) THEN
          INC=1
        ELSE IF(L.LT.999) THEN
          INC=5
        ELSE
          INC=15
        ENDIF
        ITW=L-(L/INC)*INC
        IF(ITW.EQ.0) THEN
          TST=CDABS(CFS-CFSO)/MAX(CDABS(CFS),1.0D-45)
          CFSO=CFS
        ENDIF
      ENDDO
      CF=FACT*CFS
      ERRF=TST
C
C  ****  G scattering amplitude.
C
      IF(Y.LT.1.0D-30) THEN
        CG=0.0D0
        ERRG=0.0D0
      ELSE
        P2=1.0D0
        P3=3*X
        CGS=CGM(2)
        CGSO=CGS
        CGS=CGS+CGM(3)*P3
        DO I=4,NPHM
          L=I-1
          P1=P2
          P2=P3
          P3=((L+L-1)*X*P2-L*P1)/(L-1)
          CTERM=CGM(I)*P3
          CGS=CGS+CTERM
C  ****  Convergence test.
          IF(L.LT.149) THEN
            INC=1
          ELSE IF(L.LT.999) THEN
            INC=5
          ELSE
            INC=15
          ENDIF
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST=CDABS(CGS-CGSO)/MAX(CDABS(CGS),1.0D-45)
            CGSO=CGS
          ENDIF
        ENDDO
        CG=FACT*Y*CGS
        ERRG=TST
      ENDIF
C
      IF(ERRF.LT.TOL.AND.ERRG.LT.TOL) GO TO 20
      CFO=CF
      ERRFO=ERRF
      CGO=CG
      ERRGO=ERRG
C
C  ************  TH smaller than 0.5 deg or NPH.LT.250 or ISUMP=1.
C
 10   CONTINUE
C  ****  If IZINF=1, scattering functions are calculated only for
C        TH larger than 0.5 deg.
      IF(IZINF.EQ.1.AND.TH.LT.0.008726D0) THEN
        CF=0.0D0
        CG=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        DCS=1.0D-45
        SPL=0.0D0
        RETURN
      ENDIF
C
C  ****  F scattering amplitude.
C
      P2=1.0D0
      P3=X
      CFS=CFL(1)
      CFSO=CFS
      CFS=CFS+CFL(2)*P3
      DO I=3,NPH
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CTERM=CFL(I)*P3
        CFS=CFS+CTERM
C  ****  Convergence test.
        IF(L.LT.149) THEN
          INC=1
        ELSE IF(L.LT.999) THEN
          INC=5
        ELSE
          INC=15
        ENDIF
        ITW=L-(L/INC)*INC
        IF(ITW.EQ.0) THEN
          TST=CDABS(CFS-CFSO)/MAX(CDABS(CFS),1.0D-45)
          CFSO=CFS
        ENDIF
      ENDDO
      CF=CFS
      ERRF=TST
C
C  ****  G scattering amplitude.
C
      IF(Y.LT.1.0D-30) THEN
        CG=0.0D0
        ERRG=0.0D0
      ELSE
        P2=1.0D0
        P3=3*X
        CGS=CGL(2)
        CGSO=CGS
        CGS=CGS+CGL(3)*P3
        DO I=4,NPH
          L=I-1
          P1=P2
          P2=P3
          P3=((L+L-1)*X*P2-L*P1)/(L-1)
          CTERM=CGL(I)*P3
          CGS=CGS+CTERM
C  ****  Convergence test.
          IF(L.LT.149) THEN
            INC=1
          ELSE IF(L.LT.999) THEN
            INC=5
          ELSE
            INC=15
          ENDIF
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST=CDABS(CGS-CGSO)/MAX(CDABS(CGS),1.0D-45)
            CGSO=CGS
          ENDIF
        ENDDO
        CG=Y*CGS
        ERRG=TST
      ENDIF
C  ****  The following four sentences are introduced to prevent abnormal
C        termination of the calculation when the number of (inner) phase
C        shifts is small. This solves the problem found by M. Berger.
      IF(NPH.LT.20.AND.CDABS(CTERM).LT.TOL) THEN
        ERRF=0.0D0
        ERRG=0.0D0
      ENDIF
C
C  ****  Select the most accurate method.
C
      IF(ERRFO.LT.ERRF) THEN
        CF=CFO
        ERRF=ERRFO
      ENDIF
      IF(ERRGO.LT.ERRG) THEN
        CG=CGO
        ERRG=ERRGO
      ENDIF
C
C  ****  Differential cross section (unpolarized beam).
C
 20   CONTINUE
      CF=CF*A0B
      CG=CG*A0B
      IF(IZINF.EQ.1) THEN
        XAUX=DPWAC(TH)
        CFC=CFC*A0B
        CGC=CGC*A0B
        DCSM=CDABS(CFC)**2+CDABS(CGC)**2
        CF=CF+CFC
        CG=CG+CGC
        ACF=CDABS(CF)**2
        ACG=CDABS(CG)**2
        DCS=ACF+ACG
C  ****  Scattering amplitudes that are much smaller than the Coulomb
C        ones may not be correct due to rounding off.
C        (Modified Coulomb fields only).
        IF(DCS.LT.1.0D-10*DCSM.OR.ERRFC+ERRGC.GT.1.0D0.OR.
     1    TH.LT.0.008726D0) THEN
          CF=0.0D0
          CG=0.0D0
          ERRF=1.0D0
          ERRG=1.0D0
          DCS=1.0D-45
          SPL=0.0D0
          RETURN
        ENDIF
        ERRF=ERRF+ERRFC
        ERRG=ERRG+ERRGC
      ELSE
        ACF=CDABS(CF)**2
        ACG=CDABS(CG)**2
        DCS=ACF+ACG
      ENDIF
C
      ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/MAX(DCS,1.0D-45)
      IF(ERR.GT.0.10D0) THEN
        CF=0.0D0
        CG=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        DCS=1.0D-45
        SPL=0.0D0
        RETURN
      ENDIF
C
C  ****  Asymmetry function.
C
      CSPL1=DCMPLX(0.0D0,1.0D0)*CF*DCONJG(CG)
      CSPL2=DCMPLX(0.0D0,1.0D0)*CG*DCONJG(CF)
      TST=CDABS(CSPL1-CSPL2)/MAX(CDABS(CSPL1),1.0D-45)
      IF(TST.GT.1.0D-3.AND.ERR.LT.0.01D0) THEN
        SPL=(CSPL1-CSPL2)/DCS
      ELSE
        SPL=0.0D0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE DPWAC0
C  *********************************************************************
      SUBROUTINE DPWAC0(ZZP,EV)
C
C     This subroutine computes Coulomb phase shifts and initializes the
C  calculation of the Mott differential cross section for electron or
C  positron elastic scattering by a bare point nucleus.
C
C  Input:
C     ZZP....... product of nuclear and projectile charges, that is, R
C                times the interaction energy at the distance R.
C                Negative for electrons, positive for positrons.
C     EV ....... kinetic energy of the projectile (eV).
C
C  After calling DPWAC0, the function DPWAC(TH) delivers the ratio
C  (Mott DCS / Rutherford DCS) for the scattering angle TH (rad).
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
C
      PARAMETER (SL2=SL*SL)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
C  ****  Phase shifts and partial-wave series coefficients.
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CF,CG,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
C
      CI=DCMPLX(0.0D0,1.0D0)
C
C  ************  Coulomb phase shifts.
C
      E=EV/HREV
      RUTHC=(A0B*2.0D0*ZZP*(1.0D0+E/SL2))**2
      PC=SQRT(E*(E+2.0D0*SL*SL))
      RK=PC/SL
      RK2=RK*RK
      ZETA=ZZP/SL
      W=E+SL2
      ETA=ZETA*W/PC
      RNUR=ZETA*(W+SL2)
C  ****  Negative kappa.
      DO I=1,NPC
        L=I-1
        K=-L-1
        RLAMB=SQRT(K*K-ZETA*ZETA)
        RNUI=-1.0D0*(K+RLAMB)*PC
        RNU=ATAN2(RNUI,RNUR)
        DELTAC=-CI*CLGAM(RLAMB+CI*ETA)
        DPC(I)=RNU-(RLAMB-(L+1))*PIH+DELTAC
      ENDDO
C  ****  Positive kappa.
      DMC(1)=0.0D0
      DO I=2,NPC
        L=I-1
        K=L
        RLAMB=SQRT(K*K-ZETA*ZETA)
        RNUI=-1.0D0*(K+RLAMB)*PC
        RNU=ATAN2(RNUI,RNUR)
        DELTAC=-CI*CLGAM(RLAMB+CI*ETA)
        DMC(I)=RNU-(RLAMB-(L+1))*PIH+DELTAC
      ENDDO
C
C  ****  Prints Coulomb phase shifts in file CPHASES.DAT.
C
*     OPEN(99,FILE='cphases.dat')
*     WRITE(99,1000)
*1000 FORMAT(2X,'# COULOMB PHASE SHIFTS',/2X,'#')
*     WRITE(99,1001) ZZP,E
*1001 FORMAT(2X,'# Z =',1P,E12.5,5X,'KINETIC ENERGY =',E12.5,/2X,'#')
*     WRITE(99,1002)
*1002 FORMAT(2X,'#   L',7X,'PHASE(SPIN UP)',4X,
*    1  'PHASE(SPIN DOWN)',/2X,'# ',45('-'))
*     DO I=1,NPC
*       DPI=DPC(I)
*       TT=ABS(DPI)
*       IF(TT.GT.PIH) DPI=DPI*(1.0D0-PI/TT)
*       DMI=DMC(I)
*       TT=ABS(DMI)
*       IF(TT.GT.PIH) DMI=DMI*(1.0D0-PI/TT)
*       WRITE(99,1003) I,DPI,DMI
*1003   FORMAT(3X,I5,4X,1P,E16.8,4X,E16.8)
*     ENDDO
*     CLOSE(UNIT=99)
C
C  ************  Coefficients in the partial wave expansion.
C
      CXP=CDEXP(2*CI*DPC(1))
      CFACT=1.0D0/(2.0D0*CI*RK)
      CFM(1)=(CXP-1.0D0)*CFACT
      CGM(1)=0.0D0
      DO I=2,NPC
        L=I-1
        RL=L
        CXP=CDEXP(2.0D0*CI*DPC(I))
        CXM=CDEXP(2.0D0*CI*DMC(I))
        CFM(I)=((L+1)*(CXP-1)+L*(CXM-1))*CFACT
        CGM(I)=(CXM-CXP)*CFACT
      ENDDO
C
C  ****  Reduced series.
C
      NPC1=NPC
      DO NTR=1,2
        NPC1=NPC1-1
        CFC=0.0D0
        CFP=CFM(1)
        CGC=0.0D0
        CGP=CGM(1)
        DO I=1,NPC1
          RL=I-1
          CFA=CFC
          CFC=CFP
          CFP=CFM(I+1)
          CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
          CGA=CGC
          CGC=CGP
          CGP=CGM(I+1)
          CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
        ENDDO
      ENDDO
C
C  ****  Bartlett and Watson's formula for small angles.
C
      TARG=-2.0D0*CLGAM(DCMPLX(0.5D0,ETA))*CI
      C5=CDEXP(TARG*CI)
      TARG=-2.0D0*CLGAM(DCMPLX(1.0D0,ETA))*CI
      C1=CDEXP(TARG*CI)
      BETA2=E*(E+2.0D0*SL2)/(E+SL2)**2
      WATSC=-PI*BETA2*ETA*(C5/C1)
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION DPWAC
C  *********************************************************************
      FUNCTION DPWAC(TH)
C
C     Ratio (Mott DCS / Rutherford DCS) for collisions with scattering
C  angle TH (rad). Additional information is provided through the common
C  block /CRMORU/.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),
     1              CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
C
      X=COS(TH)
      Q2=2.0D0*RK2*(1.0D0-X)
      NTEST=NPC1-5
C
C  ****  TH greater than 0.5 deg.
C
      IF(TH.LT.0.008726D0) GO TO 2
      FACT=1.0D0/(1.0D0-X)**2
C  ****  Direct scattering amplitude.
      P2=1.0D0
      P3=X
      CF=CFM(1)
      CF=CF+CFM(2)*P3
      CFA=0.0D0
      DO I=3,NPC1
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CF=CF+CFM(I)*P3
        IF(I.EQ.NTEST) CFA=CF
      ENDDO
      ERRF=CDABS(CFA-CF)/MAX(CDABS(CF),1.0D-15)
      CF=FACT*CF
C  ****  Spin-flip scattering amplitude.
      Y=SIN(TH)
      IF(Y.LT.1.0D-20) THEN
        CG=0.0D0
        ERRG=ERRF
        GO TO 1
      ENDIF
      P2=1.0D0
      P3=3*X
      CG=CGM(2)
      CG=CG+CGM(3)*P3
      CGA=0.0D0
      DO I=4,NPC1
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-L*P1)/(L-1)
        CG=CG+CGM(I)*P3
        IF(I.EQ.NTEST) CGA=CG
      ENDDO
      ERRG=CDABS(CGA-CG)/MAX(CDABS(CG),1.0D-15)
    1 CG=FACT*Y*CG
      PAV1=CDABS(CF)**2
      PAV2=CDABS(CG)**2
      ERR=2.0D0*(PAV1*ERRF+PAV2*ERRG)/(PAV1+PAV2)
      DCS=(PAV1+PAV2)*A0B2
      DPWAC=DCS*(Q2*Q2/RUTHC)
      IF(ERR.LT.1.0D-3.OR.TH.GT.0.08726D0) RETURN
C
C  ****  Bartlett and Watson's formula; used only for TH less than
C        5 deg, if needed. The computed DPWAC value may have slight
C        discontinuities, of the order of 0.01 per cent, between
C        0.5 and 5 deg.
C
 2    CONTINUE
      CF=0.0D0
      CG=0.0D0
      ERRF=1.0D10
      ERRG=1.0D10
      DPWAC=1.0D0+WATSC*SIN(0.5D0*TH)
      RETURN
      END
C  *********************************************************************
C                       FUNCTION SMOMLL
C  *********************************************************************
      FUNCTION SMOMLL(X,Y,XL,XU,NP,MOM,ILOG)
C
C     Calculates integrals of a tabulated function, Y(X), over the
C  interval (XL,XU) by using linear log-log interpolation of the input
C  table. The values of both the variable X and the function Y are
C  assumed to be non-negative.
C
C  Input arguments:
C     X(1:NP) ..... array of variable values (in increasing order).
C     Y(1:NP) ..... corresponding function values.
C     NP .......... number of points in the table.
C     XL, XU ...... limits of the integration interval.
C     MOM ......... moment order.
C     ILOG ........ optional logarithm:
C                   ILOG=1,  SMOMLL = INTEGRAL X**MOM*LOG(X)*Y(X) dX
C                   else     SMOMLL = INTEGRAL X**MOM*Y(X) dX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-12, ONEM=1.0D0-EPS, ZERO=1.0D-98)
      DIMENSION X(NP),Y(NP)
C
      IF(NP.LT.2) STOP 'SMOMLL: NP is too small.'
      IF(X(1).LT.0.0D0.OR.Y(1).LT.0.0D0) THEN
        I=1
        WRITE(6,'(A,I5,1P,2E15.7)') 'I,X(I),Y(I) =',I,X(I),Y(I)
        STOP 'SMOMLL: Negative values in the table.'
      ENDIF
      DO I=2,NP
        IF(X(I).LT.0.0D0.OR.Y(I).LT.0.0D0) THEN
          WRITE(6,'(A,I5,1P,2E15.7)') 'I,X(I),Y(I) =',I,X(I),Y(I)
          STOP 'SMOMLL: Negative values in the table.'
        ENDIF
        IF(X(I).LT.X(I-1)*ONEM) THEN
          J=I-1
          WRITE(6,'(A,I5,1P,2E15.7)') 'I,X(I),Y(I) =',J,X(J),Y(J)
          WRITE(6,'(A,I5,1P,2E15.7)') 'I,X(I),Y(I) =',I,X(I),Y(I)
          STOP 'SMOMLL: X values are in decreasing order.'
        ENDIF
      ENDDO
C
      XLOW=XL
      IF(XLOW.LT.ZERO) XLOW=ZERO
      XUP=XU
C
      IF(XLOW.GT.XUP) THEN
        WRITE(6,*) 'SMOMLL (warning): XLOW is greater than XUP.'
        WRITE(6,'(A,1P,E15.7,A,E15.7)') ' XLOW =',XLOW,', XUP =',XUP
        SMOMLL=0.0D0
        RETURN
      ENDIF
C
      IF(XLOW.GT.X(NP)) THEN
        I=NP-1
      ELSE IF(XLOW.LT.X(1)) THEN
        I=1
      ELSE
        I=1
        I1=NP
 1      IT=(I+I1)/2
        IF(XLOW.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
      ENDIF
      IL=I
C
      IF(XUP.GT.X(NP)) THEN
        I=NP-1
      ELSE IF(XUP.LT.X(1)) THEN
        I=1
      ELSE
        I=1
        I1=NP
 2      IT=(I+I1)/2
        IF(XUP.GT.X(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 2
      ENDIF
      IU=I
C
      SMOMLL=0.0D0
      IF(ILOG.EQ.1) GO TO 3
C
C  ****  SMOMLL = INTEGRAL (X**N)*Y(X) dX, MOM.GT.-100.
C
      DO I=IL,IU
        XA=MAX(XLOW,X(I))
        XB=MIN(XUP,X(I+1))
        X1L=LOG(MAX(X(I),ZERO))
        X2L=LOG(MAX(X(I+1),ZERO))
        Y1L=LOG(MAX(Y(I),ZERO))
        Y2L=LOG(MAX(Y(I+1),ZERO))
        DEN=X2L-X1L
        IF(ABS(DEN).GT.EPS) THEN  ! Interpolated values.
          YA=EXP(Y1L+(Y2L-Y1L)*(LOG(XA)-X1L)/DEN)*XA**MOM
          YB=EXP(Y1L+(Y2L-Y1L)*(LOG(XB)-X1L)/DEN)*XB**MOM
        ELSE
          YAV=EXP(0.5D0*(Y1L+Y2L))
          YA=YAV*XA**MOM
          YB=YAV*XB**MOM
        ENDIF
C
        DXL=LOG(XB)-LOG(XA)
        DYL=LOG(YB)-LOG(YA)
        IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
          AP1=1.0D0+(DYL/DXL)
          IF(ABS(AP1).GT.EPS) THEN
            DSUM=(YB*XB-YA*XA)/AP1
          ELSE
            DSUM=YA*XA*DXL
          ENDIF
        ELSE
          DSUM=0.5D0*(YA+YB)*(XB-XA)
        ENDIF
        SMOMLL=SMOMLL+DSUM
      ENDDO
      RETURN
C
C  ****  SMOMLL = INTEGRAL LOG(X)*Y(X) dX, MOM.LT.-100.
C
 3    CONTINUE
      DO I=IL,IU
        XA=MAX(XLOW,X(I))
        XB=MIN(XUP,X(I+1))
        X1L=LOG(MAX(X(I),ZERO))
        X2L=LOG(X(I+1))
        Y1L=LOG(MAX(Y(I),ZERO))
        Y2L=LOG(MAX(Y(I+1),ZERO))
        DEN=X2L-X1L
        IF(ABS(DEN).GT.ZERO) THEN
          YA=EXP(Y1L+(Y2L-Y1L)*(LOG(XA)-X1L)/DEN)*XA**MOM
          YB=EXP(Y1L+(Y2L-Y1L)*(LOG(XB)-X1L)/DEN)*XB**MOM
        ELSE
          YAV=EXP(0.5D0*(Y1L+Y2L))
          YA=YAV*XA**MOM
          YB=YAV*XB**MOM
        ENDIF
        DXL=LOG(XB)-LOG(XA)
        DYL=LOG(YB)-LOG(YA)
        IF(ABS(DXL).GT.EPS*ABS(DYL)) THEN
          AP1=1.0D0+(DYL/DXL)
          IF(ABS(AP1).GT.EPS) THEN
            APREC=1.0D0/AP1
            DSUM=(YB*XB*(LOG(XB)-APREC)-YA*XA*(LOG(XA)-APREC))*APREC
          ELSE
            DSUM=YA*XA*0.5D0*(LOG(XB)**2-LOG(XA)**2)
          ENDIF
        ELSE
          DSUM=0.5D0*(YA*LOG(XA)+YB*LOG(XB))*(XB-XA)
        ENDIF
        SMOMLL=SMOMLL+DSUM
      ENDDO
      RETURN
      END

C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C  The following subroutine performs Dirac partial-wave calculations of
C  scattering of electrons and positrons in a complex central field with
C  an imaginary (absorptive) part.
C
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C  *********************************************************************
C                      SUBROUTINE DPWAI0
C  *********************************************************************
      SUBROUTINE DPWAI0(EV,TOTCS,ABCS,NDELTA,ISCH)
C
C     This subroutine computes Dirac phase shifts, differential cross
C  sections and scattering amplitudes for elastic scattering of
C  electrons in central fields with an imaginary (absorptive) part.
C
C  Input/output arguments:
C     EV ....... effective kinetic energy of the projectile (eV).
C     TOTCS .... total cross section (cm**2).
C     ABCS ..... absorption cross section (cm**2).
C     NDELTA ... number of required phase shifts (LT.25000).
C     ISCH ..... =1: all phase shifts are computed by solving the radial
C                    equation.
C                =2: only phase shifts of selected orders are computed
C                    from the solution of the radial equation, the
C                    others are obtained by lin-log natural cubic spline
C                    interpolation. For high energies, ISCH=2 leads to a
C                    considerable reduction of the calculation time.
C
C  Input (through the common block /FIELDI/):
C     R(I) .... radial grid points (radii in increasing order). The
C               first point in the grid must be the origin, i.e. R(1)=0.
C               Repeated values are interpreted as discontinuities.
C     RV(I).... R(I) times the potential energy at R=R(I). The last
C               component, RV(NP), is assumed to be equal to the
C               asymptotic value.
C     RW(I).... R(I) times the imaginary potential (it must be negative
C               or zero).
C     IAB ..... 0 if the potential is real, 1 if it has an imaginary
C               part.
C     NP ...... number of input grid points.
C
C *** NOTE: The radii and potential values, R(I) and RV(I), are in
C           atomic units.
C
C  Output (through the common block /DCSTAB/):
C     ECS ........ total cross section (cm**2)
C                    (only for finite range fields).
C     TCS1 ....... 1st transport cross section (cm**2)
C                    (only for finite range fields).
C     TCS2 ....... 2nd transport cross section (cm**2)
C                    (only for finite range fields).
C     TH(I) ...... scattering angles (in deg)
C     XT(I) ...... values of (1-COS(TH(I)))/2.0D0.
C     DCST(I) .... differential cross section per unit solid angle at
C                    TH(I) (cm**2/sr).
C     ERROR(I) ... estimated relative uncertainty of the computed DCS
C                    value.
C     NTAB ....... number of angles in the table.
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
C  ****  Input-output.
      COMMON/FIELDI/R(NDIM),RV(NDIM),RW(NDIM),IAB,NP
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
C  ****  Link with the RADIAL package.
      COMMON/RADWF/RRR(NDIM),P(NDIM),Q(NDIM),NRT,ILAST,IER
      COMMON/RADWFI/PIM(NDIM),QIM(NDIM)
      PARAMETER (NPPG=NDIM+1)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),
     1             VD(NPPG),NVT
C  ****  Phase shifts and partial wave series coefficients.
      PARAMETER (NPC=1500,NDM=25000)
      DIMENSION CXP(NDM),CXM(NDM)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/PHASEI/DPJ(NDM),DMJ(NDM)
      DIMENSION XL(NDM),DPI(NDM),DMI(NDM),DPJI(NDM),DMJI(NDM)
      DIMENSION X(NDM),Y(NDM),SA(NDM),SB(NDM),SC(NDM),SD(NDM)
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),
     1              CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
C
      CI=DCMPLX(0.0D0,1.0D0)
      ISUMP=0
C
      OPEN(98,FILE='dpwai.dat')
      WRITE(98,2000)
 2000 FORMAT(//2X,'**** Partial wave analysis (DPWAI0) ',
     1  42('*')/)
C
      NDELT=NDELTA
      IF(NDELT.GT.NDM) THEN
        WRITE(98,2001)
 2001   FORMAT(/2X,'WARNING: NDELTA is too large')
        NDELT=NDM
      ENDIF
      IF(NDELT.LT.6) NDELT=6
      EPS=1.0D-15
      EPSCUT=1.0D-9
C
      E=EV/HREV
C
C  ****  Initialization of the RADIAL package.
C
      IF(IAB.EQ.0) THEN
        CALL VINT(R,RV,NP)
      ELSE
        CALL ZVINT(R,RV,RW,NP)
      ENDIF
      ZINF=RVG(NVT)
      IF(ABS(ZINF).GT.1.0D-10) THEN
        IZINF=1
        CALL DPWAC0(ZINF,EV)
        IF(NDELT.GT.NPC) NDELT=NPC
      ELSE
        IZINF=0
      ENDIF
C
      WRITE(98,2002) EV
 2002 FORMAT(/2X,'Kinetic energy =',1P,E12.5,' eV')
      IF(E.LT.0.0D0) THEN
        WRITE(6,2002) EV
        WRITE(98,2003)
        WRITE(6,2003)
 2003   FORMAT(//2X,'Negative energy. Stop.')
        STOP 'DPWAI0: Negative energy.'
      ENDIF
      RK=SQRT(E*(E+2.0D0*SL*SL))/SL
C
      IF(IZINF.EQ.1) WRITE(98,2004)
 2004 FORMAT(/2X,'Only inner phase shifts are tabulated')
      WRITE(98,2005)
 2005 FORMAT(/14X,'--------- Spin UP ---------',6X,
     1  '-------- Spin DOWN --------',/6X,'L',9X,
     2  'Re(phase)      Im(phase)',9X,'Re(phase)      Im(phase)',
     3  /2X,74('-'))
      ISCH0=ISCH
      IF(ISCH0.EQ.2.AND.EV.GT.1000.0D0) GO TO 1
C
C  ****  ISCH0=1, all phase shifts are computed by solving the radial
C        equation.
C
      L=0
      IF(IAB.EQ.0) THEN
        CALL DFREE(E,EPS,PHP,-1,0)
        DP(1)=PHP
        DPJ(1)=0.0D0
        DM(1)=0.0D0
        DMJ(1)=0.0D0
        CXP(1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDFREE(E,EPS,PHPR,PHPI,-1,0)
        IF(ABS(PHPI).LT.1.0D-12) PHPI=0.0D0
        DP(1)=PHPR
        DPJ(1)=PHPI
        DM(1)=0.0D0
        DMJ(1)=0.0D0
        CXP(1)=CDEXP(2.0D0*(CI*PHPR-PHPI))
      ENDIF
      IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (1).'
      WRITE(98,2006) L,DP(1),DPJ(1),DM(1),DMJ(1)
      WRITE(6,2006) L,DP(1),DPJ(1),DM(1),DMJ(1)
 2006 FORMAT(3X,I5,4X,1P,E16.8,1X,E12.5,4X,E16.8,1X,E12.5)
      NPH=1
C
      IFIRST=2
 33   CONTINUE
      ISUMP=1
      TST=0.0D0
      DO I=IFIRST,NDELT
        L=I-1
        IF(IAB.EQ.0) THEN
          CALL DFREE(E,EPS,PHP,-L-1,0)
          DP(I)=PHP
          DPJ(I)=0.0D0
          CXP(I)=CDEXP(2.0D0*CI*PHP)
        ELSE
          CALL ZDFREE(E,EPS,PHPR,PHPI,-L-1,0)
          IF(ABS(PHPI).LT.1.0D-12) PHPI=0.0D0
          DP(I)=PHPR
          DPJ(I)=PHPI
          CXP(I)=CDEXP(2.0D0*(CI*PHPR-PHPI))
        ENDIF
        IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (2).'
        IF(IAB.EQ.0) THEN
          CALL DFREE(E,EPS,PHM,L,0)
          DM(I)=PHM
          DMJ(I)=0.0D0
          CXM(I)=CDEXP(2.0D0*CI*PHM)
        ELSE
          CALL ZDFREE(E,EPS,PHMR,PHMI,L,0)
          IF(ABS(PHMI).LT.1.0D-12) PHMI=0.0D0
          DM(I)=PHMR
          DMJ(I)=PHMI
          CXM(I)=CDEXP(2.0D0*(CI*PHMR-PHMI))
        ENDIF
        IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (3).'
        NPH=I
        WRITE(98,2006) L,DP(I),DPJ(I),DM(I),DMJ(I)
        WRITE(6,2006) L,DP(I),DPJ(I),DM(I),DMJ(I)
        TST=MAX(SQRT(DP(I)**2+DPJ(I)**2),SQRT(DM(I)**2+DMJ(I)**2),
     1    SQRT(DP(I-1)**2+DM(I-1)**2))
        IF(TST.LT.EPSCUT.AND.L.GT.10) GO TO 6
C  ****  When the last phase shift (spin up) differs in more than 20 per
C  cent from the quadratic extrapolation, accumulated roundoff errors
C  may be important and the calculation of phase shifts is discontinued.
        IF(I.GT.500) THEN
          DPEXT=DP(I-3)+3.0D0*(DP(I-1)-DP(I-2))
          DPMAX=MAX(ABS(DP(I-3)),ABS(DP(I-2)),ABS(DP(I-1)),
     1              ABS(DP(I)))
          IF(ABS(DP(I)-DPEXT).GT.0.20D0*DPMAX) THEN
            NPH=I-1
            WRITE(98,2107)
            WRITE(6,2107)
 2107 FORMAT(/2X,'WARNING: Possible accumulation of round-off errors.')
            GO TO 6
          ENDIF
        ENDIF
      ENDDO
      WRITE(98,2007) TST
      WRITE(6,2007) TST
 2007 FORMAT(/2X,'WARNING: TST =',1P,E11.4,'. Check convergence.')
      GO TO 6
C
C  ****  ISCH0=2, only inner phase shifts of orders L in a given grid
C        are computed from the solution of the radial equation. Phase
C        shifts of orders not included in this grid are obtained by
C        lin-log cubic spline interpolation.
C          The adopted grid is: 0(1)100(5)300(10) ...
C
C        This is a somewhat risky procedure, which is based on the
C        observed variation of the calculated phase shifts with L for
C        atomic scattering fields. When a change of sign is found, all
C        the phases are recalculated.
C
 1    CONTINUE
      L=0
      IF(IAB.EQ.0) THEN
        CALL DFREE(E,EPS,PHP,-1,0)
        DP(1)=PHP
        DPJ(1)=0.0D0
        DM(1)=0.0D0
        DMJ(1)=0.0D0
        CXP(1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDFREE(E,EPS,PHPR,PHPI,-1,0)
        IF(ABS(PHPI).LT.1.0D-12) PHPI=0.0D0
        DP(1)=PHPR
        DPJ(1)=PHPI
        DM(1)=0.0D0
        DMJ(1)=0.0D0
        CXP(1)=CDEXP(2.0D0*(CI*PHPR-PHPI))
      ENDIF
      IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (4).'
      WRITE(98,2006) L,DP(1),DPJ(1),DM(1),DMJ(1)
      WRITE(6,2006) L,DP(1),DPJ(1),DM(1),DMJ(1)
C
      LMAX=NDELT-1
      IND=0
      IADD=1
      NADD=0
      LPP=1
 2    CONTINUE
      L=LPP
      IF(IAB.EQ.0) THEN
        CALL DFREE(E,EPS,PHP,-L-1,0)
        DP(L+1)=PHP
        DPJ(L+1)=0.0D0
        CXP(L+1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDFREE(E,EPS,PHPR,PHPI,-L-1,0)
        IF(ABS(PHPI).LT.1.0D-12) PHPI=0.0D0
        DP(L+1)=PHPR
        DPJ(L+1)=PHPI
        CXP(L+1)=CDEXP(2.0D0*(CI*PHPR-PHPI))
      ENDIF
      IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (5).'
      IF(IAB.EQ.0) THEN
        CALL DFREE(E,EPS,PHM,L,0)
        DM(L+1)=PHM
        DMJ(L+1)=0.0D0
        CXM(L+1)=CDEXP(2.0D0*CI*PHM)
      ELSE
        CALL ZDFREE(E,EPS,PHMR,PHMI,L,0)
        IF(ABS(PHMI).LT.1.0D-12) PHMI=0.0D0
        DM(L+1)=PHMR
        DMJ(L+1)=PHMI
        CXM(L+1)=CDEXP(2.0D0*(CI*PHMR-PHMI))
      ENDIF
      IF(IER.NE.0) STOP 'DPWAI0: Error in ZDPHAS (6).'
      WRITE(6,2006) L,DP(L+1),DPJ(L+1),DM(L+1),DMJ(L+1)
C
      IF(L.LT.95) THEN
        WRITE(98,2006) L,DP(L+1),DPJ(L+1),DM(L+1),DMJ(L+1)
      ELSE
        IND=IND+1
        NADD=NADD+1
        XL(IND)=L
        DPI(IND)=DP(L+1)
        DPJI(IND)=DPJ(L+1)
        DMI(IND)=DM(L+1)
        DMJI(IND)=DMJ(L+1)
        IF(IND.GT.1) THEN
          TST1=DPI(IND)*DPI(IND-1)
          TST2=DMI(IND)*DMI(IND-1)
          IF(TST1.LT.0.0D0.OR.TST2.LT.0.0D0) THEN
            IF(L.LT.600) THEN
              ISCH0=1
              IFIRST=MIN(L,94)
              GO TO 33
            ELSE
              IND=IND-1
              L=XL(IND)+0.5D0
              GO TO 3
            ENDIF
          ENDIF
C
          IF(L.GT.600.AND.NADD.GT.3) THEN
            I=IND
            DPEXT=DPI(I-3)+3.0D0*(DPI(I-1)-DPI(I-2))
            DPMAX=MAX(ABS(DPI(I-3)),ABS(DPI(I-2)),ABS(DPI(I-1)),
     1            ABS(DPI(I)))
            IF(ABS(DPI(I)-DPEXT).GT.0.20D0*DPMAX) THEN
              IND=I-1
              L=XL(IND)+0.5D0
              WRITE(98,2107)
              WRITE(6,2107)
              GO TO 3
            ENDIF
            DMEXT=DMI(I-3)+3.0D0*(DMI(I-1)-DMI(I-2))
            DMMAX=MAX(ABS(DMI(I-3)),ABS(DMI(I-2)),ABS(DMI(I-1)),
     1            ABS(DMI(I)))
            IF(ABS(DMI(I)-DMEXT).GT.0.20D0*DMMAX) THEN
              IND=I-1
              L=XL(IND)+0.5D0
              WRITE(98,2107)
              WRITE(6,2107)
              GO TO 3
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      TST=MAX(SQRT(DP(L+1)**2+DPJ(L+1)**2),SQRT(DM(L+1)**2+DMJ(L+1)**2))
      IF(TST.LT.EPSCUT) THEN
        IF(L.LT.100) THEN
          NPH=L+1
          GO TO 6
        ELSE
          GO TO 3
        ENDIF
      ENDIF
      IF(L.GE.LMAX) GO TO 3
C
      IADDO=IADD
      IF(L.GT.99) IADD=5
      IF(L.GT.299) IADD=10
      IF(L.GT.599) IADD=20
      IF(L.GT.1199) IADD=50
      IF(L.GT.2999) IADD=100
      IF(L.GT.9999) IADD=250
      IF(IADD.NE.IADDO) NADD=0
      LPP=L+IADD
      IF(LPP.GT.LMAX) LPP=LMAX
      GO TO 2
C
C  ****  Check consistency of sparsely tabulated phase shifts.
C        A discontinuity larger than 0.25*PI is considered as
C        a symptom of numerical inconsistencies.
C
 3    CONTINUE
      IF(IND.LT.5.OR.IADD.EQ.1) GO TO 6
      NPH=XL(IND)+1.5D0
      TST=0.0D0
      DO I=1,IND
        WRITE(98,2008) INT(XL(I)+0.5D0),DPI(I),DPJI(I),DMI(I),DMJI(I)
 2008   FORMAT(3X,I5,4X,1P,E16.8,1X,E12.5,4X,E16.8,1X,E12.5,'  i')
        IF(I.GT.1) THEN
          TST=MAX(TST,ABS(DPI(I)-DPI(I-1)),ABS(DMI(I)-DMI(I-1)))
        ENDIF
      ENDDO
      IF(TST.GT.0.25D0*PI) THEN
        WRITE(98,2009)
        WRITE(6,2009)
 2009   FORMAT(/2X,'ERROR: Directly computed phase shifts show',
     1    ' large discontinuities.')
        STOP 'DPWAI0: Phase shifts do not vary continuously with L.'
      ENDIF
C
C  ****  Interpolated phase shifts (lin-log cubic spline).
C
      IF(DPI(4).GT.0.0D0) THEN
        ITRAN=+1
      ELSE
        ITRAN=-1
      ENDIF
      JT=0
      DO I=1,IND
        IF(DPI(I)*ITRAN.LT.0.0D0.OR.ABS(DPI(I)).LT.1.0D-12) THEN
          GO TO 4
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=LOG(ABS(DPI(I)))
        ENDIF
      ENDDO
 4    CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DP(I)=ITRAN*EXP(SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J))))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DP(I)=0.0D0
        ENDDO
      ENDIF
C
      JT=0
      DO I=1,IND
        IF(ABS(DPJI(I)).LT.1.0D-12) THEN
          GO TO 41
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=DPJI(I)
        ENDIF
      ENDDO
 41   CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DPJ(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DPJ(I)=0.0D0
        ENDDO
      ENDIF
C
      IF(DMI(4).GT.0.0D0) THEN
        ITRAN=+1
      ELSE
        ITRAN=-1
      ENDIF
      JT=0
      DO I=1,IND
        IF(DMI(I)*ITRAN.LT.0.0D0.OR.ABS(DMI(I)).LT.1.0D-12) THEN
          GO TO 5
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=LOG(ABS(DMI(I)))
        ENDIF
      ENDDO
 5    CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DM(I)=ITRAN*EXP(SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J))))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DM(I)=0.0D0
        ENDDO
      ENDIF
C
      JT=0
      DO I=1,IND
        IF(ABS(DMJI(I)).LT.1.0D-12) THEN
          GO TO 51
        ELSE
          JT=JT+1
          X(JT)=XL(I)
          Y(JT)=DMJI(I)
        ENDIF
      ENDDO
 51   CONTINUE
      NUP=X(JT)+1.5D0
      NUP=MIN(NUP,NPH)
      CALL SPLINE(X,Y,SA,SB,SC,SD,0.0D0,0.0D0,JT)
      DO I=95+1,NUP
        RL=I-1
        CALL FINDI(RL,X,JT,J)
        DMJ(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
      ENDDO
      IF(NUP.LT.NPH) THEN
        DO I=NUP+1,NPH
          DMJ(I)=0.0D0
        ENDDO
      ENDIF
C
      TST=MAX(ABS(DP(NPH)),ABS(DM(NPH)))
      IF(TST.GT.10.0D0*EPSCUT) THEN
        WRITE(98,2007) TST
        WRITE(6,2007) TST
      ENDIF
      DO I=1,NPH
        CXP(I)=CDEXP(2.0D0*CI*DCMPLX(DP(I),DPJ(I)))
        CXM(I)=CDEXP(2.0D0*CI*DCMPLX(DM(I),DMJ(I)))
      ENDDO
C
C  ************  Coefficients in the partial-wave expansion.
C
 6    CONTINUE
      CFACT=1.0D0/(2.0D0*CI*RK)
      IF(IZINF.EQ.1) THEN
        CXPC=CDEXP(2*CI*DPC(1))
        CFL(1)=CXPC*(CXP(1)-1)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXPC=CDEXP(2.0D0*CI*DPC(I))
          CXMC=CDEXP(2.0D0*CI*DMC(I))
          CFL(I)=((L+1)*CXPC*(CXP(I)-1)+L*CXMC*(CXM(I)-1))*CFACT
          CGL(I)=(CXMC*(CXM(I)-1)-CXPC*(CXP(I)-1))*CFACT
        ENDDO
      ELSE
        CFL(1)=(CXP(1)-1.0D0)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CFL(I)=((L+1)*(CXP(I)-1)+L*(CXM(I)-1))*CFACT
          CGL(I)=(CXM(I)-CXP(I))*CFACT
        ENDDO
      ENDIF
C
C  ****  Reduced series (two iterations).
C
      IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
        DO I=1,NPH
          CFM(I)=CFL(I)
          CGM(I)=CGL(I)
        ENDDO
C
        NPHM=NPH
        DO 7 NTR=1,2
          NPHM=NPHM-1
          CFC=0.0D0
          CFP=CFM(1)
          CGC=0.0D0
          CGP=CGM(1)
          DO I=1,NPHM
            RL=I-1
            CFA=CFC
            CFC=CFP
            CFP=CFM(I+1)
            CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
            CGA=CGC
            CGC=CGP
            CGP=CGM(I+1)
            CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
          ENDDO
 7      CONTINUE
      ENDIF
C
*     OPEN(99, file='pwa-coefs0.dat')
*     WRITE(99,'(A)') '# Coefficients in the partial-wave expansions'
*     WRITE(99,'(A,I6)') '# NPH =',NPH
*     WRITE(99,'(A)') '# L, CFL(L),CGL(L) all in a.u.'
*     DO I=1,NPH
*       WRITE(99,'(I5,1P,2(2X,E14.6,E14.6))') I-1,CFL(I),CGL(I)
*     ENDDO
*     CLOSE(99)
*     IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
*       OPEN(99, file='pwa-coefs2.dat')
*       WRITE(99,'(A)') '# Coefficients in the partial-wave expansions'
*       WRITE(99,'(A)') '#   after the reduced-series transformation'
*       WRITE(99,'(A,I6)') '# NPHM =',NPHM
*       WRITE(99,'(A)') '# L, CFM(L),CGM(L) all in a.u.'
*       DO I=1,NPHM
*         WRITE(99,'(I5,1P,2(2X,E14.6,E14.6))') I-1,CFM(I),CGM(I)
*       ENDDO
*       CLOSE(99)
*     ENDIF
C
C  ****  Scattering amplitudes and DCS.
C
      WRITE(98,2010)
 2010 FORMAT(//2X,'*** Scattering amplitudes and different',
     1  'ial cross section ***')
      WRITE(98,2011)
 2011 FORMAT(/4X,'Angle',6X,'DCS',7X,'Asymmetry',4X,'Direct amplitu',
     1  'de',7X,'Spin-flip amplitude',5X,'Error',/4X,'(deg)',3X,
     2  '(cm**2/sr)',22X,'(cm)',20X,'(cm)',/2X,91('-'))
C
C  ****  Angular grid (TH in deg).
C
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
 10   CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(I.GT.NGT) STOP 'DPWAI0. The NGT parameter is too small.'
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
C
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
        CALL DPWA(THR,CF,CG,DCS,SPL,ERRF,ERRG)
        IF(MAX(ERRF,ERRG).GT.0.95D0) THEN
          ERR=1.0D0
        ELSE
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/MAX(DCS,1.0D-45)
        ENDIF
        DCST(I)=DCS
        ERROR(I)=MAX(ERR,1.0D-7)
        SPOL(I)=SPL
        WRITE(98,2012) TH(I),DCST(I),SPOL(I),CF,CG,ERROR(I)
 2012   FORMAT(1X,1P,E10.3,E12.5,1X,E10.3,2(1X,'(',E10.3,',',
     1    E10.3,')'),E10.2)
      ENDDO
C
C  ************  Total and momentum transfer cross sections.
C                Convergence test (only for finite range fields).
C
      IF(IZINF.EQ.0) THEN
        INC=5
        IF(ISUMP.EQ.1) INC=1
        TST1=0.0D0
        TST2=0.0D0
        ECS=4.0D0*PI*CFL(1)*DCONJG(CFL(1))
        TCS=0.0D0
        ECSO=ECS
        TCSO=TCS
        DO I=2,NPH
          L=I-1
          RL=L
          DECS=CFL(I)*DCONJG(CFL(I))+RL*(L+1)*CGL(I)*DCONJG(CGL(I))
          DECS=4.0D0*PI*DECS/(L+L+1)
          DTCS=CFL(L)*DCONJG(CFL(I))+DCONJG(CFL(L))*CFL(I)
     1        +(L-1)*(RL+1)*(CGL(L)*DCONJG(CGL(I))
     2        +DCONJG(CGL(L))*CGL(I))
          DTCS=4.0D0*PI*DTCS*L/((RL+L-1)*(L+L+1))
          ECS=ECS+DECS
          TCS=TCS+DTCS
C  ****  Convergence test.
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST1=ABS(ECS-ECSO)/(ABS(ECS)+1.0D-35)
            TST2=ABS(TCS-TCSO)/(ABS(TCS)+1.0D-35)
            ECSO=ECS
            TCSO=TCS
          ENDIF
        ENDDO
        TST=MAX(TST1,TST2)
        TCS=ECS-TCS
        IF(TST.GT.1.0D-5.AND.NPH.GT.40) THEN
          WRITE(98,2007) TST
          WRITE(6,2007) TST
        ENDIF
        ECS=ECS*A0B2
        TCS=TCS*A0B2
C
C  ****  ECS and TCSs are evaluated from the DCS table.
C
        ECS0=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,0,0)
        ECS1=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,1,0)
        ECS2=FOURPI*SMOMLL(XT,DCST,XT(1),XT(NTAB),NTAB,2,0)
        TST1=ABS(ECS-ECS0)/(ABS(ECS)+1.0D-35)
        WRITE(98,2013) ECS,ECS0,TST1
        WRITE(6,2013) ECS,ECS0,TST1
 2013   FORMAT(/2X,'Total elastic cross section =',1P,E13.6,' cm**2',
     1         /2X,'             From DCS table =',E13.6,
     2         '  (Rel. dif. =',E9.2,')')
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TST2=ABS(TCS-TCS1)/(ABS(TCS)+1.0D-35)
        WRITE(98,2014) TCS,TCS1,TST2
        WRITE(6,2014) TCS,TCS1,TST2
 2014   FORMAT(/2X,'1ST transport cross section =',1P,E13.6,' cm**2',
     1         /2X,'             From DCS table =',E13.6,
     2         '  (REL. DIF. =',E9.2,')')
        WRITE(98,2015) TCS2
        WRITE(6,2015) TCS2
 2015   FORMAT(/2X,'2ND transport cross section =',1P,E13.6,' cm**2')
        TST=MAX(TST1,TST2)
        IF(TST.GT.2.0D-3) THEN
          WRITE(98,2016)
          WRITE(6,2016)
        ENDIF
 2016   FORMAT(/2X,'WARNING: Relative differences are too large.',
     1         /11X,'The DCS table is not consistent.')
C
C  ****  Absorption cross section.
C
        CALL DPWA(0.0D0,CF,CG,DCS,SPL,ERRF,ERRG)
        TOTCS=(FOURPI*(A0B/RK))*(-CI*CF)
        ABCS=TOTCS-ECS
        WRITE(98,2018) TOTCS
        WRITE(6,2018) TOTCS
 2018   FORMAT(/2X,'  Grand total cross section =',1P,E13.6,' cm**2')
        WRITE(98,2019) ABCS
        WRITE(6,2019) ABCS
 2019   FORMAT(2X,'   Absorption cross section =',1P,E13.6,' cm**2')
      ENDIF
C
      WRITE(98,2017)
 2017 FORMAT(/2X,'**** DPWAI0 ended ',60('*')/)
      CLOSE(UNIT=98)
C
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C    +++++++++++++++++++++++++++++
C    +++    PROGRAM ELSCATA    +++
C    +++++++++++++++++++++++++++++
C
C
C                               F. Salvat, A. Jablonski and C.J. Powell
C                               November 30, 2003
C
C  Updated in February 2020 by F. Salvat, to work with the subroutine
C  package RADIAL.
C
C  Ref.: F. Salvat and J. M. Fernandez-Varea,
C        'RADIAL: a Fortran subroutine package for the solution of
C        radial Schrodinger and Dirac wave equations',
C        Comput. Phys. Commun. 240 (2019) 165-177.
C
C
C  ELastic SCATtering of electrons and positrons by Atoms (and ions)
C  --      ----                                     -
C
C     This program does relativistic (Dirac) partial wave calculations
C  of elastic scattering of electrons and positrons by neutral atoms and
C  positive ions with Z=1-103. The interaction between the projectile
C  and the target system is assumed to be equal to the electrostatic
C  potential energy, plus a local exchange potential in the case of
C  projectile electrons. For projectiles with relatively low energies,
C  approximate correlation-polarization and absorption corrections can
C  be introduced.
C
C     The code calculates the phase shifts, the direct and spin-flip
C  scattering amplitudes, the differential cross section (DCS) and the
C  Sherman function (spin polarization) for electrons and positrons with
C  kinetic energies larger than 10 eV. For kinetic energies larger than
C  10 MeV, the convergence of the partial-wave series is too slow and
C  the differential cross section is calculated using approximate high-
C  energy factorization methods.
C
C
C  **** The input data file.
C
C     Data are read from a formatted input file (unit 5). Each line in
C  this file consists of a 6-character keyword (columns 1-6) followed by
C  a numerical value (in free format) that is written on the columns
C  8-19. Keywords are explicitly used/verified by the program (which is
C  case sensitive!). The text after column 20 describes the input
C  quantity and its default value (in square brackets). This text is a
C  reminder for the user and is not read by the program.
C
C     Lines defining default values can be omitted from the input file.
C  The program assigns default values to input parameters that are not
C  explicitly defined, and also to those that are manifestly wrong. The
C  code stops when it finds an inconsistent input datum. The conflicting
C  quantity appears in the last line written on the screen.
C
C ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
C IZ      80         atomic number                               [none]
C MNUCL   3          rho_n (1=P, 2=U, 3=F, 4=Uu)                  [  3]
C NELEC   80         number of bound electrons                    [ IZ]
C MELEC   4          rho_e (1=TFM, 2=TFD, 3=DHFS, 4=DF, 5=file)   [  4]
C MUFFIN  0          0=free atom, 1=muffin-tin model              [  0]
C RMUF   -1.0        muffin-tin radius (cm)                  [measured]
C IELEC  -1          -1=electron, +1=positron                     [ -1]
C MEXCH   1          V_ex (0=none, 1=FM, 2=TF, 3=RT)              [  1]
C MCPOL   2          V_cp (0=none, 1=B, 2=LDA)                    [  0]
C VPOLA  -1.0        atomic polarizability (cm^3)            [measured]
C VPOLB  -1.0        b_pol parameter                          [default]
C MABS    2          W_abs (0=none, 1=LDA-I, 2=LDA-II)            [  0]
C VABSA  -1.0        absorption-potential strength, Aabs      [default]
C VABSD  -1.0        energy gap DELTA (eV)                    [default]
C IHEF    2          high-E factorization (0=no, 1=yes, 2=Born)   [  1]
C EV      1.000E2    kinetic energy (eV)                         [none]
C EV      1.000E3    optionally, more energies...
C ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
C
C  NOTE: The absorption potential W_abs is obtained from the local-
C  density approximation (LDA) by using the Lindhard dielectric function
C  of the free-electron gas. The LDA-I option only accounts for
C  electron-hole excitations, while the LDA-II option uses the full
C  dielectric function, i.e., it includes plasmon-like excitations.
C
C     The computed information is delivered in various files with the
C  extension '.dat'. The output file 'dcs_xpyyyezz.dat' contains the
C  calculated DCS for the energy x.yyyEzz (E format, in eV) in a format
C  ready for visualization with a plotting program. The output files
C  'scfield.dat' and 'scatamp.dat' contain tables of the scattering
C  potential and the scattering amplitudes of the last calculated case;
C  they are overwritten (i.e. lost) when a new case is calculated.
C
C     This program uses the subroutine packages 'elsepa2020.f' and
C  'radial.f', which are inserted into this source file by means of
C  two INCLUDE statements (see above).
C
C
      USE CONSTANTS
C
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C),
     1   INTEGER*4 (I-N)
      CHARACTER*6 KWORD
      CHARACTER*12 BUFFER,OFILE
C
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0)
C  ****  Results from the partial wave calculation.
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),
     1              ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
C
C  ****  Atomic polarizabilities of free atoms (in cm**3), from
C    Thomas M. Miller, 'Atomic and molecular polarizabilities' in
C    CRC Handbook of Chemistry and Physics, Editor-in-chief David
C    R. Linde. 79th ed., 1998-1999, pp. 10-160 to 10-174.
      DIMENSION ATPOL(103)
      DATA ATPOL/   0.666D-24, 0.205D-24, 24.30D-24, 5.600D-24,
     A   3.030D-24, 1.760D-24, 1.100D-24, 0.802D-24, 0.557D-24,
     1   3.956D-25, 24.08D-24, 10.06D-24, 6.800D-24, 5.380D-24,
     A   3.630D-24, 2.900D-24, 2.180D-24, 1.641D-24, 43.40D-24,
     2   22.80D-24, 17.80D-24, 14.60D-24, 12.40D-24, 11.60D-24,
     A   9.400D-24, 8.400D-24, 7.500D-24, 6.800D-24, 6.100D-24,
     3   7.100D-24, 8.120D-24, 6.070D-24, 4.310D-24, 3.770D-24,
     A   3.050D-24, 2.484D-24, 47.30D-24, 27.60D-24, 22.70D-24,
     4   17.90D-24, 15.70D-24, 12.80D-24, 11.40D-24, 9.600D-24,
     A   8.600D-24, 4.800D-24, 7.200D-24, 7.200D-24, 10.20D-24,
     5   7.700D-24, 6.600D-24, 5.500D-24, 5.350D-24, 4.044D-24,
     A   59.60D-24, 39.70D-24, 31.10D-24, 29.60D-24, 28.20D-24,
     6   31.40D-24, 30.10D-24, 28.80D-24, 27.70D-24, 23.50D-24,
     A   25.50D-24, 24.50D-24, 23.60D-24, 22.70D-24, 21.80D-24,
     7   21.00D-24, 21.90D-24, 16.20D-24, 13.10D-24, 11.10D-24,
     A   9.700D-24, 8.500D-24, 7.600D-24, 6.500D-24, 5.800D-24,
     8   5.100D-24, 7.600D-24, 6.800D-24, 7.400D-24, 6.800D-24,
     A   6.000D-24, 5.300D-24, 48.70D-24, 38.30D-24, 32.10D-24,
     9   32.10D-24, 25.40D-24, 24.90D-24, 24.80D-24, 24.50D-24,
     A   23.30D-24, 23.00D-24, 22.70D-24, 20.50D-24, 19.70D-24,
     1   23.80D-24, 18.20D-24, 17.50D-24, 20.00D-24/
C
C  ****  Ionization energies of neutral atoms (in eV).
C        NIST Physical Reference Data.
C        http://sed.nist.gov/PhysRefData/IonEnergy/tblNew.html
C  For astatine (Z=85), the value given below was calculated with
C  the DHFXA code (Salvat and Fernandez-Varea, UBIR-2003).
      DIMENSION EIONZ(103)
      DATA EIONZ/  13.5984D0, 24.5874D0, 5.39170D0, 9.32270D0,
     A  8.29800D0, 11.2603D0, 14.5341D0, 13.6181D0, 17.4228D0,
     1  21.5646D0, 5.13910D0, 7.64620D0, 5.98580D0, 8.15170D0,
     A  10.4867D0, 10.3600D0, 12.9676D0, 15.7596D0, 4.34070D0,
     2  6.11320D0, 6.56150D0, 6.82810D0, 6.74620D0, 6.76650D0,
     A  7.43400D0, 7.90240D0, 7.88100D0, 7.63980D0, 7.72640D0,
     3  9.39420D0, 5.99930D0, 7.89940D0, 9.78860D0, 9.75240D0,
     A  11.8138D0, 13.9996D0, 4.17710D0, 5.69490D0, 6.21710D0,
     4  6.63390D0, 6.75890D0, 7.09240D0, 7.28000D0, 7.36050D0,
     A  7.45890D0, 8.33690D0, 7.57620D0, 8.99380D0, 5.78640D0,
     5  7.34390D0, 8.60840D0, 9.00960D0, 10.4513D0, 12.1298D0,
     A  3.89390D0, 5.21170D0, 5.57690D0, 5.53870D0, 5.47300D0,
     6  5.52500D0, 5.58200D0, 5.64360D0, 5.67040D0, 6.15010D0,
     A  5.86380D0, 5.93890D0, 6.02150D0, 6.10770D0, 6.18430D0,
     7  6.25420D0, 5.42590D0, 6.82510D0, 7.54960D0, 7.86400D0,
     A  7.83350D0, 8.43820D0, 8.96700D0, 8.95870D0, 9.22550D0,
     8  10.4375D0, 6.10820D0, 7.41670D0, 7.28560D0, 8.41700D0,
     A  9.50000D0, 10.7485D0, 4.07270D0, 5.27840D0, 5.17000D0,
     9  6.30670D0, 5.89000D0, 6.19410D0, 6.26570D0, 6.02620D0,
     A  5.97380D0, 5.99150D0, 6.19790D0, 6.28170D0, 6.42000D0,
     1  6.50000D0, 6.58000D0, 6.65000D0, 4.90000D0/
C
C  ****  First excitation energies of neutral atoms (in eV).
C        NIST Physical Reference Data.
C  The value 0.0D0 indicates that the experimental value for
C  the atom was not available.
      DIMENSION EEX1Z(103)
      DATA EEX1Z/ 10.20D0, 19.82D0,  1.85D0,  2.73D0,
     A    3.58D0,  1.26D0,  2.38D0,  1.97D0, 12.70D0,
     1   16.62D0,  2.10D0,  2.71D0,  3.14D0,  0.78D0,
     A    1.41D0,  1.15D0,  8.92D0, 11.55D0,  1.61D0,
     2    1.88D0,  1.43D0,  0.81D0,  0.26D0,  0.94D0,
     A    2.11D0,  0.86D0,  0.43D0,  0.01D0,  1.38D0,
     3    4.00D0,  3.07D0,  0.88D0,  1.31D0,  1.19D0,
     A    7.87D0,  9.91D0,  0.00D0,  0.00D0,  0.00D0,
     4    0.00D0,  0.00D0,  1.34D0,  0.00D0,  0.00D0,
     A    0.00D0,  0.00D0,  0.00D0,  3.73D0,  0.00D0,
     5    0.00D0,  0.00D0,  0.00D0,  0.00D0,  8.31D0,
     A    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     6    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     A    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     7    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     A    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     8    4.67D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     A    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     9    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     A    0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     1    0.00D0,  0.00D0,  0.00D0,  0.00D0/
C
C  ****  Nearest-neighbour distances (in cm) of the elements,
C    from Ch. Kittel, 'Introduction to Solid State Physics'. 5th
C    ed. (John Wiley and Sons, New York, 1976).
C  The value -1.0D-8 indicates that the experimental value for
C  the element was not available.
      DIMENSION DNNEL(103)
      DATA DNNEL/    -1.000D-8,-1.000D-8, 3.124D-8, 2.220D-8,
     A     -1.000D-8, 1.540D-8,-1.000D-8,-1.000D-8, 1.440D-8,
     1     -1.000D-8, 3.822D-8, 3.200D-8, 2.860D-8, 2.350D-8,
     A     -1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8, 4.752D-8,
     2      3.950D-8, 3.250D-8, 2.890D-8, 2.620D-8, 2.500D-8,
     A      2.240D-8, 2.480D-8, 2.500D-8, 2.490D-8, 2.560D-8,
     3      2.660D-8, 2.440D-8, 2.450D-8, 3.160D-8, 2.320D-8,
     A     -1.000D-8,-1.000D-8, 4.837D-8, 4.300D-8, 3.550D-8,
     4      3.170D-8, 2.860D-8, 2.720D-8, 2.710D-8, 2.650D-8,
     A      2.690D-8, 2.750D-8, 2.890D-8, 2.980D-8, 3.250D-8,
     5      2.810D-8, 2.910D-8, 2.860D-8, 3.540D-8,-1.000D-8,
     A      5.235D-8, 4.350D-8, 3.730D-8, 3.650D-8, 3.630D-8,
     6      3.660D-8,-1.000D-8, 3.590D-8, 3.960D-8, 3.580D-8,
     A      3.520D-8, 3.510D-8, 3.490D-8, 3.470D-8, 3.540D-8,
     7      3.880D-8, 3.430D-8, 3.130D-8, 2.860D-8, 2.740D-8,
     A      2.740D-8, 2.680D-8, 2.710D-8, 2.770D-8, 2.880D-8,
     8     -1.000D-8, 3.460D-8, 3.500D-8, 3.070D-8, 3.340D-8,
     A     -1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8, 3.760D-8,
     9      3.600D-8, 3.210D-8, 2.750D-8, 2.620D-8, 3.100D-8,
     A      3.610D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8,
     1     -1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8/
C
C  ****  Atomic densities (number of atoms per cm**3) of the
C    elements in their natural state (from the PENELOPE database).
      DIMENSION VMOLE(103)
      DATA VMOLE/    5.0039E+19, 2.5024E+19, 4.6331E+22, 1.2349E+23,
     A   1.3202E+23, 1.0028E+23, 5.0100E+19, 5.0119E+19, 5.0093E+19,
     1   2.5024E+19, 2.5435E+22, 4.3113E+22, 6.0237E+22, 4.9959E+22,
     A   4.2774E+22, 3.7561E+22, 5.0869E+19, 2.5055E+19, 1.3277E+22,
     2   2.3290E+22, 4.0040E+22, 5.7102E+22, 7.2230E+22, 8.3158E+22,
     A   8.1555E+22, 8.4908E+22, 9.0946E+22, 9.1343E+22, 8.4912E+22,
     3   6.5692E+22, 5.0994E+22, 4.4148E+22, 4.6057E+22, 3.4321E+22,
     A   5.3301E+19, 2.4996E+19, 1.0795E+22, 1.7457E+22, 3.0271E+22,
     4   4.2949E+22, 5.5551E+22, 6.4151E+22, 7.0735E+22, 7.3944E+22,
     A   7.2621E+22, 6.8019E+22, 5.8619E+22, 4.6341E+22, 3.8340E+22,
     5   3.7084E+22, 3.3096E+22, 2.9450E+22, 2.3396E+22, 2.5161E+19,
     A   8.4865E+21, 1.5348E+22, 2.6679E+22, 2.8611E+22, 2.8677E+22,
     6   2.8808E+22, 3.0005E+22, 2.9878E+22, 2.0778E+22, 3.0256E+22,
     A   3.1181E+22, 3.1686E+22, 3.2113E+22, 3.2642E+22, 3.3228E+22,
     7   2.3422E+22, 3.3867E+22, 4.4907E+22, 5.5426E+22, 6.3219E+22,
     A   6.7980E+22, 7.1461E+22, 7.0241E+22, 6.6216E+22, 5.9069E+22,
     8   4.0668E+22, 3.4533E+22, 3.2988E+22, 2.8088E+22, 2.6857E+22,
     A   2.6728E+22, 2.4591E+19, 2.7003E+21, 1.3322E+22, 2.6711E+22,
     9   3.0417E+22, 4.0062E+22, 4.7943E+22, 5.1444E+22, 4.9981E+22,
     A   3.3869E+22, 3.2930E+22, 3.4124E+22, 3.3579E+22, 3.3446E+22,
     1  -1.0000D+00,-1.0000D+00,-1.0000D+00,-1.0000D+00/
C
C  ****  Phase shifts.
      PARAMETER (NDM=25000)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/PHASEI/DPJ(NDM),DMJ(NDM)
C
C  ************  Input data.
C
C  ****  Default model.
      IELEC =-1        ! electron
      IZ    = 0        ! no default
      NELEC = 1000     ! =Z (the present value is a flag)
      MNUCL = 3        ! Fermi nuclear charge distribution
      MELEC = 4        ! DF electron density
      MUFFIN= 0        ! free atom
      RMUF  =-1.0D0    ! free atom
      VMOL  =-1.0D0  ! free atom
      MEXCH = 1        ! FM exchange potential
      MCPOL = 0        ! no correlation-polarization
      VPOLA =-1.0D0    ! atomic polarizability
      VPOLB =-1.0D0    ! polariz. cutoff parameter
      MABS  = 0        ! no absorption
      VABSA =-1.0D0    ! absorption potential strength
      VABSD =-1.0D0    ! energy gap
      IHEF  = 1        ! high-energy factorization on
C
  100 CONTINUE
      READ(5,'(A6,1X,A12)') KWORD,BUFFER
      IF(KWORD.EQ.'IZ    ') THEN
        READ(BUFFER,*) IZ
        IF(IZ.LT.1.OR.IZ.GT.103) THEN
          WRITE(6,*) 'IZ =',IZ
          STOP 'Wrong atomic number.'
        ENDIF
      ELSE IF(KWORD.EQ.'MNUCL ') THEN
        READ(BUFFER,*) MNUCL
        IF(MNUCL.LT.1.OR.MNUCL.GT.4) MNUCL=3
      ELSE IF(KWORD.EQ.'NELEC ') THEN
        READ(BUFFER,*) NELEC
        IF(NELEC.LT.0) NELEC=IZ
        IF(NELEC.GT.IZ) THEN
          WRITE(6,*) 'NELEC =',NELEC
          STOP 'Negative ion.'
        ENDIF
      ELSE IF(KWORD.EQ.'MELEC ') THEN
        READ(BUFFER,*) MELEC
        IF(MELEC.LT.1.OR.MELEC.GT.5) MELEC=4
      ELSE IF(KWORD.EQ.'MUFFIN') THEN
        READ(BUFFER,*) MUFFIN
        IF(MUFFIN.NE.1) MUFFIN=0
      ELSE IF(KWORD.EQ.'RMUF  ') THEN
        READ(BUFFER,*) RMUF
      ELSE IF(KWORD.EQ.'IELEC ') THEN
        READ(BUFFER,*) IELEC
        IF(IELEC.NE.+1) IELEC=-1
      ELSE IF(KWORD.EQ.'MEXCH ') THEN
        READ(BUFFER,*) MEXCH
        IF(MEXCH.LT.0.OR.MEXCH.GT.3) MEXCH=1
      ELSE IF(KWORD.EQ.'MCPOL ') THEN
        READ(BUFFER,*) MCPOL
        IF(MCPOL.LT.0.OR.MCPOL.GT.2) MCPOL=0
      ELSE IF(KWORD.EQ.'VPOLA ') THEN
        READ(BUFFER,*) VPOLA
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)
      ELSE IF(KWORD.EQ.'VPOLB ') THEN
        READ(BUFFER,*) VPOLB
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
      ELSE IF(KWORD.EQ.'MABS  ') THEN
        READ(BUFFER,*) MABS
        IF(MABS.LT.0.OR.MABS.GT.2) MABS=0
      ELSE IF(KWORD.EQ.'VABSA ') THEN
        READ(BUFFER,*) VABSAI
        IF(VABSAI.GT.0.0D0) VABSA=VABSAI
      ELSE IF(KWORD.EQ.'VABSD ') THEN
        READ(BUFFER,*) VABSD
        IF(VABSD.LT.-1.0D-35) VABSD=-1.0D0
      ELSE IF(KWORD.EQ.'IHEF') THEN
        READ(BUFFER,*) IHEF
        IF(IHEF.NE.0.AND.IHEF.NE.2) IHEF=1
      ELSE IF(KWORD.EQ.'EV    ') THEN
        IF(VMOL.LT.-1.0D-35) VMOL=VMOLE(IZ)
        GO TO 200
      ELSE
        WRITE(6,*) 'Unrecognized keyword.'
        GO TO 100
      ENDIF
      GO TO 100
C
C  ****  Potential model parameters.
C
  200 CONTINUE
      IF(NELEC.EQ.1000) NELEC=IZ
C
      IF(NELEC.NE.IZ) MUFFIN=0
      IF(MUFFIN.EQ.1) THEN
        IF(RMUF.LT.1.0D-9) RMUF=0.5D0*DNNEL(IZ)
        IF(RMUF.LT.1.0D-9) THEN
          WRITE(6,*) 'RMUF = ',RMUF
          STOP 'RMUF is too small.'
        ENDIF
        VMOL=VMOLE(IZ)
        IF(VMOL.LT.1.0D0) VMOL=1.0D0/(2.0D0*RMUF)**3
C  ****  If the volume of the muffin-tin sphere is larger than that of
C  the Wigner-Seitz cell (= 1/VMOL), a cubic lattice is assumed.
        RWS=(0.75D0/(PI*VMOL))**(1.0D0/3.0D0)
        IF(RMUF.GT.RWS) THEN
          VMOL=1.0D0/(2.0D0*RMUF)**3
        ENDIF
        IHEF=0
      ELSE
        MUFFIN=0
        RMUF=200.0D-8
      ENDIF
C
      IF(IELEC.EQ.+1) MEXCH=0
C
      IF(MCPOL.EQ.1) THEN
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)  ! Default value.
        IF(VPOLA.LT.-1.0D-35) THEN
          WRITE(6,*) 'VPOLA = ',VPOLA
          STOP 'VPOLA must be positive.'
        ENDIF
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
      ELSE IF(MCPOL.EQ.2) THEN
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)
        IF(VPOLA.LT.-1.0D-35) THEN
          WRITE(6,*) 'VPOLA = ',VPOLA
          STOP 'VPOLA must be positive.'
        ENDIF
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
      ELSE
        MCPOL=0
        VPOLA=0.0D0
        VPOLB=0.0D0
      ENDIF
C
      IF(MABS.EQ.1.OR.MABS.EQ.2) THEN
        IF(VABSA.LT.-1.0D-35) THEN
          IF(MABS.EQ.1) THEN
            VABSA=2.0D0
          ELSE
            VABSA=0.75D0
          ENDIF
        ENDIF
        IF(VABSD.LT.-1.0D-35) THEN
          IF(IELEC.EQ.-1) THEN  ! Default (experimental) value.
            VABSD=EEX1Z(IZ)
            IF(VABSD.LT.-1.0D-35) VABSD=0.0D0
          ELSE
            VABSD=MAX(0.0D0,EIONZ(IZ)-6.8D0)
          ENDIF
        ENDIF
      ELSE
        MABS=0
        VABSA=0.0D0
        VABSD=1.0D0
      ENDIF
C
      OPEN(25,FILE='tcstable.dat')
      WRITE(25,1001)
 1001 FORMAT(1X,'#  Total cross section table (last run of',
     1  ' ''elscata'').', /1X,'#')
      IF(IELEC.EQ.-1) THEN
        WRITE(25,1002) IZ
 1002   FORMAT(1X,'#  Z =',I4,',   projectile: electron')
      ELSE
        WRITE(25,1003) IZ
 1003   FORMAT(1X,'#  Z =',I4,',   projectile: positron')
      ENDIF
      IF((MUFFIN.EQ.1).AND.(MABS.EQ.1.OR.MABS.EQ.2)) THEN
        WRITE(25,1004)
 1004   FORMAT(1X,'#',/1X,'#   Energy',9X,'ECS',9X,'TCS1',9X,'TCS2',
     1    9X,'ABSCS',6X,'inel-MFP',5X,'error',/1X,'#',4X,'(eV)',
     2    8X,'(cm**2)',6X,'(cm**2)',6X,'(cm**2)',6X,'(cm**2)',
     3    7X,'(cm)',/1X,'#',87('-'))
      ELSE
        WRITE(25,1005)
 1005   FORMAT(1X,'#',/1X,'#   Energy',9X,'ECS',9X,'TCS1',9X,'TCS2',
     1    9X,'ABSCS',6X,'error',/1X,'#',4X,'(eV)',8X,'(cm**2)',6X,
     2    '(cm**2)',6X,'(cm**2)',6X,'(cm**2)',/1X,'#',74('-'))
      ENDIF
C
C  ************  Partial-wave analysis.
C
  300 CONTINUE
      READ(BUFFER,*) EV
      WRITE(6,*) '   '
      WRITE(6,'(A,1P,E15.7)') 'E (eV)=',EV
C  ****  You may wish to comment off the next condition to run the
C  program for kinetic energies less that 5 eV. However, the results
C  for these energies may be highly inaccurate.
      IF(EV.LT.4.999D0) STOP 'The kinetic energy is too small.'
C
      IF(MCPOL.NE.0) THEN
        IF(EV.GT.1.0D4) THEN
          WRITE(6,*) 'WARNING: For E>10 keV the correlation-polari',
     1      'zation correction is'
          WRITE(6,*) '         switched off.'
          MCPOLC=0
          VPOLBC=0.0D0
        ELSE
          MCPOLC=MCPOL
          IF(VPOLB.LT.0.0D0) THEN
            VPOLBC=SQRT(MAX((EV-50.0D0)/16.0D0,1.0D0))
          ELSE
            VPOLBC=VPOLB
          ENDIF
        ENDIF
      ELSE
        MCPOLC=0
        VPOLBC=0.0D0
      ENDIF
      IF(MABS.NE.0) THEN
        IF(EV.GT.1.0D6) THEN
          WRITE(6,*) 'WARNING: For E>1 MeV, the absorption correc',
     1      'tion is switched off.'
          MABSC=0
        ELSE
          MABSC=MABS
        ENDIF
      ELSE
        MABSC=0
      ENDIF
C
      WRITE(BUFFER,'(1P,E12.5)') EV
      OFILE=BUFFER(2:2)//'p'//BUFFER(4:6)//'e'//BUFFER(11:12)
      OPEN(8,FILE='dcs_'//OFILE(1:8)//'.dat')
      CALL ELSEPA(IELEC,EV,IZ,NELEC,MNUCL,MELEC,MUFFIN,RMUF,VMOL,
     1  MEXCH,MCPOLC,VPOLA,VPOLBC,MABSC,VABSA,VABSD,IHEF,8)
      CLOSE(8)
C
      ERRM=0.0D0
      OPEN(8,FILE='dcs.dat')
      WRITE(8,'(1X,''#'')')
      WRITE(8,'(1X,''# Differential cross section:'',15X,
     1  ''MU=(1-COS(THETA))/2'')')
      WRITE(8,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',
     1  10X,''DCS'',8X,''Sherman''/1X,''#  (deg)'',
     2  17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',4X,''function'',
     3  /1X,''#'',62(''-''))')
      DO I=1,NTAB
        ERRM=MAX(ERRM,ERROR(I))
        WRITE(8,1006) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I)
 1006   FORMAT(1X,1P,E10.3,E13.5,3E13.5,E9.1)
      ENDDO
      CLOSE(8)
C
      IF((MUFFIN.EQ.1).AND.(MABS.EQ.1.OR.MABS.EQ.2)) THEN
        AMFP=1.0D0/(VMOL*ABCS)
        WRITE(25,1007) EV,ECS,TCS1,TCS2,ABCS,AMFP,ERRM
 1007   FORMAT(1X,1P,6E13.5,E9.1)
      ELSE
        WRITE(25,1008) EV,ECS,TCS1,TCS2,ABCS,ERRM
 1008   FORMAT(1X,1P,5E13.5,E9.1)
      ENDIF
C
      IF(MUFFIN.EQ.0.AND.MCPOL.EQ.0.AND.MABSC.EQ.0) THEN
        OPEN(8,FILE='dcsMott.dat')
        CALL MOTTSC(IELEC,IZ,EV,8)
        CLOSE(8)
      ENDIF
C
 400  CONTINUE
      READ(5,'(A6,1X,A12)',END=9999) KWORD,BUFFER
      IF(KWORD.EQ.'EV    ') GO TO 300
      GO TO 400
C
 9999 CONTINUE
      CLOSE(25)
      END
