

C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKISIG (P, T, X, XNUIK, KK, ICKWRK, RCKWRK, SIGI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKISIG (P, T, X, XNUIK, KK, ICKWRK, RCKWRK, SIGI)
C  Returns the ion species electrical conductivities given
C  collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUIK(*,*)- Real matrix, momentum-transfer collision frequencies
C              for the ions with the species;
C              dimension at least KK for the first, the total species
C              count, and at least NKKI for the second, the ion count.
C  KK        - Integer scalar, first dimension of XNUIK.
C
C  OUTPUT
C  SIGI(*)   - Real array, ion electrical conductivities (DC);
C              dimension at least NKKI, the total ion count.
C                 cgs units, GM/CM*S
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (AVAG = 6.022D23, ECHRG=1.6022D-12, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (AVAG = 6.022E23, ECHRG=1.6022E-12, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL, BIG, EXPARG
      DIMENSION T(*), X(*), XNUIK(KK,*), ICKWRK(*), RCKWRK(*),
     1          SIGI(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      DO 100 KI = 1, NKKI
         K = ICKWRK(IcKI + KI - 1)
         XDEN = TDEN * MAX(X(K),SMALLX)
         SIGI(KI) = 0.0
         DO 50 J = 1, NKK
            WTI = RCKWRK(NcWT + K - 1)
            WTJ = RCKWRK(NcWT + J - 1)
            RMASS = (WTI*WTJ)/(WTI+WTJ) / AVAG
            SIGI(KI) = SIGI(KI) + RMASS*XNUIK(J,KI)/(XDEN*ECHRG**2)
50       CONTINUE
         SIGI(KI) = 1.0/SIGI(KI)
100   CONTINUE
C
C     end of SUBROUTINE CKISIG
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKITDE  (ICKWRK, RCKWRK, ITDE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKITDE (ICKWRK, RCKWRK, ITDE)
C  Returns a set of flags indicating whether the reactions are
C  non-thermal, and if so, returns the index of the species on
C  which the reaction depends.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ITDE(*)   - Integer array, electron-impact flags for reactions;
C              dimension at least II, the total reaction count.
C              ITDE(I)= -1  reaction I is not a third-body reactions
C              ITDE(I)=  K  reaction I is a third-body reaction with
C                        temperature dependence on species # K
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), ITDE(*)
C
      DO 100 I = 1, NII
         ITDE(I) = -1
  100 CONTINUE
      DO 150 N = 0, NTDE - 1
         ITDE(ICKWRK(IcTD + N)) = ICKWRK(IcTK + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKITDE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKITR  (ICKWRK, RCKWRK, ITHB, IREV)
C
C  START PROLOGUE
C
C  SUBROUTINE CKITR  (ICKWRK, RCKWRK, ITHB, IREV)
C  Returns a set of flags indicating whether the reactions are
C  reversible or whether they contain arbitrary third bodies
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  ITHB(*)   - Integer array, third-body indices for reactions;
C              dimension at least II, the total reaction count.
C              ITHB(I)= -1  reaction I is not a third-body reactions
C              ITHB(I)=  0  reaction I is is a third-body reaction with
C                           no enhanced third body efficiencies
C              ITHB(I)=  N  reaction I is a third-body reaction with
C                        N species enhanced third-body efficiencies.
C
C  IREV(*)   - Integer array, reversibility indices and species
C              count (reactants plus products) for reactions;
C              dimension at least II, the total reaction count.
C              IREV(I)=+N, reversible reaction I has N species
C              IREV(I)=-N, irreversible reaction I has N species
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), ITHB(*), IREV(*)
C
      DO 100 I = 1, NII
         IREV(I) = ICKWRK(IcNS + I - 1)
         ITHB(I) = -1
  100 CONTINUE
      DO 150 N = 0, NTHB - 1
         ITHB(ICKWRK(IcTB + N)) = ICKWRK(IcKN + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKITR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIVIS (P, T, X, XNUIM, K, ICKWRK, RCKWRK, VISI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIVIS (P, T, X, XNUIM, K, ICKWRK, RCKWRK, VISI)
C  Returns the ion species viscosities given collision frequencies.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real workspace array; dimension at least LENRCK.
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array of mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  XNUIM     - Real scalar, momentum-transfer collision frequency
C              for an ion
C  K         - Integer scalar, species index of the ion
C  OUTPUT
C  VISI      - Real scalar, ion viscosity
C                 cgs units, GM/CM*S
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, AVAG = 6.022D23)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, AVAG = 6.022E23)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
C
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
      CALL CKMMWX (X, ICKWRK, RCKWRK, WTM)
      TDEN = (RHO/WTM)*AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      TI = T(ICKWRK(IcKTF + K - 1))
      XDEN = TDEN * MAX(X(K),SMALL)
      VISI = (4.*BOLTZ*TI*XDEN)/(PI*MAX(XNUIM,SMALL))
C
C     end of SUBROUTINE CKIVIS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKIXSM  (ICKWRK, IXSM, IXSK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKIXSM (ICKWRK, IXSM, IXSK)
C  Returns a set of flags indicating whether the reactions are ion
C  momentum-transfer cross sections.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  IXSM(*)   - Integer array, ion momentum-transfer cross-section flag;
C              dimension at least II, the total reaction count.
C              IXSM(I)= -1  reaction I is not a ion mom-transfer x-sec
C              IXSM(I)=  KI reaction I is a ion mom-trans cross-section
C                        and KI is the ion species index
C  IXSK(*)   - Integer array, species indices for the collision partner
C              of the ion momentum-transfer cross-section reactions;
C              dimension at least II, the total reaction count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), IXSM(*), IXSK(*)
C
      DO 100 I = 1, NII
         IXSM(I) = -1
         IXSK(I) = 0
  100 CONTINUE
      DO 150 N = 0, NXSM - 1
         IXSM(ICKWRK(IcXS + N)) = ICKWRK(IcXI + N)
         IXSK(ICKWRK(IcXS + N)) = ICKWRK(IcXK + N)
  150 CONTINUE
C
C     end of SUBROUTINE CKIXSM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKFKR (P, T, X, ICKWRK, RCKWRK, FWDK, REVK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKFKR (P, T, X, ICKWRK, RCKWRK, FWDK, REVK)
C  Returns the forward and reverse reaction rates for reactions
C  given pressure, temperature(s) and mole fractions.
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  FWDK(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  REVK(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), FWDK(*), REVK(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKRATX (T, RCKWRK(NcK1), ICKWRK(IcNU), ICKWRK(IcNK),
     1             RCKWRK(NcCO), ICKWRK(IcFL), ICKWRK(IcFO),
     2             ICKWRK(IcFT), ICKWRK(IcKF), RCKWRK(NcFL),
     3             ICKWRK(IcTB), ICKWRK(IcKN), RCKWRK(NcKT),
     4             ICKWRK(IcKT), RCKWRK(NcKF), RCKWRK(NcKR), FWDK,
     5             REVK, RCKWRK(NcI3), ICKWRK(IcRNU), RCKWRK(NcRNU),
     6             ICKWRK(IcORD), ICKWRK(IcKOR), RCKWRK(NcKOR),
     7             ICKWRK(IcMO), ICKWRK(IcXS))
C
C     end of SUBROUTINE CKKFKR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKFRT (P, T, ICKWRK, RCKWRK, RKFT, RKRT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKFRT (P, T, ICKWRK, RCKWRK, RKFT, RKRT)
C  Returns the forward and reverse reaction rates for reactions
C  given pressure and temperature(s).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), RKFT(*), RKRT(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RKFT, RKRT,
     9             RCKWRK(NcI1))
C
C     end of SUBROUTINE CKKFRT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKKTFL (ICKWRK, KTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKKTFL (ICKWRK, KTFL)
C  Allows the user to assign a location in the temperature array
C  to use for each gas-phase species.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  KTFL(*)   - Integer array, indices into the temperature(s) for
C              species;
C              dimension at least KK, the total species count.
C              Default value stored in ICKWRK is set to 1 in CKINIT.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), KTFL(*)
C
      DO 100 K = 1, NKK
         ICKWRK(IcKTF + K - 1) = KTFL(K)
  100 CONTINUE
C
C     end of SUBROUTINE CKKTFL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKLEN (LINC, LOUT, LI, LR, LC, IFLAG)
C
C  START PROLOGUE
C
C  SUBROUTINE CKLEN (LINC, LOUT, LENI, LENR, LENC, IFLAG)
C   Returns the lengths required for work arrays.
C
C  INPUT
C  LINC     - Integer scalar, input file unit for the linkfile.
C  LOUT     - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  LENI     - Integer scalar, minimum length required for the
C             integer work array.
C  LENR     - Integer scalar, minimum length required for the
C             real work array.
C  LENC     - Integer scalar, minimum length required for the
C             character work array.
C  IFLAG    - Integer scalar, indicates successful reading of
C             linkfile; IFLAG>0 indicates error type.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NLIST = 1)
      LOGICAL KERR, VOK, POK, LBIN
      CHARACTER*16 LIST(NLIST), FILVER, PREC, PRVERS, IFMT, CFMT,
     1             RFMT, LFMT
      PARAMETER (IFMT='(10I12)', CFMT='(8A16)', RFMT='(1P,5E24.16)',
     1           LFMT='(L8)')
C
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
      DATA LIST(1) /'1.0'/
C
      FILVER = ' '
      PRVERS   = ' '
      PREC = ' '
      LENI = 0
      LENR = 0
      LENC = 0
C
      KERR = .FALSE.
      IFLAG = 0
      REWIND LINC
      NREC = 1
C*****linkfile (gas) > binary
C      LBIN = .TRUE.
C*****END linkfile (gas) > binary
C*****linkfile (gas) > ascii
      LBIN = .FALSE.
C*****END linkfile (gas) > ascii
C
      IF (LBIN) THEN
         READ (LINC, ERR=100) FILVER
      ELSE
         READ (LINC, CFMT, ERR=100) FILVER
      ENDIF
      CALL CKCOMP (FILVER, LIST, NLIST, IND)
      IF (IND .LE. 0) THEN
         VOK = .FALSE.
      ELSE
         VOK = .TRUE.
      ENDIF
C
      IF (LBIN) THEN
         NREC = 2
         READ (LINC, ERR=100) PRVERS
         NREC = 3
         READ (LINC, ERR=100) PREC
         NREC = 4
         READ (LINC, ERR=100) KERR
      ELSE
         NREC = 2
         READ (LINC, CFMT, ERR=100) PRVERS
         NREC = 3
         READ (LINC, CFMT, ERR=100) PREC
         NREC = 4
         READ (LINC, LFMT, ERR=100) KERR
      ENDIF
C
      POK = .FALSE.
C*****precision > double
      IF (INDEX(PREC, 'DOUB') .GT. 0) POK = .TRUE.
C*****END precision > double
C*****precision > single
C      IF (INDEX(PREC, 'SING') .GT. 0) POK = .TRUE.
C*****END precision > single
C
      IF (KERR .OR. (.NOT.POK) .OR. (.NOT.VOK)) THEN
         IF (KERR) THEN
            WRITE (LOUT,'(/A,/A)')
     1      ' There is an error in the Chemkin linkfile...',
     2      ' Check CHEMKIN INTERPRETER output for error conditions.'
         ENDIF
         IF (.NOT. VOK) WRITE (LOUT, '(/A)')
     1   ' Chemkin linkfile is incompatible with Chemkin-III Library'
         IF (.NOT. POK) THEN
            WRITE (LOUT,'(/A,A)')
     1      ' Precision of Chemkin linkfile does not agree with',
     2      ' precision of Chemkin library'
         ENDIF
         IFLAG = 20
         REWIND LINC
         RETURN
      ENDIF
C
      NREC = 5
      IF (LBIN) THEN
         READ (LINC, ERR=100) LENICK, LENRCK, LENCCK
      ELSE
         READ (LINC, IFMT, ERR=100) LENICK, LENRCK, LENCCK
      ENDIF
      REWIND LINC
C
      LENI = LENICK
      LENR = LENRCK
      LENC = LENCCK
      LI   = LENI
      LR   = LENR
      LC   = LENC
      RETURN
C
  100 CONTINUE
      IFLAG = NREC
      WRITE (LOUT, '(/A,/A,I5)')
     1   ' Error reading linkfile,',
     2   ' SUBROUTINE CKLEN record index #', IFLAG
      REWIND LINC
C
C     Generic Formats - limit lines to 132 characters
C
8001  FORMAT (10I12)
8002  FORMAT (1P,5E24.16)
8003  FORMAT (8A16)
C
C     end of SUBROUTINE CKLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKDLIM (STRING, DELIM, I1, I2)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDLIM
C  returns pointers into a character string of the first and
C  second occurrences of a particular character.
C
C  Arguments:
C  STRING - Character string.
C  DELIM  - Single character.
C  I1     - Integer scalar, location in STRING of first DELIM.
C  I2     - Integer scalar, location in STRING of second DELIM.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), DELIM*1
C
      I1 = INDEX(STRING, DELIM)
      I2 = I1 + INDEX(STRING(I1+1:), DELIM)
C
C     end of SUBROUTINE CKDLIM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKDTAB (STRING)
C
C  START PROLOGUE
C
C  SUBROUTINE CKDTAB (STRING)
C  Replaces any tab character in a character string with one space.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N0
C*****END precision > single
C
      CHARACTER STRING*(*), TAB*1
C
      TAB = CHAR(9)
   10 CONTINUE
      IND = INDEX(STRING,TAB)
      IF (IND .GT. 0) THEN
         STRING(IND:IND) = ' '
         GO TO 10
      ENDIF
C
C     end of SUBROUTINE CKDTAB
      RETURN
      END
C
C----------------------------------------------------------------------C
C                                                                      C
      CHARACTER*(*) FUNCTION CKCHUP(ISTR, ILEN)
      CHARACTER*(*) ISTR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      CKCHUP = ' '
      CKCHUP = ISTR(1:ILEN)
      JJ = MIN (LEN(CKCHUP), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 05 N = 1,26
            IF (ISTR(J:J) .EQ. LCASE(N)) CKCHUP(J:J) = UCASE(N)
   05    CONTINUE
   10 CONTINUE
C
C     end of FUNCTION CKCHUP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      CHARACTER*(*) FUNCTION CKCHLO(ISTR, ILEN)
      CHARACTER*(*) ISTR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      CKCHLO = ' '
      CKCHLO = ISTR(1:ILEN)
      JJ = MIN (LEN(CKCHLO), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 05 N = 1,26
            IF (ISTR(J:J) .EQ. UCASE(N)) CKCHLO(J:J) = LCASE(N)
   05    CONTINUE
   10 CONTINUE
C
C     end of FUNCTION CKCHLO
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      INTEGER FUNCTION CKLKUP (ITEM, LIST, NLIST)
C
C  START PROLOGUE
C
C  INTEGER FUNCTION CKLKUP (ITEM, LIST, NLIST)
C
C      Looks up an item in an integer list. If an item is found,
C   it returns the first position of the item in the list. If an
C   item is not found, this routine returns the value 0.
C
C  INPUT
C
C  ITEM    - Integer scalar; Item to look up in the list
C  LIST(*) - Integer array;  List of entries
C  NLIST   - Integer scalar; Number of entries in the list
C
C  END PROLOGUE
C
      INTEGER   ITEM, LIST, NLIST, I
      DIMENSION  LIST(NLIST)
C
      CKLKUP = 0
      DO 10 I = NLIST, 1, -1
        IF (LIST(I) .EQ. ITEM) CKLKUP = I
   10 CONTINUE
C
C     end of FUNCTION CKLKUP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION CKFRCH (STR)
C
C  START PROLGUE
C
C  INTEGER FUNCTION CKFRCH (STR)
C
C  Returns the index of the first non-blank, non-tab character in
C  a string.
C
C  INPUT
C  STR   - Character string
C
C  END PROLOGUE
C
      CHARACTER STR*(*), TAB*1
      INTEGER ILEN, I
C
      ILEN = LEN(STR)
      TAB  = CHAR(9)
      CKFRCH = 0
      DO 10 I = 1, ILEN
         IF (STR(I:I).EQ.' ' .OR. STR(I:I).EQ.TAB) GO TO 10
         CKFRCH = I
         RETURN
   10 CONTINUE
C
C     end of FUNCTION CKFRCH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION CKLSCH (STR)
C
C  START PROLOGUE
C
C  INTEGER FUNCTION CKLSCH (STR)
C
C  Returns the index of the final non-blank, non-tab character in
C  a string.
C
C  INPUT
C  STR   - Character string
C
C  END PROLOGUE
C
      CHARACTER STR*(*), TAB*1, NUL*1
      INTEGER ILEN, I
C
      ILEN = LEN(STR)
      CKLSCH = 0
      TAB = CHAR(9)
      NUL = CHAR(0)
      DO 10 I = ILEN, 1, -1
         IF (STR(I:I).EQ.' ' .OR. STR(I:I).EQ.TAB .OR.
     1       STR(I:I).EQ.NUL) GO TO 10
         CKLSCH = I
         RETURN
   10 CONTINUE
C
C     end of FUNCTION CKLSCH
      RETURN
      END
C                                                                      C
      INTEGER FUNCTION CKSLEN (LINE)
C
C  BEGIN PROLOGUE
C
C  INTEGER FUNCTION CKSLEN (LINE)
C  Returns the effective length of a character string, i.e.,
C  the index of the last character before an exclamation mark (!)
C  indicating a comment.
C
C  INPUT
C  LINE     - Character string.
C
C  OUTPUT
C  CKSLEN   - Integer scalar, the effective length of LINE.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*)
      INTEGER CKLSCH, CKFRCH
      EXTERNAL CKLSCH, CKFRCH
C
      IND = CKFRCH(LINE)
      IF (IND.EQ.0 .OR. LINE(IND:IND).EQ.'!') THEN
         CKSLEN = 0
      ELSE
         IND = INDEX(LINE,'!')
         IF (IND .GT. 0) THEN
            CKSLEN = CKLSCH(LINE(1:IND-1))
         ELSE
            CKSLEN = CKLSCH(LINE)
         ENDIF
      ENDIF
C
C     end of FUNCTION CKSLEN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXMIN (X, NN, XMIN, IMIN)
C
C  START PROLOGUE
C  Returns the minimum value in an array and its location in the array.
C
C  INPUT
C  X(*)      - Real array.
C  NN        - Integer scalar; size of X.
C  OUTPUT
C  XMIN      - Real scalar.
C  IMIN      - Integer scalar; location in X of XMIN.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(NN)
C
      IMIN = 1
      XMIN = X(IMIN)
      DO 10 N = 1, NN
         IF (X(N) .LT. XMIN) THEN
            XMIN = X(N)
            IMIN = N
         ENDIF
   10 CONTINUE
C
C     end of CKXMIN
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXMAX (X, NN, XMAX, IMAX)
C
C  START PROLOGUE
C  Returns the maximum value in an array and its location in the array.
C
C  INPUT
C  X(*)      - Real array.
C  NN        - Integer scalar; size of X.
C  OUTPUT
C  XMAX      - Real scalar.
C  IMAX      - Integer scalar; location in X of XMAX.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION X(NN)
C
      IMAX = 1
      XMAX = X(IMAX)
      DO 10 N = 1, NN
         IF (X(N) .GT. XMAX) THEN
            XMAX = X(N)
            IMAX = N
         ENDIF
   10 CONTINUE
C
C     end of CKXMAX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWC (C, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWC (C, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given molar
C  concentrations;  see Eq. (5).
C
C  INPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION C(*), ICKWRK(*), RCKWRK(*)
C
      CTOT = 0.0
      WTM  = 0.0
      DO 100 K = 1, NKK
         CTOT = CTOT + C(K)
         WTM  = WTM  + C(K) * RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      WTM = WTM / CTOT
C
C     end of SUBROUTINKE CKMMWC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWX (X, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given mole
C  fractions;  see Eq. (4).
C
C  INPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION X(*), ICKWRK(*), RCKWRK(*)
C
      WTM = 0.0
      DO 100 K = 1, NKK
         WTM = WTM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKMMWX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMMWY (Y, ICKWRK, RCKWRK, WTM)
C  Returns the mean molecular weight of the gas mixture given mass
C  fractions;  see Eq. (3).
C
C  INPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WTM       - Real scalar, mean molecular weight of the mixture.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*)
C
      SUMYOW=0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      WTM = 1.0 / SUMYOW
C
C     end of SUBROUTINE CKMMWY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKMXTP (ICKWRK, MAXTP)
C
C  START PROLOGUE
C
C  SUBROUTINE CKMXTP (ICKWRK, MAXTP)
C  Returns the maximum number of temperatures used in fitting the
C  thermodynamic properties of the species.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  MXTP      - Integer scalar, maximum number of temperatures used
C              to fit the thermodynamic properties of the species.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*)
C
      MAXTP = MXTP
C
C     end of SUBROUTINE CKMXTP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNCF  (MDIM, ICKWRK, RCKWRK, NCF)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNCF  (MDIM, ICKWRK, RCKWRK, NCF)
C  Returns the elemental composition of the species
C
C  INPUT
C  MDIM      - Integer scalar, first dimension of the matrix NCF;
C              MDIM must be equal to or greater than MM, the total
C              element count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NCF(*,*)  - Real matrix, the elemental composition of the species;
C              dimension at least MM for the first, the total element
C              count, and at least KK for the second, the total species
C              count.
C              NCF(M,K) is the quantity of the element M in species K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NCF(MDIM,*)
C
      I_KNCF = IcNC - 1
      DO 150 K = 1, NKK
         DO 100 M = 1, NMM
            I_KNCF = I_KNCF + 1
            NCF(M,K) = ICKWRK(I_KNCF)
  100    CONTINUE
150   CONTINUE
C
C     end of SUBROUTINE CKNCF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNORM (ARRAY, NN)
C
C START PROLOGUE
C
C SUBROUTINE CKNORM (ARRAY, NN)
C Utility to normalize the real members of an array.
C
C INPUT
C ARRAY(*)  - Real array.
C NN        - Integer scalar; the size of ARRAY.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      SUM = CKSUM (ARRAY, NN)
      IF (SUM .NE. 0.0) THEN
         DO 20 N = 1, NN
            ARRAY(N) = ARRAY(N) / SUM
   20    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKNORM
      RETURN
      END
C                                                                      C
      SUBROUTINE CKSCAL (ARRAY, NN, SCAL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSCAL (ARRAY, NN)
C  Utility to scale the real members of an array.
C
C  INPUT
C  ARRAY(*)  - Real array.
C  NN        - Integer scalar; the size of ARRAY.
C  SCAL      - Real scalar; the multiplier for ARRAY members.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      DO 10 N = 1, NN
         ARRAY(N) = SCAL * ARRAY(N)
   10 CONTINUE
C
C     end of SUBROUTINE CKSCAL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C*****precision > double
      DOUBLE PRECISION FUNCTION CKSUM (ARRAY, NN)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      REAL FUNCTION CKSUM (ARRAY, NN)
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION ARRAY(NN)
C
      CKSUM = 0.0
      DO 10 N = 1, NN
         CKSUM = CKSUM + ARRAY(N)
   10 CONTINUE
C
C     end of FUNCTION CKSUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNPAR (LINE, NPAR, LOUT, IPAR, ISTART, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNPAR (LINE, NPAR, LOUT, IPAR, ISTART, KERR)
C  Searches a character string LINE from last to first character,
C  to create a substring IPAR containing NPAR blank-delimited numbers;
C  ISTART is the column of LINE containing IPAR. This allows format-
C  free input of combined alpha-numeric data.  For example,
C
C     input:  LINE*80   = "t1 t2 dt  300.0  3.0E3  50"
C             NPAR      = 3, the number of substrings requested
C             LOUT      = 6, a logical unit number on which to write
C                         diagnostic messages.
C     output: IPAR*80   = "300.0  3.0E3  50"
C             ISTART    = 13, the starting column in LINE of the
C                         NPAR substrings
C             KERR      = .FALSE.
C
C  INPUT
C  LINE      - Character string; length determined by calling routine.
C  NPAR      - Integer scalar, number of substrings expected.
C  LOUT      - Integer scalar, output unit for printed diagnostics.
C
C  OUTPUT
C  IPAR      - Character string, subset of LINE, containing only the
C              NPAR substrings.
C  ISTART    - Integer scalar, starting location in LINE of the NPAR
C              substrings.
C  KERR      - Logical, character length or syntax error flag.
C
C  END PROLOGUE
C
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), IPAR*(*)
      LOGICAL FOUND, KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
C----------Find Comment String (! signifies comment)
C
      ILEN = CKSLEN(LINE)
      KERR = .FALSE.
C
      IF (ILEN.GT.0) THEN
         FOUND = .FALSE.
         N = 0
         DO 40 I = ILEN, 1, -1
            IF (FOUND) THEN
               IF (LINE(I:I).EQ.' ') THEN
                  N = N+1
                  FOUND = .FALSE.
                  IF (N.EQ.NPAR) THEN
                     ISTART = I+1
                     L1 = ILEN - ISTART + 1
                     L2 = LEN(IPAR)
                     IF (L2 .GE. L1) THEN
                        IPAR = LINE(ISTART:ILEN)
                     ELSE
                        WRITE (LOUT,*)
     1               ' Error in CKNPAR...character length too small...'
                        KERR = .TRUE.
                     ENDIF
                     GO TO 100
                  ENDIF
               ENDIF
            ELSE
               IF (LINE(I:I).NE.' ') FOUND = .TRUE.
            ENDIF
   40    CONTINUE
      ENDIF
C
      WRITE (LOUT,*) ' Error in CKNPAR...',NPAR,' values not found...'
      KERR = .TRUE.
C
C     end of SUBROUTINE CKNPAR
  100 CONTINUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNU   (KDIM, ICKWRK, RCKWRK, NUKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNU   (KDIM, ICKWRK, RCKWRK, NUKI)
C  Returns the stoichiometric coefficients of the reactions;
C  see Eq. (50).
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix NUKI;
C              KDIM must be greater than or equal to KK, the total
C              species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NUKI(*,*) - Integer matrix, stoichiometric coefficients of the
C              species in the reactions;  dimension at least KK for
C              the first, the total species count, and at least II
C              for the second, the total reaction count.
C              NUKI(K,I) is the stoichiometric coefficient of
C              species K in reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NUKI(KDIM,*)
C
      DO 100 I = 1, NII
         DO 50 K = 1, NKK
            NUKI(K,I) = 0
   50    CONTINUE
  100 CONTINUE
C
      I_NK = IcNK - MXSP
      I_NU = IcNU - MXSP
      DO 200 I = 1, NII
         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
         DO 150 N = 0, MXSP - 1
            K = ICKWRK(I_NK + N)
            IF (K .NE. 0) NUKI(K,I) = NUKI(K,I) + ICKWRK(I_NU + N)
  150    CONTINUE
200   CONTINUE
C
C     end of SUBROUTINE CKNU
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKNUF   (KDIM, ICKWRK, RCKWRK, NUFKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKNUF   (KDIM, ICKWRK, RCKWRK, NUKI)
C  Returns the forward stoichiometric coefficients for reactions;
C  by definition, reactants' coefficients are negative;  see Eq. (50).
C  Contrast this subroutine with subroutine CKNU, which returns the
C  net stoichiometric coefficients for a reaction.
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix NUKI;
C              KDIM must be greater than or equal to KK, the total
C              species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  NUFKI(*,*)- Integer matrix, stoichiometric coefficients of the
C              species in the forward direction of the reactions
C              (reactants only); dimension at least KK in the first,
C              the total species count, and at least II for the
C              second, the total reaction count.
C              NUKI(K,I) is the stoichiometric coefficient of
C              species K in forward direction of reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), NUFKI(KDIM,*)
C
      DO 100 I = 1, NII
         DO 50 K = 1, NKK
            NUFKI(K,I) = 0
   50    CONTINUE
  100 CONTINUE
C
      I_NK = IcNK - MXSP
      I_NU = IcNU - MXSP
      NREAC = MXSP / 2
      DO 200 I = 1, NII
         I_NK = I_NK + MXSP
         I_NU = I_NU + MXSP
         DO 150 N = 0, NREAC - 1
            K = ICKWRK(I_NK + N)
            IF (K .EQ. 0) GO TO 200
            NUFKI(K,I) = NUFKI(K,I) + ICKWRK(I_NU + N)
  150    CONTINUE
  200 CONTINUE
C
C     end of SUBROUTINE CKNUF
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPC   (RHO, T, C, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and molar concentrations;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*)
C
      CTTOT = 0.0
      DO 100 K = 1, NKK
         CTTOT = CTTOT + C(K) * T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
      P    = RCKWRK(NcRU) * CTTOT
C
C     end of SUBROUTINE CKPC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPHAZ (ICKWRK, RCKWRK, KPHASE)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPHAZ (ICKWRK, RCKWRK, KPHASE)
C  Returns a set of flags indicating phases of the species
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  KPHASE(*) - Integer array, phases of the species;
C              dimension at least KK, the total species count.
C              KPHASE(K)=-1, species K is solid
C              KPHASE(K)= 0, species K is gaseous
C              KPHASE(K)=+1, species K is liquid
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), KPHASE(*)
C
      DO 100 K = 1, NKK
         KPHASE(K) = ICKWRK(IcPH + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKPHAZ
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPNT (LSAVE, LOUT, NPOINT, V, P, LI, LR, LC, IERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPNT (LSAVE, LOUT, NPOINT, VERS, PREC, LENI, LENR,
C                    LENC, KERR)
C  Reads from a file information about a Chemkin linkfile, and
C  pointers for work arrays.
C
C  INPUT
C  LSAVE     - Integer scalar, input unit for binary data file.
C  LOUT      - Integer scalar, formatted output file unit number.
C
C  OUTPUT
C  NPOINT    - Integer scalar, total pointers count.
C  VERS      - Real scalar, version number of the Chemkin linkfile.
C  PREC      - Character string, machine precision of the linkfile.
C  LENI      - Integer scalar, length required for integer work array.
C  LENR      - Integer scalar, length required for real work array.
C  LENC      - Integer scalar, length required for character work array.
C  KERR      - Logical, error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
C
C     Data about the machine dependent constants is carried in
C
      COMMON/MACH/SMALL,BIG,EXPARG
      LOGICAL KERR, IERR
      CHARACTER*16 PREC, FILVER, PRVERS, P, V
C
C      THIS STATEMENT WILL NOT COMPILE, MACHINE-DEPENDENT CONSTANTS
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
      EXPARG = LOG(BIG)
C
      KERR = .FALSE.
      READ (LSAVE, ERR=100)
     *                FILVER,   PREC,   LENI,   LENR,   LENC,
C
C     include file for CHEMKIN-III cklib.f, dated: March 1, 1966
C
C     Integer constants
C
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     4   MXORD, KEL, NKKI,
C
C     Integer pointers to character arrays in CCKWRK
C
     5   IcMM, IcKK,
C
C     Integer pointers to integer arrays in ICKWRK
C
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL,
     7   IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, IcEI,
     8   IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, IcTD,
     9   IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
C
C     Integer pointers to real variables and arrays in RCKWRK
C
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,
     2   NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
C
C     END include file for cklib.f
C
      NPOINT = 102
      V = FILVER
      P = PREC
      LI = LENI
      LR = LENR
      LC = LENC
      IERR = KERR
      RETURN
C
  100 CONTINUE
      WRITE (LOUT, *) ' Error reading Chemkin linkfile data...'
      KERR   = .TRUE.
      IERR   = KERR
      NPOINT = 0
      FILVER   = ' '
      PRVERS     = ' '
      V      = FILVER
      PREC   = ' '
      P      = PREC
C
C     end of SUBROUTINE CKPNT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPX   (RHO, T, X, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPX   (RHO, T, X, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and mole fractions;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      SUMXW = 0.0
      SUMXT = 0.0
      DO 100 K = 1, NKK
         SUMXW = SUMXW + X(K)*RCKWRK(NcWT + K - 1)
         SUMXT = SUMXT + X(K)*T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
      P = RHO * RCKWRK(NcRU) * SUMXT / SUMXW
C
C     end of SUBROUTINE CKPX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
C
C  START PROLOGUE
C
C  SUBROUTINE CKPY   (RHO, T, Y, ICKWRK, RCKWRK, P)
C  Returns the pressure of the gas mixture given mass density,
C  temperature(s) and mass fractions;  see Eq. (1).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      SUMYOW = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUMYOW = SUMYOW + Y(K+1) * T(ICKWRK(IcKTF + K))
     1                            /RCKWRK(NcWT + K)
150   CONTINUE
      P = RHO * RCKWRK(NcRU) * SUMYOW
C
C     end of SUBROUTINE CKPY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQC   (T, C, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQC   (T, C, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given temperature(s)
C  and molar concentrations;  see Eqs. (51) and (58).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKRATX (T, C, ICKWRK(IcNU), ICKWRK(IcNK), RCKWRK(NcCO),
     2             ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcFT),
     3             ICKWRK(IcKF), RCKWRK(NcFL), ICKWRK(IcTB),
     3             ICKWRK(IcKN), RCKWRK(NcKT), ICKWRK(IcKT),
     4             RCKWRK(NcKF), RCKWRK(NcKR), RCKWRK(NcI1),
     5             RCKWRK(NcI2), RCKWRK(NcI3), ICKWRK(IcRNU),
     6             RCKWRK(NcRNU), ICKWRK(IcORD), ICKWRK(IcKOR),
     7             RCKWRK(NcKOR), ICKWRK(IcMO), ICKWRK(IcXS))
C
      NIM1 = NII - 1
      DO 100 I = 0, NIM1
         Q(I+1) = RCKWRK(NcI1 + I) - RCKWRK(NcI2 + I)
  100 CONTINUE
C
C     end of SUBROUTINE CKQC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQNUE (Q, T, P, C, ICKWRK, RCKWRK, XNUES, XNUEH)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQNUE (Q, T, P, C, ICKWRK, RCKWRK, XNUES, XNUEH)
C  Returns the electron momentum-transfer collision frequencies for
C  all species, using averages for unspecified species
C
C  INPUT
C  Q(*)      - Real array, rates of production of the reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm2
C  C(*)      - Real array, concentrations of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  OUTPUT
C  XNUES(*)  - Real array, momentum-transfer collision frequencies for
C              the electrons with the species;
C              dimension at least KK, the total species count.
C  XNUEH     - Real scalar, total momentum-transfer collision frequency
C              for the electrons
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, PERM0 = 8.854D-7, ECHRG=1.6022D-12,
     1           AVAG = 6.022D23, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, PERM0 = 8.854E-7, ECHRG=1.6022E-12,
C     1           AVAG = 6.022E23, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION Q(*), T(*), C(*), ICKWRK(*), RCKWRK(*), XNUES(*)
C
      DO 40 K = 1, NKK
         XNUES(K) = 0.0
 40   CONTINUE
C
      XNUEH = 0.0
      SUMQN = 0.0
      DO 100 I1 = 1, NMOM
         I = ICKWRK(IcMO + I1 - 1)
         K = ICKWRK(IcMK + I1 -1)
         XNUES(K) = Q(I) * C(K)*AVAG
         XNUEH = XNUEH + XNUES(K)
         SUMQN = SUMQN + Q(I)
 100  CONTINUE
      IF (NMOM.NE.0) THEN
         QAV = SUMQN / NMOM
      ELSE
         QAV = 0.0
      ENDIF
C
C    CALCULATE NUEI AND NUEE:
C
      EDEN = C(KEL)*AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      EMASS = RCKWRK(NcWT + KEL - 1)/AVAG
      BOLTZ = RCKWRK(NcRU)/AVAG
      DO 200 KI = 1, NKKI
         K = ICKWRK(IcKI + KI - 1)
         IF (EDEN.GT.1.0) THEN
            TK = T(ICKWRK(IcKTF + K - 1))
            XDEN = C(K)*AVAG
            Z = ABS(FLOAT(ICKWRK(IcCH + K - 1)))
            XLAM = (PERM0 * BOLTZ * TE / ECHRG**2)**1.5
     1           * 12.0 * PI / (Z*SQRT(EDEN))
            XNUES(K) = XDEN * LOG(XLAM) * (4.*SQRT(2.*PI)/3.)
     1           * (EMASS/(BOLTZ*TK))**1.5
     2           * ((ECHRG**2)/(4.*PI*PERM0*EMASS))**2
         ELSE
            XNUES(K) = SMALL
         ENDIF
200   CONTINUE
      IF (EDEN.GT.1.0) THEN
         XNUES(KEL) = EDEN * (8./3.) * SQRT(PI) *
     1                (EMASS/(BOLTZ*TE))**1.5 *
     2            ((ECHRG**2)/(4*PI*PERM0*EMASS))**2 * LOG(XLAM)
      ELSE
         XNUES(KEL) = SMALL
      ENDIF
C
C    FILL IN ANY GAPS IN DATA PROVIDED BY AVERAGE NEUTRAL FREQ
C
      DO 110 K = 1, NKK
         IF (XNUES(K).LE.0.0) THEN
            XNUES(K) = QAV * C(K)*AVAG
            XNUEH = XNUEH + XNUES(K)
         ENDIF
110   CONTINUE
C
C     end of SUBROUTINE CKQNUE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQNUI (Q, T, P, C, KDIM, ICKWRK, RCKWRK, XNUIK, XNUIM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQNUI (Q, T, P, C, KDIM, ICKWRK, RCKWRK, XNUIK, XNUIM)
C  Returns the ion momentum-transfer collision frequencies from
C  collision cross-sections for all species, using averages for
C  unspecified species.
C
C  INPUT
C  Q(*)      - Real array, rates of production for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm2.
C  C(*)      - Real array, concentrations of the mixtures;
C              dimension at least KK, the total species count.
C  KDIM      - Integer scalar, first dimension of the matrix XNUIK;
C              KDIM must be at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  XNUIK(*,*)- Real matrix, momentum-transfer collision frequencies
C              for the ions with all other species;
C              dimension at least KK for the first, the total
C              species count, and at least NKKI for the second, the
C              total ion count.
C  XNUIM(*)  - Real array, momentim-averaged momentum-transfer
C              collision frequencies for the ions;
C              dimension at least NKKI, the total ion count.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (PI = 3.14159265, PERM0 = 8.854D-7, ECHRG=1.6022D-12,
     1           AVAG = 6.022D23, SMALLX = 1.D-50)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (PI = 3.14159265, PERM0 = 8.854E-7, ECHRG=1.6022E-12,
C     1           AVAG = 6.022E23, SMALLX = 1.E-30)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION Q(*), T(*), C(*), ICKWRK(*), RCKWRK(*),
     1          XNUIK(KDIM,*), XNUIM(*)
C
      DO 10 KI = 1, NKKI
         DO 5 K = 1, NKK
            XNUIK(K,KI) = 0.0
 5       CONTINUE
         XNUIM(KI) = 0.0
 10   CONTINUE
C
      EDEN = C(KEL)*AVAG
      TE = T(ICKWRK(IcKTF+KEL-1))
      BOLTZ = RCKWRK(NcRU)/AVAG
      DEBYE = SQRT (PERM0 * BOLTZ * TE / (EDEN * ECHRG**2))
      EMASS = RCKWRK(NcWT+KEL-1)/AVAG
C
      SUMQN = 0.0
      NEUT = 0
      DO 100 N = 1, NXSM
         I = ICKWRK(IcXS + N - 1)
         QXS = Q(I)
         K = ICKWRK(IcXK + N - 1)
         KIK = ICKWRK(IcXI + N - 1)
         DO 50 KI1 = 1, NKKI
            IF (KIK .EQ. ICKWRK(IcKI + KI1 - 1)) KI = KI1
 50      CONTINUE
         XMSS = RCKWRK(NcWT + KIK - 1)/AVAG
         TI = T(ICKWRK(IcKTF + KIK - 1))
         TK = T(ICKWRK(IcKTF + K - 1))
C         RMASS = (RCKWRK(NcWT + KIK - 1)*RCKWRK(NcWT + K - 1))
C     1              / (RCKWRK(NcWT+KIK-1) + RCKWRK(NcWT+K-1)) / AVAG
         XDEN = C(K)*AVAG
         XNUIK(K,KI) = XDEN * SQRT(3.*BOLTZ*TI/XMSS) * QXS
         SUMQN = SUMQN + QXS
         NEUT = NEUT + 1
         WTFAC = 2.*RCKWRK(NcWT + K - 1)
     1                / (RCKWRK(NcWT + KIK - 1)+RCKWRK(NcWT + K - 1))
         XNUIM(KI) = XNUIM(KI) + XNUIK(K,KI)*WTFAC
 100  CONTINUE
C
C  USE AVERAGE X-SECTION OF NEUTRAL SPECIES TO DETERMINE MISSING DATA
C
      QNAVG = SUMQN / NEUT
      DO 200 KI = 1, NKKI
         KIK = ICKWRK(IcKI + KI - 1)
         XMSS = RCKWRK(NcWT + KIK - 1)/AVAG
         TI = T(ICKWRK(IcKTF + KIK - 1))
         DO 150 K = 1, NKK
            IF (ICKWRK(IcCH+K-1).EQ.0 .AND. XNUIK(K,KI) .LE. 0.0) THEN
               XDEN = C(K)*AVAG
               XNUIK(K,KI) = XDEN * SQRT(3.*BOLTZ*TI/XMSS) * QNAVG
               WTFAC = 2.*RCKWRK(NcWT + K - 1)
     1                  / (RCKWRK(NcWT+KI-1)+RCKWRK(NcWT + K - 1))
               XNUIM(KI) = XNUIM(KI) + XNUIK(K,KI)*WTFAC
            ELSEIF (ICKWRK(IcCH + K - 1).NE.0 .AND. K.NE.KEL) THEN
               XDEN = C(K)*AVAG
               TK = T(ICKWRK(IcKTF + K - 1))
               B0FAC = ECHRG**2 / (12. * PI * PERM0 * BOLTZ * TK)
               B0BAR = ABS(ICKWRK(IcCH + K - 1))
     1                 * ABS(ICKWRK(IcCH + KIK - 1)) * B0FAC
               XLAM = DEBYE / B0BAR
               XNUIK(K,KI) = XDEN * (8./3.) * SQRT(PI) *
     1                      (EMASS/(BOLTZ*TK))**1.5 *
     2                      (ECHRG**2/(4*PI*PERM0*EMASS))**2 * LOG(XLAM)
            ELSEIF (K.EQ.KEL) THEN
               IF (EDEN.GT.1.) THEN
                  XIDEN = C(KIK)*AVAG
                  Z = ABS(ICKWRK(IcCH + KIK - 1))
                  XLAM = (PERM0 * BOLTZ * TE / (ECHRG**2))**1.5
     1                * 12.0 * PI / (Z * SQRT(EDEN))
                  XNUIK(KEL,KI) = XIDEN * LOG(XLAM)
     1                         * (4.0*SQRT(2.*PI)/3.0)
     2                         * (EMASS/(BOLTZ*TI))**1.5
     3                         * ((ECHRG**2)/(4.*PI*PERM0*EMASS))**2
               ELSE
                  XNUIK(KEL,KI) = SMALL
               ENDIF
            ENDIF
 150     CONTINUE
 200  CONTINUE
C
C     end of SUBROUTINE CKQNUI
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQXP  (P, T, X, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQXP  (P, T, X, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given pressure,
C  temperature(s) and mole fractions;  see Eqs. (51) and (58).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQYP  (P, T, Y, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQYP  (P, T, Y, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given pressure,
C  temperature(s) and mass fractions;  see Eqs. (51) and (58).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, Mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKQYR  (RHO, T, Y, ICKWRK, RCKWRK, Q)
C
C  START PROLOGUE
C
C  SUBROUTINE CKQYR  (RHO, T, Y, ICKWRK, RCKWRK, Q)
C  Returns the rates of progress for reactions given mass density,
C  temperature(s) and mass fractions; see Eqs. (51) and (58).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Q(*)      - Real array, rates of progress for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), Q(*)
C
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKQC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, Q)
C
C     end of SUBROUTINE CKQYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKR2CH (RNUM, STR, I, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKR2CH (RNUM, STR, I, KERR)
C  Returns a character string representation of a real number
C  and the effective length of the string.
C
C  INPUT
C  RNUM      - Real scalar, to be converted to a string;
C              the maximum magnitude of RNUM is machine-dependent.
C
C  OUTPUT
C  STR      - Character string, left-justified representation of RNUM;
C             i.e., RNUM=  0.0      returns STR=" 0.00"
C                   RNUM= -10.5     returns STR="-1.05E+01"
C                   RNUM= 1.86E-100 returns in STR=" 1.86E-100"
C                   the minimum length of STR required is 5
C  I        - Integer scalar, total non-blank characters in RNUM.
C  KERR     - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STR*(*)
      LOGICAL KERR, IERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C*****exponent range > +/-30
C      SMALL = 1.0E-30
C      BIG   = 1.0E+30
C*****END exponent range > +/-30
C*****exponent range > +/-300
      SMALL = 10.0D0**(-300)
      BIG   = 10.0D0**(+300)
C*****END exponent range > +/-300
C
      ILEN = LEN(STR)
      STR = ' '
      KERR = .FALSE.
      I = 0
C
      RDUM = ABS(RNUM)
      IF (ILEN .LT.10) THEN
         IF (RDUM.GE.9.995 .OR. RDUM.LT.0.995 .OR. ILEN.LT.5) THEN
            KERR = .TRUE.
            RETURN
         ENDIF
      ENDIF
C
      STR = ' 0.00'
      I = 5
      IF (RDUM .EQ. 0.0) RETURN
C
C     convert RDUM to a value between 1.0 and 10.0
C
      IF (RDUM.GT.BIG .OR. RDUM.LT. SMALL) THEN
         KERR = .TRUE.
         RETURN
      ENDIF
C
      IE = 0
   10 CONTINUE
C
C     find a real number 10 < rnum <= 1
C
      IF (RDUM .GE. 10.0) THEN
         IE = IE - 1
         RDUM = RDUM / 10.0
         GO TO 10
C
      ELSEIF (RDUM .LT. 1.0) THEN
         IE = IE + 1
         RDUM = RDUM * 10.0
         GO TO 10
      ENDIF
C
      IF (NINT(100.0*RDUM) .GT. 100.0*RDUM) THEN
         RDUM = RDUM + .005
         IF (RDUM .GE. 10.0) GO TO 10
      ENDIF
      IVAL = INT(RDUM)
      IF (RNUM .LT. 0.0) THEN
         STR(1:) = '- .0'
      ELSE
         STR(1:) = '  .0'
      ENDIF
      CALL CKI2CH (IVAL, STR(2:2), LT, IERR)
      NREM = 100.0*RDUM - 100*IVAL
      IF (NREM .GE. 10) THEN
         CALL CKI2CH (NREM, STR(4:5), LT, IERR)
      ELSE
         CALL CKI2CH (NREM, STR(5:5), LT, IERR)
      ENDIF
C
      IF (ABS(IE) .NE. 0) THEN
         IF (IE .LT. 0) THEN
            STR(6:8) = 'E+0'
         ELSE
            STR(6:8) = 'E-0'
         ENDIF
         I = 9
         IF (ABS(IE) .GE. 10) THEN
            CALL CKI2CH (ABS(IE), STR(8:), LT, IERR)
            IF (ABS(IE) .GT. 99) I = 10
         ELSE
            CALL CKI2CH (ABS(IE), STR(9:9), LT, IERR)
         ENDIF
      ENDIF
      I = CKLSCH(STR)
C
C     end of SUBROUTINE CKR2CH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRAEX (I, RCKWRK, RA)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRAEX (I, RCKWRK, RA)*
C  Get/put the Pre-exponential coefficient of the Ith reaction
C
C  INPUT
C  I         - Integer scalar, reaction index;
C              I > 0 gets RA(I) from RCKWRK
C              I < 0 puts RA(I) into RCKWRK
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  If I < 1:
C  RA        - Real scalar, pre-exponential coefficient for reaction I.
C                 cgs units, mole-cm-sec-K
C
C  OUTPUT
C  If I > 1:
C  RA        - Real scalar, pre-exponential coefficient for reaction I.
C                 cgs units, mole-cm-sec-K
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*)
C
      NI = NcCO + (IABS(I)-1)*(NPAR+1)
      IF (I .GT. 0) THEN
         RA = RCKWRK(NI)
      ELSE
         RCKWRK(NI) = RA
      ENDIF
C
C     end of SUBROUTINE CKRAEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRAT  (RCKWRK, ICKWRK, II, KK, MAXSP, MAXTB, RU,
     1                   PATM, T, C, NSPEC, NU, NUNK, NPAR, PAR, NREV,
     2                   IREV, RPAR, NFAL, IFAL, IFOP, IFLO, KFAL, NFAR,
     3                   FPAR, NLAN, NLAR, ILAN, PLT, NRLT, IRLT, RPLT,
     4                   NTHB, ITHB, NTBS, AIK, NKTB, SMH, RKFT, RKRT,
     5                   RKF, RKR, EQK, CTB, NRNU, IRNU, RNU, NORD,
     6                   IORD, MXORD, KORD, RORD, NEIM, IEIM, IEIMT,
     7                   NJAN, NJAR, IJAN, PJAN, NFT1, NF1R, IFT1,
     8                   PF1, NMOM, IMOM, NXSM, IXSM, NTDE, ITDE, ITDK,
     9                   KTFL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRAT  (RCKWRK, ICKWRK, II, KK, MAXSP, MAXTB, RU,
C                     PATM, T, C, NSPEC, NU, NUNK, NPAR, PAR, NREV,
C                     IREV, RPAR, NFAL, IFAL, IFOP, IFLO, KFAL, NFAR,
C                     FPAR, NLAN, NLAR, ILAN, PLT, NRLT, IRLT, RPLT,
C                     NTHB, ITHB, NTBS, AIK, NKTB, SMH, RKFT, RKRT,
C                     RKF, RKR, EQK, CTB, NRNU, IRNU, RNU, NORD,
C                     IORD, MXORD, KORD, RORD, NEIM, IEIM, IEIMT,
C                     NJAN, NJAR, IJAN, PJAN, NFT1, NF1R, IFT1,
C                     PF1, NMOM, IMOM, NXSM, IXSM, NTDE, ITDE, ITDK,
C                     KTFL)**
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  II        - Integer scalar, total reaction count.
C  KK        - Integer scalar, total species count.
C  MAXSP     - Integer scalar, maximum number of species allowed in a
C              reaction; in the current formulation MAXSP=12.
C  MAXTB     - Integer scalar, maximum number of third bodies allowed
C              in a reaction.
C  RU        - Real scalar, universal gas constant.
C                 cgs units, 8.314510E7 ergs/(mole*K)
C  PATM      - Real scalar, pressure of one standard atmosphere.
C                 cgs units, 1.01325E6 dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                   cgs units, mole/cm**3
C  NSPEC(*)  - Integer array, count of participant species in reactions,
C              and the flag for reversibility of reactions;
C              dimension at least II, the total reaction count.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  NPAR      - Integer scalar, total number of parameters in the
C              Arrhenius rate expression for reactions;
C              in the current formulation NPAR=4.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  NREV      - Integer scalar, total number of reactions with
C              explicit reverse parameters.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C              IREV(N) is the reaction index for the Nth reaction
C              with explicit reverse parameters.
C  RPAR(*,*) - Real matrix,  reverse Arrhenius rate coefficients for
C              the NREV reactions; dimension at least NPAR for the
C              first, and at least NREV for the second.
C              RPAR(L,N) is the Lth coefficient for the Nth reaction
C              with explicit reverse rate coefficients.
C                 cgs units same as PAR(*,*)
C  NFAL      - Integer scalar, total number of pressure-dependent
C              reactions.
C  IFAL(*)   - Integer array,  reaction indices for the NFAL reactions;
C              dimension at least NFAL.
C              IFAL(N) = I, reaction I is a pressure-dependent reaction.
C  IFOP(*)   - Integer array, formulation type for the NFAL
C              reactions; dimension at least NFAL.
C  IFLO(*)   - Integer array, pressure-dependence option (unimolecular,
C              or chemically activated) for the NFAL reactions;
C              dimension at least NFAL.
C  KFAL(*)   - Integer array, Array of
C  NFAR      - Integer scalar, total number of additional parameters
C              allowed for reaction formulations.
C  FPAR(*,*) - Real matrix of additional parameters for the NFAL
C              reactions; dimension at least NFAR for the first, and
C              at least NFAL for the second.
C              FPAR(L,N) is the Lth parameter for the Nth pressure-
C              dependent reaction.
C  NLAN      - Integer scalar, total number of Landau-Teller reactions.
C  NLAR      - Integer scalar, number of additional parameters required
C              in a Landau-Teller rate expression;
C              in the current formulation NLAR=2.
C  ILAN(*)   - Integer array, reaction indices for the NLAN reactions;
C              dimension at least NLAN.
C              ILAN(N) is the reaction index for the Nth Landau-
C              Teller reaction.
C  PLT(*,*)  - Real matrix, the additional parameters for the NLAN
C              reactions; dimension at least NLAR for the first and at
C              least NLAN for the second.
C              PLAN(L,N) is the Lth parameter for the Nth
C              Landau-Teller reaction.
C  NRLT      - Integer scalar, total number of Landau-Teller reactions
C              with explicit reverse parameters.
C  IRLT(*)   - Integer array, reaction indices for the NRLT reactions;
C              dimension at least NRLT.
C              IRLT(N) is the reaction index for the Nth reaction
C              with Landau-Teller reverse parameters.
C  RPLT(*,*) - Real matrix, the additional rate parameters for the
C              NRLT reactions; dimension at least NLAR for the first
C              and at least NRLT for the second.
C              RPLT(L,N) is the Lth reverse parameter for the Nth
C              Landau-Teller reaction with reverse parameters.
C  NTHB      - Integer scalar, total number of third-body reactions.
C  ITHB(*)   - Integer array, reaction indices for the NTHB reactions;
C              dimension at least NTHB.
C              ITHB(N) is the reaction index for the Nth third-body
C              reaction.
C  NTBS(*)   - Integer array, total enhanced third-body count for a
C              third-body reaction; dimension at least NTHB.
C              NTBS(N) is the number of enhanced third bodies for
C              the Nth third-body reaction.
C  AIK(*,*)  - Real matrix, enhancement factors of third bodies
C              for the NTHB reactions; dimension MAXTB for the first,
C              the maximum number of enhancement factors, and NTHB
C              for the second.
C              AIK(L,N) is the enhancement factor for the Lth
C              enhanced third body in the Nth third-body reqction.
C  NKTB(*,*) - Integer matrix, species indices for the enhanced
C              third bodies in the NTHB reactions; dimension MAXTB
C              for the first and NTHB for the second.
C
C  NRNU      - Integer scalar, total number of reactions with real
C              stoichiometric coefficients.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  CTB(*)   -  Real array, concentration of third bodies for the
C              reactions; dimension at least II, the total reaction
C              count.
C              CTB(I) is the third-body concentration for reaction I.
C  NORD      - Integer scalar, total number of species changed-order
C              reaction.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  MXORD     - Integer scalar, maximum number of species change-orders
C              allowed in a reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  NEIM      - Integer scalar, total number of electron-impact
C              reactions.
C  IEIM(*)   - Integer array, reaction indices for the NEIM reactions;
C              dimension at least NEIM.
C              IEIM(N) is the reaction index for the Nth electron-
C              impact reaction.
C  IEIMT(*)  - Integer array, temperature-dependence indices for the
C              NEIM reactions; dimension at least NEIM.
C              IEIMT(N) is a pointer into a temperature array for
C              the Nth electron-impact reaction.
C  NJAN      - Integer scalar, total number of Janev, Langer, et al.
C              reactions.
C  NJAR      - Integer scalar, the number of additional rate
C              parameters required for Janev reactions.
C  IJAN(*)   - Integer array, reaction indices of the NJAN reactions;
C              dimension at least NJAN.
C              IJAN(N) is the reaction index for the Nth Janev et al.
C              reaction.
C  PJAN(*,*) - Real matrix, rate parameters for the NJAN reactions;
C              dimension at least NJAR for the first and NJAN for the
C              second.
C              PJAN(L,N) is the Lth parameter for the Nth Janev et al.
C              reaction.
C  NFT1      - Integer scalar, total number of fit-type reactions.
C  IFT1(*)   - Integer array, reaction indices for the NFT1 reactions;
C              dimension at least NFT1.
C              IFT1(N) is the reaction index for the Nth fit-type
C              reaction.
C  NF1R      - Integer scalar, number of additional rate parameters
C              required in a fit-type reactions.
C  PF1(*,*)  - Real matrix, the additional rate parameters for the
C              NFT1 reactions; dimension at least NF1R for the first
C              and at least NFT1 for the second.
C              PF1(L,N) is the Lth fit parameter for the Nth fit-type
C              reaction.
C  NMOM      - Integer scalar, total number of electron momentum-
C              transfer reactions.
C  IMOM(*)   - Integer array, reaction indices for the NMOM reactions;
C              dimension at least NMOM.
C              IMOM(N) is the reaction index for the Nth electron
C              momentum-transfer reaction.
C  NXSM      - Integer scalar, total number of ion momemtum-
C              transfer cross-section reactions.
C  IXSM(*)   - Integer array, reaction indices for the NXSM reactions;
C              dimension at least NXSM.
C  NTDE      - Integer scalar, total number of non-thermal-
C              equilibrium reactions.
C  ITDE(*)   - Integer array, reaction indices for the NTDE reactions;
C              dimension at least NTDE.
C              ITDE(N) is the reaction index for the Nth non-
C              thermal-equilibrium reaction.
C  KTFL(*)   - Integer array, indices into the temperature array,
C              for the KK species, as required by the NTDE reactions;
C              dimension at least KK.
C              KTFL(K) is the temperature array index for species K.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT; dimension at least KK.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units ,depend on the reaction
C  EQKC(*)   - Real array, equilibrium constants in concentration
C              units for reactions; dimension at least II.
C              EQKC(I) is the equilibrium constant for reaction I.
C                 cgs units, (mole/cm**3)**some power,
C                 depends on reaction
C  RKF(*)    - Real array
C              dimension at least II.
C  RKR(*)    - Real array
C              dimension at least II.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RCKWRK(*), ICKWRK(*), C(*), NSPEC(*), NU(MAXSP,*),
     1          NUNK(MAXSP,*), PAR(NPAR,*), IREV(*), RPAR(NPAR,*),
     2          ILAN(*), IRLT(*), PLT(NLAR,*), RPLT(NLAR,*),
     3          IFAL(*), IFOP(*), IFLO(*), KFAL(*), FPAR(NFAR,*),
     4          ITHB(*), NTBS(*), AIK(MAXTB,*), NKTB(MAXTB,*),
     5          SMH(*), RKFT(*), RKRT(*), RKF(*), RKR(*), EQK(*),
     6          CTB(*), IRNU(*), RNU(MAXSP,*), IORD(*), KORD(MXORD,*),
     7          RORD(MXORD,*), T(*), IEIM(*), IEIMT(*), IJAN(*),
     8          PJAN(NJAR,*), IFT1(*), PF1(NF1R,*), IMOM(*), IXSM(*),
     9          ITDE(*), ITDK(*), KTFL(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
     1             RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU, IEIM,
     2             IEIMT, IJAN, PJAN, IFT1, PF1, ITDE, ITDK, KTFL,
     3             RKFT, RKRT, EQK)
      CALL CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL, FPAR,
     1             ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF, RKR, CTB,
     2             IRNU, RNU, IORD, KORD, RORD, IMOM, IXSM)
C
C     end of SUBROUTINE CKRAT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
     1                   RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU,
     2                   IEIM, IEIMT, IJAN, PJAN, IFT1, PF1, ITDE,
     3                   ITDK, KTFL, RKFT, RKRT, EQK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATT (RCKWRK, ICKWRK, T, NSPEC, NU, NUNK, PAR, IREV,
C                     RPAR, ILAN, PLT, IRLT, RPLT, SMH, IRNU, RNU,
C                     IEIM, IEIMT, IJAN, PJAN, IFT1, PF1, ITDE,
C                     ITDK, KTFL, RKFT, RKRT, EQK)
C
C  INPUT
C
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  NSPEC(*)  - Integer array, total number of participant species
C              for reactions, and the flag for reversibility of
C              reactions; dimension at least II.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  IREV(*)   - Integer array, reaction indices for the NREV reactions;
C              dimension at least NREV.
C              IREV(N) is the reaction index for the Nth reaction
C              with explicit reverse parameters.
C  RPAR(*,*) - Real matrix,  reverse Arrhenius rate coefficients for
C              the NREV reactions; dimension at least NPAR for the
C              first, and at least NREV for the second.
C              RPAR(L,N) is the Lth coefficient for the Nth reaction
C              reaction with explicit reverse rate coefficients.
C                 cgs units same as PAR(*,*)
C  ILAN(*)   - Integer array, reaction indices for the NLAN reactions;
C              dimension at least NLAN.
C              ILAN(N) is the reaction index for the Nth Landau-
C              Teller reaction.
C  PLT(*,*)  - Real matrix, the additional parameters for the NLAN
C              reactions; dimension at least NLAR for the first and at
C              least NLAN for the second.
C              PLAN(L,N) is the Lth parameter for the Nth
C              Landau-Teller reaction.
C  IRLT(*)   - Integer array, reaction indices for the NRLT reactions;
C              dimension at least NRLT.
C              IRLT(N) is the reaction index for the Nth reaction
C              with Landau-Teller reverse parameters.
C  RPLT(*,*) - Real matrix, the additional rate parameters for the
C              NRLT reactions; dimension at least NLAR for the first
C              and at least NRLT for the second.
C              RPLT(L,N) is the Lth reverse parameter for the Nth
C              Landau-Teller reaction with reverse parameters.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  IEIM(*)   - Integer array, reaction indices for the NEIM reactions;
C              dimension at least NEIM.
C              IEIM(N) is the reaction index for the Nth electron-
C              impact reaction.
C  IEIMT(*)  - Integer array, temperature-dependence indices for the
C              NEIM reactions; dimension at least NEIM.
C              IEIMT(N) is a pointer into a temperature array for
C              the Nth electron-impact reaction.
C  IJAN(*)   - Integer array, reaction indices of the NJAN reactions;
C              dimension at least NJAN.
C              IJAN(N) is the reaction index for the Nth Janev et al.
C              reaction.
C  PJAN(*,*) - Real matrix, rate parameters for the NJAN reactions;
C              dimension at least NJAR for the first and NJAN for the
C              second.
C              PJAN(L,N) is the Lth parameter for the Nth Janev et al.
C              reaction.
C  IFT1(*)   - Integer array, reaction indices for the NFT1 reactions;
C              dimension at least NFT1.
C              IFT1(N) is the reaction index for the Nth fit-type
C              reaction.
C  PF1(*,*)  - Real matrix, the additional rate parameters for the
C              NFT1 reactions; dimension at least NF1R for the first
C              and at least NFT1 for the second.
C              PF1(L,N) is the Lth fit parameter for the Nth fit-type
C              reaction.
C  ITDE(*)   - Integer array, reaction indices for the NTDE reactions;
C              dimension at least NTDE.
C              ITDE(N) is the reaction index for the Nth non-
C              thermal-equilibrium reaction.
C  ITDK(*)   - Integer array, special species for the NTDE reactions;
C              dimension at least NTDE.
C              ITDK(N) is the indentifying species, K, whose associated
C              temperature (KTFL(K)) should be used in calculating
C              quantities associated with the reaction.
C  KTFL(*)   - Integer array, indices into the temperature array,
C              for the KK species, as required by the NTDE reactions;
C              dimension at least KK.
C              KTFL(K) is the temperature array index for species K.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for the species,
C              SMH(K) = S(K)/R - H(K)/RT; dimension at least KK.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units depend on the reaction
C  EQK(*)    - Real array, equilibrium constants in concentration
C              units for reactions; dimension at least II.
C              EQKC(I) is the equilibrium constant for reaction I.
C                 cgs units (mole/cm**3)**some power,
C                 depends on reaction
C               NOTE: EQK(I) as returned from this routine may be
C                     in error. Use CKEQ instead.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      include 'ckstrt.h'
C     Integer arrays
      DIMENSION ICKWRK(*), NSPEC(NII), NU(MXSP,NII), NUNK(MXSP,NII),
     1          IREV(NREV), ILAN(NLAN), IRLT(NRLT), IRNU(NRNU),
     2          IEIM(NEIM), IEIMT(NEIM), IJAN(NJAN), IFT1(NFT1),
     3          ITDE(NTDE), ITDK(NTDE), KTFL(*),
C     Real arrays
     4          RCKWRK(*), PAR(NPAR+1,NII), RPAR(NPAR+1,NREV),
     5          PLT(NLAR,NII), RPLT(NLAR,NRLT), SMH(NKK), RKFT(NII),
     6          RKRT(NII), RNU(MXSP,NRNU), T(*), EQK(NII),
     7          PJAN(NJAR,NJAN), PF1(NF1R,NFT1)
C
      COMMON /MACH/ SMALL, BIG, EXPARG
      INTEGER CKLKUP
      EXTERNAL CKLKUP, CKSMH
C
C     Find Gibbs/ RT for all species in the mechanism
C     Note: CKSMH takes a vector of temperatures, and the G/RT
C           value for each species will be evaluated at the
C           temperature associated with that species.
C
      CALL CKSMH (T, ICKWRK, RCKWRK, SMH)
      TEMP = T(1)
      ALOGT = LOG(TEMP)
      TINV = 1.0/TEMP
      RU = RCKWRK(NcRU)
      PATM = RCKWRK(NcPA)
      PFAC = PATM / (RU * TEMP)
C
C.....Default calculation of the forward rate constant..................
C     Put it in a loop that vectorizes
C
      DO 30 I = 1, NII
        RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*ALOGT - PAR(3,I)*TINV)
   30 CONTINUE
C
C.....Default way to calculate the equilibrium constant.................
C     Loop over NII -> largest expense in the routine, so optimize)
C
      DO 50 I = 1, NII
C
C     First decide on whether the equilibrium constant needs to be
C     calculated. -> irreversible rxns don't need their equilibrium
C     constants calculated. REV rxns don't either, but a lookup
C     function on a generic loop over all gas reactions would be too
C     expensive.
C
       IF (NSPEC(I) .GT. 0) THEN
C
C         Initialize the net mole change number and the
C         the net DeltaG value with the contributions from
C         the first reactant and product
C
         NUSUMK = NU(1,I) + NU(7,I)
         SUMSMH = NU(1,I)*SMH(NUNK(1,I)) + NU(7,I)*SMH(NUNK(7,I))
         IF (NUNK(2,I) .NE. 0) THEN
           NUSUMK = NUSUMK + NU(2,I)
           SUMSMH = SUMSMH + NU(2,I)*SMH(NUNK(2,I))
           IF (NUNK(3,I) .NE. 0) THEN
             NUSUMK = NUSUMK + NU(3,I)
             SUMSMH = SUMSMH + NU(3,I)*SMH(NUNK(3,I))
             IF (NUNK(4,I) .NE. 0) THEN
               NUSUMK = NUSUMK + NU(4,I)
               SUMSMH = SUMSMH + NU(4,I)*SMH(NUNK(4,I))
               IF (NUNK(5,I) .NE. 0) THEN
                 NUSUMK = NUSUMK + NU(5,I)
                 SUMSMH = SUMSMH + NU(5,I)*SMH(NUNK(5,I))
                 IF (NUNK(6,I) .NE. 0) THEN
                   NUSUMK = NUSUMK + NU(6,I)
                   SUMSMH = SUMSMH + NU(6,I)*SMH(NUNK(6,I))
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
           NUSUMK = NUSUMK + NU(8,I)
           SUMSMH = SUMSMH + NU(8,I)*SMH(NUNK(8,I))
           IF (NUNK(9,I) .NE. 0) THEN
             NUSUMK = NUSUMK + NU(9,I)
             SUMSMH = SUMSMH + NU(9,I)*SMH(NUNK(9,I))
             IF (NUNK(10,I) .NE. 0) THEN
               NUSUMK = NUSUMK + NU(10,I)
               SUMSMH = SUMSMH + NU(10,I)*SMH(NUNK(10,I))
               IF (NUNK(11,I) .NE. 0) THEN
                 NUSUMK = NUSUMK + NU(11,I)
                 SUMSMH = SUMSMH + NU(11,I)*SMH(NUNK(11,I))
                 IF (NUNK(12,I) .NE. 0) THEN
                   NUSUMK = NUSUMK + NU(12,I)
                   SUMSMH = SUMSMH + NU(12,I)*SMH(NUNK(12,I))
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
C
C        Calculate the concentration equilibrium constant,
C        Protecting against overflow in the exponential
C
         IF (NUSUMK .NE. 0) THEN
           EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**NUSUMK)
         ELSE
           EQK(I) = EXP(MIN(SUMSMH,EXPARG))
         ENDIF
C
C        Calculate the reverse rate constant from the forward
C        rate constant and the equilibrium constant
C
          RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
        ELSE
C
C        Case or irreversible reactions
C
          RKRT(I) = 0.0
        ENDIF
   50 CONTINUE
C
C.....Fix-up look for rxn's with real-stoichiometry.....................
C
C     Must completely redo the calculation of the Gibbs free energy
C     of reaction because the stoichiometric coefficients have been
C     overridden.
C     There is also no need to do this loop if this real stoich
C     reaction is also a REV reaction.
C
      DO 70 N = 1, NRNU
         I = IRNU(N)
         ISREV = CKLKUP(I, IREV, NREV)
         IF (NSPEC(I).GT.0 .AND. ISREV.EQ.0) THEN
           RNUSUM = RNU(1,N) + RNU(7,N)
           SUMSMH = RNU(1,N)*SMH(NUNK(1,I)) + RNU(7,N)*SMH(NUNK(7,I))
           DO 60 L = 2, 6
              IF (NUNK(L,I) .EQ. 0) GO TO 61
                 SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
                 RNUSUM = RNUSUM + RNU(L,N)
   60      CONTINUE
   61      CONTINUE
           DO 62 L = 8, 12
             IF (NUNK(L,I) .EQ. 0) GO TO 63
                SUMSMH = SUMSMH + RNU(L,N)*SMH(NUNK(L,I))
                RNUSUM = RNUSUM + RNU(L,N)
   62      CONTINUE
   63      CONTINUE
           IF (RNUSUM .NE. 0.0) THEN
              EQK(I) = EXP(MIN(SUMSMH,EXPARG)) * (PFAC**RNUSUM)
           ELSE
              EQK(I) = EXP(MIN(SUMSMH,EXPARG))
           ENDIF
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         END IF
   70 CONTINUE
C
C.....Rxns with explicitly defined reverse rate constants...............
C
C     Fix-up loop for reactions which have a defined reverse
C     rate constant. We will recalculate the reverse rate
C     constant in this case
C
      DO 90 N = 1, NREV
         I = IREV(N)
         RKRT(I) = RPAR(1,N) * EXP(RPAR(2,N)*ALOGT - RPAR(3,N)*TINV)
   90 CONTINUE
C
C.....Landau-Teller reactions...........................................
C
      DO 200 N = 1, NLAN
         I = ILAN(N)
C
C        modify K_f for reaction I
         TFAC = PLT(1,N)/TEMP**(1.0/3.0) + PLT(2,N)/TEMP**(2.0/3.0)
         RKFT(I) = RKFT(I) * EXP(TFAC)
C
C        Bail out for irreversible reactions
         IF (NSPEC(I) .LE. 0) GO TO 200
C
C        Lookup whether this is an Landau-Teller Rxn with an
C        explicitly given reverse reaction rate
C        - Note: If this is such a rxn, need to modify calculation
C                to add more parameters.
C                Even if it isn't such a reaction, we just changed
C                RKFT(I) above, so need to recalculate RKRT(I)
C
         ISRLT = CKLKUP(I, IRLT, NRLT)
         IF (ISRLT .EQ. 0) THEN
C           new K_r = K_eq / new K_f if no explicit reverse LT
C           parameters
            RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ELSE
C           new K_r = old K_r, modified by explicit parameters
            TFAC = RPLT(1,ISRLT)/TEMP**(1.0/3.0)
     1           + RPLT(2,ISRLT)/TEMP**(2.0/3.0)
            RKRT(I) = RKRT(I) * EXP(TFAC)
         ENDIF
  200 CONTINUE
C
C.....Electron-impact reactions.........................................
C
      DO 300 N = 1, NEIM
C
C        No change in rate constant unless temperature other
C        than first temperature is associated with this reaction.
C
         IF (IEIMT(N) .EQ. 1) GO TO 300
C
C        new K_f for reaction I is Arrhenius expression with new
C        temperature. Still check to see if TEMP is different than
C        T(1)
C
         TEMP = T(IEIMT(N))
         IF (TEMP .EQ. T(1)) GO TO 300
         I = IEIM(N)
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I) * LOG(TEMP) - PAR(3,I)/TEMP)
C
C          If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 300
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C             new K_r with explicit parameters and new temperature
            RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV)*LOG(TEMP)
     1              - RPAR(3,ISREV)/TEMP)
         ELSE
C
C            Fix concentration term in equilibrium constant, using the
C            previously calculated Gibbs Free energy part, which is
C            still good.
C
           PFAC2 = PATM / (RU*TEMP)
C          Does reaction have real stoichiometry?
           ISREAL = CKLKUP(I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
             L = ISREAL
             RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1              + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2              + RNU(11,L)+RNU(12,L)
             IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
             NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1              + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
             IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF
C
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
C
  300 CONTINUE
C
C.....Non-thermal-equilibrium, species-temperature-dependent reactions:.
C
      DO 400 N = 1, NTDE
C
C        Get the temperature index for the species which controls
C        which temperature to use in this reaction
C
         KTEMP = KTFL(ITDK(N))
C
C        No change in rate constant unless temperature other
C        than first temperature is associated with this reaction.
C
         IF (KTEMP .EQ. 1) GO TO 400
C
         TEMP = T(KTEMP)
         IF (TEMP .EQ. T(1)) GOTO 400
         I = ITDE(N)
C
C        Calculate a new K_f for reaction I, Arrhenius expression
C        with the new temperature
C
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 400
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C          Fix concentration term in equilibrium constant, using the
C          previously calculated Gibbs Free energy part, which is
C          still good.
C
           PFAC2 = PATM / (RU*TEMP)
C            real coefficients?
           ISREAL = CKLKUP(I, IRNU, NRNU)
           IF (ISREAL .GT. 0) THEN
             L = ISREAL
             RNUSUM= RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1             + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2             + RNU(11,L) +RNU(12,L)
             IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
           ELSE
             NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1              + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)+NU(12,I)
             IF (NUSUMK.NE.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** NUSUMK
           ENDIF
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
  400 CONTINUE
C
C.....reactions using fit#1:  k = A * T^B * exp(v1/T+v2/T^2+v3/T^3...)..
C
      DO 500 N = 1, NFT1
         IF (N .EQ. 1) BLOG = LOG(BIG)
         I = IFT1(N)
C
C        Species temperature-dependent reaction?
C
         N2 = CKLKUP(I, ITDE, NTDE)
         IF (N2 .EQ. 0) THEN
           TEMP = T(1)
         ELSE
           TEMP = T(KTFL(ITDK(N2)))
         ENDIF
C
C        Electron-impact reaction?
C
         N2 = CKLKUP(I, IEIM, NEIM)
         IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C        Calculate a new K_f based on fit#1, and possibly a
C        different temperature
C
         SUMJ = 0.0
         DO 470 J = 1, NF1R
            ROOTJ = 1.0/J
            IF (TEMP.GE.BIG**ROOTJ) THEN
               TEMPPW = BIG
            ELSE
               TEMPPW = TEMP**J
            ENDIF
            IF (SUMJ .GE. BLOG) THEN
               SUMJ = BLOG
            ELSE
               SUMJ = SUMJ + PF1(J,N)/TEMPPW
            ENDIF
 470     CONTINUE
         SUMJ = MIN (SUMJ, BLOG)
         RKFT(I) = MIN(BIG, (PAR(1,I) * TEMP**PAR(2,I) * EXP(SUMJ)))
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 500
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C         Fix concentration term in equilibrium constant, using the
C         previously calculated Gibbs Free energy part, which is
C         still good. Only needs to be done, if we have a new
C         temperature
C
           IF (TEMP .NE. T(1)) THEN
C
              PFAC2 = PATM / (RU * TEMP)
C             real coefficients?
              ISREAL = CKLKUP(I, IRNU, NRNU)
              IF (ISREAL .GT. 0) THEN
                 L = ISREAL
                 RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1                  + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2                  + RNU(11,L)+RNU(12,L)
                 IF (RNUSUM .NE. 0.0) EQK(I)=EQK(I)*(PFAC2/PFAC)**RNUSUM
              ELSE
                 NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)
     1                  + NU(6,I)+NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)
     2                  + NU(11,I)+NU(12,I)
                 IF (NUSUMK .NE. 0) EQK(I)=EQK(I)* (PFAC2/PFAC)**NUSUMK
              ENDIF
           ENDIF
C
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
  500 CONTINUE
C
C.....jannev, langer, evans & post - type reactions.....................
C
      DO 600 N = 1, NJAN
         I = IJAN(N)
C
C        Specify the temperature to be used below by looking up
C        whether this rxn has a special temperature variable
C        associated with it. If it doesn't, use the default first
C        temperature variable.
C
         N2 = CKLKUP(I, ITDE, NTDE)
         IF (N2 .EQ. 0) THEN
           TEMP = T(1)
         ELSE
           TEMP = T(KTFL(ITDK(N2)))
         ENDIF
C
C        Lookup up whether this rxn is an electron impact rxn
C        If it is, override the temperature variable determined
C        in the previous loop.
C
         N2 = CKLKUP (I, IEIM, NEIM)
         IF (N2 .NE. 0) TEMP = T(IEIMT(N2))
C
C        Re-evaluate Arrhenius K_f, possibly different temperature,
C        then modify for jannev expression
C
         RKFT(I) = PAR(1,I) * EXP(PAR(2,I)*LOG(TEMP) - PAR(3,I)/TEMP)
C        convert E- temperature to eV's
         TEV = TEMP / 11595.0
         SUMJ = 0.0
         DO 530 J = 1, NJAR
            SUMJ = SUMJ + PJAN(J,N) * (LOG(TEV))**(J-1)
  530    CONTINUE
         RKFT(I) =  MIN(BIG, RKFT(I) * EXP(SUMJ))
C
C        If reaction is irreversible, we can stop here
C
         IF (NSPEC(I) .LT. 0) GO TO 600
C
C        Lookup whether the rxn has explicit parameters for the
C        reverse rate constant, and branch accordingly
C
         ISREV = CKLKUP(I, IREV, NREV)
         IF (ISREV .GT. 0) THEN
C          new K_r with explicit parameters and new temperature
           RKRT(I) = RPAR(1,ISREV) * EXP(RPAR(2,ISREV) *LOG(TEMP)
     1             - RPAR(3,ISREV)/TEMP)
         ELSE
C
C          Fix concentration term in equilibrium constant, using the
C          previously calculated Gibbs Free energy part, which is
C          still good. Only needs to be done, if we have a new
C          temperature
C
           IF (TEMP .NE. T(1)) THEN
C
             PFAC2 = PATM / (RU * TEMP)
C            real coefficients?
             ISREAL = CKLKUP(I, IRNU, NRNU)
             IF (ISREAL .GT. 0) THEN
               L = ISREAL
               RNUSUM = RNU(1,L)+RNU(2,L)+RNU(3,L)+RNU(4,L)+RNU(5,L)
     1                + RNU(6,L)+RNU(7,L)+RNU(8,L)+RNU(9,L)+RNU(10,L)
     2                + RNU(11,L)+RNU(12,L)
               IF (RNUSUM.NE.0.0) EQK(I)=EQK(I)* (PFAC2/PFAC) ** RNUSUM
             ELSE
               NUSUMK = NU(1,I)+NU(2,I)+NU(3,I)+NU(4,I)+NU(5,I)+NU(6,I)
     1                + NU(7,I)+NU(8,I)+NU(9,I)+NU(10,I)+NU(11,I)
     2                + NU(12,I)
               IF (NUSUMK.NE.0) EQK(I) = EQK(I) * (PFAC2/PFAC) ** NUSUMK
             ENDIF
           ENDIF
C          new K_r = new K_f / new K_eq
           RKRT(I) = RKFT(I) / MAX(EQK(I), SMALL)
         ENDIF
 600  CONTINUE
C
C     end of SUBROUTINE CKRATT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL,
     1                   FPAR, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,
     2                   RKR, CTB, IRNU, RNU, IORD, KORD, RORD, IMOM,
     3                   IXSM)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRATX (T, C, NU, NUNK, PAR, IFAL, IFOP, IFLO, KFAL,
C                     FPAR, ITHB, NTBS, AIK, NKTB, RKFT, RKRT, RKF,
C                     RKR, CTB, IRNU, RNU, IORD, KORD, RORD, IMOM,
C                     IXSM)
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  NSPEC(*)  - Integer array, total number of participant species
C              for reactions, and the flag for reversibility of
C              reactions; dimension at least II.
C              NSPEC(I) = +N, reaction I is reversible and has
C                          N participant species (reactants + products)
C                       = -N, reaction I is irreversible, etc.
C  NU(*,*)   - Integer matrix, stoichiometric coefficients for
C              species in reactions;
C              dimension at least MAXSP for the first and at least II
C              for the second.
C              NU(N,I) is the stoichiometric coefficient of the Nth
C              species in reaction I, and
C              NU < 0 if the Nth species is a reactant,
C              NU > 0 if the Nth species is a product.
C  NUNK(*,*) - Integer matrix, indices of species in reactions;
C              dimension at least MAXSP for the first, and at least
C              II for the second.
C              NUNK(N,I) is the species index for the Nth species in
C              reaction I.
C  PAR(*,*)  - Real matrix, Arrhenius coefficients for reactions;
C              dimension at least NPAR for the first, and at least
C              II for the second.  For any reaction I,
C              PAR(1,I) is the pre-exponential constant
C                 cgs units, mole-cm-sec-K
C              PAR(2,I) is the temperature dependent exponent
C                 cgs units, none
C              PAR(3,I) is the activation energy
C                 cgs units, K
C              PAR(4,I) is used as a perturbation factor in
C              sensitivity analyses.
C  IFAL(*)   - Integer array, reaction indices for the NFAL reactions;
C              dimension at least NFAL.
C              IFAL(N) = I, reaction I is a pressure-dependent reaction.
C  IFOP(*)   - Integer array, formulation type for the NFAL reactions;
C              dimension at least NFAL.
C  IFLO(*)   - Integer array, pressure-dependence type for the NFAL
C              reactions (unimolecular vs chemically activated);
C              dimension at least NFAL.
C  KFAL(*)   - Integer array, array of 3rd-body species indices for the
C              NFAL reactions; 0 to use total concentration of mixture,
C              K for concentration of species K.
C  FPAR(*,*) - Real matrix of parameters for the NFAL reactions;
C              dimension at least NFAR for the first, and at least NFAL
C              for the second.
C              FPAR(L,N) is the Lth parameter for the Nth pressure-
C              dependent reaction.
C  ITHB(*)   - Integer array, reaction indices for the NTHB reactions;
C              dimension at least NTHB.
C              ITHB(N) is the reaction index for the Nth third-body
C              reaction.
C  NTBS(*)   - Integer array, total number of enhanced third-bodies
C              in a third-body reaction; dimension at least NTHB.
C              NTBS(N) is the total enhanced third-body count for
C              the Nth third-body reaction.
C  AIK(*,*)  - Real matrix, enhancement factors of third bodies
C              for the NTHB reactions; dimension MAXTB for the first,
C              the maximum number of enhancement factors, and NTHB
C              for the second.
C              AIK(L,N) is the enhancement factor for the Lth
C              enhanced third body in the Nth third-body reqction.
C  NKTB(*,*) - Integer matrix, species indices for the enhanced
C              third bodies in the NTHB reactions; dimension MAXTB
C              for the first and NTHB for the second.
C  RKFT(*)   - Real array, temperature-dependent portion of the
C              forward reaction rates for reactions; dimension at
C              least II.
C              RKFT(I) is the temperature-dependent portion of the
C              forward reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  RKRT(*)   - Real array, temperature-dependent portion of reverse
C              reaction rates for reactions; dimension at least II.
C              RKRT(I) is the temperature-dependent portion of the
C              reverse reaction rate for reaction I.
C                 cgs units, depend on the reaction
C  CTB(*)   -  Real array, concentration of third bodies for the
C              reactions;
C              dimension at least II, the total reaction count.
C              CTB(I) is the third-body concentration for reaction I.
C  IRNU(*)   - Integer array, reaction indices for the NRNU reactions;
C              dimension at least NRNU.
C              IRNU(N) is the reaction index for the Nth reaction
C              with real stoichiometric coefficients.
C  RNU(*,*)  - Real matrix, stoichiometric coefficients for the NRNU
C              reactions; dimension at least MAXSP for the first and
C              at least NRNU for the second.
C              RNU(L,N) is the Lth stoichiometric coefficient for
C              the Nth reaction with real stoichiometry.
C  IORD(*)   - Integer array, reaction indices for the NORD reactions;
C              dimension at least NORD.
C              IORD(N) is the index of the Nth change-order reaction.
C  KORD(*,*) - Integer matrix, species indices for the order changes in
C              the NORD reactions; dimension at least MXORD for the
C              first and at least NORD for the second.
C              KORD(L,N) is the species index for the Lth order change
C              in the Nth change-order reaction.
C              KORD < 0 indicates change in forward order;
C              KORD > 0 indicates change in reverse order.
C  RORD(*,*) - Real matrix, order values for the NORD reactions;
C              dimension at least MXORD for the first and at least NORD
C              for the second.
C              RORD(L,N) is the order for the Lth order change in the
C              Nth change-order reaction.
C  IMOM(*)   - Integer array, reaction indices for the NMOM reactions;
C              dimension at least NMOM.
C              IMOM(N) is the reaction index for the Nth electron
C              momentum-transfer reaction.
C  IXSM(*)   - Integer array, reaction indices for the NXSM reactions;
C              dimension at least NXSM.
C  OUTPUT
C  RKF(*)-
C              dimension at least II, the total reaction count.
C  RKR(*)-
C              dimension at least II, the total reaction count.
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ZERO=0.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ZERO=0.0)
C*****END precision > single
C
      include 'ckstrt.h'
C     Integer arrays
      DIMENSION NU(MXSP,NII), NUNK(MXSP,NII),
     1          IFAL(NFAL), IFOP(NFAL), IFLO(NFAL), KFAL(NFAL),
     2          ITHB(NTHB), NTBS(NTHB), NKTB(MXTB,NTHB),
     3          IRNU(NRNU), IORD(NORD), KORD(MXORD, NORD),
     4          IMOM(NMOM), IXSM(NXSM),
C     Real arrays
     5          T(*), C(NKK), PAR(NPAR+1,NII), FPAR(NFAR,NII),
     6          AIK(MXTB,NTHB), RKFT(NII), RKRT(NII), RKF(NII),
     7          RKR(NII), CTB(NII), RNU(MXSP,NRNU), RORD(MXORD,NORD)
      COMMON /MACH/ SMALL,BIG,EXPARG
C
C     third-body reactions
C
      IF (NTHB .GT. 0) THEN
         CTOT = 0.0
         DO 10 K = 1, NKK
            CTOT = CTOT + C(K)
   10    CONTINUE
         DO 80 N = 1, NTHB
            I = ITHB(N)
            CTB(I) = CTOT
            DO 70 L = 1, NTBS(N)
               CTB(I) = CTB(I) + (AIK(L,N)-1.0)*C(NKTB(L,N))
   70       CONTINUE
   80    CONTINUE
      ENDIF
C
C     If pressure correction:
C
      IF (NFAL .GT. 0) THEN
         TK = T(1)
         ALOGT = LOG(TK)
C
         DO 90 N = 1, NFAL
C
            I = IFAL(N)
            IF (KFAL(N) .GT. 0) THEN
C              third-body species named
               CONC = C(KFAL(N))
            ELSE
C              no species named
               CONC = CTB(I)
               CTB(I) = 1.0
            ENDIF
C
            IF (IFLO(N) .EQ. 0) THEN
C              unimolecular reaction (RKFT is K_inf)
               PR = (FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/TK))
     1              * CONC / RKFT(I)
               PCOR = PR / (1.0 + PR)
            ELSE
C              chemically activated reaction (RKFT is K_zero)
               PR = RKFT(I) * CONC /
     1              (FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)/TK))
               PCOR = 1.0 / (1.0 + PR)
            ENDIF
C
            IF (IFOP(N) .GT. 1) THEN
C
C              Fcentering (F^x)
C
               PRLOG = LOG10(MAX(PR,SMALL))
C
               IF (IFOP(N) .EQ. 2) THEN
C
C              8-PARAMETER SRI FORM
C
                  XP = 1.0/(1.0 + PRLOG**2)
                  FC = ((FPAR(4,N)*EXP(-FPAR(5,N)/TK)
     1                   + EXP(-TK/FPAR(6,N))) **XP)
     2                  * FPAR(7,N) * TK**FPAR(8,N)
C
               ELSE
C
C              6-PARAMETER TROE FORM
C
                  FCENT = (1.0-FPAR(4,N)) * EXP(-TK/FPAR(5,N))
     1                  +       FPAR(4,N) * EXP(-TK/FPAR(6,N))
C
C              7-PARAMETER TROE FORM
C
                  IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)/TK)
C
                  FCLOG = LOG10(MAX(FCENT,SMALL))
                  XN    = 0.75 - 1.27*FCLOG
                  CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
                  FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
                  FC = 10.0**FLOG
               ENDIF
               PCOR = FC * PCOR
            ENDIF
C
            RKFT(I) = RKFT(I) * PCOR
            RKRT(I) = RKRT(I) * PCOR
   90    CONTINUE
      ENDIF
C
C       Modify rate constants for third body collision reactions
C       now that we have taken care of the unimolecular reaction
C       special case
C
      DO 100 N = 1, NTHB
        I = ITHB(N)
        RKFT(I) = RKFT(I) * CTB(I)
        RKRT(I) = RKRT(I) * CTB(I)
  100 CONTINUE
C
C     Calculate the rate of each reaction from the rate
C     constant and product of the concentration raized
C     to their stoichiometric powers.
C      - NU(1,I) check is used, because real stoich. coefficient
C        rxns set NU(*,I)=0. They are fixed up in a later
C        loop so there is no need to initialize RKF and RKR.
C
      DO 150 I = 1, NII
         IF (NU(1,I) .EQ. 0) GO TO 150
C        4th parameter may be perturbation factor
         RKF(I) = RKFT(I)*C(NUNK(1,I))**IABS(NU(1,I)) * PAR(4,I)
         RKR(I) = RKRT(I)*C(NUNK(7,I))**NU(7,I) * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I)= RKF(I) * C(NUNK(2,I))**IABS(NU(2,I))
            IF (NUNK(3,I).NE.0) THEN
               RKF(I) = RKF(I) * C(NUNK(3,I))**IABS(NU(3,I))
               IF (NUNK(4,I).NE.0) THEN
                  RKF(I) = RKF(I) * C(NUNK(4,I))**IABS(NU(4,I))
                  IF (NUNK(5,I).NE.0) THEN
                     RKF(I) = RKF(I) * C(NUNK(5,I))**IABS(NU(5,I))
                     IF (NUNK(6,I) .NE. 0)
     1                  RKF(I) = RKF(I) * C(NUNK(6,I))**IABS(NU(6,I))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I)*C(NUNK(8,I))**NU(8,I)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I)*C(NUNK(9,I))**NU(9,I)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I)*C(NUNK(10,I))**NU(10,I)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I)*C(NUNK(11,I))**NU(11,I)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I)*C(NUNK(12,I))**NU(12,I)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  150 CONTINUE
C
      DO 160 N = 1, NRNU
C        For real coefficients, must ensure non-negative concentrations
         I = IRNU(N)
         RKF(I) = RKFT(I) * MAX(ZERO,C(NUNK(1,I))) ** ABS(RNU(1,N)) 
     1                    * PAR(4,I)
         RKR(I) = RKRT(I) * MAX(ZERO,C(NUNK(7,I))) ** RNU(7,N) 
     1            * PAR(4,I)
         IF (NUNK(2,I) .NE. 0) THEN
            RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(2,I))) ** ABS(RNU(2,N))
            IF (NUNK(3,I) .NE. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(3,I))) ** ABS(RNU(3,N))
               IF (NUNK(4,I) .NE. 0) THEN
                  RKF(I) = RKF(I) * 
     1                     MAX(ZERO,C(NUNK(4,I))) ** ABS(RNU(4,N))
                  IF (NUNK(5,I) .NE. 0) THEN
                     RKF(I) = RKF(I) * 
     1                        MAX(ZERO,C(NUNK(5,I))) ** ABS(RNU(5,N))
                     IF (NUNK(6,I) .NE. 0) THEN
                       RKF(I) = RKF(I) * MAX(ZERO,C(NUNK(6,I))) 
     1                          ** ABS(RNU(6,N))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (NUNK(8,I) .NE. 0) THEN
            RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(8,I))) ** RNU(8,N)
            IF (NUNK(9,I) .NE. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(9,I))) ** RNU(9,N)
               IF (NUNK(10,I) .NE. 0) THEN
                  RKR(I) = RKR(I) * MAX(ZERO,C(NUNK(10,I))) ** RNU(10,N)
                  IF (NUNK(11,I) .NE. 0) THEN
                     RKR(I) = RKR(I) * 
     1                        MAX(ZERO,C(NUNK(11,I))) ** RNU(11,N)
                     IF (NUNK(12,I) .NE. 0) THEN
                        RKR(I) = RKR(I) * 
     1                           MAX(ZERO,C(NUNK(12,I))) ** RNU(12,N)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  160 CONTINUE
C
C       Handle reactions with different rxn orders than their
C       stoichiometry
C
      DO 200 N = 1, NORD
         I = IORD(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = RKRT(I) * PAR(4,I)
         DO 190 L = 1, MXORD
            NK = KORD(L,N)
            IF (NK .LT. 0) THEN
               RKF(I) = RKF(I) * MAX(ZERO,C(-NK)) ** RORD(L,N)
            ELSEIF (NK .GT. 0) THEN
               RKR(I) = RKR(I) * MAX(ZERO,C(NK)) ** RORD(L,N)
            ENDIF
  190    CONTINUE
  200 CONTINUE
C
C RETURN PURE COLLISION FREQUENCIES FOR E MOMENTUM-TRANSFER (NOT RATES)
C
      DO 300 N = 1, NMOM
         I = IMOM(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = 0.0
 300  CONTINUE
C
C RETURN PURE COLLISION X-SECTIONS FOR ION MOMENTUM-TRANSFER (NOT RATES)
C
      DO 350 N = 1, NXSM
         I = IXSM(N)
         RKF(I) = RKFT(I) * PAR(4,I)
         RKR(I) = 0.0
 350  CONTINUE
C
C     end of SUBROUTINE CKRATX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRCXP  (P, T, X, ICKWRK, RCKWRK, RCFT, RCRT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRCXP  (P, T, X, ICKWRK, RCKWRK, RCFT, RCRT)
C     Returns the forward and reverse rate constants for all reactions
C     given pressure, temperature and mole fractions;
C     see Eqs. (51) and (58).  Note this subroutine will calculate
C     a value for the reverse rate constant irrespective of
C     whether the reaction was deemed reversible in the interpretor
C     file.  Also note that the concentration of third bodies
C     for third body reactions is included in the returned rate
C     constant.  The units for the rate constant will depend
C     on the number of reactants.
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T(*)   - Temperature.
C                   cgs units - K
C                   Data type - real array
C     X      - Mole fractions of the species.
C                   cgs units - none
C                   Data type - real array
C                   Dimension X(*) at least KK, the total number of
C                   species.
C     ICKWRK - Array of integer workspace.
C                   Data type - integer array
C                   Dimension ICKWRK(*) at least LENIWK.
C     RCKWRK - Array of real work space.
C                   Data type - real array
C                   Dimension RCKWRK(*) at least LENRWK.
C
C  OUTPUT
C     RCFT      - Rate constant for the forward reaction.
C                   cgs units - mole/(cm**3*sec) *
C                           (cm**3/mole)**(sum of forward stoich. coeff)
C                   Data type - real array
C                   Dimension RCFT(*) at least II, the total number
C                   of reactions.
C     RCRT      - Rate constant for the forward reaction.
C                   cgs units - mole/(cm**3*sec) *
C                           (cm**3/mole)**(sum of reverse stoich. coeff)
C                   Data type - real array
C                   Dimension RCRT(*) at least II, the total number
C                   of reactions.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), X(NKK), RCFT(NII), RCRT(NII), T(*)
C      CHARACTER*64 VERSN
C      DATA VERSN/'@(#) ckrcxp.f 1.5 12/22/95'/
C
      EXTERNAL CKRATT, CKRTCN
C
C  Replace the normal NSPEC vector, located at ICKWRK(IcNS),
C  with an all positive one. Set the NREAC vector to be negative
C  to flag this. The NREAC vector is never used!
C    This way CKRATT will calculate the reverse rate constant for
C  reactions that are defined as irreversible in the mechanism
C
      NIM1 = NII - 1
      DO 10 I = 0, NIM1
        IF (ICKWRK(IcNS + I) .LT. 0) THEN
           ICKWRK(IcNS + I) = - ICKWRK(IcNS + I)
           ICKWRK(IcNR + I) = - ICKWRK(IcNR + I)
        ENDIF
 10   CONTINUE
C
C  Call the normal temperature dependent part of the rate constant
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     1             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     2             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     3             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     4             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     5             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     6             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     7             ICKWRK(IcTK), ICKWRK(IcKTF), RCFT, RCRT,
     8             RCKWRK(NcI1))
C
C  Flip the normal NSPEC vector, located at ICKWRK(IcNS),
C  back to its previous state
C
      DO 20 I = 0, NIM1
        IF (ICKWRK(IcNR + I) .LT. 0) THEN
           ICKWRK(IcNS + I) = - ICKWRK(IcNS + I)
           ICKWRK(IcNR + I) = - ICKWRK(IcNR + I)
        ENDIF
 20   CONTINUE
C
C  Calculate the concentrations of all species
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
C  Fix up the rate constants due to fall-off and third-body
C  effects
C
      CALL CKRTCN(RCKWRK, ICKWRK, NII, NKK, MXTB, RCKWRK(NcRU),
     1            RCKWRK(NcPA), T, RCKWRK(NcK1), NPAR+1, RCKWRK(NcCO),
     2            NFAL, ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcKF),
     3            NFAR, RCKWRK(NcFL), NTHB, ICKWRK(IcTB), ICKWRK(IcKN),
     4            RCKWRK(NcKT), ICKWRK(IcKT), RCFT, RCRT, RCKWRK(NcI4))
C
C     end of SUBROUTINE CKRCXP
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKRTCN (RCKWRK, ICKWRK, II, KK, MAXTB, RU, PATM, T, C,
     1                   NPAR, PAR, NFAL, IFAL, IFOP, KFAL, NFAR, FPAR,
     2                   NTHB, ITHB, NTBS, AIK, NKTB, RCFT, RCRT, CTB)
C
C  START PROLOGUE
C
C   This subroutine modifies the forward and reverse rate constants
C obtained from CKRATT to account for those parts of CKRATX that
C don't involve multiplications of the concentrations of reactants
C or products. This specifically includes third body effects and
C and fall-off effects.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      COMMON /MACH/ SMALL, BIG, EXPARG
      DIMENSION ICKWRK(*), RCKWRK(*), C(KK), PAR(NPAR,II),
     1          FPAR(NFAR,NFAL), AIK(MAXTB,NTHB), RCFT(II), RCRT(II),
     2          CTB(II), IFAL(NFAL), IFOP(NFAL), KFAL(NFAL),
     3          NTBS(NTHB), ITHB(NTHB), NKTB(MAXTB,NTHB)
C
C         PRELIMINARIES
C
      ALOGT = LOG(T)
      AINVT = 1.0/T
C
C  Third-body reactions -
C       Note that fall-off reactions with no special fall-off species
C       are alwys declared as third body reactions. That way CTB(I)
C       will be defined correctly in the fall-off reaction section
C       below.
C
      IF (NTHB .GT. 0) THEN
C
C         Find the total concentration
C
         CTOT = 0.0
         DO 30 K = 1, KK
           CTOT = CTOT + C(K)
   30    CONTINUE
         DO 80 N = 1, NTHB
           CTB(ITHB(N)) = CTOT
           DO 40 L = 1, NTBS(N)
             CTB(ITHB(N)) = CTB(ITHB(N)) + (AIK(L,N)-1.0)*C(NKTB(L,N))
   40      CONTINUE
   80   CONTINUE
      ENDIF
C
C  Corrections for fall-off reactions
C
      DO 90 N = 1, NFAL
C
C         Store the reaction number
C
        I = IFAL(N)
C
C         Store the special species, if there is one
C
        K = KFAL(N)
C
C         Calculate the low pressure reaction rate
C
        RKLOW = FPAR(1,N) * EXP(FPAR(2,N)*ALOGT - FPAR(3,N)*AINVT)
C
C         Find the concentration of the third body
C
        IF (K .EQ. 0) THEN
           PR = RKLOW * CTB(I) / RCFT(I)
        ELSE
           PR = RKLOW * C(K)   / RCFT(I)
        ENDIF
C
C         Unitize the effective concentration for fall-off
C         reactions. We don't want to multiply the rxn rate
C         constant by an extra concentration factor
C
        CTB(I) = 1.0
C
C         This is the Lindemann form , i.e., IFOP(N) = 1
C
        PCOR = PR / (1.0 + PR)
C
        IF (IFOP(N) .GT. 1) THEN
            PRLOG = LOG10(MAX(PR,SMALL))
            IF (IFOP(N) .EQ. 2) THEN
C
C              SRI FORM
C
               XP = 1.0/(1.0 + PRLOG**2)
               FC = ( ( FPAR(4,N)*EXP(-FPAR(5,N)*AINVT)
     $            + EXP(-T/FPAR(6,N)) )**XP ) * FPAR(7,N) * T**FPAR(8,N)
C
            ELSE
C
C              6-PARAMETER TROE FORM
C
               FCENT = (1.0-FPAR(4,N)) * EXP(-T/FPAR(5,N))
     $               + FPAR(4,N) * EXP(-T/FPAR(6,N))
C
C              7-PARAMETER TROE FORM
C
               IF (IFOP(N) .EQ. 4) FCENT = FCENT + EXP(-FPAR(7,N)*AINVT)
C
               FCLOG = LOG10(MAX(FCENT,SMALL))
               XN    = 0.75 - 1.27*FCLOG
               CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
               FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
               FC = 10.0**FLOG
            ENDIF
            PCOR = PCOR * FC
         ENDIF
C
C           Correct both the forward and reverse rate constant
C           by the k/k_inf, PCOR
C
         RCFT(I) = RCFT(I) * PCOR
         RCRT(I) = RCRT(I) * PCOR
C
   90 CONTINUE
C
C Multiply the rate constant by the third body factor
C
      DO 140 N = 1, NTHB
         I = ITHB(N)
         RCFT(I) = RCFT(I)*CTB(I)
         RCRT(I) = RCRT(I)*CTB(I)
  140 CONTINUE
C
C Multiply the rate constant by the perturbation factor.
C
      DO 150 I = 1, II
         RCFT(I) = RCFT(I)*PAR(4,I)
         RCRT(I) = RCRT(I)*PAR(4,I)
  150 CONTINUE
C
C     end of SUBROUTINE CKRTCN
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE CKRDEX (I, RCKWRK, RD)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRDEX (I, RCKWRK, RD)*
C  Get/put the perturbation factor of the Ith reaction
C
C  INPUT
C  I         - Integer scalar, reaction index;
C              I > 0 gets RD(I) from RCKWRK
C              I < 0 puts RD(I) into RCKWRK
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  If I < 1:
C  RD        - Real scalar, perturbation factor for reaction I;
C              cgs units, mole-cm-sec-K.
C
C  OUTPUT
C  If I > 1:
C  RD        - Real scalar, perturbation factor for reaction I;
C              cgs units, mole-cm-sec-K.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (ONE = 1.0D0)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C      PARAMETER (ONE = 1.0)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*)
C
      NI = NcCO + (IABS(I)-1)*(NPAR+1) + NPAR
      IF (I .GT. 0) THEN
         RD = RCKWRK(NI)
      ELSE
C
C          Assign the perturbation factor
C
         RCKWRK(NI) = RD
      ENDIF
C
C     end of SUBROUTINE CKRDEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHEX (K, RCKWRK, A6)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHEX (K, RCKWRK, A6)
C
C  Returns an array of the sixth thermodynamic polynomial
C  coefficients for a species, or changes their value,
C  depending on the sign of K.
C
C  INPUT
C  K         - Integer scalar, species index;
C              K>0 gets A6(*) from RCKWRK,
C              K<0 puts A6(*) into RCKWRK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  A6(*)     - Real array, the 6th thermodynamic polynomial
C              coefficients for species K, over the number
C              of fit temperature ranges; dimension at least (MXTP-1),
C              where MXTP is the maximum number of temperatures used
C              to divide the thermodynamic fits.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION RCKWRK(*), A6(*)
C
      DO 100 L = 1, MXTP-1
         NA6 = NCAA + (L-1)*NCP2 + (IABS(K)-1)*NCP2T + NCP
         IF (K .GT. 0) THEN
            A6(L) = RCKWRK(NA6)
         ELSE
            RCKWRK(NA6) = A6(L)
         ENDIF
  100 CONTINUE
C
C     end of SUBROUTINE CKRHEX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOC (P, T, C, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOC (P, T, C, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and molar concentrations;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*)
C
      RHO  = 0.0
      DO 100 K = 1, NKK
         RHO = RHO + C(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKRHOC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOX (P, T, X, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and mole fractions;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      SUMT = 0.0
      SUMW = 0.0
      DO 100 K = 1, NKK
         SUMW = SUMW + X(K)*RCKWRK(NcWT + K - 1)
         SUMT = SUMT + X(K)*T(ICKWRK(IcKTF + K - 1))
  100 CONTINUE
C
      RHO = P * SUMW / (SUMT * RCKWRK(NcRU))
C
C     end of SUBROUTINE CKRHOX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRHOY (P, T, Y, ICKWRK, RCKWRK, RHO)
C  Returns the mass density of the gas mixture given pressure,
C  temperature(s) and mass fractions;  see Eq. (2).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      SUM = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUM = SUM + Y(K+1) * T(ICKWRK(IcKTF + K)) / RCKWRK(NcWT + K)
150   CONTINUE
      RHO = P / (SUM * RCKWRK(NcRU))
C
C     end of SUBROUTINE CKRHOY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
C
C  START PROLOGUE
C
C  SUBROUTINE CKRP   (ICKWRK, RCKWRK, RU, RUC, PA)
C  Returns universal gas constants and the pressure of one standard
C  atmosphere
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  RU        - Real scalar, universal gas constant.
C                 cgs units, 8.314510E7 ergs/(mole*K)
C  RUC       - Real scalar, universal gas constant used only in
C              conjuction with activation energy.
C                 preferred units, RU / 4.184 cal/(mole*K)
C  PA        - Real scalar, pressure of one standard atmosphere.
C                 cgs units, 1.01325E6 dynes/cm**2
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*)
C
      RU  = RCKWRK(NcRU)
      RUC = RCKWRK(NcRC)
      PA  = RCKWRK(NcPA)
C
C     end of SUBROUTINE CKRP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSAVE (LOUT, LSAVE, ICKWRK, RCKWRK, CCKWRK)
C  Writes to a binary file information about a Chemkin linkfile,
C  pointers for the Chemkin Library, and Chemkin work arrays.
C
C  INPUT
C  LOUT      - Integer scalar, formatted output file unit number.
C  LSAVE     - Integer scalar, binary output file unit number.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
      COMMON /CKCONS/ PREC, FILVER, PRVERS, KERR, LENI, LENR, LENC
C
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER*(*) CCKWRK(*)
      CHARACTER*16 FILVER, PRVERS, PREC
      LOGICAL KERR
C
      NPOINT = 85
      WRITE (LSAVE, ERR=999)
     *                FILVER, PREC, LENI, LENR, LENC,
C
C     include file for CHEMKIN-III cklib.f, dated: March 1, 1966
C
C     Integer constants
C
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD,
     4   MXORD, KEL, NKKI,
C
C     Integer pointers to character arrays in CCKWRK
C
     5   IcMM, IcKK,
C
C     Integer pointers to integer arrays in ICKWRK
C
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, IcLT, IcRL,
     7   IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, IcEI,
     8   IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, IcTD,
     9   IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
C
C     Integer pointers to real variables and arrays in RCKWRK
C
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,
     2   NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
C
C     END include file for cklib.f
C
C
      WRITE (LSAVE, ERR=999) (ICKWRK(L), L = 1, LENI)
      WRITE (LSAVE, ERR=999) (RCKWRK(L), L = 1, LENR)
      WRITE (LSAVE, ERR=999) (CCKWRK(L), L = 1, LENC)
      RETURN
C
  999 CONTINUE
      WRITE (LOUT,*)' Error writing Chemkin binary file information...'
      KERR = .TRUE.
C
C     end of SUBROUTINE CKSAVE
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSBML (P, T, X, ICKWRK, RCKWRK, SBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSBML (P, T, X, ICKWRK, RCKWRK, SBML)*
C  Returns the mean entropy of the mixture in molar units given
C  pressure, temperature(s) and mole fractions; see Eq. (42).
C
C  INPUT
C  P         - Real scalar, pressure.
C                cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SBML      - Real scalar, mean entropy.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      RLNP = RCKWRK(NcRU) * LOG(P / RCKWRK(NcPA))
      SBML = 0.0
      DO 100 K = 1, NKK
         SBML = SBML + X(K) * ( RCKWRK(NcK1 + K - 1) -
     1          RCKWRK(NcRU)*LOG(MAX(X(K),SMALL)) - RLNP )
  100 CONTINUE
C
C     end of SUBROUTINE CKSBML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSBMS (P, T, Y, ICKWRK, RCKWRK, SBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSBMS (P, T, Y, ICKWRK, RCKWRK, SBMS)*
C  Returns the mean entropy of the mixture in mass units given pressure,
C  temperature(s) and mass fractions; see Eq. (43).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SBMS      - Real scalar, mean entropy.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /MACH/ SMALL,BIG,EXPARG
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKSML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKYTX (Y, ICKWRK, RCKWRK, RCKWRK(NcK2))
      CALL CKMMWY(Y, ICKWRK, RCKWRK, WTM)
      RLNP = RCKWRK(NcRU) * LOG (P / RCKWRK(NcPA))
C
      SUM = 0.0
      RU = RCKWRK(NcRU)
      NKM1 = NKK - 1
      DO 100 K = 0, NKM1
         SUM = SUM + RCKWRK(NcK2 + K) *
     1             ( RCKWRK(NcK1 + K) - RU *
     3               LOG(MAX(RCKWRK(NcK2 + K),SMALL)) - RLNP)
  100 CONTINUE
      SBMS = SUM / WTM
C
C     end of SUBROUTINE CKSBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMH  (T, ICKWRK, RCKWRK, SMH)*
C  Returns the array of entropies minus enthalpies for species.
C  It is normally not called directly by the user.
C
C  INPUT
C  T(*)      - Real array, temepratures;
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SMH(*)    - Real array, entropy minus enthalpy for species,
C              SMH(K) = S(K)/R - H(K)/RT;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SMH(*)
      SAVE TN1, TNLOG, TNHALF, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TNHALF = TN1 / 2
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 6
         TN3 = TN3 / 12
         TN4 = TN4 / 20
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TKHALF = TK1 / 2
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TKHALF = TNHALF
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMH(K+1) = RCKWRK(NA1) * (TKLOG - 1.0)
     1            + RCKWRK(NA1+1) * TKHALF + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            - RCKWRK(NA1 + NCP1 - 1) / TK1
     4            + RCKWRK(NA1 + NCP2 - 1)
 250  CONTINUE
C
C     end of SUBROUTINE CKSMH
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSML  (T, ICKWRK, RCKWRK, SML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSML  (T, ICKWRK, RCKWRK, SML)
C  Returns the standard state entropies in molar units.
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SML(*)    - Real array, standard state entropies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(mole*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SML(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SML(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            + RCKWRK(NA1 + NCP2 - 1))
250   CONTINUE
C
C     end of SUBROUTINE CKSML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSMS  (T, ICKWRK, RCKWRK, SMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSMS  (T, ICKWRK, RCKWRK, SMS)
C  Returns the standard state entropies in mass units; see Eq. (28).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SMS(*)    - Real array, standard state entropies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/(gm*K)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SMS(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SMS(K+1) = RCKWRK(NcRU) * (RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     1            + RCKWRK(NA1 + NCP2 - 1) ) / RCKWRK(NcWT + K)
  250 CONTINUE
C
C     end of SUBROUTINE CKSMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,
     1                   RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSNUM (LINE, NEXP, LOUT, KRAY, NN, KNUM, NVAL,
C                     RVAL, KERR)
C  Search a character string, LINE, for (1) a character substring which
C  may also appear in an array of character substrings KRAY, and
C  (2) some number of character substrings representing numbers.
C  In the case of (1), if the character substring appears in KRAY,
C  KNUM is its index position.
C  In the case of (2), the character substrings are converted to
C  NVAL real numbers and stored in RVAL, until NEXP are converted.
C
C  This allows format-free input of combined alpha-numeric data.
C  For example, the subroutine might be called to find a Chemkin
C  species index and convert the other substrings to real values:
C
C     input:  LINE    = "N2  1.2"
C             NEXP    = 1, the number of values expected
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C             KRAY(*) = "H2" "O2" "N2" "H" "O" "N" "OH" "H2O" "NO"
C             NN      = 9, the number of entries in KRAY(*)
C     output: KNUM    = 3, the index number of the substring in
C                       KRAY(*) which corresponds to the first
C                       substring in LINE
C             NVAL    = 1, the number of values found in LINE
C                       following the first substring
C             RVAL(*) = 1.200E+00, the substring converted to a number
C             KERR    = .FALSE.
C  INPUT
C  LINE      - Character string; length depends on calling routine.
C  NEXP      - Integer scalar, number of values to be found in LINE.
C              If NEXP < 0, then IABS(NEXP) values are expected, but
C              it is not an error condition if less values are found.
C  LOUT      - Integer scalar, formatted output file unit.
C  KRAY(*)   - Character string array.
C  NN        - Integer scalar, total number of character strings
C              in KRAY.
C
C  OUTPUT
C  KNUM      - Integer scalar, index of character string in KRAY
C              which corresponds to the first substring in LINE.
C  NVAL      - Integer scalar, count of real values found in LINE.
C  RVAL(*)   - Real array, real values found in LINE; dimension at least
C              NEXP.
C  KERR      - Logical, syntax or dimensioning error flag;
C              corresponding string not found, or total of
C              values found is not the number of values expected,
C              will result in KERR = .TRUE.
C
C  END PROLOGUE
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), KRAY(*)*(*), ISTR*80
      DIMENSION RVAL(*)
      LOGICAL KERR, IERR
      INTEGER CKFRCH, CKLSCH, CKSLEN
      EXTERNAL CKFRCH, CKLSCH, CKSLEN
C
      NVAL = 0
      KERR = .FALSE.
      ILEN = CKSLEN(LINE)
      IF (ILEN .LE. 0) RETURN
C
      I1 = CKFRCH(LINE(1:ILEN))
      I3 = INDEX(LINE(I1:ILEN),' ')
      IF (I3 .EQ. 0) I3 = ILEN - I1 + 1
      I2 = I1 + I3
      ISTR = ' '
      ISTR = LINE(I1:I2-1)
C
      CALL CKCOMP (ISTR, KRAY, NN, KNUM)
      IF (KNUM.EQ.0) THEN
         LT = MAX (CKLSCH(ISTR), 1)
         WRITE (LOUT,'(A)')
     1   ' Error in CKSNUM...'//ISTR(1:LT)//' not found...'
         KERR = .TRUE.
      ENDIF
C
      ISTR = ' '
      ISTR = LINE(I2:ILEN)
      IF (NEXP .NE. 0)
     1      CALL CKXNUM (ISTR, NEXP, LOUT, NVAL, RVAL, IERR)
C
C     end of SUBROUTINE CKSNUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSOR  (T, ICKWRK, RCKWRK, SOR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSOR  (T, ICKWRK, RCKWRK, SOR)
C  Returns the nondimensional entropies;  see Eq. (21).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  SOR(*)    - Real array, nondimensional entropies for species;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), SOR(*)
      SAVE TN1, TNLOG, TN2, TN3, TN4
      DATA TN1/1.0/
C
      IF (T(1) .NE. TN1) THEN
C        FIRST of the species-specific temperature array (default)
         TN1 = T(1)
         TNLOG = LOG(TN1)
         TN2 = TN1*TN1
         TN3 = TN1*TN2
         TN4 = TN1*TN3
         TN2 = TN2 / 2
         TN3 = TN3 / 3
         TN4 = TN4 / 4
      ENDIF
      NKM1 = NKK - 1
      DO 250 K = 0, NKM1
         IF (ICKWRK(IcKTF + K) .NE. 1) THEN
C           different temperature required by this species
            TK1 = T(ICKWRK(IcKTF + K))
            TKLOG = LOG(TK1)
            TK2 = TK1*TK1
            TK3 = TK1*TK2
            TK4 = TK1*TK3
            TK2 = TK2 / 2
            TK3 = TK3 / 3
            TK4 = TK4 / 4
         ELSE
            TK1 = TN1
            TKLOG = TNLOG
            TK2 = TN2
            TK3 = TN3
            TK4 = TN4
         ENDIF
C
C        number of temperature ranges for this species
         NTR = ICKWRK(IcNT + K) - 1
C        location of FIRST set of thermodynamic fit coefficients
         NA1 = NcAA + K*NCP2T
C        location of upper limit of FIRST temperature range
         KTEMP = NcTT + K*MXTP + 1
C
  200    CONTINUE
         IF (NTR.GT.1 .AND. TK1.GT.RCKWRK(KTEMP)) THEN
C           Remaining number of temperature ranges
            NTR = NTR - 1
C           Location of next set of fit coefficients
            NA1 = NA1 + NCP2
            KTEMP = KTEMP + 1
C           Check against next temperature, unless last
            IF (NTR .GT. 1) GO TO 200
         ENDIF
C
         SOR(K+1) = RCKWRK(NA1) * TKLOG
     1            + RCKWRK(NA1+1) * TK1 + RCKWRK(NA1+2) * TK2
     2            + RCKWRK(NA1+3) * TK3 + RCKWRK(NA1+4) * TK4
     3            + RCKWRK(NA1 + NCP2 - 1)
250   CONTINUE
C
C     end of SUBROUTINE CKSOR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
C  Returns an array of substrings in a character string with blanks
C  or tabs as delimiters
C
C  INPUT
C  LINE      - Character string; length determined by calling routine.
C  LOUT      - Integer scalar, formatted output file unit.
C  NDIM      - Integer scalar, dimension of a character string array.
C
C  OUTPUT
C  SUB(*)    - Character string array, the character substrings of
C              LINE; dimension SUB at least NDIM.
C  NFOUND    - Integer scalar, count of substrings found in LINE.
C  KERR      - Logical, error flag; dimensioning errors will result in
C              KERR = .TRUE.
C
C  END PROLOGUE
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER SUB(*)*(*), LINE*(*)
      LOGICAL KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
      NFOUND = 0
      ILEN = LEN(SUB(1))
C
      IEND = 0
      KERR = .FALSE.
   25 CONTINUE
C
      ISTART = IEND + 1
      DO 100 L = ISTART, CKSLEN(LINE)
C
         IF (LINE(L:L) .NE. ' ' .AND. LINE(L:L).NE.CHAR(9)) THEN
            IEND   = INDEX(LINE(L:), ' ')
            IF (IEND .EQ. 0) THEN
               IEND = CKSLEN(LINE)
            ELSE
               IEND = L + IEND - 1
            ENDIF
            IF (IEND-L+1 .GT. ILEN) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...substring too long'
               KERR = .TRUE.
            ELSEIF (NFOUND+1 .GT. NDIM) THEN
               WRITE (LOUT,*) ' Error in CKSUBS...NDIM too small'
               KERR = .TRUE.
            ELSE
               NFOUND = NFOUND + 1
               SUB(NFOUND) = LINE(L:IEND)
            ENDIF
            GO TO 25
         ENDIF
C
  100 CONTINUE
C
C     end of SUBROUTINE CKSUBS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYME (CCKWRK, LOUT, ENAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYME (CCKWRK, LOUT, ENAME, KERR)*
C  Returns the character strings of element names.
C
C  INPUT
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C  LOUT      - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  ENAME(*)  - Character string array, element names; dimension at
C              least MM, the total element count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      CHARACTER CCKWRK(*)*(*), ENAME(*)*(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(ENAME(1))
      DO 150 M = 1, NMM
         LT = CKLSCH(CCKWRK(IcMM+M-1))
         ENAME(M) = ' '
         IF (LT .LE. ILEN) THEN
            ENAME(M) = CCKWRK(IcMM+M-1)
         ELSE
            WRITE (LOUT,'(A)')
     1      ' Error in CKSYME...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
C
C     end of SUBROUTINE CKSYME
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYMR (I, LOUT, ICKWRK, RCKWRK, CCKWRK, LT, ISTR,
     1                   KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYMR (I, ICKWRK, RCKWRK, CCKWRK, LT, ISTR, KERR)*
C  Returns a character string which describes the Ith reaction,
C  and the effective length of the character string.
C
C  INPUT
C  I         - Integer scalar, reaction index.
C  LOUT      - Integer scalar, formatted output file unit.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENCCK.
C
C  OUTPUT
C  ISTR      - Character string, description of reaction I.
C  LT        - Integer scalar, number of non-blank characters in ISTR.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      COMMON /CMIN/ CKMIN
C
      DIMENSION ICKWRK(*), RCKWRK(*)
      CHARACTER CCKWRK(*)*(*), ISTR*(*), IDUM*80
      LOGICAL KERR, IERR
      INTEGER CKLKUP, CKLSCH
      EXTERNAL CKLKUP, CKLSCH
C
      ISTR = ' '
      ILEN = LEN(ISTR)
      KERR = .FALSE.
C
      ISRNU = CKLKUP (I, ICKWRK(IcRNU), NRNU)
      ISFAL = CKLKUP (I, ICKWRK(IcFL),  NFAL)
      ISTHB = CKLKUP (I, ICKWRK(IcTB),  NTHB)
      ISWL  = CKLKUP (I, ICKWRK(IcWL),  NWL)
C
      KFAL = 0
      IF (ISFAL.GT.0) KFAL = ICKWRK(IcKF + ISFAL - 1)
      IF (ISWL.GT.0)  WL   = RCKWRK(NcWL + ISWL  - 1)
C
      NSTART = 1
      NEND = MXSP/2
   50 CONTINUE
      N1 = NSTART
      N2 = NEND
      DO 100 N = N1, N2
         K = ICKWRK(IcNK + (I-1)*MXSP + N - 1)
C
         IF (K .GT. 0) THEN
C           append this species
            IF (N .GT. N1) THEN
C              need '+'
               NEXT = CKLSCH(ISTR)+1
               IF (NEXT .GT. ILEN) GO TO 100
               ISTR(NEXT:NEXT) = '+'
            ENDIF
            IDUM = ' '
            IF (ISRNU .EQ. 0) THEN
C              integer stoichiometric coefficient
               NU = ABS(ICKWRK(IcNU + (I-1)*MXSP + N - 1))
            ELSE
C              real stoichiometric coefficient
               RNU = ABS(RCKWRK(NcRNU + (ISRNU-1)*MXSP + N - 1))
C              integer part of coefficient
               NU = RNU
               IF (RNU-NU .NE. 0.0) THEN
C                 RNU has decimal part > 0, else can treat
C                 as whole number
                  NU = 0
C                 convert RNU to character string
                  CALL CKR2CH (RNU, IDUM, L, IERR)
                  IF (IERR) GO TO 200
                  IND = INDEX(IDUM,'.')
C                 restrict decimal part to 3 digits
                  L = MIN (L, IND+3)
               ENDIF
            ENDIF
            IF (NU .GT. 1) THEN
C              convert NU to character string
               CALL CKI2CH (NU, IDUM, L, IERR)
               IF (IERR) GO TO 200
            ENDIF
C
            IF (IDUM .NE. ' ') THEN
C              coefficient > 1 added to string
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+L .GT. ILEN) GO TO 400
               ISTR(NEXT:) = IDUM
            ENDIF
C
C           species name
            NEXT = CKLSCH(ISTR) + 1
            L = CKLSCH(CCKWRK(IcKK + K - 1))
            IF (NEXT+L .GT. ILEN) GO TO 400
            ISTR(NEXT:) = CCKWRK(IcKK + K - 1)
         ENDIF
C
         IF (K.EQ.0 .OR. N.EQ.N2) THEN
C           last product or reactant - supplemental info
            IF (ISWL .GT. 0) THEN
C              radiation enhancement factor
               IF (N .EQ. MXSP/2 .AND. WL.LT.0 .OR.
     1             N .EQ. MXSP   .AND. WL.GT.0) THEN
                   NEXT = CKLSCH(ISTR) + 1
                   IF (NEXT+3 .GT. ILEN) GO TO 400
                   ISTR(NEXT:) = '+HV'
               ENDIF
            ENDIF
C
            IF (ISFAL .GT. 0) THEN
C              pressure-dependence third-body species
               IDUM = '(+'
               IF (KFAL .GT. 0) THEN
                  IDUM(3:) = CCKWRK(IcKK + KFAL - 1)
               ELSE
                  IDUM(3:) = 'M'
               ENDIF
               NEXT = CKLSCH(IDUM) + 1
               IDUM(NEXT:) = ')'
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+CKLSCH(IDUM)-1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = IDUM
            ELSEIF (ISTHB .GT. 0) THEN
C              third-body reaction
               NEXT = CKLSCH(ISTR) + 1
               IF (NEXT+1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '+M'
            ENDIF
C
C           quit at end of products
            IF (N .GT. MXSP/2) GO TO 150
C
C           delimeter at end of reactants
            NEXT = CKLSCH(ISTR) + 1
            IF (ICKWRK(IcNS+I-1) .LT. 0) THEN
               IF (NEXT+1 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '=>'
            ELSE
               IF (NEXT+2 .GT. ILEN) GO TO 400
               ISTR(NEXT:) = '<=>'
            ENDIF
            NSTART = NEND + 1
            NEND = MXSP
            GO TO 50
         ENDIF
  100 CONTINUE
C
  150 CONTINUE
      LT = CKLSCH(ISTR)
      RETURN
C
  200 CONTINUE
      WRITE (LOUT, 300)
  300 FORMAT (' Error in CKSYMR...character string length too small')
      ISTR = ' '
      LT = 0
      KERR = .TRUE.
      RETURN
C
  400 CONTINUE
      WRITE (LOUT, 500)
  500 FORMAT (' Syntax error in CKSYMR...')
      ISTR = ' '
      LT = 0
      KERR = .TRUE.
C
C     end of SUBROUTINE CKSYMR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKSYMS (CCKWRK, LOUT, KNAME, KERR)*
C  Returns the character strings of species names
C
C  INPUT
C  CCKWRK(*) - Character string workspace array;
C              dimension at least LENRCK.
C  LOUT      - Integer scalar, formatted output file unit.
C
C  OUTPUT
C  KNAME(*)  - Character string array, species names;
C              dimension at least KK, the total species count.
C  KERR      - Logical, character length error flag.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      CHARACTER CCKWRK(*)*(*), KNAME(*)*(*)
      LOGICAL KERR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
      KERR = .FALSE.
      ILEN = LEN(KNAME(1))
      DO 150 K = 1, NKK
         LT = CKLSCH(CCKWRK(IcKK + K - 1))
         KNAME(K) = ' '
         IF (LT .LE. ILEN) THEN
            KNAME(K) = CCKWRK(IcKK + K - 1)
         ELSE
            WRITE (LOUT,*)
     1      ' Error in CKSYM...character string length too small '
            KERR = .TRUE.
         ENDIF
150   CONTINUE
C
C     end of SUBROUTINE CKSYMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKTHB  (KDIM, ICKWRK, RCKWRK, AKI)
C
C  START PROLOGUE
C
C  SUBROUTINE CKTHB  (KDIM, ICKWRK, RCKWRK, AKI)
C  Returns matrix of enhanced third body coefficients; see Eq. (58).
C
C  INPUT
C  KDIM      - Integer scalar, first dimension of the matrix AKI;
C              KDIM must be at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  AKI(*,*)  - Real matrix, enhanced third body efficiencies of the
C              species in the reactions;
C              dimension at least KK for first, the total species count,
C              and at least II for the second, the total reaction count.
C              AKI(K,I) is the enhanced efficiency of species K in
C              reaction I.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), AKI(KDIM,*)
C
      DO 150 I = 1, NII
         DO 100 K = 1, NKK
            AKI(K,I) = 1.0
  100    CONTINUE
  150 CONTINUE
C
      I_KTHB = IcKT
      I_KFAC = NcKT
      DO 250 N = 0, NTHB - 1
         I = ICKWRK(IcTB + N)
         DO 200 L = 0, ICKWRK(IcKN + N) - 1
            K        = ICKWRK(I_KTHB + L)
            AKI(K,I) = RCKWRK(I_KFAC + L)
  200    CONTINUE
         I_KTHB = I_KTHB + MXTB
         I_KFAC = I_KFAC + MXTB
  250 CONTINUE
C
C     end of SUBROUTINE CKTHB
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUBML (T, X, ICKWRK, RCKWRK, UBML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUBML (T, X, ICKWRK, RCKWRK, UBML)
C  Returns the mean internal energy of the mixture in molar units;
C  see Eq. (39).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UBML      - Real scalar, mean internal energy.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKUML (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      UBML = 0.0
      DO 100 K = 1, NKK
         UBML = UBML + X(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKUBM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUBMS (T, Y, ICKWRK, RCKWRK, UBMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUBMS (T, Y, ICKWRK, RCKWRK, UBMS)
C  Returns the mean internal energy of the mixture in mass units;
C  see Eq. (40).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UBMS      - Real scalar, mean internal energy.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*)
C
      CALL CKUMS (T, ICKWRK, RCKWRK, RCKWRK(NcK1))
C
      UBMS = 0.0
      DO 100 K = 1, NKK
         UBMS = UBMS + Y(K)*RCKWRK(NcK1 + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKUBMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUML  (T, ICKWRK, RCKWRK, UML)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUML  (T, ICKWRK, RCKWRK, UML)
C  Returns the internal energies in molar units;  see Eq. (23).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UML(*)    - Real array, internal energies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), UML(*)
C
      CALL CKHML (T, ICKWRK, RCKWRK, UML)
      RU = RCKWRK(NcRU)
      DO 150 K = 1, NKK
         UML(K) = UML(K) - RU * T(ICKWRK(IcKTF + K - 1))
150   CONTINUE
C
C     end of SUBROUTINE CKUML
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKUMS  (T, ICKWRK, RCKWRK, UMS)
C
C  START PROLOGUE
C
C  SUBROUTINE CKUMS  (T, ICKWRK, RCKWRK, UMS)
C  Returns the internal energies in mass units;  see Eq. (30).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  UMS(*)    - Real array, internal energies for species;
C              dimension at least KK, the total species count.
C                 cgs units, ergs/gm
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), ICKWRK(*), RCKWRK(*), UMS(*)
C
      CALL CKHMS (T, ICKWRK, RCKWRK, UMS)
      RU = RCKWRK(NcRU)
      DO 150 K = 1, NKK
         TEMP = T(ICKWRK(IcKTF + K - 1))
         UMS(K) = UMS(K) - TEMP*RU/RCKWRK(NcWT + K - 1)
150   CONTINUE
C
C     end of SUBROUTINE CKUMS
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWC   (T, C, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWC   (T, C, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given
C  temperature(s) and molar concentrations;  see Eq. (49).
C
C  INPUT
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), C(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKRATT (RCKWRK, ICKWRK, T, ICKWRK(IcNS), ICKWRK(IcNU),
     2             ICKWRK(IcNK), RCKWRK(NcCO), ICKWRK(IcRV),
     3             RCKWRK(NcRV), ICKWRK(IcLT), RCKWRK(NcLT),
     4             ICKWRK(IcRL), RCKWRK(NcRL), RCKWRK(NcK1),
     5             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcEI),
     6             ICKWRK(IcET), ICKWRK(IcJN), RCKWRK(NcJN),
     7             ICKWRK(IcF1), RCKWRK(NcF1), ICKWRK(IcTD),
     8             ICKWRK(IcTK), ICKWRK(IcKTF), RCKWRK(NcKF),
     9             RCKWRK(NcKR), RCKWRK(NcI1))
      CALL CKRATX (T, C, ICKWRK(IcNU), ICKWRK(IcNK), RCKWRK(NcCO),
     2             ICKWRK(IcFL), ICKWRK(IcFO), ICKWRK(IcFT),
     3             ICKWRK(IcKF), RCKWRK(NcFL), ICKWRK(IcTB),
     3             ICKWRK(IcKN), RCKWRK(NcKT), ICKWRK(IcKT),
     4             RCKWRK(NcKF), RCKWRK(NcKR), RCKWRK(NcI1),
     5             RCKWRK(NcI2), RCKWRK(NcI3), ICKWRK(IcRNU),
     6             RCKWRK(NcRNU), ICKWRK(IcORD), ICKWRK(IcKOR),
     7             RCKWRK(NcKOR), ICKWRK(IcMO), ICKWRK(IcXS))
      CALL CKDOT (RCKWRK(NcI1), RCKWRK(NcI2), ICKWRK, RCKWRK,
     1            RCKWRK(NcK1), RCKWRK(NcK2))
C
      NKM1 = NKK - 1
      DO 25 K = 0, NKM1
         WDOT(K+1) = RCKWRK(NcK1 + K) - RCKWRK(NcK2 + K)
   25 CONTINUE
C
C     end of SUBROUTINE CKWC
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWL   (ICKWRK, RCKWRK, WL)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWL   (ICKWRK, RCKWRK, WL)
C  Returns a set of flags providing information on the wave length
C  of photon radiation
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WL(*)     - Real array, radiation wavelengths for reactions;
C              dimension at least II, total reaction count.
C                 cgs units, angstrom.
C              WL(I)= 0.  reaction I does not have radiation as
C                         either a reactant or product
C              WL(I)=-A   reaction I has radiation of wavelength A
C                         as a reactant
C              WL(I)=+A   reaction I has radiation of wavelength A
C                         as a product
C              If A = 1.0 then no wavelength information was given;
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), WL(*)
C
      DO 100 I = 1, NII
         WL(I) = 0.0
  100 CONTINUE
      DO 150 N = 1, NWL
         WL(ICKWRK(IcWL+N-1)) = RCKWRK(NcWL+N-1)
  150 CONTINUE
C
C     end of SUBROUTINE CKWL
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWT   (ICKWRK, RCKWRK, WT)
C  Returns the molecular weights of the species
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WT(*)     - Real array, molecular weights of the species;
C              dimension at least KK, the total species count.
C                 cgs units, gm/mole
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*), RCKWRK(*), WT(*)
C
      DO 100 K = 1, NKK
         WT(K) = RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
C     end of SUBROUTINE CKWT
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWXP  (P, T, X, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWXP  (P, T, X, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mole fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKXTCP (P, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWXP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWXR  (RHO, T, X, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWXR  (RHO, T, X, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given mass
C  density, temperature(s) and mole fractions;  see Eq. (49).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKXTCR (RHO, T, X, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWXR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mass fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWYP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYPK  (P, T, Y, RKFT, RKRT, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYPK  (P, T, Y, RKFT, RKRT, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given pressure,
C  temperature(s) and mass fractions;  see Eq. (49).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  RKFT(*)   - Real array, forward reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  RKRT(*)   - Real array, reverse reaction rates for reactions;
C              dimension at least II, the total reaction count.
C                 cgs units, depends on the reaction
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), RKFT(*), RKRT(*), ICKWRK(*), RCKWRK(*),
     1          WDOT(*)
C
      DO 25 I = 1, NII
         RCKWRK(NcKF + I - 1) = RKFT(I)
         RCKWRK(NcKR + I - 1) = RKRT(I)
   25 CONTINUE
      CALL CKYTCP (P, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK1))
      CALL CKRATX (T, RCKWRK(NcK1), ICKWRK(IcNU), ICKWRK(IcNK),
     2             RCKWRK(NcCO), ICKWRK(IcFL), ICKWRK(IcFO),
     3             ICKWRK(IcFT), ICKWRK(IcKF), RCKWRK(NcFL),
     3             ICKWRK(IcTB), ICKWRK(IcKN), RCKWRK(NcKT),
     4             ICKWRK(IcKT), RCKWRK(NcKF), RCKWRK(NcKR),
     5             RCKWRK(NcI1), RCKWRK(NcI2), RCKWRK(NcI3),
     6             ICKWRK(IcRNU), RCKWRK(NcRNU), ICKWRK(IcORD),
     7             ICKWRK(IcKOR), RCKWRK(NcKOR), ICKWRK(IcMO),
     8             ICKWRK(IcXS))
      CALL CKDOT (RCKWRK(NcI1), RCKWRK(NcI2), ICKWRK, RCKWRK,
     1            RCKWRK(NcK1), RCKWRK(NcK2))
C
      NKM1 = NKK - 1
      DO 50 K = 0, NKM1
         WDOT(K+1) = RCKWRK(NcK1 + K) - RCKWRK(NcK2 + K)
   50 CONTINUE
C
C     end of SUBROUTINE CKWYPK
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C
C  START PROLOGUE
C
C  SUBROUTINE CKWYR  (RHO, T, Y, ICKWRK, RCKWRK, WDOT)
C  Returns the molar production rates of the species given mass
C  density, temperature and mass fractions;  see Eq. (49).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature;
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  WDOT(*)   - Real array, chemical production rates of the species;
C              dimension at least KK, the total species count.
C                 cgs units, moles/(cm**3*sec)
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), WDOT(*)
C
      CALL CKYTCR (RHO, T, Y, ICKWRK, RCKWRK, RCKWRK(NcK4))
      CALL CKWC   (T, RCKWRK(NcK4), ICKWRK, RCKWRK, WDOT)
C
C     end of SUBROUTINE CKWYR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
C  Searches a character string, LINE, for NEXP space-delimited
C  substrings representing numbers, until NVAL real values are
C  converted and stored in the array, RVAL.
C  This allows format-free input of numerical data.  For example:
C
C     input:  LINE    = " 0.170E+14 0 47780.0"
C             NEXP    = 3, the number of values requested
C             LOUT    = 6, a logical unit number on which to write
C                       diagnostic messages.
C     output: NVAL    = 3, the number of values found
C             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
C             KERR    = .FALSE.
C
C  INPUT
C  LINE   - Character string, length established by calling program.
C  NEXP   - Integer scalar, number of real values to be found in LINE;
C           If NEXP < 0 then IABS(NEXP) values are expected, but
C           it is not an error condition if fewer values are found.
C  LOUT   - Integer scalar, output unit for printed diagnostics.
C
C  OUTPUT
C  NVAL   - Integer scalar, count of real values found in LINE.
C  RVAL   - Real array, values converted from characters in LINE;
C           dimension at least NEXP.
C  KERR   - Logical, syntax or dimensioning error flag.
C
C  END PROLOGUE
C
C     A '!' will comment out a line, or remainder of the line.
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*), ITEMP*80
      DIMENSION RVAL(*), RTEMP(80)
      LOGICAL KERR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
C
C----------Find Comment String (! signifies comment)
C
      ILEN = CKSLEN(LINE)
      NVAL = 0
      KERR = .FALSE.
C
      IF (ILEN .LE. 0) RETURN
      IF (ILEN .GT. 80) THEN
         WRITE (LOUT,*)     ' Error in CKXNUM...line length > 80 '
         WRITE (LOUT,'(A)') LINE
         KERR = .TRUE.
         RETURN
      ENDIF
C
      ITEMP = LINE(1:ILEN)
      CALL CKDTAB (ITEMP)
      IF (NEXP .LT. 0) THEN
         CALL CKPARR (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
      ELSE
         CALL CKPARR (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
         IF (IERR .EQ. 1) THEN
            WRITE (LOUT, *)    ' Syntax errors in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
         ELSEIF (NVAL .NE. NEXP) THEN
            WRITE (LOUT,*) ' Error in CKXNUM...'
            WRITE (LOUT,'(A)') LINE
            KERR = .TRUE.
            WRITE (LOUT,*) NEXP,' values expected, ',
     1                     NVAL,' values found.'
         ENDIF
      ENDIF
      IF (NVAL .LE. IABS(NEXP)) THEN
         DO 20 N = 1, NVAL
            RVAL(N) = RTEMP(N)
   20    CONTINUE
      ENDIF
C
C     end of SUBROUTINE CKXNUM
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTCP (P, T, X, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given pressure, temperature(s)
C  and mole fractions;  see Eq. (10).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUMXT = 0.0
      DO 100 K = 1, NKK
         SUMXT = SUMXT + X(K)*T(ICKWRK(IcKTF + K - 1))
 100  CONTINUE
      PRUT = P/(RCKWRK(NcRU)*SUMXT)
      DO 150 K = 1, NKK
         C(K) = X(K)*PRUT
150   CONTINUE
C
C     end of SUBROUTINE CKXTCP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTCR (RHO, T, X, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTCR (RHO, T, X, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given mass density, temperature(s),
C  and mole fractions;  see Eq. (11).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), X(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K)*RCKWRK(NcWT + K - 1)
  100 CONTINUE
      RHOW = RHO / SUM
      DO 200 K = 1, NKK
         C(K) = X(K)*RHOW
200   CONTINUE
C
C     end of SUBROUTINE CKXTCR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
C
C  START PROLOGUE
C
C  SUBROUTINE CKXTY  (X, ICKWRK, RCKWRK, Y)
C  Returns the mass fractions given mole fractions; see Eq. (9).
C
C  INPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION X(*), ICKWRK(*), RCKWRK(*), Y(*)
C
      SUM = 0.0
      DO 100 K = 1, NKK
         SUM = SUM + X(K) * RCKWRK(NcWT + K - 1)
  100 CONTINUE
C
      DO 200 K = 1, NKK
         Y(K) = X(K) * RCKWRK(NcWT + K - 1) / SUM
200   CONTINUE
C
C     end of SUBROUTINE CKXTY
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCP (P, T, Y, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given pressure, temperature(s)
C  and mass fractions;  see Eq. (7).
C
C  INPUT
C  P         - Real scalar, pressure.
C                 cgs units, dynes/cm**2
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), C(*)
C
      SUMYOW = 0.0
      NKM1 = NKK - 1
      DO 150 K = 0, NKM1
         SUMYOW = SUMYOW + Y(K+1)*T(ICKWRK(IcKTF + K))
     1                          /RCKWRK(NcWT + K)
150   CONTINUE
      SUMYOW = SUMYOW*RCKWRK(NcRU)
      DO 200 K = 1, NKK
         C(K) = P*Y(K)/(SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
C
C     end of SUBROUTINE CKYTCP
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTCR (RHO,T, Y, ICKWRK, RCKWRK, C)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTCR (RHO,T, Y, ICKWRK, RCKWRK, C)
C  Returns the molar concentrations given mass density, temperature(s),
C  and mass fractions;  see Eq. (8).
C
C  INPUT
C  RHO       - Real scalar, mass density.
C                 cgs units, gm/cm**3
C  T(*)      - Real array, temperature(s); dimension is determined by
C              the application program to be the total number of
C              species temperatures, nominally 1.
C                 cgs units, K
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  C(*)      - Real array, concentrations of the species;
C              dimension at least KK, the total species count.
C                 cgs units, mole/cm**3
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION T(*), Y(*), ICKWRK(*), RCKWRK(*), C(*)
C
      DO 150 K = 1, NKK
         C(K) = RHO * Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
C
C     end of SUBROUTINE CKYTCR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
C
C  START PROLOGUE
C
C  SUBROUTINE CKYTX  (Y, ICKWRK, RCKWRK, X)
C  Returns the mole fractions given mass fractions;  see Eq. (6).
C
C  INPUT
C  Y(*)      - Real array, mass fractions of the mixture;
C              dimension at least KK, the total species count.
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C  RCKWRK(*) - Real    workspace array; dimension at least LENRCK.
C
C  OUTPUT
C  X(*)      - Real array, mole fractions of the mixture;
C              dimension at least KK, the total species count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*), X(*)
C
      SUMYOW = 0.0
      DO 150 K = 1, NKK
         SUMYOW = SUMYOW + Y(K)/RCKWRK(NcWT + K - 1)
150   CONTINUE
      DO 200 K = 1, NKK
         X(K) = Y(K) / (SUMYOW*RCKWRK(NcWT + K - 1))
200   CONTINUE
C
C     end of SUBROUTINE CKYTX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE PKINDX (ICKWRK, KELECT, KKION)
C
C  START PROLOGUE
C
C  SUBROUTINE PKINDX (ICKWRK, KELECT, KKION)
C  Returns plasma indices for the particular reaction mechanism.
C
C  INPUT
C  ICKWRK(*) - Integer workspace array; dimension at least LENICK.
C
C  OUTPUT
C  KELECT    - Integer scalar, species array index for the electron.
C  KKION     - Integer scalar, total ion count.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END precision > single
C
      INCLUDE 'ckstrt.h'
      DIMENSION ICKWRK(*)
C
      KKION  = NKKI
      KELECT = KEL
C
C     end of SUBROUTINE PKINDX
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
C   BEGIN PROLOGUE  CKPARI
C   REFER TO  IPGETI
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851725   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses integer variables from a character variable.  Called
C            by IPGETI, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  CKPARI may be used for parsing an input record that contains integer
C  values, but was read into a character variable instead of directly
C  into integer variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by CKPARI to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of IVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1,   2,,40000   , ,60'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set IVAL(1) = 1
C     (2) set IVAL(2) = 2
C     (3) leave IVAL(3) unchanged
C     (4) set IVAL(4) = 40000
C     (5) leave IVAL(5) unchanged
C     (6) set IVAL(6) = 60
C
C   CKPARI will print diagnostics on the default output device, if
C   desired.
C
C   CKPARI is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume IVAL = (0, 0, 0) and NEXPEC = 3 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2 ,   3 45 '         (2, 3, 45)                0       3
C  '2.15,,3'               (2, 0, 3)                 1       0
C  '3X, 25, 2'             (0, 0, 0)                 1       0
C  '10000'                 (10000, 0, 0)             2       1
C
C      Assume IVAL = (0, 0, 0, 0) and NEXPEC = -4 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1, 2'                  (1, 2)                    0       2
C  ',,37  400'             (0, 0, 37, 400)           0       4
C  ' 1,,-3,,5'             (1, 0, -3, 0)             3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  IVAL (I,O) - the integer value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to IVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  CKFRCH,CKLSCH
C   END PROLOGUE  CKPARI
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION IVAL(*)
      CHARACTER *8 FMT(14)
      LOGICAL OKINCR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C   FIRST EXECUTABLE STATEMENT  CKPARI
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = CKLSCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set false when a space follows
C--- an integer value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR .OR. NC .EQ. IE) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into integer array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      IES = NC - 1
      NCH = IES - IBS + 1
      DATA FMT/' (I1)', ' (I2)', ' (I3)', ' (I4)', ' (I5)',
     1   ' (I6)', ' (I7)', ' (I8)', ' (I9)', '(I10)',
     2   '(I11)', '(I12)', '(I13)', '(I14)'/
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) IVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .NE. 0 .AND. ICARD .GE. 0) THEN
         IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
         IF (IERR .EQ. 1)
     1   WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
         IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
         IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      ENDIF
C
C     end of SUBROUTINE CKPARI
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
C   BEGIN PROLOGUE  CKPARR
C   REFER TO  IPGETR
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851625   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses real variables from a character variable.  Called
C            by IPGETR, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  CKPARR may be used for parsing an input record that contains real
C  values, but was read into a character variable instead of directly
C  into real variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by CKPARR to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of RVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1.,   2,,4.e-5   , ,6.e-6'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set RVAL(1) = 1.0
C     (2) set RVAL(2) = 2.0
C     (3) leave RVAL(3) unchanged
C     (4) set RVAL(4) = 4.0E-05
C     (5) leave RVAL(5) unchanged
C     (6) set RVAL(6) = 6.0E-06
C
C   CKPARR will print diagnostics on the default output device, if
C   desired.
C
C   CKPARR is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
C  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
C  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
C  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
C  '1.0'                   (1.0, 0.0, 0.0)           2       1
C
C      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1.,2.'                 (1.0, 2.0)                0       2
C  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
C  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  RVAL (I,O) - the real value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to RVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  CKLSCH
C   END PROLOGUE  CKPARR
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION RVAL(*)
      CHARACTER *8 FMT(22)
      LOGICAL OKINCR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
C
C   FIRST EXECUTABLE STATEMENT  CKPARR
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = CKLSCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set negative when a space follows
C--- a real value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into real array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      DATA FMT/     ' (E1.0)', ' (E2.0)', ' (E3.0)', ' (E4.0)',
     1   ' (E5.0)', ' (E6.0)', ' (E7.0)', ' (E8.0)', ' (E9.0)',
     2   '(E10.0)', '(E11.0)', '(E12.0)', '(E13.0)', '(E14.0)',
     3   '(E15.0)', '(E16.0)', '(E17.0)', '(E18.0)', '(E19.0)',
     4   '(E20.0)', '(E21.0)', '(E22.0)'/
      IES = NC - 1
      NCH = IES - IBS + 1
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) RVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .NE. 0 .AND. ICARD .GE. 0) THEN
         IF (ICARD .NE. 0) WRITE(LOUT,'(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
         IF (IERR .EQ. 1)
     1   WRITE(LOUT,'(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
         IF (IERR .EQ. 2) WRITE(LOUT,'(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
         IF (IERR .EQ. 3) WRITE(LOUT,'(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      ENDIF
C
C     end of SUBROUTINE CKPARR
      RETURN
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C     OBSOLETE ROUTINES...
C
      CHARACTER*(*) FUNCTION UPCASE (STR, ILEN)
      CHARACTER STR*(*), CKCHUP*128
      INTEGER ILEN
      EXTERNAL CKCHUP
      UPCASE = CKCHUP (STR, ILEN)
      RETURN
      END
C
      CHARACTER*(*) FUNCTION LOCASE (STR, ILEN)
      CHARACTER STR*(*), CKCHLO*128
      INTEGER ILEN
      EXTERNAL CKCHLO
      LOCASE = CKCHLO (STR, ILEN)
      RETURN
      END
C
      INTEGER FUNCTION IPPLEN (STR)
      CHARACTER*(*) STR
      INTEGER CKSLEN
      EXTERNAL CKSLEN
      IPPLEN = CKSLEN(STR)
      RETURN
      END
C
      INTEGER FUNCTION ILASCH (STR)
      CHARACTER*(*) STR
      INTEGER CKLSCH
      EXTERNAL CKLSCH
      ILASCH = CKLSCH(STR)
      RETURN
      END
C
      INTEGER FUNCTION IFIRCH (STR)
      CHARACTER*(*) STR
      INTEGER CKFRCH
      EXTERNAL CKFRCH
      IFIRCH = CKFRCH(STR)
      RETURN
      END
C
      SUBROUTINE IPPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
C*****precision > double
      DOUBLE PRECISION RVAL(*)
C*****END precision > double
C*****precision > single
C      REAL RVAL(*)
C*****END precision > single
      CHARACTER*(*) STRING
      INTEGER ICARD, NEXPEC, NFOUND, IERR, LOUT
      CALL CKPARR (STRING,ICARD,NEXPEC,RVAL,NFOUND,IERR,LOUT)
      RETURN
      END
C
      SUBROUTINE IPPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
      INTEGER ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT
      CHARACTER*(*) STRING
      DIMENSION IVAL(*)
      CALL CKPARI (STRING,ICARD,NEXPEC,IVAL,NFOUND,IERR,LOUT)
      RETURN
      END
