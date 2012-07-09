      MODULE OBSCHDMODULE
         INTEGER, SAVE, POINTER    ::NQCH,NQCCH,NQTCH,IUCHOBSV
         INTEGER, SAVE, DIMENSION(:),   POINTER ::NQOBCH
         INTEGER, SAVE, DIMENSION(:),   POINTER ::NQCLCH
         INTEGER, SAVE, DIMENSION(:),   POINTER ::IOBTS
         REAL,    SAVE, DIMENSION(:),   POINTER ::FLWSIM
         REAL,    SAVE, DIMENSION(:),   POINTER ::FLWOBS
         REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
         REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
         REAL,    SAVE, DIMENSION(:,:), POINTER ::QCELL
         CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
      TYPE OBSCHDTYPE
         INTEGER,            POINTER ::NQCH,NQCCH,NQTCH,IUCHOBSV
         INTEGER,     DIMENSION(:),  POINTER ::NQOBCH
         INTEGER,     DIMENSION(:),  POINTER ::NQCLCH
         INTEGER,     DIMENSION(:),  POINTER ::IOBTS
         REAL,        DIMENSION(:),  POINTER ::FLWSIM
         REAL,        DIMENSION(:),  POINTER ::FLWOBS
         REAL,        DIMENSION(:),  POINTER ::TOFF
         REAL,        DIMENSION(:),  POINTER ::OTIME
         REAL,        DIMENSION(:,:),POINTER ::QCELL
         CHARACTER*12,DIMENSION(:),  POINTER ::OBSNAM
      END TYPE
      TYPE(OBSCHDTYPE), SAVE  ::OBSCHDDAT(10)
      END MODULE OBSCHDMODULE



      SUBROUTINE OBS2CHD7AR(IUCHOB,IGRID)
C     ******************************************************************
C     ALLOCATE AND READ DATA FOR FLOW OBSERVATIONS AT CONSTANT-HEAD
C     BOUNDARY CELLS
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY: NCOL,NROW,NLAY,NPER,NSTP,PERLEN,TSMULT,ISSFLG,
     1                  IOUT,ITRSS
      USE OBSCHDMODULE
C
      CHARACTER*200 LINE
C     ------------------------------------------------------------------
      ALLOCATE(NQCH,NQTCH,NQCCH,IUCHOBSV)
C
C1------INITIALIZE VARIABLEA.
      ZERO=0.0
      IERR=0
      NT=0
      NC=0
C
C2------IDENTIFY PROCESS
      WRITE(IOUT,14) IUCHOB
   14 FORMAT(/,' OBS2CHD7 -- CONSTANT-HEAD BOUNDARY FLOW OBSERVATIONS',
     &    /,' VERSION 2.0, 02/28/2006       INPUT READ FROM UNIT ',I3)
C
C3------ITEM 1
      CALL URDCOM(IUCHOB,IOUT,LINE)
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQCH,DUM,IOUT,IUCHOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQCCH,DUM,IOUT,IUCHOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NQTCH,DUM,IOUT,IUCHOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUCHOBSV,DUM,IOUT,IUCHOB)
      WRITE (IOUT,17) NQCH, NQCCH, NQTCH
   17 FORMAT (/,
     &    ' NUMBER OF FLOW-OBSERVATION CONSTANT-HEAD-CELL GROUPS:',I5,/,
     &    '   NUMBER OF CELLS IN CONSTANT-HEAD-CELL GROUPS......:',I5,/,
     &    '   NUMBER OF CONSTANT-HEAD-CELL FLOWS................:',I5)
      IF(NQTCH.LE.0) THEN
         WRITE(IOUT,*) ' NQTCH LESS THAN OR EQUAL TO 0'
         CALL USTOP(' ')
      END IF
      IF(IUCHOBSV.GT.0) THEN
         WRITE(IOUT,21) IUCHOBSV
   21    FORMAT(1X,
     1      'CH OBSERVATIONS WILL BE SAVED ON UNIT...............:',I5)
      ELSE
         WRITE(IOUT,22)
   22    FORMAT(1X,'CH OBSERVATIONS WILL NOT BE SAVED IN A FILE')
      END IF
C
C4------ALLOCATE ARRAYS
      ALLOCATE (NQOBCH(NQCH))
      ALLOCATE (NQCLCH(NQCH))
      ALLOCATE (IOBTS(NQTCH))
      ALLOCATE (FLWSIM(NQTCH))
      ALLOCATE (FLWOBS(NQTCH))
      ALLOCATE (TOFF(NQTCH))
      ALLOCATE (OTIME(NQTCH))
      ALLOCATE (QCELL(4,NQCCH))
      ALLOCATE (OBSNAM(NQTCH))
      DO 19 N=1,NQTCH
      OTIME(N)=ZERO
      FLWSIM(N)=ZERO
   19 CONTINUE
C
C5------READ AND WRITE TIME-OFFSET MULTIPLIER FOR FLOW-OBSERVATION TIMES
      READ(IUCHOB,*) TOMULTCH
      WRITE (IOUT,520) TOMULTCH
  520 FORMAT (/,' OBSERVED CONSTANT-HEAD-CELL FLOW DATA',/,
     &' -- TIME OFFSETS ARE MULTIPLIED BY: ',G12.5)
C
C6------LOOP THROUGH CELL GROUPS.
      DO 120 IQ = 1,NQCH
C
C7------READ ITEM 3
        READ (IUCHOB,*) NQOBCH(IQ), NQCLCH(IQ)
        WRITE (IOUT,525) IQ, 'CHD', NQCLCH(IQ), NQOBCH(IQ)
  525   FORMAT (/,'   GROUP NUMBER: ',I3,'   BOUNDARY TYPE: ',A,
     &         '   NUMBER OF CELLS IN GROUP: ',I5,/,
     &         '   NUMBER OF FLOW OBSERVATIONS: ',I5,//,
     &         40X,'OBSERVED',/,
     &         20X,'REFER.',12X,'BOUNDARY FLOW',/,
     &      7X,'OBSERVATION',2X,'STRESS',4X,'TIME',5X,'GAIN (-) OR',/,
     &         2X,'OBS#    NAME',6X,'PERIOD   OFFSET',5X,'LOSS (+)')
C
C8------SET FLAG FOR SETTING ALL FACTORS TO 1
        IFCTFLG = 0
        IF (NQCLCH(IQ).LT.0) THEN
          IFCTFLG = 1
          NQCLCH(IQ) = -NQCLCH(IQ)
        ENDIF
C
C9------READ TIME STEPS, MEASURED FLOWS, AND WEIGHTS.
        NT1 = 1 + NT
        NT2 = NQOBCH(IQ) + NT
        DO 30 N = NT1, NT2
C
C10-----READ ITEM 4
          READ (IUCHOB,*) OBSNAM(N), IREFSP, TOFFSET, FLWOBS(N)
          WRITE (IOUT,535) N, OBSNAM(N), IREFSP, TOFFSET, FLWOBS(N)
  535     FORMAT(1X,I5,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
          CALL UOBSTI(OBSNAM(N),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,
     &                IOBTS(N),PERLEN,TOFF(N),TOFFSET,TOMULTCH,TSMULT,1,
     &                OTIME(N))
   30   CONTINUE
C
C11-----READ LAYER, ROW, COLUMN, AND FACTOR (ITEM 5)
        NC1 = NC + 1
        NC2 = NC + NQCLCH(IQ)
        WRITE (IOUT,540)
  540   FORMAT (/,'       LAYER  ROW  COLUMN    FACTOR')
        DO 40 L = NC1, NC2
          READ (IUCHOB,*) (QCELL(I,L),I=1,4)
          IF(QCELL(4,L).EQ.0. .OR. IFCTFLG.EQ.1) QCELL(4,L) = 1.
          WRITE (IOUT,550) (QCELL(I,L),I=1,4)
  550     FORMAT (4X,F8.0,F6.0,F7.0,F9.2)
          K = QCELL(1,L)
          I = QCELL(2,L)
          J = QCELL(3,L)
          IF (K.LE.0 .OR. K.GT.NLAY .OR .J.LE.0 .OR. J.GT.NCOL .OR.
     &        I.LE.0 .OR. I.GT.NROW) THEN
            WRITE (IOUT,590)
  590       FORMAT (/,' ROW OR COLUMN NUMBER INVALID',
     &        ' -- STOP EXECUTION (OBS2CHD7AR)',/)
            IERR = 1
          ENDIF
   40   CONTINUE
C
C12-----UPDATE COUNTERS.
        NC = NC2
        NT = NT2
  120 CONTINUE
C
C13-----STOP IF THERE WERE ANY ERRORS WHILE READING.
      IF (IERR.GT.0) THEN
        WRITE(IOUT,620)
  620   FORMAT (/,' ERROR:  SEE ABOVE FOR ERROR MESSAGE AND "STOP',
     &        ' EXECUTION" (OBS2CHD7AR)')
        CALL USTOP(' ')
      ENDIF
C
C14-----RETURN.
      CALL SOBS2CHD7PSV(IGRID)
      RETURN
      END
      SUBROUTINE OBS2CHD7SE(KKPER,IGRID)
C     ******************************************************************
C     CALCULATE SIMULATED EQUIVALENTS TO OBSERVED CONSTANT-HEAD FLOWS
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,    ONLY:IBOUND,IOUT
      USE OBSBASMODULE,ONLY:ITS
      USE OBSCHDMODULE
C
      DOUBLE PRECISION RATE
C     ------------------------------------------------------------------
      CALL SOBS2CHD7PNT(IGRID)
C
C1------INITIALIZE VARIABLES
      ZERO = 0.0
      NC = 0
      NT1 = 1
C
C2------LOOP THROUGH BOUNDARY FLOW CELL GROUPS
      DO 60 IQ = 1, NQCH
        NT2 = NT1 + NQOBCH(IQ) - 1
C
C3--------LOOP THROUGH THE OBSERVATION TIMES FOR THIS CELL GROUP.
        DO 40 NT = NT1, NT2
C
C4--------WAS THERE A MEASUREMENT AT THIS BOUNDARY THIS TIME STEP?
          IF (IOBTS(NT).EQ.ITS .OR.
     &        (IOBTS(NT).EQ.ITS-1.AND.TOFF(NT).GT.ZERO)) THEN
C
C5------YES -- LOOP THROUGH CELLS.
            NC1 = NC + 1
            NC2 = NC + NQCLCH(IQ)
            DO 30 N = NC1, NC2
              K = QCELL(1,N)
              I = QCELL(2,N)
              J = QCELL(3,N)
              IF (IBOUND(J,I,K).GE.0) THEN
                WRITE(IOUT,500) K,I,J,KKPER
  500           FORMAT(/,
     &' *** ERROR: CONSTANT-HEAD FLOW OBSERVATION SPECIFIED FOR CELL (',
     &I3,',',I5,',',I5,'),',/,
     &12X,'BUT THIS CELL IS NOT CONSTANT-HEAD IN STRESS PERIOD ',I4,/
     &12X,'-- STOP EXECUTION (OBS2CHD7SE)')
                CALL USTOP(' ')
              ENDIF
C
C6------CALL SUBROUTINE TO CALCULATE CONSTANT-HEAD FLOW FOR CELL
              CALL SOBS2CHD7FFLW(J,I,K,RATE)
C
C7------SUM VALUES FROM INDIVIDUAL CELLS.
C7------CALCULATE FACTOR FOR TEMPORAL INTERPOLATION
   20         FACT = 1.0
              IF (TOFF(NT).GT.ZERO) THEN
                IF (IOBTS(NT).EQ.ITS) FACT = 1. - TOFF(NT)
                IF (IOBTS(NT).EQ.ITS-1) FACT = TOFF(NT)
              ENDIF
C
C8------ACCUMULATE FLOWS FOR THE SIMULATED OBSERVATION.
              FLWSIM(NT) = FLWSIM(NT) + RATE*FACT*QCELL(4,N)
C
Cx------END OF LOOP FOR CELLS IN ONE GROUP.
   30       CONTINUE
C
          ENDIF
C
Cx------END OF LOOP FOR OBSERVATION TIMES FOR ONE GROUP
   40   CONTINUE
C
C10-----UPDATE CELL AND TIME COUNTERS 
        NC = NC + NQCLCH(IQ)
        NT1 = NT2 + 1
C
Cx------END OF LOOP FOR CELL GROUPS.
   60 CONTINUE
C
C11-----RETURN.
      RETURN
      END
      SUBROUTINE SOBS2CHD7FFLW(J,I,K,RATE)
C     ******************************************************************
C     CALCULATE CONSTANT-HEAD BOUNDARY FLOW FOR A GIVEN CELL
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,HNEW,CR,CC,CV,BOTM,NBOTM,
     1                       NCOL,NROW,NLAY,LAYHDT,LBOTM
      USE GWFBASMODULE, ONLY:ICHFLG
C
      DOUBLE PRECISION HD,X1,X2,X3,X4,X5,X6,RATE
C     ------------------------------------------------------------------
C
C6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
      ZERO=0.
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      X6=ZERO
C
C7------CALCULATE FLOW THROUGH THE LEFT FACE.
C7------COMMENTS A-C APPEAR ONLY IN THE SECTION HEADED BY COMMENT 7,
C7------BUT THEY APPLY IN A SIMILAR MANNER TO SECTIONS 8-12.
C
C7A-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C7A-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
C7A-----ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL
C7A-----WHEN ICHFLG IS 0.
      IF(J.EQ.1) GO TO 30
      IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
      IF(ICHFLG.EQ.0 .AND. IBOUND(J-1,I,K).LT.0) GO TO 30
C
C7B-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
      HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
      X1=HDIFF*CR(J-1,I,K)
C
C8------CALCULATE FLOW THROUGH THE RIGHT FACE.
   30 IF(J.EQ.NCOL) GO TO 60
      IF(IBOUND(J+1,I,K).EQ.0) GO TO 60
      IF(ICHFLG.EQ.0 .AND. IBOUND(J+1,I,K).LT.0) GO TO 60
      HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
      X2=HDIFF*CR(J,I,K)
C
C9------CALCULATE FLOW THROUGH THE BACK FACE.
   60 IF(I.EQ.1) GO TO 90
      IF (IBOUND(J,I-1,K).EQ.0) GO TO 90
      IF(ICHFLG.EQ.0 .AND. IBOUND(J,I-1,K).LT.0) GO TO 90
      HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
      X3=HDIFF*CC(J,I-1,K)
C
C10-----CALCULATE FLOW THROUGH THE FRONT FACE.
   90 IF(I.EQ.NROW) GO TO 120
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 120
      IF(ICHFLG.EQ.0 .AND. IBOUND(J,I+1,K).LT.0) GO TO 120
      HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
      X4=HDIFF*CC(J,I,K)
C
C11-----CALCULATE FLOW THROUGH THE UPPER FACE.
  120 IF(K.EQ.1) GO TO 150
      IF (IBOUND(J,I,K-1).EQ.0) GO TO 150
      IF(ICHFLG.EQ.0 .AND. IBOUND(J,I,K-1).LT.0) GO TO 150
      HD=HNEW(J,I,K)
      IF(LAYHDT(K).EQ.0) GO TO 122
      TMP=HD
      TOP=BOTM(J,I,LBOTM(K)-1)
      IF(TMP.LT.TOP) HD=TOP
  122 HDIFF=HD-HNEW(J,I,K-1)
      X5=HDIFF*CV(J,I,K-1)
C
C12-----CALCULATE FLOW THROUGH THE LOWER FACE.
  150 IF(K.EQ.NLAY) GO TO 180
      IF(IBOUND(J,I,K+1).EQ.0) GO TO 180
      IF(ICHFLG.EQ.0 .AND. IBOUND(J,I,K+1).LT.0) GO TO 180
      HD=HNEW(J,I,K+1)
      IF(LAYHDT(K+1).EQ.0) GO TO 152
      TMP=HD
      TOP=BOTM(J,I,LBOTM(K+1)-1)
      IF(TMP.LT.TOP) HD=TOP
  152 HDIFF=HNEW(J,I,K)-HD
      X6=HDIFF*CV(J,I,K)
C
C13-----SUM THE FLOWS THROUGH SIX FACES OF CONSTANT HEAD CELL
 180  RATE=X1+X2+X3+X4+X5+X6
C
C-----RETURN
      RETURN
      END
      SUBROUTINE OBS2CHD7OT(IGRID)
C     ******************************************************************
C     WRITE ALL OBSERVATIONS TO LISTING FILE.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY: IOUT
      USE OBSCHDMODULE
      DOUBLE PRECISION SQ,SUMSQ
C     ------------------------------------------------------------------
      CALL SOBS2CHD7PNT(IGRID)
C
C1------WRITE OBSERVATIONS TO LISTING FILE.
      WRITE(IOUT,17)
   17 FORMAT(1X,/,1X,'CONSTANT HEAD FLOW OBSERVATIONS',/,
     1  1X,'OBSERVATION     OBSERVED      SIMULATED',/
     2  1X,'  NAME            VALUE         VALUE      DIFFERENCE',/
     3  1X,'-------------------------------------------------------')
      SUMSQ=0.
      DO 100 N=1,NQTCH
      DIFF=FLWOBS(N)-FLWSIM(N)
      SQ=DIFF*DIFF
      SUMSQ=SUMSQ+SQ
      WRITE(IOUT,27) OBSNAM(N),FLWOBS(N),FLWSIM(N),DIFF
   27 FORMAT(1X,A,1P,3G14.6)
  100 CONTINUE
      WRITE(IOUT,28) SUMSQ
   28 FORMAT(1X,/,1X,'SUM OF SQUARED DIFFERENCE:',1P,E15.5)
C
C2------WRITE OBSERVATIONS TO SEPARATE FILE.
      IF(IUCHOBSV.GT.0) CALL UOBSSV(IUCHOBSV,NQTCH,FLWSIM,FLWOBS,
     1                              OBSNAM,0)
C
C3------RETURN.
      RETURN
      END
      SUBROUTINE OBS2CHD7DA(IGRID)
C  Deallocate OBSCHD memory
      USE OBSCHDMODULE
C
      DEALLOCATE(NQCH)
      DEALLOCATE(NQTCH)
      DEALLOCATE(NQCCH)
      DEALLOCATE(IUCHOBSV)
      DEALLOCATE(NQOBCH)
      DEALLOCATE(NQCLCH)
      DEALLOCATE(IOBTS)
      DEALLOCATE(FLWSIM)
      DEALLOCATE(FLWOBS)
      DEALLOCATE(TOFF)
      DEALLOCATE(OTIME)
      DEALLOCATE(QCELL)
      DEALLOCATE(OBSNAM)
C
      RETURN
      END
      SUBROUTINE SOBS2CHD7PNT(IGRID)
C  Change OBSCHD data to a different grid.
      USE OBSCHDMODULE
C
      NQCH=>OBSCHDDAT(IGRID)%NQCH
      NQTCH=>OBSCHDDAT(IGRID)%NQTCH
      NQCCH=>OBSCHDDAT(IGRID)%NQCCH
      IUCHOBSV=>OBSCHDDAT(IGRID)%IUCHOBSV
      NQOBCH=>OBSCHDDAT(IGRID)%NQOBCH
      NQCLCH=>OBSCHDDAT(IGRID)%NQCLCH
      IOBTS=>OBSCHDDAT(IGRID)%IOBTS
      FLWSIM=>OBSCHDDAT(IGRID)%FLWSIM
      FLWOBS=>OBSCHDDAT(IGRID)%FLWOBS
      TOFF=>OBSCHDDAT(IGRID)%TOFF
      OTIME=>OBSCHDDAT(IGRID)%OTIME
      QCELL=>OBSCHDDAT(IGRID)%QCELL
      OBSNAM=>OBSCHDDAT(IGRID)%OBSNAM
C
      RETURN
      END
      SUBROUTINE SOBS2CHD7PSV(IGRID)
C  Save OBSCHD data for a grid.
      USE OBSCHDMODULE
C
      OBSCHDDAT(IGRID)%NQCH=>NQCH
      OBSCHDDAT(IGRID)%NQTCH=>NQTCH
      OBSCHDDAT(IGRID)%NQCCH=>NQCCH
      OBSCHDDAT(IGRID)%IUCHOBSV=>IUCHOBSV
      OBSCHDDAT(IGRID)%NQOBCH=>NQOBCH
      OBSCHDDAT(IGRID)%NQCLCH=>NQCLCH
      OBSCHDDAT(IGRID)%IOBTS=>IOBTS
      OBSCHDDAT(IGRID)%FLWSIM=>FLWSIM
      OBSCHDDAT(IGRID)%FLWOBS=>FLWOBS
      OBSCHDDAT(IGRID)%TOFF=>TOFF
      OBSCHDDAT(IGRID)%OTIME=>OTIME
      OBSCHDDAT(IGRID)%QCELL=>QCELL
      OBSCHDDAT(IGRID)%OBSNAM=>OBSNAM
C
      RETURN
      END
