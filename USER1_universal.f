C=======================================================================
C
C     USER1 : CFRP THERMOCHEMICAL DECOMPOSITION FOR NX / TMG
C             universal version (any geometry, any material orientation)
C
C     WHAT IT DOES
C       Two-step decomposition kinetics advance a degree of
C       decomposition alpha in every element of one TMG group. alpha
C       changes the heat capacity C, the decomposition heat load Q and
C       the conductances GG of the group, every time step (KODE = 1).
C
C     WHAT IS READ FROM NX / TMG (nothing is hardcoded here)
C       group name          : GROUPNAME below
C       base rho,cp,k1,k2,k3: the NX material of the group (solver
C                             input XML next to the solver log)
C       material directions : first/second direction vectors per element
C                             in the same XML (third = cross product)
C       base C and GG       : TMG arrays at the first call
C       reference temperature: initial element temperature
C       through-thickness   : the material axis with the lowest base k
C
C     HOW THE PROPERTIES ARE SCALED
C       C  = C0  * [rho(alpha)*cp(T,alpha)] / [rho0*cp(Tref,0)]
C       GG = GG0 / f      (GG is a RESISTANCE in the TMG API)
C            f = k(T,alpha)/k(Tref,0), with the through-thickness law
C            for the link between the two thickness faces of a hex and
C            the in-plane law for every other link.
C       Q  = decomposition heat - gas enthalpy, added to the boundary
C            elements (TMG moves solid heat loads there).
C     The laws and kinetics live in CFRP_MODEL (edit it for another
C     material). Units: SI, temperatures in degC.
C
C     LIMITS (checked and reported at start-up)
C       - hexahedra (brick-like) aligned with the material directions
C       - NX material must be ORTHO or isotropic with CONSTANT values
C       - the first USER1 call with KODE=1 must see the original GG, C
C
C=======================================================================

      SUBROUTINE USER1(GG,T,C,Q,QD,R,TIME,DT,IT,KODE,NOCON,MAXNO,
     + ICONV,DUM1,DUM2,DTP,TF)

      IMPLICIT NONE

      DOUBLE PRECISION T(*)
      DOUBLE PRECISION TIME
      DOUBLE PRECISION DT
      REAL GG(*)
      REAL C(*)
      REAL Q(*)
      REAL QD(*)
      REAL R(*)
      REAL DUM1
      REAL DUM2
      REAL DTP
      REAL TF
      INTEGER IT
      INTEGER KODE
      INTEGER NOCON
      INTEGER MAXNO
      INTEGER ICONV(*)

C     ---- user settings ----------------------------------------------
C     GROUPNAME : NX group (long name) of the decomposing elements
C     XMLCFG    : solver input XML. Leave blank: the name is derived
C                 from the solver log. Only needed if that fails.
C     DIAGINT   : seconds between one-line summaries
C     PRTINT    : seconds between through-thickness temperature tables
C     PROBE     : depth [m] of the interpolated probe temperature
C     THKCOS    : |cos| above which a face counts as a thickness face
C     ZRATIO    : conductance-pattern test (see CFRP_TOPO)
      CHARACTER*80 GROUPNAME
      CHARACTER*100 XMLCFG
      DOUBLE PRECISION DIAGINT
      DOUBLE PRECISION PRTINT
      DOUBLE PRECISION PROBE
      DOUBLE PRECISION THKCOS
      DOUBLE PRECISION ZRATIO
      PARAMETER (GROUPNAME='Plate_3D')
      PARAMETER (XMLCFG=' ')
      PARAMETER (DIAGINT=10.0D0)
      PARAMETER (PRTINT=50.0D0)
      PARAMETER (PROBE=3.0D-3)
      PARAMETER (THKCOS=0.7071D0)
      PARAMETER (ZRATIO=0.05D0)

C     ---- array sizes -----------------------------------------------
      INTEGER MAXELEM
      INTEGER MAXCOND
      INTEGER MAXEXT
      INTEGER MAXCOL
      INTEGER MAXSAMP
      PARAMETER (MAXELEM=50000)
      PARAMETER (MAXCOND=100000)
      PARAMETER (MAXEXT=1000000)
      PARAMETER (MAXCOL=2000)
      PARAMETER (MAXSAMP=15)

C     ---- group, element data, connectivity --------------------------
      CHARACTER*7 SNAME
      CHARACTER*80 LNAME
      CHARACTER*300 XMLNAME
      INTEGER GLIST(MAXELEM)
      INTEGER GNUMELEM
      INTEGER GID
      INTEGER GLENGTH
      REAL ELPROP(10,MAXELEM)
      REAL ELGEOM(3,4,MAXELEM)
      INTEGER NODES(8,1)
      REAL NODELIST(3,1)
      INTEGER ICNDGR(4,MAXCOND)
      INTEGER NGRCOND
      INTEGER ICNDTMP(4,200)
      INTEGER BOUNDARYLIST(200)
      INTEGER NCTMP
      INTEGER NB_TMP
      LOGICAL ISBOUND(MAXELEM)
      LOGICAL ISCFRP(MAXELEM)
      INTEGER EXTMAP(MAXEXT)

C     ---- saved TMG baseline -----------------------------------------
      REAL C0(MAXELEM)
      REAL C0B(MAXELEM)
      REAL GG0(MAXCOND)
      LOGICAL UPDATED(MAXCOND)

C     ---- per-element state ------------------------------------------
      DOUBLE PRECISION ALPHA1(MAXELEM)
      DOUBLE PRECISION ALPHA2(MAXELEM)
      DOUBLE PRECISION QSTORE(MAXELEM)
      DOUBLE PRECISION QWATT(MAXELEM)
      DOUBLE PRECISION VOLE(MAXELEM)
      DOUBLE PRECISION RR(MAXELEM)
      DOUBLE PRECISION TREF(MAXELEM)
      DOUBLE PRECISION FPAR(MAXELEM)
      DOUBLE PRECISION FPERP(MAXELEM)
      DOUBLE PRECISION SUMV(MAXELEM)
      DOUBLE PRECISION SUMVR(MAXELEM)
      DOUBLE PRECISION PITCH(MAXELEM)
      DOUBLE PRECISION WZ(MAXELEM)
      DOUBLE PRECISION WL(MAXELEM)
      DOUBLE PRECISION E1(3,MAXELEM)
      DOUBLE PRECISION E2(3,MAXELEM)
      DOUBLE PRECISION ETA(3,MAXELEM)
      LOGICAL HAVEOR(MAXELEM)

C     ---- topology ---------------------------------------------------
      INTEGER NBND(MAXELEM)
      INTEGER BND(8,MAXELEM)
      INTEGER NADJ(MAXELEM)
      INTEGER ADJN(4,MAXELEM)
      INTEGER BLIST(MAXELEM)
      INTEGER NBLIST
      INTEGER ZB1(MAXELEM)
      INTEGER ZB2(MAXELEM)
      INTEGER NBP(MAXELEM)
      INTEGER NBM(MAXELEM)
      INTEGER METH(MAXELEM)
      INTEGER NMETH(0:3)
      INTEGER LCLS(MAXCOND)
      INTEGER LIDX(MAXCOND)
      INTEGER LFIRST(MAXELEM)
      INTEGER LNUM(MAXELEM)
      INTEGER COLN(MAXCOL)
      DOUBLE PRECISION COLS(MAXCOL)
      INTEGER NCOL
      INTEGER SAMPJ(MAXSAMP)
      INTEGER NSAMP
      REAL GGEXP(MAXSAMP)

C     ---- material model and base properties -------------------------
      DOUBLE PRECISION KIN(9)
      DOUBLE PRECISION HDEC(2)
      DOUBLE PRECISION DENSR
      DOUBLE PRECISION CPV(4)
      DOUBLE PRECISION CPE(4)
      DOUBLE PRECISION KPV(4)
      DOUBLE PRECISION KPE(4)
      DOUBLE PRECISION KNV(4)
      DOUBLE PRECISION KNE(4)
      DOUBLE PRECISION CPG(4)
      DOUBLE PRECISION RHO0
      DOUBLE PRECISION CP0
      DOUBLE PRECISION BK(3)
      INTEGER TAX
      INTEGER MATUID
      INTEGER NDIFF
      INTEGER NNOOR
      INTEGER XSTAT
      INTEGER IU
      LOGICAL XOK

C     ---- time and bookkeeping ---------------------------------------
      DOUBLE PRECISION LASTTIME
      DOUBLE PRECISION TPRINT
      DOUBLE PRECISION NEXTPRT
      DOUBLE PRECISION T2NEXT
      DOUBLE PRECISION SAMPT(3)
      INTEGER ISAMP
      INTEGER N2PERS
      INTEGER CNT(8)
      INTEGER NK(0:9)
      INTEGER NCALL
      LOGICAL FIRST
      LOGICAL DISABLED
      LOGICAL NEWTIME
      LOGICAL OPN

C     ---- locals ------------------------------------------------------
      INTEGER N,J,K,B,L,NMOD,NEQ,NREV,NOTH,IC
      INTEGER NHI,NLO
      DOUBLE PRECISION TC,TK,A1OLD,A2OLD,A1NEW,A2NEW,R1,R2
      DOUBLE PRECISION ALPHA,STEPDT,QCHEM,RHOGDOT,DH,DRHO
      DOUBLE PRECISION CPM,KPM,KNM,CP00,KP00,KN00
      DOUBLE PRECISION SUMC0,SUMCB,SUMVOL,RCPM,RCP0,QW,WLE
      DOUBLE PRECISION XN,V1(3),V2(3),V3(3),AZ,AL,DL
      DOUBLE PRECISION FPMIN,FPMAX,FNMIN,FNMAX,AMAX,TMINC,TMAXC
      DOUBLE PRECISION CFRP_HGAS
      DOUBLE PRECISION F1,F2

      SAVE
      DATA FIRST /.TRUE./
      DATA DISABLED /.FALSE./
      DATA LASTTIME /-1.0D30/
      DATA TPRINT /-1.0D30/
      DATA NEXTPRT /0.0D0/
      DATA T2NEXT /0.0D0/
      DATA NCOL /0/
      DATA NSAMP /0/
      DATA NCALL /0/
      DATA ISAMP /1/
      DATA N2PERS /0/
      DATA SAMPT /1.0D0,150.0D0,1.0D30/
      DATA NK /10*0/
      DATA CNT /8*0/

C=======================================================================
C     CALL LOG: which KODE values does TMG really send?
C=======================================================================

      NCALL=NCALL+1
      IF (KODE .GE. 0 .AND. KODE .LE. 9) NK(KODE)=NK(KODE)+1
      IF (NCALL .LE. 6) THEN
         PRINT *,' USER1 CALL',NCALL,' KODE=',KODE,' TIME=',TIME
      ENDIF

C=======================================================================
C     INITIALIZATION  (first call with KODE = 1: GG, C still original)
C=======================================================================

      IF (FIRST .AND. KODE .EQ. 1) THEN

         FIRST=.FALSE.

         PRINT *,' '
         PRINT *,'=============================================='
         PRINT *,' CFRP USER1 : thermochemical decomposition'
         PRINT *,'=============================================='

         IF (MAXNO .GT. MAXELEM .OR. NOCON .GT. MAXCOND) THEN
            PRINT *,' ERROR: model too large: elements',MAXNO,
     +              ' (max',MAXELEM,') conductances',NOCON,
     +              ' (max',MAXCOND,')'
            DISABLED=.TRUE.
         ENDIF

C        ---- group ---------------------------------------------------
         IF (.NOT. DISABLED) THEN
            LNAME=GROUPNAME
            SNAME='       '
            GID=0
            GLENGTH=0
            CALL LONGNAME(SNAME,LNAME,GID,GLENGTH,2)
            GNUMELEM=0
            CALL NAMAR(SNAME,GLIST,GNUMELEM)
            PRINT *,' group "',GROUPNAME(1:LEN_TRIM(GROUPNAME)),
     +              '" : ',GNUMELEM,' elements'
            IF (GNUMELEM .LE. 0 .OR. GNUMELEM .GT. MAXELEM) THEN
               PRINT *,' ERROR: group not found or too large'
               DISABLED=.TRUE.
            ENDIF
         ENDIF

C        ---- element data and TMG baseline ----------------------------
         IF (.NOT. DISABLED) THEN
            CALL READPROP2(ELPROP,ELGEOM,NODES,NODELIST,
     +                     'Y','Y','N','N')
            DO J=1,MAXELEM
               ISCFRP(J)=.FALSE.
               ISBOUND(J)=.FALSE.
               HAVEOR(J)=.FALSE.
               NBND(J)=0
               NADJ(J)=0
               C0B(J)=0.0
               E1(1,J)=1.0D0
               E1(2,J)=0.0D0
               E1(3,J)=0.0D0
               E2(1,J)=0.0D0
               E2(2,J)=1.0D0
               E2(3,J)=0.0D0
            END DO
            DO J=1,MAXEXT
               EXTMAP(J)=0
            END DO
            SUMC0=0.0D0
            SUMCB=0.0D0
            SUMVOL=0.0D0
            NBLIST=0
            DO J=1,GNUMELEM
               N=GLIST(J)
               IF (N .LE. 0 .OR. N .GT. MAXELEM) THEN
                  PRINT *,' ERROR: internal element out of range',N
                  DISABLED=.TRUE.
               ELSE
                  ISCFRP(N)=.TRUE.
                  ALPHA1(N)=0.0D0
                  ALPHA2(N)=0.0D0
                  QSTORE(N)=0.0D0
                  QWATT(N)=0.0D0
                  C0(N)=C(N)
                  VOLE(N)=DBLE(ELPROP(1,N))
                  TREF(N)=T(N)
                  FPAR(N)=1.0D0
                  FPERP(N)=1.0D0
                  RR(N)=1.0D0
                  SUMC0=SUMC0+DBLE(C0(N))
                  SUMVOL=SUMVOL+VOLE(N)
                  IF (ICONV(N) .GE. 1 .AND. ICONV(N) .LE. MAXEXT) THEN
                     EXTMAP(ICONV(N))=N
                  ELSE
                     PRINT *,' ERROR: element label too large:',
     +                       ICONV(N),' (increase MAXEXT)'
                     DISABLED=.TRUE.
                  ENDIF
               ENDIF
            END DO
         ENDIF

C        ---- boundary elements of every hexahedron --------------------
         IF (.NOT. DISABLED) THEN
            DO J=1,GNUMELEM
               N=GLIST(J)
               NCTMP=0
               NB_TMP=0
               CALL CONDNOS(N,ICNDTMP,NCTMP,BOUNDARYLIST,NB_TMP)
               IF (NB_TMP .GT. 8) THEN
                  PRINT *,' ERROR: element',ICONV(N),' has',NB_TMP,
     +                    ' boundary elements (max 8)'
                  DISABLED=.TRUE.
                  NB_TMP=8
               ENDIF
               NBND(N)=NB_TMP
               DO L=1,NB_TMP
                  B=BOUNDARYLIST(L)
                  BND(L,N)=B
                  IF (B .GT. 0 .AND. B .LE. MAXELEM) THEN
                     IF (.NOT. ISBOUND(B)) THEN
                        ISBOUND(B)=.TRUE.
                        NBLIST=NBLIST+1
                        BLIST(NBLIST)=B
                        C0B(B)=C(B)
                        SUMCB=SUMCB+DBLE(C0B(B))
                     ENDIF
                     IF (NADJ(B) .LT. 4) THEN
                        NADJ(B)=NADJ(B)+1
                        ADJN(NADJ(B),B)=N
                     ENDIF
                  ENDIF
               END DO
            END DO
            DO J=1,NOCON
               GG0(J)=GG(J)
               UPDATED(J)=.FALSE.
            END DO
            NGRCOND=0
            CALL CONDNOS2(SNAME,ICNDGR,NGRCOND)
            IF (NGRCOND .LE. 0 .OR. NGRCOND .GT. MAXCOND) THEN
               PRINT *,' ERROR: CONDNOS2 returned',NGRCOND,
     +                 ' conductances (max',MAXCOND,')'
               DISABLED=.TRUE.
            ENDIF
         ENDIF

C        ---- NX material and orientation (solver input XML) ----------
         IF (.NOT. DISABLED) THEN
            IU=31
   10       OPN=.FALSE.
            INQUIRE(UNIT=IU,OPENED=OPN)
            IF (OPN .AND. IU .LT. 99) THEN
               IU=IU+1
               GOTO 10
            ENDIF
            CALL CFRP_OPENXML(XMLCFG,IU,XMLNAME,XOK)
            IF (.NOT. XOK) THEN
               PRINT *,' ERROR: solver input XML not found. Put the',
     +                 ' NX .xml in the run folder or set XMLCFG.'
               DISABLED=.TRUE.
            ELSE
               CALL CFRP_READXML(IU,EXTMAP,MAXEXT,E1,E2,HAVEOR,RHO0,
     +              CP0,BK,MATUID,NDIFF,NNOOR,GNUMELEM,GLIST,XSTAT)
               CLOSE(IU)
               PRINT *,' NX input: ',XMLNAME(1:LEN_TRIM(XMLNAME))
               IF (XSTAT .EQ. 1 .OR. MATUID .LE. 0) THEN
                  PRINT *,' ERROR: material of the group not found'
                  DISABLED=.TRUE.
               ELSE IF (XSTAT .EQ. 2) THEN
                  PRINT *,' ERROR: the NX material properties must be',
     +                    ' CONSTANT (not tables) for this routine'
                  DISABLED=.TRUE.
               ENDIF
               IF (NDIFF .GT. 0) THEN
                  PRINT *,' WARNING: group uses more than one',
     +                    ' material; first one used'
               ENDIF
            ENDIF
         ENDIF

         IF (.NOT. DISABLED) THEN

C           through-thickness axis = material axis with the lowest k
            TAX=3
            XN=BK(3)
            IF (BK(2) .LT. XN*(1.0D0-1.0D-6)) THEN
               TAX=2
               XN=BK(2)
            ENDIF
            IF (BK(1) .LT. XN*(1.0D0-1.0D-6)) THEN
               TAX=1
               XN=BK(1)
            ENDIF

C           material axes of every element: e1, e2, e3=e1 x e2
            DO J=1,GNUMELEM
               N=GLIST(J)
               XN=DSQRT(E1(1,N)**2+E1(2,N)**2+E1(3,N)**2)
               IF (XN .LT. 1.0D-12) THEN
                  E1(1,N)=1.0D0
                  E1(2,N)=0.0D0
                  E1(3,N)=0.0D0
                  XN=1.0D0
               ENDIF
               DO K=1,3
                  V1(K)=E1(K,N)/XN
               END DO
               XN=DSQRT(E2(1,N)**2+E2(2,N)**2+E2(3,N)**2)
               IF (XN .LT. 1.0D-12) THEN
                  E2(1,N)=0.0D0
                  E2(2,N)=1.0D0
                  E2(3,N)=0.0D0
                  XN=1.0D0
               ENDIF
               DO K=1,3
                  V2(K)=E2(K,N)/XN
               END DO
               V3(1)=V1(2)*V2(3)-V1(3)*V2(2)
               V3(2)=V1(3)*V2(1)-V1(1)*V2(3)
               V3(3)=V1(1)*V2(2)-V1(2)*V2(1)
               XN=DSQRT(V3(1)**2+V3(2)**2+V3(3)**2)
               IF (XN .LT. 1.0D-12) THEN
                  V3(1)=0.0D0
                  V3(2)=0.0D0
                  V3(3)=1.0D0
                  XN=1.0D0
               ENDIF
               DO K=1,3
                  V3(K)=V3(K)/XN
                  IF (TAX .EQ. 1) ETA(K,N)=V1(K)
                  IF (TAX .EQ. 2) ETA(K,N)=V2(K)
                  IF (TAX .EQ. 3) ETA(K,N)=V3(K)
               END DO
            END DO

C           ---- model laws and the NX base --------------------------
            CALL CFRP_MODEL(KIN,HDEC,DENSR,CPV,CPE,KPV,KPE,KNV,KNE,
     +                      CPG)
            N=GLIST(1)
            CALL CFRP_LAWS(TREF(N),0.0D0,CPV,CPE,KPV,KPE,KNV,KNE,
     +                     CPM,KPM,KNM)
            PRINT *,' base material (NX, uid',MATUID,'): rho =',
     +              REAL(RHO0),' cp =',REAL(CP0)
            PRINT *,'   k1,k2,k3 =',REAL(BK(1)),REAL(BK(2)),
     +              REAL(BK(3)),' -> through-thickness axis =',TAX
            PRINT *,'   model law at Tref=',REAL(TREF(N)),': cp =',
     +              REAL(CPM),' k(in-plane) =',REAL(KPM),
     +              ' k(thickness) =',REAL(KNM)
            PRINT *,'   orientation vectors read for',
     +              GNUMELEM-NNOOR,' of',GNUMELEM,' elements',
     +              ' (others: global X,Y)'
            RCPM=0.0D0
            IF (SUMVOL .GT. 0.0D0) RCPM=(SUMC0+SUMCB)/SUMVOL
            PRINT *,'   TMG capacitance / volume =',REAL(RCPM),
     +              '  vs rho*cp(NX) =',REAL(RHO0*CP0)
            IF (SUMC0 .LT. 0.2D0*(SUMC0+SUMCB)) THEN
               PRINT *,'   (capacitance sits on the boundary elements)'
            ENDIF

C           ---- link index, thickness faces, link classes -----------
            CALL CFRP_LINKINDEX(NGRCOND,ICNDGR,ISCFRP,GLIST,GNUMELEM,
     +                          LNUM,LFIRST,LIDX)
            CALL CFRP_TOPO(GLIST,GNUMELEM,ISCFRP,NBND,BND,NADJ,ADJN,
     +                     ELGEOM,ELPROP,ETA,VOLE,LFIRST,LNUM,LIDX,
     +                     ICNDGR,GG0,MAXELEM,THKCOS,ZRATIO,ZB1,ZB2,
     +                     NBP,NBM,PITCH,METH,LCLS,NMETH)
            NHI=0
            NLO=0
            DO J=1,NGRCOND
               IF (GG0(ICNDGR(1,J)) .LT. 0.0) NLO=NLO+1
               IF (LCLS(J) .EQ. 1) NHI=NHI+1
            END DO
            PRINT *,' conductance map =',NGRCOND,' links (',
     +              REAL(NGRCOND)/REAL(GNUMELEM),' per element),',
     +              NLO,' with negative resistance'
            PRINT *,' thickness faces found: topology',NMETH(1),
     +              ' elimination',NMETH(2),' conductance pattern',
     +              NMETH(3),' UNRESOLVED',NMETH(0)
            PRINT *,' thickness links =',NHI
            IF (NMETH(0) .GT. 0) THEN
               PRINT *,' WARNING: unresolved elements keep their',
     +                 ' original conductances'
            ENDIF

C           heat-load split between thickness faces and side faces
            DO J=1,GNUMELEM
               N=GLIST(J)
               WZ(N)=0.5D0
               WL(N)=0.0D0
               IF (METH(N) .GT. 0 .AND. PITCH(N) .GT. 0.0D0) THEN
                  AZ=VOLE(N)/PITCH(N)
                  DL=DSQRT(VOLE(N)/PITCH(N))
                  IF (DL .GT. 0.0D0) THEN
                     AL=VOLE(N)/DL
                     WZ(N)=AZ/(2.0D0*AZ+4.0D0*AL)
                     WL(N)=AL/(2.0D0*AZ+4.0D0*AL)
                  ENDIF
               ENDIF
            END DO

C           thickness column for the temperature table and link audit
            CALL CFRP_COLUMN(GLIST(1),NBP,NBM,ETA,ELGEOM,MAXCOL,NCOL,
     +                       COLN,COLS)
            PRINT *,' thickness column: ',NCOL,' layers'
            NSAMP=0
            N=GLIST(1)
            DO L=1,LNUM(N)
               IF (NSAMP .LT. MAXSAMP) THEN
                  NSAMP=NSAMP+1
                  SAMPJ(NSAMP)=LIDX(LFIRST(N)+L-1)
               ENDIF
            END DO
            PRINT *,'=============================================='
         ENDIF

      ENDIF

      IF (DISABLED) RETURN
      IF (FIRST) RETURN

C=======================================================================
C     KODE = 2 : after the temperatures are calculated.
C     Did the GG values written at KODE = 1 survive the solve?
C=======================================================================

      IF (KODE .EQ. 2) THEN
         IF (NSAMP .GT. 0 .AND. TIME .GT. 1.0D0 .AND.
     +       TIME .GE. T2NEXT .AND. N2PERS .LT. 3) THEN
            NEQ=0
            NREV=0
            NOTH=0
            DO K=1,NSAMP
               IC=ICNDGR(1,SAMPJ(K))
               IF (ABS(GG(IC)-GGEXP(K)) .LE.
     +             1.0E-5*ABS(GGEXP(K))) THEN
                  NEQ=NEQ+1
               ELSE IF (ABS(GG(IC)-GG0(IC)) .LE.
     +                  1.0E-5*ABS(GG0(IC))) THEN
                  NREV=NREV+1
               ELSE
                  NOTH=NOTH+1
               ENDIF
            END DO
            PRINT *,' GG persistence at KODE=2, t=',TIME,
     +              ': kept',NEQ,' reverted',NREV,' other',NOTH,
     +              ' of',NSAMP
            T2NEXT=TIME+DIAGINT
            N2PERS=N2PERS+1
         ENDIF
         RETURN
      ENDIF

C=======================================================================
C     KODE = 3 : end of run
C=======================================================================

      IF (KODE .EQ. 3) THEN
         PRINT *,' CFRP FINAL SUMMARY t=',TIME,' calls per KODE 0..9:',
     +           (NK(K),K=0,9)
         CALL CFRP_SAMPLES(TIME,NSAMP,SAMPJ,ICNDGR,LCLS,GG,GG0,ICONV)
         IF (NCOL .GT. 1) CALL CFRP_PROFILE(TIME,T,NCOL,COLN,COLS,
     +                                      PROBE,ICONV)
         RETURN
      ENDIF

      IF (KODE .NE. 1) RETURN

C=======================================================================
C     KODE = 1 : before the temperatures are calculated
C=======================================================================

      NEWTIME=.FALSE.
      IF (TIME .GT. LASTTIME+1.0D-12) NEWTIME=.TRUE.

      FPMIN=1.0D30
      FPMAX=-1.0D30
      FNMIN=1.0D30
      FNMAX=-1.0D30
      AMAX=0.0D0
      TMINC=1.0D30
      TMAXC=-1.0D30
      F1=KIN(4)
      F2=KIN(8)
      DRHO=RHO0*(1.0D0-DENSR)

      DO J=1,GNUMELEM
         N=GLIST(J)
         TC=T(N)
         TK=TC+273.15D0
         IF (TK .LT. 1.0D0) TK=1.0D0

C        ---- chemistry: advance once per new time value ------------
         IF (NEWTIME .AND. TIME .GT. 1.0D-12) THEN
            A1OLD=ALPHA1(N)
            A2OLD=ALPHA2(N)
            IF (LASTTIME .LT. -1.0D20) THEN
               STEPDT=DT
            ELSE
               STEPDT=TIME-LASTTIME
            ENDIF
            IF (STEPDT .LE. 0.0D0) STEPDT=DT
            IF (STEPDT .GT. 0.0D0) THEN
               CALL CFRP_RK4(KIN,A1OLD,A2OLD,TK,STEPDT,A1NEW,A2NEW)
               R1=(A1NEW-A1OLD)/STEPDT
               R2=(A2NEW-A2OLD)/STEPDT
            ELSE
               A1NEW=A1OLD
               A2NEW=A2OLD
               R1=0.0D0
               R2=0.0D0
            ENDIF
            ALPHA1(N)=A1NEW
            ALPHA2(N)=A2NEW
C           decomposition heat [W/m3] minus the enthalpy of the gas
            QCHEM=DRHO*(F1*HDEC(1)*R1+F2*HDEC(2)*R2)
            RHOGDOT=DRHO*(F1*R1+F2*R2)
            DH=CFRP_HGAS(CPG,TC)-CFRP_HGAS(CPG,TREF(N))
            IF (DH .LT. 0.0D0) DH=0.0D0
            QSTORE(N)=QCHEM-RHOGDOT*DH
         ENDIF

C        ---- degree of decomposition and property factors ----------
         ALPHA=F1*ALPHA1(N)+F2*ALPHA2(N)
         IF (ALPHA .LT. 0.0D0) ALPHA=0.0D0
         IF (ALPHA .GT. 1.0D0) ALPHA=1.0D0

         CALL CFRP_LAWS(TC,ALPHA,CPV,CPE,KPV,KPE,KNV,KNE,CPM,KPM,KNM)
         CALL CFRP_LAWS(TREF(N),0.0D0,CPV,CPE,KPV,KPE,KNV,KNE,
     +                  CP00,KP00,KN00)
         IF (KP00 .GT. 0.0D0 .AND. KPM .GT. 0.0D0) THEN
            FPAR(N)=KPM/KP00
         ELSE
            FPAR(N)=1.0D0
         ENDIF
         IF (KN00 .GT. 0.0D0 .AND. KNM .GT. 0.0D0) THEN
            FPERP(N)=KNM/KN00
         ELSE
            FPERP(N)=1.0D0
         ENDIF
         IF (CP00 .GT. 0.0D0) THEN
            RR(N)=(1.0D0-(1.0D0-DENSR)*ALPHA)*CPM/CP00
         ELSE
            RR(N)=1.0D0
         ENDIF

C        ---- capacitance of the CG element ---------------------------
         C(N)=REAL(DBLE(C0(N))*RR(N))

C        ---- decomposition heat of this element [W] -----------------
         QWATT(N)=QSTORE(N)*VOLE(N)

         IF (FPAR(N) .LT. FPMIN) FPMIN=FPAR(N)
         IF (FPAR(N) .GT. FPMAX) FPMAX=FPAR(N)
         IF (FPERP(N) .LT. FNMIN) FNMIN=FPERP(N)
         IF (FPERP(N) .GT. FNMAX) FNMAX=FPERP(N)
         IF (ALPHA .GT. AMAX) AMAX=ALPHA
         IF (TC .LT. TMINC) TMINC=TC
         IF (TC .GT. TMAXC) TMAXC=TC
      END DO

C=======================================================================
C     CAPACITANCE ON THE BOUNDARY ELEMENTS (only if TMG keeps it there)
C=======================================================================

      DO J=1,NBLIST
         B=BLIST(J)
         SUMV(B)=0.0D0
         SUMVR(B)=0.0D0
      END DO
      DO J=1,GNUMELEM
         N=GLIST(J)
         DO K=1,NBND(N)
            B=BND(K,N)
            IF (B .GT. 0 .AND. B .LE. MAXELEM) THEN
               SUMV(B)=SUMV(B)+VOLE(N)
               SUMVR(B)=SUMVR(B)+VOLE(N)*RR(N)
            ENDIF
         END DO
      END DO
      DO J=1,NBLIST
         B=BLIST(J)
         IF (SUMV(B) .GT. 0.0D0 .AND. C0B(B) .NE. 0.0) THEN
            C(B)=REAL(DBLE(C0B(B))*SUMVR(B)/SUMV(B))
         ENDIF
      END DO

C=======================================================================
C     DECOMPOSITION HEAT LOAD: added (Q is the total load) to the two
C     thickness faces and the side faces of each hexahedron
C=======================================================================

      DO J=1,GNUMELEM
         N=GLIST(J)
         QW=QWATT(N)
         IF (QW .NE. 0.0D0) THEN
            IF (METH(N) .GT. 0 .AND. NBND(N) .GT. 2) THEN
               Q(ZB1(N))=Q(ZB1(N))+REAL(QW*WZ(N))
               Q(ZB2(N))=Q(ZB2(N))+REAL(QW*WZ(N))
               WLE=(1.0D0-2.0D0*WZ(N))/DBLE(NBND(N)-2)
               DO K=1,NBND(N)
                  B=BND(K,N)
                  IF (B .NE. ZB1(N) .AND. B .NE. ZB2(N)) THEN
                     Q(B)=Q(B)+REAL(QW*WLE)
                  ENDIF
               END DO
            ELSE
               Q(N)=Q(N)+REAL(QW)
            ENDIF
         ENDIF
      END DO

C=======================================================================
C     CONDUCTANCES (GG is a resistance)
C=======================================================================

      DO J=1,NOCON
         UPDATED(J)=.FALSE.
      END DO
      DO J=1,8
         CNT(J)=0
      END DO
      NMOD=0
      CALL CFRP_UPDATE_GROUP(NGRCOND,ICNDGR,GG,GG0,UPDATED,NOCON,
     +     ISCFRP,FPAR,FPERP,LCLS,NMOD,CNT)
      DO K=1,NSAMP
         GGEXP(K)=GG(ICNDGR(1,SAMPJ(K)))
      END DO

C=======================================================================
C     DIAGNOSTICS (bounded)
C=======================================================================

      IF (NEWTIME) THEN
         IF (TIME .GE. TPRINT+DIAGINT-1.0D-9) THEN
            PRINT *,' CFRP t=',REAL(TIME),' alpha max=',REAL(AMAX),
     +              ' T(CG) min/max=',REAL(TMINC),REAL(TMAXC)
            PRINT *,'      factor in-plane',REAL(FPMIN),REAL(FPMAX),
     +              ' thickness',REAL(FNMIN),REAL(FNMAX),
     +              ' | links updated',NMOD,
     +              ' (thickness',CNT(1),', negative',CNT(3),
     +              ') skipped',CNT(4)+CNT(5)+CNT(6)
            TPRINT=TIME
         ENDIF
         IF (ISAMP .LE. 2) THEN
            IF (TIME .GE. SAMPT(ISAMP)-1.0D-9) THEN
               CALL CFRP_SAMPLES(TIME,NSAMP,SAMPJ,ICNDGR,LCLS,GG,GG0,
     +                           ICONV)
               ISAMP=ISAMP+1
            ENDIF
         ENDIF
         IF (NCOL .GT. 1 .AND. TIME .GE. NEXTPRT-1.0D-9) THEN
            CALL CFRP_PROFILE(TIME,T,NCOL,COLN,COLS,PROBE,ICONV)
            NEXTPRT=NEXTPRT+PRTINT
         ENDIF
         LASTTIME=TIME
      ENDIF

      RETURN
      END

C=======================================================================
C     MATERIAL MODEL  -  the ONLY place with decomposition-model data
C
C     Edit this routine to use another material / another kinetic
C     model. Nothing else in the file depends on these numbers.
C     (Values below: T700/M21, Tranchard et al., Materials 2017.)
C     Base properties (rho, cp, k1,k2,k3) are NOT here: they are read
C     from the NX material (solver input XML) and from TMG itself.
C
C     Polynomial laws:  P(T) = c(1) + c(2)*T + c(3)*T**2 + c(4)*T**3
C     with T in degC and results in SI units.
C=======================================================================

      SUBROUTINE CFRP_MODEL(KIN,HDEC,DENSR,CPV,CPE,KPV,KPE,KNV,KNE,
     +                      CPG)

      IMPLICIT NONE

C     KIN(1..4)  step 1: A [1/s], Ea [J/mol], order n, fraction
C     KIN(5..8)  step 2: A [1/s], Ea [J/mol], order n, fraction
C     KIN(9)     autocatalysis constant of step 2
C     HDEC(1..2) heat of decomposition per kg lost [J/kg], + = exothermic
C     DENSR      final density / virgin density
C     CPV,CPE    specific heat  virgin / degraded  [J/(kg K)]
C     KPV,KPE    in-plane conductivity  virgin / degraded  [W/(m K)]
C     KNV,KNE    through-thickness conductivity virgin / degraded
C     CPG        specific heat of the released gas [J/(kg K)]
      DOUBLE PRECISION KIN(9),HDEC(2),DENSR
      DOUBLE PRECISION CPV(4),CPE(4),KPV(4),KPE(4),KNV(4),KNE(4)
      DOUBLE PRECISION CPG(4)
      INTEGER I

      KIN(1)=7.219381629371594D0
      KIN(2)=5.87107D4
      KIN(3)=1.1125D0
      KIN(4)=0.5540D0
      KIN(5)=5.330893422367024D8
      KIN(6)=1.46674D5
      KIN(7)=2.0825D0
      KIN(8)=0.9479D0
      KIN(9)=7.488589075216437D0

      HDEC(1)=2.5954D5
      HDEC(2)=-1.5222D5

      DENSR=1165.0D0/1575.0D0

      DO I=1,4
         CPV(I)=0.0D0
         CPE(I)=0.0D0
         KPV(I)=0.0D0
         KPE(I)=0.0D0
         KNV(I)=0.0D0
         KNE(I)=0.0D0
         CPG(I)=0.0D0
      END DO

      CPV(1)=687.31D0
      CPV(2)=2.8773D0

      CPE(1)=662.53D0
      CPE(2)=2.599D0
      CPE(3)=-2.0761D-3
      CPE(4)=5.1327D-7

      KPV(1)=2.7811D0
      KPV(2)=7.4675D-3

      KPE(1)=2.3005D0
      KPE(2)=6.4898D-3
      KPE(3)=-3.5481D-6

      KNV(1)=0.61391D0
      KNV(2)=1.1113D-3

      KNE(1)=7.7211D-2
      KNE(2)=2.3397D-4
      KNE(3)=-4.1903D-7
      KNE(4)=7.4228D-10

      CPG(1)=1256.6D0
      CPG(2)=1.0610D0
      CPG(3)=-9.2485D-4
      CPG(4)=3.5977D-7

      RETURN
      END


C=======================================================================
C     polynomial and law helpers
C=======================================================================

      DOUBLE PRECISION FUNCTION CFRP_POLY(C,T)
      IMPLICIT NONE
      DOUBLE PRECISION C(4),T
      CFRP_POLY=C(1)+T*(C(2)+T*(C(3)+T*C(4)))
      RETURN
      END

C     mixture laws at temperature TC [degC] and decomposition degree AL
      SUBROUTINE CFRP_LAWS(TC,AL,CPV,CPE,KPV,KPE,KNV,KNE,CPM,KPM,KNM)
      IMPLICIT NONE
      DOUBLE PRECISION TC,AL,CPM,KPM,KNM
      DOUBLE PRECISION CPV(4),CPE(4),KPV(4),KPE(4),KNV(4),KNE(4)
      DOUBLE PRECISION CFRP_POLY
      CPM=(1.0D0-AL)*CFRP_POLY(CPV,TC)+AL*CFRP_POLY(CPE,TC)
      KPM=(1.0D0-AL)*CFRP_POLY(KPV,TC)+AL*CFRP_POLY(KPE,TC)
      KNM=(1.0D0-AL)*CFRP_POLY(KNV,TC)+AL*CFRP_POLY(KNE,TC)
      RETURN
      END

C     enthalpy of the gas, integral of CPG dT from 0 degC to TC
      DOUBLE PRECISION FUNCTION CFRP_HGAS(CPG,TC)
      IMPLICIT NONE
      DOUBLE PRECISION CPG(4),TC
      CFRP_HGAS=TC*(CPG(1)+TC*(CPG(2)/2.0D0+TC*(CPG(3)/3.0D0
     +          +TC*CPG(4)/4.0D0)))
      RETURN
      END


C=======================================================================
C     RK4 STEP OF THE TWO-REACTION KINETICS (T held constant)
C     adaptive sub-stepping, same scheme as the validated version
C=======================================================================

      SUBROUTINE CFRP_RK4(KIN,A1,A2,TK,DT,AN1,AN2)

      IMPLICIT NONE
      DOUBLE PRECISION KIN(9)
      DOUBLE PRECISION A1,A2,TK,DT,AN1,AN2
      DOUBLE PRECISION K1A,K1B,K2A,K2B
      DOUBLE PRECISION K3A,K3B,K4A,K4B
      DOUBLE PRECISION X1,X2,H,TLEFT,RTOT,F1R,F2R
      DOUBLE PRECISION SCALE
      INTEGER NSUB
      INTEGER MAXSUB
      PARAMETER (MAXSUB=50000)

      AN1=MAX(0.0D0,MIN(1.0D0,A1))
      AN2=MAX(0.0D0,MIN(1.0D0,A2))
      IF (AN1+AN2 .GT. 1.0D0) THEN
         SCALE=1.0D0/(AN1+AN2)
         AN1=AN1*SCALE
         AN2=AN2*SCALE
      ENDIF
      IF (DT .LE. 0.0D0) RETURN

      TLEFT=DT
      NSUB=0
   10 CONTINUE
      IF (TLEFT .LE. MAX(1.0D-14,DT*1.0D-12)) GOTO 90
      CALL CFRP_RATE(KIN,AN1,AN2,TK,F1R,F2R)
      RTOT=F1R+F2R
      IF (RTOT .GT. 1.0D-20) THEN
         H=MIN(TLEFT,0.05D0/RTOT)
      ELSE
         H=TLEFT
      ENDIF
      IF (H .LE. 0.0D0) GOTO 90
      IF (NSUB .GE. MAXSUB) THEN
         IF (RTOT .GT. 1.0D-20) THEN
            SCALE=F1R/RTOT
            X1=1.0D0-AN1-AN2
            IF (X1 .LT. 0.0D0) X1=0.0D0
            AN1=AN1+X1*SCALE
            AN2=AN2+X1*(1.0D0-SCALE)
         ENDIF
         GOTO 90
      ENDIF
      CALL CFRP_RATE(KIN,AN1,AN2,TK,K1A,K1B)
      X1=AN1+0.5D0*H*K1A
      X2=AN2+0.5D0*H*K1B
      CALL CFRP_RATE(KIN,X1,X2,TK,K2A,K2B)
      X1=AN1+0.5D0*H*K2A
      X2=AN2+0.5D0*H*K2B
      CALL CFRP_RATE(KIN,X1,X2,TK,K3A,K3B)
      X1=AN1+H*K3A
      X2=AN2+H*K3B
      CALL CFRP_RATE(KIN,X1,X2,TK,K4A,K4B)
      X1=AN1+(H/6.0D0)*
     +   (K1A+2.0D0*K2A+2.0D0*K3A+K4A)
      X2=AN2+(H/6.0D0)*
     +   (K1B+2.0D0*K2B+2.0D0*K3B+K4B)
      AN1=MAX(0.0D0,MIN(1.0D0,X1))
      AN2=MAX(0.0D0,MIN(1.0D0,X2))
      IF (AN1+AN2 .GT. 1.0D0) THEN
         SCALE=1.0D0/(AN1+AN2)
         AN1=AN1*SCALE
         AN2=AN2*SCALE
      ENDIF
      TLEFT=TLEFT-H
      NSUB=NSUB+1
      GOTO 10
   90 CONTINUE
      RETURN
      END


C=======================================================================
C     REACTION RATES (step 1: n-th order, step 2: n-th order with
C     autocatalysis)
C=======================================================================

      SUBROUTINE CFRP_RATE(KIN,A1,A2,TK,F1R,F2R)

      IMPLICIT NONE
      DOUBLE PRECISION KIN(9)
      DOUBLE PRECISION A1,A2,TK,F1R,F2R
      DOUBLE PRECISION GASR
      DOUBLE PRECISION REM,EX1,EX2
      PARAMETER (GASR=8.31446261815324D0)

      REM=1.0D0-A1-A2
      IF (REM .LT. 0.0D0) REM=0.0D0
      IF (REM .GT. 1.0D0) REM=1.0D0
      IF (TK .LE. 1.0D0) THEN
         F1R=0.0D0
         F2R=0.0D0
         RETURN
      ENDIF
      EX1=DEXP(-KIN(2)/(GASR*TK))
      F1R=KIN(1)*REM**KIN(3)*EX1
      EX2=DEXP(-KIN(6)/(GASR*TK))
      F2R=KIN(5)*REM**KIN(7)*EX2*(1.0D0+KIN(9)*A2)
      IF (F1R .LT. 0.0D0) F1R=0.0D0
      IF (F2R .LT. 0.0D0) F2R=0.0D0
      RETURN
      END


C=======================================================================
C     LOCATE AND OPEN THE SOLVER INPUT XML
C     1) name derived from the solver log (units 17/18 are
C        <name>_verbose.log / <name>_report.log, XML = <name>.xml)
C     2) the file name given in XMLCFG
C=======================================================================

      SUBROUTINE CFRP_OPENXML(XMLCFG,IU,FNAME,OK)

      IMPLICIT NONE
      CHARACTER*(*) XMLCFG
      CHARACTER*(*) FNAME
      INTEGER IU
      LOGICAL OK

      CHARACTER*300 NM
      CHARACTER*300 CAND
      LOGICAL OPN,EX
      INTEGER U,P,ISTAT

      OK=.FALSE.
      FNAME=' '

      DO U=17,18
         NM=' '
         OPN=.FALSE.
         INQUIRE(UNIT=U,OPENED=OPN,IOSTAT=ISTAT)
         IF (OPN .AND. ISTAT .EQ. 0) THEN
            INQUIRE(UNIT=U,NAME=NM,IOSTAT=ISTAT)
            IF (ISTAT .EQ. 0) THEN
               P=INDEX(NM,'_verbose.log')
               IF (P .EQ. 0) P=INDEX(NM,'_report.log')
               IF (P .GT. 1 .AND. .NOT. OK) THEN
                  CAND=NM(1:P-1)//'.xml'
                  EX=.FALSE.
                  INQUIRE(FILE=CAND,EXIST=EX)
                  IF (EX) THEN
                     OPEN(UNIT=IU,FILE=CAND,STATUS='OLD',
     +                    IOSTAT=ISTAT)
                     IF (ISTAT .EQ. 0) THEN
                        FNAME=CAND
                        OK=.TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      END DO

      IF (.NOT. OK .AND. XMLCFG .NE. ' ') THEN
         EX=.FALSE.
         INQUIRE(FILE=XMLCFG,EXIST=EX)
         IF (EX) THEN
            OPEN(UNIT=IU,FILE=XMLCFG,STATUS='OLD',IOSTAT=ISTAT)
            IF (ISTAT .EQ. 0) THEN
               FNAME=XMLCFG
               OK=.TRUE.
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


C=======================================================================
C     small text helpers
C=======================================================================

C     integer starting exactly at position START of STR (digits only)
      SUBROUTINE CFRP_GETINT(STR,START,VAL,OK)
      IMPLICIT NONE
      CHARACTER*(*) STR
      INTEGER START,VAL
      LOGICAL OK
      INTEGER I,D
      VAL=0
      OK=.FALSE.
      I=START
   10 IF (I .LE. LEN(STR)) THEN
         D=ICHAR(STR(I:I))-ICHAR('0')
         IF (D .GE. 0 .AND. D .LE. 9) THEN
            VAL=VAL*10+D
            OK=.TRUE.
            I=I+1
            GOTO 10
         ENDIF
      ENDIF
      RETURN
      END

C     first integer after "<E>"  (element label of a <E> line)
      SUBROUTINE CFRP_EID(STR,ID)
      IMPLICIT NONE
      CHARACTER*(*) STR
      INTEGER ID,P
      LOGICAL OK
      ID=0
      P=INDEX(STR,'<E>')
      IF (P .EQ. 0) RETURN
      P=P+3
   10 IF (P .LE. LEN(STR)) THEN
         IF (STR(P:P) .EQ. ' ') THEN
            P=P+1
            GOTO 10
         ENDIF
      ENDIF
      CALL CFRP_GETINT(STR,P,ID,OK)
      IF (.NOT. OK) ID=0
      RETURN
      END


C=======================================================================
C     READ THE NX MODEL (solver input XML)
C
C     - material of the group's elements: density, specific heat and
C       the three conductivities (must be "Constant" in the NX UI)
C     - orientation of every group element: first and second material
C       direction vectors (the third one is their cross product)
C
C     EXTMAP(external label) = internal element number (0 = not in group)
C     STAT: 0 ok, 1 no material found, 2 base property is not constant
C=======================================================================

      SUBROUTINE CFRP_READXML(IU,EXTMAP,MAXEXT,E1,E2,HAVEOR,
     +     RHO0,CP0,BK,MATUID,NDIFF,NNOOR,GNUM,GLIST,STAT)

      IMPLICIT NONE

      INTEGER IU,MAXEXT,MATUID,NDIFF,NNOOR,GNUM,STAT
      INTEGER EXTMAP(*)
      INTEGER GLIST(*)
      DOUBLE PRECISION E1(3,*),E2(3,*)
      LOGICAL HAVEOR(*)
      DOUBLE PRECISION RHO0,CP0,BK(3)

      INTEGER MAXMAT
      PARAMETER (MAXMAT=40)
      INTEGER MUID(MAXMAT)
      DOUBLE PRECISION MRHO(MAXMAT),MCP(MAXMAT),MK(3,MAXMAT)
      LOGICAL MNC(MAXMAT)
      CHARACTER*600 LINE
      LOGICAL INMAT,INEL,INOR,OK
      INTEGER NMAT,CURM,CURP,ELMAT,ID,N,P,P1,P2,K,I,UIDV
      DOUBLE PRECISION X,V(6)

      STAT=0
      MATUID=0
      NDIFF=0
      NMAT=0
      CURM=0
      CURP=0
      ELMAT=-1
      INMAT=.FALSE.
      INEL=.FALSE.
      INOR=.FALSE.

   10 CONTINUE
      LINE=' '
      READ(IU,'(A)',END=800,ERR=800) LINE
      P=INDEX(LINE,CHAR(13))
      IF (P .GT. 0) LINE(P:P)=' '

      IF (INDEX(LINE,'<MaterialList>') .GT. 0) INMAT=.TRUE.
      IF (INDEX(LINE,'</MaterialList>') .GT. 0) INMAT=.FALSE.
      IF (INDEX(LINE,'<ElementList>') .GT. 0) INEL=.TRUE.
      IF (INDEX(LINE,'</ElementList>') .GT. 0) INEL=.FALSE.
      IF (INDEX(LINE,'<MaterialOrientationVector3D>') .GT. 0)
     +   INOR=.TRUE.
      IF (INDEX(LINE,'</MaterialOrientationVector3D>') .GT. 0)
     +   INOR=.FALSE.
      IF (INDEX(LINE,'</ElementAssociatedDataList>') .GT. 0) GOTO 800

C     ---- materials -------------------------------------------------
      IF (INMAT) THEN
         P=INDEX(LINE,'<Material uid="')
         IF (P .GT. 0) THEN
            CALL CFRP_GETINT(LINE,P+15,UIDV,OK)
            IF (OK .AND. NMAT .LT. MAXMAT) THEN
               NMAT=NMAT+1
               CURM=NMAT
               MUID(NMAT)=UIDV
               MRHO(NMAT)=0.0D0
               MCP(NMAT)=0.0D0
               MNC(NMAT)=.FALSE.
               DO K=1,3
                  MK(K,NMAT)=0.0D0
               END DO
            ELSE
               CURM=0
            ENDIF
            CURP=0
            GOTO 10
         ENDIF
         IF (INDEX(LINE,'<Property name="') .GT. 0) THEN
            CURP=0
            IF (INDEX(LINE,'name="Mass Density"') .GT. 0) CURP=1
            IF (INDEX(LINE,'name="Specific Heat"') .GT. 0) CURP=2
            IF (INDEX(LINE,'name="Thermal Conductivity X"') .GT. 0)
     +         CURP=3
            IF (INDEX(LINE,'name="Thermal Conductivity Y"') .GT. 0)
     +         CURP=4
            IF (INDEX(LINE,'name="Thermal Conductivity Z"') .GT. 0)
     +         CURP=5
            IF (INDEX(LINE,'name="Thermal Conductivity"') .GT. 0)
     +         CURP=6
            IF (CURP .GT. 0 .AND. CURM .GT. 0) THEN
               IF (INDEX(LINE,'type="Constant"') .EQ. 0)
     +            MNC(CURM)=.TRUE.
            ENDIF
            GOTO 10
         ENDIF
         IF (INDEX(LINE,'<Value>') .GT. 0 .AND. CURP .GT. 0 .AND.
     +       CURM .GT. 0) THEN
            P1=INDEX(LINE,'<Value>')+7
            P2=INDEX(LINE,'</Value>')-1
            IF (P2 .GE. P1) THEN
               READ(LINE(P1:P2),*,ERR=20) X
               IF (CURP .EQ. 1) MRHO(CURM)=X
               IF (CURP .EQ. 2) MCP(CURM)=X
               IF (CURP .EQ. 3) MK(1,CURM)=X
               IF (CURP .EQ. 4) MK(2,CURM)=X
               IF (CURP .EQ. 5) MK(3,CURM)=X
               IF (CURP .EQ. 6) THEN
                  MK(1,CURM)=X
                  MK(2,CURM)=X
                  MK(3,CURM)=X
               ENDIF
   20          CONTINUE
            ENDIF
            CURP=0
            GOTO 10
         ENDIF
      ENDIF

C     ---- element list: which material do the group elements use? ---
      IF (INEL) THEN
         IF (INDEX(LINE,'<Set ') .GT. 0) THEN
            ELMAT=-1
            P=INDEX(LINE,' material="')
            IF (P .GT. 0) THEN
               CALL CFRP_GETINT(LINE,P+11,ELMAT,OK)
               IF (.NOT. OK) ELMAT=-1
            ENDIF
         ELSE IF (INDEX(LINE,'<E>') .GT. 0) THEN
            CALL CFRP_EID(LINE,ID)
            N=0
            IF (ID .GE. 1 .AND. ID .LE. MAXEXT) N=EXTMAP(ID)
            IF (N .GT. 0) THEN
               IF (MATUID .EQ. 0) THEN
                  MATUID=ELMAT
               ELSE IF (MATUID .NE. ELMAT) THEN
                  NDIFF=NDIFF+1
               ENDIF
            ENDIF
         ENDIF
         GOTO 10
      ENDIF

C     ---- material orientation vectors ------------------------------
      IF (INOR .AND. INDEX(LINE,'<E>') .GT. 0) THEN
         P1=INDEX(LINE,'<E>')+3
         P2=INDEX(LINE,'</E>')-1
         IF (P2 .GT. P1) THEN
            READ(LINE(P1:P2),*,ERR=10) ID,(V(K),K=1,6)
            N=0
            IF (ID .GE. 1 .AND. ID .LE. MAXEXT) N=EXTMAP(ID)
            IF (N .GT. 0) THEN
               DO K=1,3
                  E1(K,N)=V(K)
                  E2(K,N)=V(K+3)
               END DO
               HAVEOR(N)=.TRUE.
            ENDIF
         ENDIF
      ENDIF
      GOTO 10

  800 CONTINUE

C     ---- pick the material of the group ----------------------------
      STAT=1
      DO I=1,NMAT
         IF (MUID(I) .EQ. MATUID) THEN
            RHO0=MRHO(I)
            CP0=MCP(I)
            DO K=1,3
               BK(K)=MK(K,I)
            END DO
            STAT=0
            IF (MNC(I)) STAT=2
         ENDIF
      END DO

      NNOOR=0
      DO I=1,GNUM
         N=GLIST(I)
         IF (.NOT. HAVEOR(N)) NNOOR=NNOOR+1
      END DO

      RETURN
      END


C=======================================================================
C     LINK INDEX: map positions grouped by parent hexahedron
C=======================================================================

      SUBROUTINE CFRP_LINKINDEX(NGRCOND,ICNDGR,ISCFRP,GLIST,GNUM,
     +     LNUM,LFIRST,LIDX)

      IMPLICIT NONE
      INTEGER NGRCOND,GNUM
      INTEGER ICNDGR(4,*),GLIST(*),LNUM(*),LFIRST(*),LIDX(*)
      LOGICAL ISCFRP(*)
      INTEGER I,J,N,PAR,POS

      DO I=1,GNUM
         LNUM(GLIST(I))=0
      END DO
      DO J=1,NGRCOND
         PAR=ICNDGR(4,J)
         IF (PAR .GT. 0 .AND. PAR .LE. 50000) THEN
            IF (ISCFRP(PAR)) LNUM(PAR)=LNUM(PAR)+1
         ENDIF
      END DO
      POS=1
      DO I=1,GNUM
         N=GLIST(I)
         LFIRST(N)=POS
         POS=POS+LNUM(N)
         LNUM(N)=0
      END DO
      DO J=1,NGRCOND
         PAR=ICNDGR(4,J)
         IF (PAR .GT. 0 .AND. PAR .LE. 50000) THEN
            IF (ISCFRP(PAR)) THEN
               LIDX(LFIRST(PAR)+LNUM(PAR))=J
               LNUM(PAR)=LNUM(PAR)+1
            ENDIF
         ENDIF
      END DO

      RETURN
      END


C=======================================================================
C     THICKNESS FACES AND LINK CLASSES (no global x/y/z is used)
C
C     For every hexahedron N and each of its boundary elements B the
C     direction  d = (centroid of the neighbouring hexahedron sharing B)
C                    - (centroid of N)
C     is compared with the through-thickness material direction ETA(N)
C     (from the NX orientation vectors). Faces with |cos(d,ETA)| >=
C     THKCOS are the two "thickness faces".
C       METH = 1  both thickness faces found from the mesh topology
C       METH = 2  found by elimination (the only boundary elements left)
C       METH = 3  found from the conductance pattern (the top<->bottom
C                 link has by far the smallest |GG0|, flat bricks)
C       METH = 0  not resolved -> conductances of this element are left
C                 unchanged (counted and reported)
C     LCLS(link) = 1 : link between the two thickness faces (scaled with
C                      the through-thickness conductivity)
C                  0 : every other link (in-plane conductivity)
C                 -1 : element unresolved, link not modified
C=======================================================================

      SUBROUTINE CFRP_TOPO(GLIST,GNUM,ISCFRP,NBND,BND,NADJ,ADJN,
     +     ELGEOM,ELPROP,ETA,VOLE,LFIRST,LNUM,LIDX,ICNDGR,GG0,MAXEL,
     +     THKCOS,ZRATIO,ZB1,ZB2,NBP,NBM,PITCH,METH,LCLS,NMETH)

      IMPLICIT NONE

      INTEGER GNUM,MAXEL
      INTEGER GLIST(*),NBND(*),BND(8,*),NADJ(*),ADJN(4,*)
      INTEGER LFIRST(*),LNUM(*),LIDX(*),ICNDGR(4,*)
      INTEGER ZB1(*),ZB2(*),NBP(*),NBM(*),METH(*),LCLS(*),NMETH(0:3)
      LOGICAL ISCFRP(*)
      REAL ELGEOM(3,4,*),ELPROP(10,*),GG0(*)
      DOUBLE PRECISION ETA(3,*),VOLE(*),PITCH(*)
      DOUBLE PRECISION THKCOS,ZRATIO

      INTEGER I,K,L,N,B,M,J,ID,E1,E2,NZ,NU,KA,KB,KF,KU,KV,MM
      INTEGER AX(8),NBR(8)
      DOUBLE PRECISION DOT(8),CSA(8),DX,DY,DZ,DL,CS,G,GMIN,GMIN2,PP
      LOGICAL KNOWN

      DO K=0,3
         NMETH(K)=0
      END DO

      DO I=1,GNUM
         N=GLIST(I)
         ZB1(N)=0
         ZB2(N)=0
         NBP(N)=0
         NBM(N)=0
         METH(N)=0
         PITCH(N)=ABS(VOLE(N))**(1.0D0/3.0D0)
         NZ=0
         NU=0

C        ---- direction of every boundary element of N --------------
         DO K=1,NBND(N)
            AX(K)=0
            NBR(K)=0
            DOT(K)=0.0D0
            CSA(K)=0.0D0
            KNOWN=.FALSE.
            DX=0.0D0
            DY=0.0D0
            DZ=0.0D0
            B=BND(K,N)
            IF (B .GT. 0 .AND. B .LE. MAXEL) THEN
               DO L=1,NADJ(B)
                  M=ADJN(L,B)
                  IF (M .NE. N .AND. NBR(K) .EQ. 0) THEN
                     IF (ISCFRP(M)) NBR(K)=M
                  ENDIF
               END DO
               IF (NBR(K) .GT. 0) THEN
                  M=NBR(K)
                  DX=DBLE(ELGEOM(1,1,M))-DBLE(ELGEOM(1,1,N))
                  DY=DBLE(ELGEOM(2,1,M))-DBLE(ELGEOM(2,1,N))
                  DZ=DBLE(ELGEOM(3,1,M))-DBLE(ELGEOM(3,1,N))
                  KNOWN=.TRUE.
               ELSE IF (ELPROP(1,B) .GT. 0.0) THEN
                  DX=DBLE(ELGEOM(1,1,B))-DBLE(ELGEOM(1,1,N))
                  DY=DBLE(ELGEOM(2,1,B))-DBLE(ELGEOM(2,1,N))
                  DZ=DBLE(ELGEOM(3,1,B))-DBLE(ELGEOM(3,1,N))
                  KNOWN=.TRUE.
               ENDIF
            ENDIF
            IF (KNOWN) THEN
               DL=DSQRT(DX*DX+DY*DY+DZ*DZ)
               IF (DL .GT. 1.0D-30) THEN
                  DOT(K)=DX*ETA(1,N)+DY*ETA(2,N)+DZ*ETA(3,N)
                  CS=DOT(K)/DL
                  CSA(K)=ABS(CS)
                  IF (CSA(K) .GE. THKCOS) THEN
                     AX(K)=3
                  ELSE
                     AX(K)=1
                  ENDIF
               ENDIF
            ENDIF
            IF (AX(K) .EQ. 3) NZ=NZ+1
            IF (AX(K) .EQ. 0) NU=NU+1
         END DO

C        ---- choose the two thickness faces ------------------------
         KA=0
         KB=0
         IF (NZ .GE. 2) THEN
            DO K=1,NBND(N)
               IF (AX(K) .EQ. 3) THEN
                  IF (KA .EQ. 0) THEN
                     KA=K
                  ELSE IF (CSA(K) .GT. CSA(KA)) THEN
                     KB=KA
                     KA=K
                  ELSE IF (KB .EQ. 0) THEN
                     KB=K
                  ELSE IF (CSA(K) .GT. CSA(KB)) THEN
                     KB=K
                  ENDIF
               ENDIF
            END DO
            METH(N)=1
         ELSE IF (NZ+NU .EQ. 2) THEN
            DO K=1,NBND(N)
               IF (AX(K) .EQ. 3 .OR. AX(K) .EQ. 0) THEN
                  IF (KA .EQ. 0) THEN
                     KA=K
                  ELSE
                     KB=K
                  ENDIF
               ENDIF
            END DO
            METH(N)=2
         ELSE IF (NZ+NU .GT. 2) THEN
C           conductance pattern: smallest |GG0| link
            GMIN=1.0D30
            GMIN2=1.0D30
            KF=0
            DO K=1,NBND(N)
               IF (AX(K) .EQ. 3) KF=K
            END DO
            KU=0
            KV=0
            DO L=1,LNUM(N)
               J=LIDX(LFIRST(N)+L-1)
               ID=ICNDGR(1,J)
               E1=ICNDGR(2,J)
               E2=ICNDGR(3,J)
               G=ABS(DBLE(GG0(ID)))
               KA=0
               KB=0
               DO K=1,NBND(N)
                  IF (BND(K,N) .EQ. E1) KA=K
                  IF (BND(K,N) .EQ. E2) KB=K
               END DO
               IF (KA .GT. 0 .AND. KB .GT. 0) THEN
                  MM=0
                  IF (NZ .EQ. 1) THEN
                     IF (KA .EQ. KF .AND. AX(KB) .EQ. 0) MM=1
                     IF (KB .EQ. KF .AND. AX(KA) .EQ. 0) MM=1
                  ELSE
                     IF (AX(KA) .EQ. 0 .AND. AX(KB) .EQ. 0) MM=1
                  ENDIF
                  IF (MM .EQ. 1) THEN
                     IF (G .LT. GMIN) THEN
                        GMIN2=GMIN
                        GMIN=G
                        KU=KA
                        KV=KB
                     ELSE IF (G .LT. GMIN2) THEN
                        GMIN2=G
                     ENDIF
                  ENDIF
               ENDIF
            END DO
            KA=0
            KB=0
            IF (KU .GT. 0 .AND. GMIN .LT. ZRATIO*GMIN2) THEN
               IF (NZ .EQ. 1) THEN
                  KA=KF
                  IF (KU .EQ. KF) THEN
                     KB=KV
                  ELSE
                     KB=KU
                  ENDIF
               ELSE
                  KA=KU
                  KB=KV
               ENDIF
               METH(N)=3
            ENDIF
         ENDIF

         IF (METH(N) .GT. 0 .AND. KA .GT. 0 .AND. KB .GT. 0) THEN
            ZB1(N)=BND(KA,N)
            ZB2(N)=BND(KB,N)

C           neighbours across the thickness faces and layer pitch
            PP=0.0D0
            DO K=1,2
               IF (K .EQ. 1) THEN
                  KF=KA
               ELSE
                  KF=KB
               ENDIF
               IF (NBR(KF) .GT. 0 .AND. DOT(KF) .NE. 0.0D0) THEN
                  IF (DOT(KF) .GT. 0.0D0) THEN
                     NBP(N)=NBR(KF)
                  ELSE
                     NBM(N)=NBR(KF)
                  ENDIF
                  IF (PP .EQ. 0.0D0) PP=ABS(DOT(KF))
               ENDIF
            END DO
            IF (PP .EQ. 0.0D0) THEN
               DO K=1,2
                  IF (K .EQ. 1) THEN
                     KF=KA
                  ELSE
                     KF=KB
                  ENDIF
                  IF (NBR(KF) .EQ. 0 .AND. DOT(KF) .NE. 0.0D0 .AND.
     +                PP .EQ. 0.0D0) PP=2.0D0*ABS(DOT(KF))
               END DO
            ENDIF
            IF (PP .GT. 0.0D0) PITCH(N)=PP
         ELSE
            METH(N)=0
         ENDIF
         NMETH(METH(N))=NMETH(METH(N))+1

C        ---- class of every link of this hexahedron ----------------
         DO L=1,LNUM(N)
            J=LIDX(LFIRST(N)+L-1)
            E1=ICNDGR(2,J)
            E2=ICNDGR(3,J)
            IF (METH(N) .EQ. 0) THEN
               LCLS(J)=-1
            ELSE IF ((E1 .EQ. ZB1(N) .AND. E2 .EQ. ZB2(N)) .OR.
     +               (E1 .EQ. ZB2(N) .AND. E2 .EQ. ZB1(N))) THEN
               LCLS(J)=1
            ELSE
               LCLS(J)=0
            ENDIF
         END DO
      END DO

      RETURN
      END


C=======================================================================
C     THICKNESS COLUMN: chain of hexahedra through the thickness,
C     starting from element N0 (used only for the temperature printout)
C=======================================================================

      SUBROUTINE CFRP_COLUMN(N0,NBP,NBM,ETA,ELGEOM,MAXCOL,NCOL,COLN,
     +     COLS)

      IMPLICIT NONE
      INTEGER N0,MAXCOL,NCOL
      INTEGER NBP(*),NBM(*),COLN(*)
      DOUBLE PRECISION ETA(3,*),COLS(*)
      REAL ELGEOM(3,4,*)
      INTEGER N,NSTEP,START,NXT,I

      N=N0
      NSTEP=0
   10 IF (NBM(N) .GT. 0 .AND. NSTEP .LT. MAXCOL) THEN
         N=NBM(N)
         NSTEP=NSTEP+1
         GOTO 10
      ENDIF
      START=N
      NCOL=1
      COLN(1)=N
      COLS(1)=0.0D0
   20 IF (NBP(N) .GT. 0 .AND. NCOL .LT. MAXCOL) THEN
         NXT=NBP(N)
         NCOL=NCOL+1
         COLN(NCOL)=NXT
         COLS(NCOL)=0.0D0
         DO I=1,3
            COLS(NCOL)=COLS(NCOL)+(DBLE(ELGEOM(I,1,NXT))
     +                 -DBLE(ELGEOM(I,1,START)))*ETA(I,START)
         END DO
         N=NXT
         GOTO 20
      ENDIF

      RETURN
      END


C=======================================================================
C     CONDUCTANCE UPDATE FOR THE WHOLE GROUP (GG is a RESISTANCE)
C
C     R_new = R0 / f      f = K(T,alpha)/K(base)  of the link's direction
C     LCLS = 1 : thickness link -> through-thickness factor FPERP
C     LCLS = 0 : any other link  -> in-plane factor FPAR
C     The sign of R0 is kept (negative-resistance links are scaled too).
C=======================================================================

      SUBROUTINE CFRP_UPDATE_GROUP(NGRCOND,ICNDGR,GG,GG0,UPDATED,NOCON,
     +     ISCFRP,FPAR,FPERP,LCLS,NMOD,CNT)

      IMPLICIT NONE

      INTEGER NGRCOND,NOCON,NMOD
      INTEGER ICNDGR(4,*),LCLS(*),CNT(*)
      REAL GG(*),GG0(*)
      LOGICAL UPDATED(*),ISCFRP(*)
      DOUBLE PRECISION FPAR(*),FPERP(*)

      INTEGER J,ID,PAR
      DOUBLE PRECISION G0,F

      DO J=1,NGRCOND
         ID=ICNDGR(1,J)
         PAR=ICNDGR(4,J)
         IF (ID .LE. 0 .OR. ID .GT. NOCON .OR. PAR .LE. 0 .OR.
     +       PAR .GT. 50000) THEN
            CNT(5)=CNT(5)+1
            GOTO 100
         ENDIF
         IF (.NOT. ISCFRP(PAR)) THEN
            CNT(5)=CNT(5)+1
            GOTO 100
         ENDIF
         IF (UPDATED(ID)) THEN
            CNT(6)=CNT(6)+1
            GOTO 100
         ENDIF
         IF (LCLS(J) .LT. 0) THEN
            CNT(4)=CNT(4)+1
            GOTO 100
         ENDIF
         G0=DBLE(GG0(ID))
         IF (G0 .EQ. 0.0D0) THEN
            CNT(5)=CNT(5)+1
            GOTO 100
         ENDIF
         IF (LCLS(J) .EQ. 1) THEN
            F=FPERP(PAR)
         ELSE
            F=FPAR(PAR)
         ENDIF
         IF (F .LE. 0.0D0) THEN
            CNT(5)=CNT(5)+1
            GOTO 100
         ENDIF
         GG(ID)=REAL(G0/F)
         UPDATED(ID)=.TRUE.
         NMOD=NMOD+1
         IF (LCLS(J) .EQ. 1) THEN
            CNT(1)=CNT(1)+1
         ELSE
            CNT(2)=CNT(2)+1
         ENDIF
         IF (G0 .LT. 0.0D0) CNT(3)=CNT(3)+1
  100    CONTINUE
      END DO

      RETURN
      END


C=======================================================================
C     LINK AUDIT: the links of one hexahedron, baseline and current GG
C=======================================================================

      SUBROUTINE CFRP_SAMPLES(TIME,NSAMP,SAMPJ,ICNDGR,LCLS,GG,GG0,
     +     ICONV)

      IMPLICIT NONE
      DOUBLE PRECISION TIME
      INTEGER NSAMP
      INTEGER SAMPJ(*),ICNDGR(4,*),LCLS(*),ICONV(*)
      REAL GG(*),GG0(*)
      INTEGER J,K,ID
      REAL RATIO

      IF (NSAMP .LE. 0) RETURN
      PRINT *,' '
      PRINT *,' LINK AUDIT t=',TIME,' (class 1 = thickness link)'
      DO K=1,NSAMP
         J=SAMPJ(K)
         ID=ICNDGR(1,J)
         RATIO=0.0
         IF (GG0(ID) .NE. 0.0) RATIO=GG(ID)/GG0(ID)
         PRINT *,'   id',ID,' ext',ICONV(ICNDGR(2,J)),
     +           ICONV(ICNDGR(3,J)),' class',LCLS(J),' GG0=',GG0(ID),
     +           ' GG=',GG(ID),' GG/GG0=',RATIO
      END DO
      PRINT *,' '

      RETURN
      END


C=======================================================================
C     THROUGH-THICKNESS TEMPERATURE PRINTOUT
C     CG temperatures of the thickness column; the exposed face is the
C     hotter end. Depth is measured from that face along the column.
C=======================================================================

      SUBROUTINE CFRP_PROFILE(TIME,T,NCOL,COLN,COLS,PROBE,ICONV)

      IMPLICIT NONE
      DOUBLE PRECISION TIME,PROBE
      DOUBLE PRECISION T(*),COLS(*)
      INTEGER NCOL
      INTEGER COLN(*),ICONV(*)

      INTEGER I,IE
      DOUBLE PRECISION SLO,SHI,DE,DEP,D1,D2,W,TP
      LOGICAL HOTHI,FOUND

      IF (NCOL .LT. 2) RETURN
      SLO=COLS(1)-0.5D0*(COLS(2)-COLS(1))
      SHI=COLS(NCOL)+0.5D0*(COLS(NCOL)-COLS(NCOL-1))
      HOTHI=(T(COLN(NCOL)) .GT. T(COLN(1)))

      PRINT *,' '
      PRINT *,' CFRP PROFILE t=',REAL(TIME),'  thickness =',
     +        REAL(SHI-SLO),' layers =',NCOL
      IF (TIME .GT. 1.0D0) THEN
         IF (HOTHI) THEN
            PRINT *,'   exposed (hotter) face = end of the column at',
     +              ' +direction (last layer)'
         ELSE
            PRINT *,'   exposed (hotter) face = end of the column at',
     +              ' -direction (first layer)'
         ENDIF
      ENDIF
      PRINT *,'   depth from exposed face [m], T [C], element label'
      DO I=1,NCOL
         IF (HOTHI) THEN
            IE=NCOL+1-I
         ELSE
            IE=I
         ENDIF
         IF (HOTHI) THEN
            DEP=SHI-COLS(IE)
         ELSE
            DEP=COLS(IE)-SLO
         ENDIF
         PRINT *,'   ',REAL(DEP),REAL(T(COLN(IE))),ICONV(COLN(IE))
      END DO

      FOUND=.FALSE.
      TP=0.0D0
      DO I=1,NCOL-1
         IF (HOTHI) THEN
            D1=SHI-COLS(NCOL+1-I)
            D2=SHI-COLS(NCOL-I)
            IF (.NOT. FOUND .AND. PROBE .GE. D1 .AND. PROBE .LE. D2)
     +      THEN
               W=(PROBE-D1)/(D2-D1)
               TP=T(COLN(NCOL+1-I))+W*(T(COLN(NCOL-I))
     +            -T(COLN(NCOL+1-I)))
               FOUND=.TRUE.
            ENDIF
         ELSE
            D1=COLS(I)-SLO
            D2=COLS(I+1)-SLO
            IF (.NOT. FOUND .AND. PROBE .GE. D1 .AND. PROBE .LE. D2)
     +      THEN
               W=(PROBE-D1)/(D2-D1)
               TP=T(COLN(I))+W*(T(COLN(I+1))-T(COLN(I)))
               FOUND=.TRUE.
            ENDIF
         ENDIF
      END DO
      IF (FOUND .AND. TIME .GT. 1.0D0) THEN
         PRINT *,'   >>> T at depth',REAL(PROBE),
     +           ' from the exposed face =',REAL(TP)
      ENDIF
      PRINT *,' '

      RETURN
      END
