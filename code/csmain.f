CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  AGB CIRCUMSTELLAR ENVELOPE CHEMICAL MODEL                           C
C     VERSION: 30 OCTOBER 2023                                         C
C  BY M. A. CORDINER(1), T. J. MILLAR(2) AND A. J. MARKWICK(3)         C
C     (1) NASA GODDARD SPACE FLIGHT CENTER                             C
C     (2) QUEEN'S UNIVERSITY BELFAST                                   C
C     (3) THE UNIVERSITY OF MANCHESTER                                 C
C                                                                      C
C  SOLVES FOR CHEMICAL ABUNDANCES AS A FUNCTION OF RADIUS IN A UNIFORM C
C  SPHERICALLY-SYMMETRIC OUTFLOW WITH CONSTANT VELOCITY AND MASS-LOSS  C
C  RATE. HANDLES UP TO 995 SPECIES (INCLUDING UP TO 10 CONSERVED       C
C  SPECIES) LINKED BY UP TO 10000 BINARY REACTIONS                     C
C                                                                      C
C  ADAPTED BY CATHERINE WALSH(4) FOR EASE OF RUNNING VIA SCRIPTING     C
C                                                                      C
C     (4) UNIVERSITY OF LEEDS                                          C
C                                                                      C
C  ADAPTED BY MARIE VAN DE SANDE(5) TO INCLUDE A CLUMPY OUTFLOW        C
C  (VAN DE SANDE ET AL. 2018, A&A, 616, A106)                          C                                                                C
C                                                                      C
C     (5) LEIDEN OBSERVATORY                                           C
C                                                                      C
C  ADAPTED BY MARIE VAN DE SANDE AND TOM MILLAR TO INCLUDE STELLAR     C
C  AND COMPANION PHOTONS                                               C
C                                                                      C
C  (VAN DE SANDE & MILLAR 2019, APJ, 873, 36,                          C                                           C
C  (VAN DE SANDE & MILLAR 2022, MNRAS, 510, 1204)                      C                                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PROGRAM CSE
      CHARACTER*12 SPECI(1000),PARENT(200)
      CHARACTER*11 JAC
      CHARACTER*500 FRATES,FSPECS,FOUTF,FOUTN,INFILE,FPP,FCD,
     *   OUTFOLDER,VELMODE,CLUMPMODE,FPPCD
      DOUBLE PRECISION Y(1000),K(10000),X(10),TOTAL(10),COL(1000),
     *   B(1000,0:200),F(1000,0:200),PABUND(200),RADII(0:200),
     *   RWORK(2000000)
      DOUBLE PRECISION HNR,ACCR,TEMP,V,A_G,X_G,H2COL,AUV,MLOSS,MH,MU,
     *   KB,ALBEDO,ZETA,GAMCO,RINIT,RFINAL,RAD,R,ROUT,RTOL,ATOL,PI,
     *   CORATE,RANA,AUV_AV,R_STAR,R_START,T_STAR,EPSIL,BETA,
     *   RESOLUTION,SENSITIVITY,FVOL,L,FMASS,FNOEXT,FIC,H,HNRMEAN,
     *   H2COLEFF,RSCALE,RDUST,GSTAR,DELTA_AUV,DISTANCE,GETGSTAR,
     *   H2COLDUSTEFF,DELTA_H2COLEFF,GE0,GBIN,HNRDUST,H2COLDUST,TBIN,
     *   AUVDUST,HNRMEANDUST,RSCALE_INPUT,RDUST_INPUT,BINSCALE,RBIN
      DOUBLE PRECISION ALF(10000,5),BET(10000,5),
     *   GAM(10000,5),TLOWER(10000,5),TUPPER(10000,5)
      DOUBLE PRECISION GETHNR,GETTEM,GETH2,GETAUV,GETRAD,GETCOR
      INTEGER IWORK(1025),I,J,USPEC,URATES,UFRAC,UNUM,NSPEC,NCONS,
     *   NPAR,ICO,IN2,IRUN,IS,IF,LI,LTAG,ITASK,ISTATE,IOPT,LRW,LIW,MF,
     *   IANA,UANA,UIN,UCO,UN2,NSTEPS,UPP,UCD,FULLANA,IANA_ORIG,CLUMP,
     *   CLUMPING,ISTELLAR,IBIN,ICLUMP,ITOL,IPAR,NTOT,NEXTRA,IH,IH2
      INTEGER RTYPE(10000),NTR(10000),NREAC,ICOR,TR,ICOIP,ICOAP
      CHARACTER*10 EXT
      CHARACTER*20 OUTF,OUTN,PP,CD
      CHARACTER*100 OUT
      CHARACTER DUMMY


      EXTERNAL DIFFUN
      
      COMMON/BL1/ K,X,TOTAL,ACCR,HNR
      COMMON/BL2/ MLOSS,V
      COMMON/BL3/ Y,X_G,A_G,ZETA,ALBEDO,GAMCO,AUV_AV,ICO,GSTAR,GBIN,
     *   RSCALE,RDUST,DELTA_AUV,ISTELLAR,IBIN,RBIN,TBIN
      COMMON/BL4/ ALF,BET,GAM,TLOWER,TUPPER,RTYPE,NTR,NREAC,ICOR,ICOIP,
     *  ICOAP
      COMMON/BLC/ PI,MH,MU,KB
      COMMON/BLTEMP/ R_STAR,T_STAR,EPSIL
      COMMON/BLSTEPS/ NSTEPS,RESOLUTION
      COMMON/ANA/ OUTFOLDER,SENSITIVITY
      COMMON/CLM/ CLUMPMODE,ICLUMP,FVOL,L,FIC
      
      
C Data I/O unit numbers
      DATA USPEC,URATES,UFRAC,UNUM,UANA,NEXTRA,UIN/1,2,3,4,5,8,9/
      DATA UPP,UCD/7,8/
C Physical constants
      DATA PI,MH,KB/3.1415927,1.6605E-24,1.3807E-16/
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C                 FIRST STEP: READ INPUTFILE                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

C  READ IN INPUT PARAMETERS 
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
            
      CALL GETARG(1,INFILE)
      
      INFILE = TRIM(INFILE)
      
      IF(INFILE.EQ.'') THEN
          WRITE(*,*) 'ERROR: No input file specified!'
          WRITE(*,*) 'Useage: % ./csmodel input_parameters.txt'
          STOP
      END IF
      
      WRITE(*,*) 'Input file = ', INFILE
      
      OPEN(UNIT=UIN,FILE=INFILE)
      DO I=1,3
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,R_STAR
      READ(UIN,*) DUMMY,DUMMY,T_STAR
      READ(UIN,*) DUMMY,DUMMY,MLOSS
      READ(UIN,*) DUMMY,DUMMY,V
      READ(UIN,*) DUMMY,DUMMY,EPSIL
      DO I=1,2
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,A_G
      READ(UIN,*) DUMMY,DUMMY,X_G
      READ(UIN,*) DUMMY,DUMMY,ALBEDO
      READ(UIN,*) DUMMY,DUMMY,ZETA
      READ(UIN,*) DUMMY,DUMMY,GAMCO
      READ(UIN,*) DUMMY,DUMMY,AUV_AV
      DO I=1,2
         READ(UIN,*)
      END DO      
      READ(UIN,*) DUMMY,DUMMY,RINIT
      READ(UIN,*) DUMMY,DUMMY,RFINAL
      READ(UIN,*) DUMMY,DUMMY,RESOLUTION
      DO I=1,2
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,IANA_ORIG
      READ(UIN,*) DUMMY,DUMMY,RANA
      READ(UIN,*) DUMMY,DUMMY,FULLANA
      READ(UIN,*) DUMMY,DUMMY,SENSITIVITY
      DO I=1,2
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,RTOL
      READ(UIN,*) DUMMY,DUMMY,ATOL
      DO I=1,2
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,CLUMPMODE
      READ(UIN,*) DUMMY,DUMMY,FVOL
      READ(UIN,*) DUMMY,DUMMY,L
      READ(UIN,*) DUMMY,DUMMY,FIC
      DO I=1,2
         READ(UIN,*)
      END DO      
      READ(UIN,*) DUMMY,DUMMY,ISTELLAR
      READ(UIN,*) DUMMY,DUMMY,IBIN
      READ(UIN,*) DUMMY,DUMMY,RSCALE_INPUT
      READ(UIN,*) DUMMY,DUMMY,RDUST_INPUT
      READ(UIN,*) DUMMY,DUMMY,TBIN
      READ(UIN,*) DUMMY,DUMMY,RBIN
      DO I=1,2
         READ(UIN,*)
      END DO
      READ(UIN,*) DUMMY,DUMMY,FRATES
      READ(UIN,*) DUMMY,DUMMY,FSPECS
      READ(UIN,*)
      READ(UIN,*) DUMMY,DUMMY,OUTFOLDER
      CLOSE(UNIT=UIN)
      
      
C  OPEN RATE, SPECIES AND OUTPUT FILES
      OPEN(UNIT=USPEC, FILE=FSPECS, STATUS='OLD')

      
C  INTERNAL PHOTON PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  RADIUS AT WHICH PHOTORATES DUE TO BINARY PHOTONS ARE CALCULATED
      RSCALE = 50*R_STAR

C  CONVERT RDUST FROM UNITS OF R_STAR TO CM
      RDUST = RDUST_INPUT*R_STAR 
      
C  SET RATIO OF INTEGRATED UV FLUX AT RSCALE = 50 * R_STAR IN IRC+20126
C  TO INTEGRATED DRAINE ISM FIELD
C  INTEGRATED DRAINE UV FIELD = 5.89e11 PHOTONS M-2 S-1 SR-1 OVER 912-2050 A
      GSTAR = GETGSTAR(T_STAR,RSCALE/R_STAR)
      write(*,*) 'GSTAR IN MAIN CODE = ', GSTAR
      
C  SET RADIUS AT WHICH PHOTORATES DUE TO BINARY PHOTONS ARE CALCULATED
C  THIS RADIUS IS TAKEN TO BE THE SAME AS FOR STELLAR PHOTONS,
C  ENABLING TO DIRECTLY COMPARE PHOTORATES
      BINSCALE = RSCALE
C SET RATIO OF INTEGRATED BINARY UV FLUX AT RSCALE TO INTEGRATED DRAINE FIELD AS ABOVE
      GBIN = GETGSTAR(TBIN,BINSCALE/RBIN)
      write(*,*) 'GBIN IN MAIN CODE = ', GBIN
      
      
      
      
C  READ IN CHEMICAL SPECIES : NAMES AND ABUNDANCES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      WRITE(*,*) 'Reading species file...'
      READ(USPEC,1)
 1    FORMAT(:)

C  INITIALISE NUMBER OF SPECIES
      NSPEC = 0
C  INITIALISE NUMBER OF CONSERVED SPECIES
      NCONS = 0
C  INITIALISE NUMBER OF PARENT SPECIES
      NPAR = 0

      I = 1
C  LOOP TO READ SPECIES NAMES
 100  READ(USPEC,2,END=101) J,SPECI(I)
 2       FORMAT(1X,I4,1X,A12)
C  DETERMINE INDEX OF CO
         IF(SPECI(I).EQ.'CO') ICO = I
C  DETERMINE INDEX OF H
         IF(SPECI(I).EQ.'H') IH = I
         I=I+1
         IF (J.EQ.9999) THEN
            I=I-1
            NSPEC = I - 1
            GO TO 102
         ENDIF
      GO TO 100

C  LOOP TO READ CONSERVED SPECIES
 102  IF (SPECI(I).EQ.'Conserved:') THEN
         NCONS = 1
 103     READ(USPEC,3,END=101) J,SPECI(I),TOTAL(NCONS)
 3          FORMAT(1X,I4,1X,A12,1X,0PE8.2) 
C  DETERMINE INDEX OF H2
         IF(SPECI(I).EQ.'H2') IH2 = I
            I=I+1
            NCONS = NCONS + 1
            IF (J.EQ.9999) THEN
               I=I-1
               NCONS = NCONS - 2
               GO TO 104
            ENDIF
         GO TO 103
      ENDIF

C  LOOP TO READ PARENT ABUNDANCES
 104  IF (SPECI(I).EQ.'Parents:') THEN
         NPAR = 1
 105     READ(USPEC,*,END=101) PARENT(NPAR),PABUND(NPAR)
C CALCULATE THE APPROPRIATE VALUE OF MU (average mass per H molecule) TAKING INTO ACCOUNT ABUNDANCE OF He
            IF (PARENT(NPAR).EQ.'He') THEN 
C  WHEN USING FRACTIONAL ABUNDANCE WRT H2
C               MU = 2.0 + 4.0*PABUND(NPAR)*2.
                MU = 1.0 + 4.0*PABUND(NPAR)
                WRITE(*,*) 'AVERAGE MASS PER H ATOM = ',MU
            ENDIF
            NPAR = NPAR + 1
         GO TO 105
      ENDIF

 101  CONTINUE

      NPAR = NPAR - 1
      NTOT = NSPEC+NCONS
      
      WRITE(*,*) 'Species Read: ',NTOT
      WRITE(*,*) 'Conserved species: ',NCONS
      WRITE(*,*) 'Parent species: ',NPAR

C  CATCH ERRORS IN INPUT SPECIES
      IF(NSPEC.LT.1) CALL DIE
     *   ("No species defined in model -- check species input file    ")
      IF(NPAR.LT.1) CALL DIE
     *   ("No parents defined in model -- check species input file    ")
     
C  READ RATE FILE AND EXTRACT RATE DATA FOR EACH REACTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(*,*) 'Reading rate file...'
      CALL READR(FRATES,URATES)     
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C                   SECOND STEP: INITIALISATION                       C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC CLUMPING LOOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  IF A TWO-COMPONENT MODEL IS REQUESTED (CLUMPMODE=POROSITY),
C  ICLUMP IS SET TO 1 SO THAT BOTH COMPONENTS ARE CALCULATED
      WRITE(*,*) 'CLUMPMODE ', CLUMPMODE
      IF (CLUMPMODE.EQ.'POROSITY') THEN
            CLUMPING=1
            WRITE(*,*) 'CLUMPING = ',ICLUMP
      ELSE
            CLUMPING = 0
      END IF
      WRITE(*,*) 'ICLUMP = ',ICLUMP


C  CLUMPING LOOP
      ICLUMP = -1
 900  ICLUMP = ICLUMP+1 
      IANA = IANA_ORIG

C  INITIALISE OUTPUT FILES   
      IF (ICLUMP.EQ.0.AND.CLUMPMODE.EQ.'POROSITY') THEN
            IF (FIC.EQ.0 .OR. FIC.EQ.1) THEN   
C                 ONE-COMPONENT OUTFLOW (EVERYTHING INSIDE OF CLUMPS)
C                 OR NOTHING INSIDE OF CLUMPS (SMOOTH CASE)
C                 GO DIRECTLY TO FINAL LOOP
                  GO TO 900
            ELSE
                  OUTF='csfrac_interclump'
                  OUTN='csnum_interclump'
                  PP='csphyspar_interclump'
                  CD='cscoldens_interclump'
      WRITE(*,*) 'Start calculation of interclump component...'	
            END IF
      ELSE IF (ICLUMP.EQ.1) THEN
            IF (CLUMPMODE.EQ.'POROSITY'.AND.FIC.EQ.0) THEN
C                 EVERYTHING INSIDE OF CLUMPS: JUST ONE, MAJOR COMPONENT 
                  OUTF='csfrac_clump'
                  OUTN='csnum_clump'
                  PP='csphyspar_clump'
                  CD='cscoldens_clump'
      WRITE(*,*) 'Starting calculation of only clump component...'
            ELSE IF (CLUMPMODE.EQ.'POROSITY'.AND.FIC.EQ.1) THEN
C                 NOTHING INSIDE OF CLUMPS: SMOOTH CASE 
C                 Nothing is in clumps, smooth case
                  OUTF='csfrac_interclump'
                  OUTN='csnum_interclump'
                  PP='csphyspar_interclump'
                  CD='cscoldens_interclump'
      WRITE(*,*) 'Starting calculation of only interclump component...'
            ELSE
C                 MINOR (CLUMPED) COMPONENT OF POROUS OUTFLOW
                  OUTF='csfrac_clump'
                  OUTN='csnum_clump'
                  PP='csphyspar_clump'
                  CD='cscoldens_clump'
      WRITE(*,*) 'Starting calculation of clump component...'
            END IF
            
      ELSE
C           IN CASE OF SMOOTH: JUST ONE COMPONENT
            OUTF='csfrac_smooth'
            OUTN='csnum_smooth'
            PP='csphyspar_smooth'
            CD='cscoldens_smooth'
            WRITE(*,*) 'Starting calculation of smooth outflow...'
      END IF

      
C  INITIALISE OUTPUT FOLDER AND OUTPUT FILES      
      EXT='.out'
      FOUTF= TRIM(OUTFOLDER)//TRIM(OUTF)//TRIM(EXT)
      FOUTN= TRIM(OUTFOLDER)//TRIM(OUTN)//TRIM(EXT)
      FPP= TRIM(OUTFOLDER)//TRIM(PP)//TRIM(EXT)
      FCD= TRIM(OUTFOLDER)//TRIM(CD)//TRIM(EXT)
      
      OPEN(UNIT=UFRAC, FILE=FOUTF)
      OPEN(UNIT=UNUM, FILE=FOUTN)
      OPEN(UNIT=UPP, FILE=FPP)
      OPEN(UNIT=UCD, FILE=FCD)    
 
      
      
C  SET INITIAL PHYSICAL PARAMETERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  H2 NUMBER DENSITY - ACCOUNT FOR POROSITY
      HNRMEAN = GETHNR(RINIT)
      IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.0) THEN
            HNR = FIC*HNRMEAN
      ELSE IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.1) THEN	
            HNR = (HNRMEAN/FVOL)*(1-(1-FVOL)*FIC)
      ELSE 
            HNR = HNRMEAN
      END IF
      
C  H2 NUMBER DENSITY AT DUST CONDENSATION RADIUS
      HNRMEANDUST = GETHNR(RDUST)
      IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.0) THEN
            HNRDUST = FIC*HNRMEANDUST
      ELSE IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.1) THEN	
            HNRDUST = (HNRMEANDUST/FVOL)*(1-(1-FVOL)*FIC)
      ELSE 
            HNRDUST = HNRMEANDUST
      END IF 
      
C  KINETIC TEMPERATURE
      TEMP = GETTEM(RINIT)
      
C  H2 COLUMN DENSITY OF GAS 
C  (EFF = EFFECTIVE COL DENS IN A CLUMPY OUTFLOW)
      H2COL = GETH2(RINIT,HNRMEAN)
      H2COLEFF = GETH2(RINIT,HNR)
C  H2 COLUMN DENSITY AT DUST CONDENSATION RADIUS
      H2COLDUST = GETH2(RDUST,HNRMEANDUST)
      H2COLDUSTEFF = GETH2(RDUST,HNRDUST)

C  RADIAL UV EXTINCTION BETWEEN R AND INFINITY
      AUV = GETAUV(H2COL,RINIT)
C  RADIAL UV EXTINCTION FROM RDUST ONWARDS
      AUVDUST = GETAUV(H2COLDUST,RDUST)
C  UV EXTINCTION BETWEEN RDUST AND R
      DELTA_AUV = AUVDUST - AUV
      
C  RADIATION FIELD STRENGTH
      RAD = GETRAD(AUV)
      
      
C  INITIALIZE ABUNDANCE ARRAY
      DO I = 1,NSPEC
         Y(I)= 0.0
      END DO

C  INITIALIZE 0TH ELEMENT OF OUTPUT DATA STORAGE ARRAYS
      DO I = 1,NTOT+NEXTRA
         B(I,0) = 0.0
         F(I,0) = 0.0
      END DO

C  SET INITIAL ABUNDANCES
      DO I = 1,NPAR
         DO J = 1,NTOT
            IF(SPECI(J).EQ.PARENT(I)) THEN
            Y(J) = PABUND(I)
C  LOAD FIRST ELEMENT OF OUTPUT ARRAYS
            B(J,0) = PABUND(I)
            F(J,0) = PABUND(I) * HNR
            WRITE(*,106) SPECI(J),PABUND(I)
 106        FORMAT(3X,A12,1PE8.2)
            ENDIF
         END DO
      END DO

C  LOAD 0TH ELEMENT OF CONSERVED SPECIES ARRAYS
      J = 1
      DO I = NSPEC+1,NTOT
         B(I,0) = TOTAL(J)
         F(I,0) = TOTAL(J) * HNR
         J = J+1
      END DO

C  INITIAL CO PHOTODISSOCIATION RATE
      CORATE = GETCOR(H2COLEFF,Y(ICO),V,AUV)
      WRITE(*,*) 'CORATE PH RATE = ', CORATE

C  PHYSICAL PARAMETER NAME STRINGS
      SPECI(NTOT+1) = 'n(Htot)'
      SPECI(NTOT+2) = 'TEMP'
      SPECI(NTOT+3) = 'A_UV EFF'
      SPECI(NTOT+4) = 'RAD. EFF'
      SPECI(NTOT+5) = 'CO K(PH)'
      SPECI(NTOT+6) = 'CO K(IP)'
      SPECI(NTOT+7) = 'CO K(AP)'
      SPECI(NTOT+8) = 'DELTA_AUV'
      
C  LOAD 0TH ELEMENTS OF PHYSICAL PARAMETER ARRAYS
      B(NTOT+1,0) = HNR
      B(NTOT+2,0) = TEMP
      B(NTOT+3,0) = AUV
      B(NTOT+4,0) = RAD
C  CO PHOTODISSCIATION, ISM RATE
      B(NTOT+5,0) = K(ICOR)
C  CO PHOTODISSCIATION, STELLAR PHOTON RATE
      IF(ISTELLAR.EQ.1) THEN
        B(NTOT+6,0) = K(ICOIP)
      ELSE
        B(NTOT+6,0) = 0.0
      ENDIF 
C  CO PHOTODISSCIATION, COMPANION PHOTON RATE
      IF(IBIN.EQ.1) THEN
        B(NTOT+7,0) = K(ICOAP)
      ELSE
        B(NTOT+7,0) = 0.0
      ENDIF
      B(NTOT+8,0) = DELTA_AUV

      
C START DVODE SOLVER
CCCCCCCCCCCCCCCCCCCCC
      
C  INITIALISE INTEGRATION STEP VARIABLES
      R = RINIT
      ROUT = RINIT
      IRUN = 0
      RADII(IRUN) = ROUT


C  INITIALISE DVODE SOLVER VARIABLES
      LIW = NTOT + 30
      LRW = 22 + (9*NTOT) + (2*(NTOT**2))
      ITOL   = 1
      ITASK  = 1
      ISTATE = 1
      IOPT   = 0
      JAC    = 'DUMMYMATRIX'
      MF     = 22

      WRITE(*,*) 'Starting ODE solver DVODE...'


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C                 THIRD STEP: MAIN SOLVER LOOP                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C INCREASE COUNTER, SELECT RADIUS AND VELOCITY 
 300  ROUT = 10.0**(LOG10(ROUT) + RESOLUTION)
C  STEP NUMBER INCREMENT
      IRUN = IRUN + 1
      
C  CALL VARIABLE-COEFFICIENT ODE SOLVER (DVODE)
C  TO CALCULATE CHEMICAL ABUNDANCES (Y) WITH RESPECT TO RADIUS (ROUT)
C  BY INTEGRATION OF DY/DR VALUES CALCULATED BY SUBROUTINE 'DIFFUN'  
      CALL DVODE (DIFFUN,NSPEC,Y,R,ROUT,ITOL,RTOL,ATOL,ITASK,
     &            ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)

      IF(ISTATE.LT.0) THEN
         ISTATE = 2
      ENDIF       
      
      WRITE(*,301) ROUT
 301  FORMAT(1X,"Radius: ",1PE8.2," cm")
    

C  LOAD OUTPUT ARRAYS
CCCCCCCCCCCCCCCCCCCCCCC

C  PHYSICAL PARAMETERS
      HNRMEAN = GETHNR(ROUT)
      IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.0) THEN
            HNR = FIC*HNRMEAN
      ELSE IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.1) THEN	
            HNR = (HNRMEAN/FVOL)*(1-(1-FVOL)*FIC)
      ELSE 
            HNR = HNRMEAN
      END IF 
            
      HNRMEANDUST = GETHNR(RDUST)
      IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.0) THEN
            HNRDUST = FIC*HNRMEANDUST
      ELSE IF (CLUMPMODE.EQ.'POROSITY'.AND.ICLUMP.EQ.1) THEN	
            HNRDUST = (HNRMEANDUST/FVOL)*(1-(1-FVOL)*FIC)
      ELSE 
            HNRDUST = HNRMEANDUST
      END IF 

C  KINETIC TEMPERATURE
      TEMP = GETTEM(ROUT)
      
C  H2 COLUMN DENSITY OF GAS 
C  (EFF = EFFECTIVE COL DENS IN A CLUMPY OUTFLOW)
      H2COL = GETH2(ROUT,HNRMEAN)
      H2COLDUST = GETH2(RDUST,HNRMEANDUST)
C  H2 COLUMN DENSITY AT DUST CONDENSATION RADIUS
      H2COLEFF = GETH2(ROUT,HNR)
      H2COLDUSTEFF = GETH2(RDUST,HNRDUST)

C  RADIAL UV EXTINCTION
      AUV = GETAUV(H2COL,ROUT)
C  RADIAL UV EXTINCTION FROM RDUST ONWARDS
      AUVDUST = GETAUV(H2COLDUST,RDUST)
C  UV EXTINCTION BETWEEN RDUST AND R
      DELTA_AUV = AUVDUST - AUV
      
C  RADIATION FIELD STRENGTH
      RAD = GETRAD(AUV)
      

C  WRITE PHYSICAL PARAMETER ARRAYS
      B(NTOT+1,IRUN) = HNR
      B(NTOT+2,IRUN) = TEMP
      B(NTOT+3,IRUN) = AUV
      B(NTOT+4,IRUN) = RAD
      B(NTOT+5,IRUN) = K(ICOR)
      IF(ISTELLAR.EQ.1) THEN
        B(NTOT+6,IRUN) = K(ICOIP)
      ELSE
        B(NTOT+6,IRUN) = 0.0
      ENDIF 
      IF(IBIN.EQ.1) THEN
        B(NTOT+7,IRUN) = K(ICOAP)
      ELSE
        B(NTOT+7,IRUN) = 0.0
      ENDIF
      B(NTOT+8,IRUN) = DELTA_AUV
      
C  CHEMICAL ABUNDANCES AND NUMBER DENSITIES
      DO I=1,NSPEC
         B(I,IRUN) = Y(I)
         F(I,IRUN) = Y(I)*HNR
      END DO

      DO I=1,NCONS
         B(NSPEC+I,IRUN) = X(I)
         F(NSPEC+I,IRUN) = X(I)*HNR
      END DO
      
      RADII(IRUN) = ROUT

C  CALL THE SUBROUTINE TO ANALYSE THE OUTPUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ANALYSE ROUTINE AT A SINGLE RADIUS
      IF(IANA.EQ.1.AND.ROUT.GE.RANA) THEN
      WRITE(*,*) 'Analysing chemistry...'
      CALL ANALYSE(SPECI,B,K,HNR,ROUT,IRUN,UANA,NTOT,FRATES,URATES)
      WRITE(*,*) 'Resuming model...'
      IANA = 0
      ENDIF

C  ANALYSE ROUTINE AT ALL RADII
      IF(FULLANA.EQ.1) THEN
      WRITE(*,*) 'Analysing chemistry...'
      CALL ANALYSE(SPECI,B,K,HNR,ROUT,IRUN,UANA,NTOT,FRATES,URATES)
      WRITE(*,*) 'Resuming model...'
      IANA = 0
      ENDIF


C  TEST FOR END CONDITION
CCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(ROUT.LT.RFINAL) GO TO 300

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C                        STEP FOUR: OUTPUT                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  

C  CALCULATE COLUMN DENSITIES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 400 I = 1,NTOT
         DO 401 J = 0,IRUN-1
            COL(I)=COL(I)+((RADII(J+1)-RADII(J))*(F(I,J)+
     *       (0.5*(F(I,J+1)-F(I,J)))))
 401     CONTINUE
 400  CONTINUE
      WRITE(*,310) FOUTF, FOUTN, FPPCD 
 310  FORMAT(1X,'Writing results to',1X,A25,1X,',',1X,A25,'and',1X,A35)
 
  
C  OUTPUT ROUTINE
CCCCCCCCCCCCCCCCCC	  
C  FRACTIONAL ABUNDANCES AND NUMBER DENSITIES
      LI=NTOT+1
      IS=1
      IF=10
      WRITE(UFRAC,500)
 500  FORMAT(1X,' ** CALCULATED FRACTIONAL ABUNDANCES (/H2) ** '//)
      WRITE(UNUM,520)
 520  FORMAT(1X,' ** CALCULATED NUMBER DENSITIES (/CM^3) ** '//)
 501  WRITE(UFRAC,502)(SPECI(I),I=IS,IF)
 502  FORMAT(3X,'RADIUS',6X,10(1A12))
      WRITE(UNUM,502)(SPECI(I),I=IS,IF)
      LTAG=11
      DO 504 I=0,IRUN
         IF(I.EQ.LTAG) LTAG=LTAG+10
         WRITE(UFRAC,503) RADII(I),(B(J,I),J=IS,IF)
         WRITE(UNUM,503) RADII(I),(F(J,I),J=IS,IF)
 503  FORMAT(1X,11(1PE12.3E3))
 504  CONTINUE
      WRITE(UFRAC,502)(SPECI(I),I=IS,IF)
      WRITE(UNUM,502)(SPECI(I),I=IS,IF)
      WRITE(UFRAC,510)
      WRITE(UNUM,510)
      IS=IS+10
      IF=IF+10
      IF(IF.GT.(NTOT)) IF=NTOT
      IF(IS.GE.LI) GO TO 505
      GO TO 501
 
C OUTPUT PHYSICAL PARAMETERS
 505  WRITE(UPP,506)
 514  FORMAT(3X,'RADIUS',6X,11(1A12))
 513  FORMAT(1X,12(1PE12.3E3))
 506  FORMAT(1X,' **  PHYSICAL PARAMETERS  **'//)
      WRITE(UPP,514)(SPECI(I),I=NTOT+1,NTOT+NEXTRA)
      DO 507 I=0,IRUN
         WRITE(UPP,513) RADII(I),(B(J,I),J=NTOT+1,NTOT+NEXTRA)
 507   CONTINUE

C  OUTPUT COLUMN DENSITIES
      WRITE(UCD,510)
 510  FORMAT(//)
      WRITE(UCD,511)
 511  FORMAT(19X,'**  CALCULATED COLUMN DENSITIES  **'//)
      WRITE(UCD,512)(SPECI(I),COL(I),I=1,NTOT)
 512  FORMAT(3(1X,A12,'= ',1PE11.3,2X,:))

      CLOSE(UNIT=UFRAC)
      CLOSE(UNIT=UNUM)
      CLOSE(UNIT=UPP)
      CLOSE(UNIT=UCD)
            
      WRITE(*,*) "DONE!"
      
      IF (ICLUMP.LT.CLUMPING) GO TO 900

      END
