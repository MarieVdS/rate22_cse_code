CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO CALCULATE GSTAR, THE RATIO OF THE STELLAR UV FLUX      C
C  AT 50 STELLAR RADII TO THE DRAINE FLUX                             C
C  RATIO OF INTEGRATED FLUXES                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION GETGSTAR(T,RADIUS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMDA, RADIUS
      DOUBLE PRECISION T
      CHARACTER*10 STAR,FIELD
      CHARACTER*6 MOLECULE
      CHARACTER*3 PROCESS
      CHARACTER*20 PPFILE
C T = EFFECTIVE TEMPERATURE OF STAR
      DIMENSION LAMDA(29999),FLUX(29999)
      DATA C/2.9979E10/,pie/3.14159/

       
C  OPEN OUTPUT FILE TO WRITE UNSHIELDED PHOTORATES
      STAR = 'BB'
      OPEN(UNIT=3,FILE='inner_unshielded_ratio_new_'//TRIM(STAR)//
     * '.out') 
      
      WRITE(3,999) RADIUS
      WRITE (3,296)
      
      SUMDRAINE = 0.0
      DO 2 I=1,1139
      LAMDA(I) = 911.0 + I

C INTEGRATED DRAINE INTERSTELLAR UV FIELD FLUX 
      IF (LAMDA(I).LT.2000.0)  THEN
        FLUX(I) = 3.2028E+15*(LAMDA(I)**(-3.0))
     *  -5.1542E+18*(LAMDA(I)**(-4.0))+ 2.0546E+21*(LAMDA(I)**(-5.0)) 
      END IF
      IF (LAMDA(I).GE.2000.0) THEN
        FLUX(I) = 732.0*(LAMDA(I)**0.7)
      END IF
 2    CONTINUE

      LAMDA(1) = 2000.0
      FLUXUPPER = -0.5*(3.2028E+15*(LAMDA(1)**(-2.0))) +(1./3.)*
     *   (5.1542E+18*(LAMDA(1)**(-3.0))) -
     *   0.25*(2.0546E+21*(LAMDA(1)**(-4.0))) 
      LAMDA(2) = 912.0
      FLUXLOWER = -0.5*(3.2028E+15*(LAMDA(2)**(-2.0))) +(1./3.)*
     *   (5.1542E+18*(LAMDA(2)**(-3.0))) -
     *   0.25*(2.0546E+21*(LAMDA(2)**(-4.0))) 
      SUMDRAINE1 = (FLUXUPPER - FLUXLOWER)
      LAMDA(1) = 2050.0
      FLUXUPPER = (732.0/1.7)*(LAMDA(1)**1.7)
      LAMDA(2) = 2000.0
      FLUXLOWER = (732.0/1.7)*(LAMDA(2)**1.7)
      SUMDRAINE2 = FLUXUPPER - FLUXLOWER
      SUMDRAINETOT = SUMDRAINE1 + SUMDRAINE2

C  PHYSICAL CONSTANTS IN SI UNITS
C  CALCULATE BAND RADIANCE FROM BB AT TEMPERATURE T
C  FORMULAE BELOW FOR BAND RADIANCE IN PHOTONS PER UNIT AREA PER SECOND PER STERADIAN 
C  TAKEN FROM DOCUMENT ON SPECTRAL CALCULATOR WEBSITE
      NITER = 5
      CSI = C/100.0
      BOLTZ = 1.3807E-23
      HPLANCK = 6.6261E-34
C SET WAVELENGTHS IN CM-1 from angstroms
      SIGMA1 = 1.0e8/912.0
      SIGMA2 = 1.0e8/2050.0
C SET X1 AND X2
      X1 = 100.0*(HPLANCK*CSI*SIGMA1)/(BOLTZ*T)
      X2 = 100.0*(HPLANCK*CSI*SIGMA2)/(BOLTZ*T) 
C NOW DO SUM 
      SUM1 = 0.0
      SUM2 = 0.0
      DO 996 N = 1,NITER
      SUM1 = SUM1+((X1*X1/N)+2.0*X1/(N*N)+2.0/(N*N*N))*EXP(-N*X1)
      SUM2 = SUM2+((X2*X2/N)+2.0*X2/(N*N)+2.0/(N*N*N))*EXP(-N*X2) 
 996  END DO
C NOW CALCULATE DIFFERNCE TO GET BAND RADIANCE
      SUMBB = SUM2 - SUM1
      SUMBB = SUMBB*(2.0*CSI)*(BOLTZ*T/(HPLANCK*CSI))**3
C NOW DIVIDE BY 10000 IN ORDER TO GET PHOTONS CM-2 S-1 STER-1
      SUMBB = SUMBB*1.0E-4

C*** THIS HAS BEEN CHECKED AGAINST THE SPECTRAL CALCULATOR WEBSITE 

C NOW SCALE TO RADIUS
      SUMBB = PIE*SUMBB/(RADIUS*RADIUS)

      GETGSTAR = SUMBB/SUMDRAINETOT
      WRITE (3,297) T, GETGSTAR
 
 296  FORMAT(2X,'TEMP_BB',5X,'BB/DRAINE')    
 297  FORMAT(2X,F6.1,5X,1PE11.4)        
 298  FORMAT(2X,//,'** PHOTO RATES (S-1) AT 50 R_* IN ',A,' **' 
     *  //)
 998  FORMAT (2X,'DRAINE INTERSTELLAR UV FIELD (VAN DISHOECK 1988)')
 999  FORMAT (2X, 'RATIO OF STELLAR BB TO DRAINE FIELD - INTEGRATED
     * 912-2050A'//'BLACKBODY RADIATION FIELD AT TEMPERATURE T_BB
     * AND RADIUS =',F5.1,1X,'STELLAR RADII'//)

      RETURN
      END

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO CALCULATE SELF-SHIELDED CO PHOTODISSOCIATION RATE   C 
C  USING ONE-BAND APPROXIMATION                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION GETCOR(H2COL,XCO,V,AUV)
C  H2COL = RADIAL H2 COLUMN DENSITY (OUTWARDS)
C  XCO = CO/H2 FRACTIONAL ABUNDANCE
C  V = OUTFLOW VELOCITY
C  AUV = RADIAL EXTINCTION AT 1000 ANGSTROMS

      DOUBLE PRECISION FRACE,LAMDAE,FOSCE,BANDS,GE0,TAUE,GAMMAD,
     *   BETAE,H2COL,XCO,V,AUV,EARG

C  FRACTIONAL POPULATION OF LOWER LEVEL
      FRACE = 1.0/3.0
C  EFFECTIVE WAVELENGTH (IN CM)
      LAMDAE = 1000.0*1.0E-08
C  EFFECTIVE DISSOCIATIVE OSCILLATOR STRENGTH
      FOSCE = 0.017
C  EFFECTIVE NUMBER OF BANDS
      BANDS = 1.0
C  UNSHIELDED PHOTODISSOCIATION RATE OF CO
      GE0 = 2.4E-10
C  CALCULATE EFFECTIVE OPTICAL DEPTH OF CO AT RADIUS
      TAUE = 0.0265*FRACE*FOSCE*LAMDAE*H2COL*XCO*(1.0/V)
C  CALCULATE CONTINUUM SHIELDING BY DUST (MORRIS AND JURA)
      GAMMAD = EARG(-1.644*(AUV**0.86))
C  THE FOLLOWING IS THE MORRIS/JURA APPROXIMATION TO THE FULL INTEGRAL
      BETAE = (1.0-EARG(-1.5*TAUE))/(1.5*TAUE)
C  CALCULATE CO PHOTODISSOCIATION RATE
      GETCOR = GE0*BETAE*GAMMAD*BANDS
      

      
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION FOR THE TOTAL H NUMBER DENSITY AS A FUNCTION OF RADIUS  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION GETHNR(RADIUS)
      DOUBLE PRECISION RADIUS,PI,MLOSS,MSUN,SEC,MU,MH,V,KB
      COMMON/BL2/ MLOSS,V
      COMMON/BLC/ PI,MH,MU,KB

      MSUN = 1.989E33
      SEC = 31557600.0

     
      GETHNR =(MLOSS*MSUN/SEC)/(4.0*PI*(RADIUS**2.0)*MU*MH*V)
      
      RETURN
      END
      
      

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION FOR THE TEMPERATURE AS A FUNCTION OF RADIUS  C
C  (FIT TO CROSAS & MENTEN 1997), MIN = 10 K             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION GETTEM(R)
      DOUBLE PRECISION R,R_STAR,T_STAR,EPSIL,R_EPSIL
      CHARACTER*500 TEMPMODE
      COMMON/BLTEMP/ R_STAR,T_STAR,EPSIL,TEMPMODE,R_EPSIL
      
C ORIGINAL TEMPERATURE PROFILE      
      IF (TEMPMODE.EQ.'EPSILON') THEN
      GETTEM = T_STAR*(R/R_EPSIL)**(-EPSIL)
      ENDIF 
      
      IF (GETTEM.LT.10.0) GETTEM = 10.0

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION FOR TOTAL H RADIAL COLUMN DENSITY TO OUTER EDGE OF ENVELOPE  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION GETH2(RADIUS,HNR)
      DOUBLE PRECISION RADIUS,HNR

      GETH2 = HNR*RADIUS

      RETURN
      END


      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION FOR EXTINCTION AT 1000 ANGSTROMS (NEJAD AND MILLAR)   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION GETAUV(H2COL,RADIUS)
C  H2COL = RADIAL TOTAL H COLUMN DENSITY (OUTWARDS)
      CHARACTER*500 CLUMPMODE,TEMPMODE
      DOUBLE PRECISION H2COL,AUV_AV,GAMCO,ALBEDO,ZETA,A_G,X_G,Y(1000),
     *    GSTAR,RSCALE,RDUST,DELTA_AUV,AUV,K,H,A,B,C,RADIUS,TEFF,
     *    TEFF_R,TEFF_INF,R_STAR,PI,MH,MU,KB,T_STAR,EPSIL,
     *    FVOL,L,FIC,BFILL,CFILL,R,R_EPSIL,GBIN,RBIN,TBIN
      INTEGER ICO,ICLUMP,ISTELLAR,IBIN
      COMMON/BL3/ Y,X_G,A_G,ZETA,ALBEDO,GAMCO,AUV_AV,ICO,GSTAR,GBIN,
     *   RSCALE,RDUST,DELTA_AUV,ISTELLAR,IBIN,RBIN,TBIN
      COMMON/BLC/ PI,MH,MU,KB
      COMMON/CLM/ CLUMPMODE,ICLUMP,FVOL,L,FIC
      COMMON/BLTEMP/ R_STAR,T_STAR,EPSIL,TEMPMODE,R_EPSIL

C  COLUMN DENSITY H2/A(V) = 1.87E21 ATOMS CM-2 MAG-1 (ISM VALUE)
C  http://adsabs.harvard.edu/abs/1995A&A...293..889P
C  POROSITY DESCRIPTION 
C  https://ui.adsabs.harvard.edu/abs/2018A%26A...616A.106V/abstract

C  WHEN USING FRACTIONAL ABUNDANCE WRT H2
C     AUV = (AUV_AV*2.0*H2COL)/1.87E21
      AUV = (AUV_AV*H2COL)/1.87E21

      IF (CLUMPMODE.EQ.'POROSITY') THEN            
      R = RADIUS      
      H = L*(R_STAR**(-2.0/3.0))/FVOL 
      BFILL = 1.0-(1.0-FVOL)*FIC  
C DEFINE FOUR TERMS FOR USE IN CALCULATING TAU_EFF/AUV_EFF
      CFILL = AUV*H*BFILL
      F1 = DATAN(1.0 + (SQRT(2.0)*(R**(1.0/12.0))*((CFILL)**(-0.25))))
      F2 = DATAN(1.0 - (SQRT(2.0)*(R**(1.0/12.0))*((CFILL)**(-0.25))))
      F3 = 0.5*DLOG(DSQRT(CFILL*R)+SQRT(2.0)*((CFILL)**0.25)*
     *  (R**(7.0/12.0)) + R**(2.0/3.0))
      F4 = 0.5*DLOG(DSQRT(CFILL*R)-SQRT(2.0)*((CFILL)**0.25)*
     *  (R**(7.0/12.0)) + R**(2.0/3.0))
      F5 = (3.0/(2.0*SQRT(2.0)))*((AUV*R)**(0.25))*
     *   ((H*BFILL)**(-0.75))*(1.0-FIC)
C  INTERCLUMP
      GETAUV = AUV*FIC + F5*(PI-F1+F2-F3+F4)
      
      ELSE
         GETAUV = AUV
      END IF
      
      
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO CALCULATE RADIATION FIELD STRENGTH AT 1000 A  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION GETRAD(AUV)
C  AUV = RADIAL EXTINCTION AT 1000 ANGSTROMS
C  
      CHARACTER*500 CLUMPMODE,TEMPMODE
      DOUBLE PRECISION AUV,CETA(6),GA(6),W(6),PI,SUM,GAMMA,TAU,
     * TCL,TEFF,EARG,RADIUS,FVOL,L,H2COL,HNR,KAPPA,GETH2,GETHNR
      DOUBLE PRECISION MH,MU,KB,CORR,H,SUM2,B,FIC,
     * T_STAR,EPSIL,R_STAR,K,A,C,TEFF_INF,TEFF_PAR,TEFF_R,
     * AUV_AV,GAMCO,ALBEDO,ZETA,A_G,X_G,Y(1000),MLOSS,V,GSTAR,
     *   RSCALE,RDUST,DELTA_AUV,R_EPSIL,GBIN
      INTEGER I,ICLUMP,ICO,ISTELLAR,IBIN

C  AVERAGES OVER A SPHERICAL OUTFLOW FOLLOWING MORRIS & JURA 1983
C  https://articles.adsabs.harvard.edu/pdf/1983ApJ...264..546M

      PI = 4.*ATAN(1.)

C   SET COEFFICIENTS FOR CALCULATION OF RADIATION FIELD
      W(1) = 0.17132449
      W(2) = 0.36076157
      W(3) = 0.46791393
      W(4) = W(1)
      W(5) = W(2)
      W(6) = W(3)
      GA(1) = 0.93246951
      GA(2) = 0.66120939
      GA(3) = 0.23861919
      GA(4) = -GA(1)
      GA(5) = -GA(2)
      GA(6) = -GA(3)

      SUM = 0.0
      DO I=1,6
      CETA(I) = (PI*GA(I)+PI)/2.0
      SUM=SUM+(W(I)*(SIN(CETA(I))*EXP((-AUV*CETA(I))/SIN(CETA(I)))))
      END DO
      SUM = (PI/4.0)*SUM
      
      GETRAD = SUM
                        
      RETURN
      END


      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ALTERNATIVE EXP FUNCTION TO AVOID EXPONENTIAL ERROR  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION FUNCTION EARG(TZZ)
      DOUBLE PRECISION TZZ

      IF(TZZ.LE.-100.0) EARG = 0.0
      IF(TZZ.GT.-20.0.AND.TZZ.LT.710.) EARG = EXP(TZZ)
      IF(TZZ.GE.710.) EARG = 8D307

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   DUMMY ROUTINE TO SATISFY LOADER   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PEDERV(N,T,Y,PD,N0)

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   ERROR HANDLING ROUTINE  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DIE(ERROR)
      CHARACTER*59 ERROR

      WRITE(*,*) "ERROR: ",ERROR
      WRITE(*,*) "...Aborting......................."

      STOP
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  BEGIN NUMERICAL SUBROUTINES REQUIRED BY DVODE  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
C
C---
C
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
C
C---
C
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
C
C---
C
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
C
C---
C
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
C
C---
C
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
C
C---
C
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
C
C---
C
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
C
C---
C
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  FUNCTION TO READ RATE FILE AND STORE REQUIRED DATA  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE READR(FRATES,URATES)
C  FRATES = RATE FILE NAME
C  URATES = RATE FILE UNIT NUMBER

      CHARACTER*500 FRATES,SDATA,TOKEN(100)
      CHARACTER*12 RE1,RE2,PR1,PR2,PR3,PR4,STAR
      CHARACTER*1 SEPC,SEP
      DOUBLE PRECISION ALF(10000,5),BET(10000,5),GAM(10000,5),
     *   TLOWER(10000,5),TUPPER(10000,5)
      INTEGER RTYPE(10000),NTR(10000),NREAC,I,J,K,L,URATES,ICOR,
     *   ISTRIN,ICOIP,ICOAP

C  COMMON BLOCK FOR REACTION DATA STORAGE
      COMMON/BL4/ ALF,BET,GAM,TLOWER,TUPPER,RTYPE,NTR,NREAC,ICOR,
     *   ICOIP,ICOAP

C  TOKEN (DATA SUBSTRING SEPARATOR)
      SEPC = ':'
      SEP = SEPC

C  OPEN RATE FILE
      OPEN(UNIT=URATES, FILE=FRATES, STATUS='OLD')

C  OPEN A NEW TEXT FILE TO WRITE THE BASIC REACTION DATA INTO
      OPEN(UNIT=10, FILE='rates.dat')

C  DATA FILE READ LOOP
      J = 0
 100  J = J+1
      READ(URATES,101,END=300) SDATA
 101  FORMAT(A)

C  SEARCH THROUGH DATA STRING (SDATA) AND FIND N SUB-STRINGS (TOKEN(I))

C  INITIALISE TOKEN COUNTER
      I = 0

C  BEGIN TOKENISATION OF THE DATA STRING 
 200  I = I + 1
C  TOKEN IS A STRING (DELIMITED BY ""), WHICH MAY CONTAIN THE SEPARATOR
         IF (SDATA(1:1).EQ.'"') THEN
            SDATA = SDATA(2:)
            SEP = '"'
            ISTRIN = 1
         END IF
C  FIND TOKEN POSITION
         L = INDEX(SDATA,SEP) - 1
C  POSITIVE LENGTH TOKEN
         IF (L.GT.0) THEN
            TOKEN(I) = SDATA(:L)
            L = L + 1
C  ZERO-LENGTH TOKEN
         ELSE IF (L.EQ.0) THEN
            TOKEN(I) = ""
            L = 1
C  NO MORE TOKENS FOUND         
         ELSE IF (L.EQ.-1) THEN
            TOKEN(I) = SDATA(1:)
            L = 0
         END IF
C  REMOVE PREVIOUS TOKEN FROM DATA STRING
         SDATA = SDATA(L+1:)
C  IF TOKEN WAS A STRING THEN REMOVE THE FOLLOWING CHARACTER AS WELL
         IF (ISTRIN.EQ.1) THEN
            ISTRIN = 0
            SEP = SEPC
            SDATA = SDATA(2:)
         END IF        
C  CONTINUE LOOPING UNTIL NO MORE TOKENS REMAIN
      IF (I.LE.100.AND.L.GT.0) GO TO 200 

C  IF LESS THAN 14 TOKENS WERE READ, REACTION DATA IS INCOMPLETE
      IF(I.LT.14) CALL DIE
     *   ("Unable to interpret reaction string -- check rate file     ")


C  STORE INDEX, REACTANTS AND PRODUCTS
      READ(TOKEN(1),*) I
      RE1 = TOKEN(3)
      RE2 = TOKEN(4)
      PR1 = TOKEN(5)
      PR2 = TOKEN(6)
      PR3 = TOKEN(7)
      PR4 = TOKEN(8)
      STAR = TOKEN(18)
C  NUMBER OF TEMPERATURE RANGES
      READ(TOKEN(9),*) NTR(J)

C  LOOP THROUGH AND STORE REACTION DATA FOR DIFFERENT TEMPERATURE RANGES
      DO K = 1,NTR(J)
         L = K - 1
         READ(TOKEN(10+9*L),*) ALF(J,K)
         READ(TOKEN(11+9*L),*) BET(J,K)
         READ(TOKEN(12+9*L),*) GAM(J,K)
         READ(TOKEN(13+9*L),*) TLOWER(J,K)
         READ(TOKEN(14+9*L),*) TUPPER(J,K)
      END DO

C  DETERMINE REACTION TYPE:

C  DEFAULT BINARY REACTION
      RTYPE(J) = 0
C  COSMIC RAY REACTION
      IF(RE2.EQ.'CRP') RTYPE(J) = 1
C  COSMIC RAY-INDUCED PHOTON REACTION
      IF(RE2.EQ.'CRPHOT') RTYPE(J) = 2
C  PHOTON REACTION
      IF(RE2.EQ.'PHOTON') RTYPE(J) = 3
C  INTERNAL PHOTONS - CROSS-SECTIONS NOT AVAILABLE
      IF(RE2.EQ.'INPHOTON') RTYPE(J) = 4
C  INTERNAL PHOTONS - CROSS-SECTIONS AVAILABLE
      IF(RE2.EQ.'INPHOTON'.AND.STAR.EQ.'CWLeo') RTYPE(J) = 5
C  INTERNAL PHOTONS - CROSS-SECTIONS AVAILABLE - IONISATION RATES
      IF(RE2.EQ.'INPHOTON'.AND.STAR.EQ.'Millar 2018') RTYPE(J) = 5
C  INTERNAL DISK PHOTONS - CROSS-SECTIONS NOT AVAILABLE
      IF(RE2.EQ.'ACPHOTON') RTYPE(J) = 6
C  INTERNAL PHOTONS - CROSS-SECTIONS AVAILABLE
      IF(RE2.EQ.'ACPHOTON'.AND.STAR.EQ.'CWLeo') RTYPE(J) = 7
C  INTERNAL PHOTONS - CROSS-SECTIONS NOT AVAILABLE - IONISATION RATES
      IF(RE2.EQ.'ACPHOTON'.AND.STAR.EQ.'Millar 2018') RTYPE(J) = 8

C  DETERMINE INDEX OF CO PHOTODISSOCIATION REACTION
      IF(RE1.EQ.'CO'.AND.RE2.EQ.'PHOTON') ICOR = J
C
C IDENTIFY REACTION NUMBER OF CO + INPHOTON
      IF(RE1.EQ.'CO'.AND.RE2.EQ.'INPHOTON') ICOIP = J
C
C IDENTIFY REACTION NUMBER OF CO + ACPHOTON
      IF(RE1.EQ.'CO'.AND.RE2.EQ.'ACPHOTON') ICOAP = J

C  WRITE FORMATTED OUTPUT (ONLY WRITES DATA FOR FIRST TEMP. RANGE)
      WRITE(10,1) I,RE1,RE2,PR1,PR2,PR3,PR4,ALF(J,1),BET(J,1),GAM(J,1)
 1    FORMAT(I4,1X,A12,A12,A12,A12,A5,A5,1PE8.2,0PF7.1,0PF10.1)

      GO TO 100
 300  CONTINUE

      NREAC = J - 1

      WRITE(*,*) 'Reactions read: ',NREAC

C  CATCH BASIC RATE FILE ERROR
      IF(NREAC.LT.1) CALL DIE
     *   ("No reactions defined -- check rate file                    ")

      CLOSE(10)
      CLOSE(URATES)

      RETURN
      END
