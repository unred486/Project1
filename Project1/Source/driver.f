      PROGRAM DRIVER                                                            
      use var
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)                                  
                                                                               
C          LENLWK allocates the logical working space                           
C          LENIWK allocates the integer working space                           
C          LENRWK allocates the real working space                              
C          LENCWK allocates the character working space                         
C          LENSYM is the length of a character string                    

C      PARAMETER (LENLWK=4500,LENIWK=1000000,LENRWK=2000000000,
C     1           LENCWK=400, LENSYM=16)                                                     
      PARAMETER (LENLWK=4500,LENIWK=7000000,LENRWK=200000000,
     1           LENCWK=4000, LENSYM=16)                                                     

C                                                                               
C          NMAX is the total number of grid points allowed                      
      PARAMETER ( NMXROW=366, NMXCOL=10, NORD = 20)
      INTEGER :: NMAX =1000
C                                                                               
C          All storage needed by the flame code is allocated in the             
C          following three arrays.  LWORK is for logical variables,             
C          IWORK is for integer variables, RWORK is of real variables,          
C          CWORK is for character variables. WEG if of real variables,
C          IWEG is of integer variables                          
C                                                                               
      LOGICAL :: LWORK(LENLWK)
      DIMENSION IWORK(LENIWK)
      DIMENSION RWORK(LENRWK)                     
      CHARACTER CWORK(LENCWK)*(LENSYM)
c       COMMON /RI/ RWORK(LENRWK), IWORK(LENIWK)
C                                                                               
C          LIN is the unit from which user input is read                        
C          LOUT is the unit to which printed output is written                  
C          LRIN is the unit from which the restart file is read                 
C          LROUT is the unit to which the save file is written                  
C          LRCRVR is the unit to which the scratch save file is written         
C          LINKCK is unit from which the Chemkin Linking file is read           
C          LINTP is unit from which the Transport Linking file is read          
      DATA LIN/7/, LOUT/6/, LRIN/14/, LROUT/15/, LRCRVR/16/,                    
     1     LINKCK/25/, LINKTP/35/, LH2O/17/, LCO2/18/, LCO/19/,
     2     LRATE/20/, LTESTT/40/ 

      character(255), parameter :: sphdifid =
     1 '$Id: driver.f 10 2019-06-05 20:17:05Z hlee $'
      character(255), parameter :: sphdifrev=
     1 '$Revision: 10 $'
      character(255), parameter :: sphdifdate=
     1 '$Date: 2019-06-05 15:17:05 $'

 
C
      WRITE(*,10) 
10    FORMAT (
     1/'         ##################################################',
	2/'',
	3/'           SPHERICAL DIFFUSION FLAME CODE IS RUNNING',
	4/'               PLEASE DO NOT CLOSE THIS WINDOW !!',
     5/'',
	6/'         ##################################################')

C   
      WRITE(*,*) ""
      WRITE(*,*) "Version number, version #:", trim(sphdifrev)
      WRITE(*,*) "Version commited on :", trim(sphdifdate)
      WRITE(*,*) ""      
      WRITE(*,*) ""
C            open the user input file
      OPEN(UNIT=LIN, FORM='FORMATTED', FILE='sphdif.in')
C            open the printed output file
      OPEN(UNIT=LOUT,FORM='FORMATTED', FILE='sphdif.out')  
C            open the restart file
      OPEN(UNIT=LRIN, FORM='UNFORMATTED', 
     +     CONVERT = 'BIG_ENDIAN',
     +     FILE='restart')
C            open the save output file
      OPEN(UNIT=LROUT, FORM='UNFORMATTED',     
     +     CONVERT='BIG_ENDIAN',
     +     FILE='save')
C            open the recover output file
      OPEN(UNIT=LRCRVR, FORM='UNFORMATTED',
     +     CONVERT='BIG_ENDIAN',
     +     FILE='recover')
C            open the Chemkin Link file
      OPEN(UNIT=LINKCK, FORM='UNFORMATTED', 
     +     CONVERT='BIG_ENDIAN',
     +     FILE='cklink')
C            open the Transport Link file
      OPEN(UNIT=LINKTP, FORM='UNFORMATTED', 
     +     CONVERT='BIG_ENDIAN',
     +     FILE='tplink')
      OPEN(UNIT=LTESTT,FORM='FORMATTED', ACTION= 'WRITE', 
     +     STATUS='REPLACE', FILE='test.out') 

      CALL SPHDIFB (NMAX, LIN, LOUT, LINKCK, LINKTP, LRIN, LROUT,                
     1             LRCRVR, LENLWK, LWORK, LENIWK, IWORK, LENRWK,                
     2             RWORK,LENCWK, CWORK, LRATE, NMXROW, NMXCOL, NORD,
     3             LH2O, LCO2, LCO,LTESTT)   
                                    
      STOP                                                                      
      END
