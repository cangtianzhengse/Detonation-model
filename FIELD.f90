!  FIELD.f90 
!
!  FUNCTIONS:
!  FIELD - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: FIELD
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    PROGRAM FIELD

    IMPLICIT NONE
    
    
    INTEGER         ::  M, PARTICLE, Field_X, Field_T, READSTATE, L, Ns, K
    
    REAL            ::  Rho, P, r, U, T, S, dx
    
    REAL            ::  PRPR, PUPR, PPPR, PTPR, PSPR, PP2PR, PU2PR, Ind, GAMA
    
    REAL,ALLOCATABLE::  CVMS(:), CPMS(:), HMS(:), UMS(:), PRODUCTION(:), Mlw(:), PRHOPR(:)
    
    REAL            ::  V(1,16), W1(1,16), W2(1,16), E1(1,16), E2(1,16), X(1,33)
    
    CHARACTER(100)  ::  OPath1, OPath2, FILENAME1, FNAME1, FILENAME2, FNAME2, FMT1, FMT2
    
    CHARACTER(500)  ::  VARs
    
    CHARACTER(100),ALLOCATABLE  ::  SPNAME(:)
    
    CHARACTER(1)    ::  Slash
    
    Ns = 10
    
    dx = 5E-4
    
    ALLOCATE( CVMS(Ns), CPMS(Ns), HMS(Ns), UMS(Ns), PRODUCTION(Ns), Mlw(Ns), PRHOPR(Ns), SPNAME(Ns) )
    
    SPNAME(1) = 'H2'
    SPNAME(2) = 'H'
    SPNAME(3) = 'O'
    SPNAME(4) = 'O2'
    SPNAME(5) = 'OH'
    SPNAME(6) = 'H2O'
    SPNAME(7) = 'HO2'
    SPNAME(8) = 'H2O2'
    SPNAME(9) = 'N2'
    SPNAME(10) = 'AR'
    
    M = 1
    
    Slash = '\'
    
    Opath1 = 'E:\ ÛƒÍ∫ÆºŸ\C1_1.82_Track\C1_1.82_Track\C1_1.82_Track\CONDENSED'
    
    OPath2 = 'E:\ ÛƒÍ∫ÆºŸ\C1_1.82_Track\C1_1.82_Track\C1_1.82_Track\FIELD(3 point)'
    
    WRITE(FMT1,'(A,I3,A)')'(', 16,'(1P,E22.15,2X) )'
    
    WRITE(FMT2,'(A,I3,A)')'(', 33,'(1P,E22.15,2X) )' 
    
    DO Field_T = 1, 141
        
        WRITE(FILENAME1,'(A,I4.4,A)')'Field_T', Field_T, '.plt'        
        FNAME1 = trim(OPath1)//Slash//trim(FILENAME1)        
        OPEN(1, FILE=trim(FNAME1), ACTION='READ')
        
        WRITE(FILENAME2,'(A,I4.4,A)')'FIELD', Field_T, '.plt'        
        FNAME2 = trim(OPath2)//Slash//trim(FILENAME2)       
        OPEN(2, FILE=trim(FNAME2), ACTION='WRITE')
        
        WRITE(VARs,'(A)') 'Variables = "r" "Rho" "U" "P" "T"'
        DO K=1,Ns
            WRITE(VARs,'(A)')trim(VARs)//' '//'"'//trim(SPNAME(K))//'"'
        END DO
        WRITE(VARs,'(A)')trim(VARs)//' "PRhoPr" "PUPr" "PPPr" "PTPr" "PSPr"'
        DO K=1,Ns
            WRITE(VARs,'(A)')trim(VARs)//' '//'"P'//trim(SPNAME(K))//'Pr"'
        END DO
        WRITE(VARs,'(A)')trim(VARs)//' "PP2Pr" "PU2Pr" "S"'
        WRITE(2,'(A)')trim(VARs)
        
        READ(1,*)
        
        Ind = 1
        
        L = 0
        
        READSTATE = 0
        
        DO WHILE ((Ind .NE. 2.0265E+05))
            
            L = L + 1
            
            IF (L .EQ. 1) THEN             
                !READ(1,FMT1)    E2  
                !READ(1,FMT1)    E1
                !READ(1,FMT1)    V
                !READ(1,FMT1)    W1
                !READ(1,FMT1)    W2
                READ(1,FMT1)    E1
                READ(1,FMT1)    V
                READ(1,FMT1)    W1
            END IF
            
            IF (L .GT. 1) THEN              
                !E2  = E1
                !E1  = V
                !V   = W1
                !W1  = W2
                !READ(1,FMT1)    W2
                E1  = V
                V   = W1
                READ(1,FMT1)    W1
                
                !DO WHILE (W2(1,1) .EQ. W1(1,1))
                !    
                !    READ(1,FMT1,IOSTAT=READSTATE)    W2
                !    
                !END DO
                
                DO WHILE (W1(1,1) .EQ. V(1,1))
                    
                    READ(1,FMT1,IOSTAT=READSTATE)    W1
                    
                END DO
                
            END IF          
            
            !Ind     = W2(1,4)
            Ind     = W1(1,4)
            
            X(1,1:5) = V(1,1:5)
            ! X, Rho, U, P, T
            X(1,6:15) = V(1,7:16)
            ! Species Concentrention
            !X(16:20),X(21:30),X(31:32)
            !1st derivative: Rho, U, P, T, S
            !ist derivative: Rhok
            !2nd derivative: P, U
            
            r   = V(1,1)
            Rho = V(1,2)
            U   = V(1,3)
            P   = V(1,4)
            T   = V(1,5)
            
            !PRPR    = (W2(1,2)-8*W1(1,2)+8*E1(1,2)-E2(1,2))/(12*dx)
            !PUPR    = (W2(1,3)-8*W1(1,3)+8*E1(1,3)-E2(1,3))/(12*dx)
            !PPPR    = (W2(1,4)-8*W1(1,4)+8*E1(1,4)-E2(1,4))/(12*dx)
            !PTPR    = (W2(1,5)-8*W1(1,5)+8*E1(1,5)-E2(1,5))/(12*dx)
            PRPR    = (E1(1,2)-W1(1,2))/(2*dx)
            PUPR    = (E1(1,3)-W1(1,3))/(2*dx)
            PPPR    = (E1(1,4)-W1(1,4))/(2*dx)
            PTPR    = (E1(1,5)-W1(1,5))/(2*dx)
            S       = PUPR + M*U/r
            !PP2PR   = (-W2(1,4)+16*W1(1,4)-30*V(1,4)+16*E1(1,4)-E2(1,4))/(12*dx**2)
            !PU2PR   = (-W2(1,3)+16*W1(1,3)-30*V(1,3)+16*E1(1,3)-E2(1,3))/(12*dx**2)
            PP2PR   =(W1(1,4)-2*V(1,4)+E1(1,4))/dx**2
            PU2PR   =(W1(1,3)-2*V(1,3)+E1(1,3))/dx**2
            PSPR    = PU2PR + M*PUPR/r - M*U/r**2
            
            DO K = 1,Ns
                ! PRHOPR data
                !PRHOPR(K) = (W2(1,6+K)-8*W1(1,6+K)+8*E1(1,6+K)-E2(1,6+K))/(12*dx)  
                PRHOPR(K) = (E1(1,6+K)-W1(1,6+K))/(2*dx)  
            ENDDO
            
            X(1,16)   = PRPR
            X(1,17)   = PUPR
            X(1,18)   = PPPR
            X(1,19)   = PTPR
            X(1,20)   = PSPR
            X(1,21:30)= PRHOPR
            X(1,31)   = PP2PR
            X(1,32)   = PU2PR
            X(1,33)   = S
            
            WRITE(2,FMT2)   X(1,:)
            

            
        END DO
        
        CLOSE(1)
        CLOSE(2)
        
    END DO  
                  
    END PROGRAM FIELD

