    PROGRAM BUDGET
    
    USE MOD_CKLIBs
    USE MOD_AJLIBs
    IMPLICIT NONE
    
        
    INTEGER               ::  M, PARTICLE, Field_X, Field_T, READSTATE, L, Ns, K
    
    REAL(KDPC)            ::  Rho, P, r, U, T, S, dx
    
    REAL(KDPC)            ::  PRPR, PUPR, PPPR, PTPR, PSPR, PP2PR, PU2PR, GAMA, Ct
    
    REAL(KDPC),ALLOCATABLE::  CVMS(:), CPMS(:), HMS(:), UMS(:), RHOK(:), PRHOKPR(:), dCVMS(:)
    
    REAL(KDPC),ALLOCATABLE::  C(:), Cdot(:), Ddot(:), Wdot(:), dCdot(:,:), dDdot(:,:), dWdot(:,:)
    
    REAL(KDPC)            ::  V(1,33), X(1,100)
    
    CHARACTER(100)  ::  OPath1,OPath2, FILENAME1, FNAME1, FILENAME2, FNAME2, FMT1, FMT2
    
    CHARACTER(5000)  ::  VARs
    
    CHARACTER(100),ALLOCATABLE  ::  SPNAME(:)
    
    CHARACTER(1)    ::  Slash
    
    LOGICAL         ::  KERR
    
    REAL            ::  Null1, Null2
    
    Ns = 10
    
    dx = 5E-4
    
    ALLOCATE( CVMS(Ns), dCVMS(Ns), CPMS(Ns), HMS(Ns), UMS(Ns),RHOK(Ns), PRHOKPR(Ns), SPNAME(Ns) )
    ALLOCATE( C(Ns), Cdot(Ns), Ddot(Ns), Wdot(Ns), dCdot(Ns,Ns+2), dDdot(Ns,Ns+2), dWdot(Ns,Ns+2))
    
    SPNAME(1)       = 'H2'
    SPNAME(2)       = 'H'
    SPNAME(3)       = 'O'
    SPNAME(4)       = 'O2'
    SPNAME(5)       = 'OH'
    SPNAME(6)       = 'H2O'
    SPNAME(7)       = 'HO2'
    SPNAME(8)       = 'H2O2'
    SPNAME(9)       = 'N2'
    SPNAME(10)      = 'AR'
    
    CALL CKINIT('chem.ASC')
    CALL AJINIT
    
    M = 1
    
    Slash = '\'
    
    Opath1 = 'E:\ ÛƒÍ∫ÆºŸ\C1_1.82_Track\C1_1.82_Track\C1_1.82_Track\FIELD'
    
    OPath2 = 'E:\ ÛƒÍ∫ÆºŸ\C1_1.82_Track\C1_1.82_Track\C1_1.82_Track\BUDGETDATA'
    
    WRITE(FMT1,'(A,I3,A)')'(', 33,'(1P,E22.15,2X) )' 
    
    WRITE(FMT2,'(A,I3,A)')'(', 100,'(1P,E22.15,2X) )' 
    
    DO Field_T = 1, 141
        
        WRITE(FILENAME1,'(A,I4.4,A)')'FIELD', Field_T, '.plt'        
        FNAME1 = trim(OPath1)//Slash//trim(FILENAME1)       
        OPEN(1, FILE=trim(FNAME1), ACTION='READ')
        READ(1,*)
        
        WRITE(FILENAME2,'(A,I4.4,A)')'BUDGET', Field_T, '.plt'        
        FNAME2 = trim(OPath2)//Slash//trim(FILENAME2)       
        OPEN(2, FILE=trim(FNAME2), ACTION='WRITE')
        
        WRITE(VARs,'(A)') 'Variables = "DRhoDt" "DUDt" "DTDt1" "DTDt2" "DTDt" "DPDt1" "DPDt2" "DPDt"'
        
        DO K=1,Ns
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DRho'//trim(SPNAME(K))//'Dt1"'
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DRho'//trim(SPNAME(K))//'Dt2"'
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DRho'//trim(SPNAME(K))//'Dt"'
        END DO
        ! 1 and 2 means the reaction contribution and divergence contribution
        
        WRITE(VARs,'(A)')trim(VARs)//' "DPUPrDt1" "DPUPrDt2" "DPUPrDt3" "DPUPrDt"'
        WRITE(VARs,'(A)')trim(VARs)//' "DPRhoPrDt1" "DPRhoPrDt2" "DPRhoPrDt"'
        
        DO K=1,Ns
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DPRho'//trim(SPNAME(K))//'PrDt1"'
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DPRho'//trim(SPNAME(K))//'PrDt2"'
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DPRho'//trim(SPNAME(K))//'PrDt3"'
            WRITE(VARs,'(A)')trim(VARs)//' '//'"DPRho'//trim(SPNAME(K))//'PrDt"'
        END DO
        
        WRITE(VARs,'(A)')trim(VARs)//' "DPTPrDt1" "DPTPrDt2" "DPTPrDt3" "DPTPrDt4" "DPTPrDt5" "DPTPrDt"'
        
        WRITE(VARs,'(A)')trim(VARs)//' "DPPPrDt1" "DPPPrDt2" "DPPPrDt3" "DPPPrDt4" "DPPPrDt"'
        
        WRITE(VARs,'(A)')trim(VARs)//' "P" "U" "T" "r"'
        
        WRITE(2,'(A)')trim(VARs)
        
        READSTATE   = 0
        
        DO WHILE (READSTATE .EQ. 0)
            
            READ(1,*,IOSTAT=READSTATE)   V
            
            r           = V(1,1)
            Rho         = V(1,2)
            U           = V(1,3)
            P           = V(1,4)
            T           = V(1,5)
            RHOK        = V(1,6:15)
            PRPR        = V(1,16)
            PUPR        = V(1,17)
            PPPR        = V(1,18)
            PTPR        = V(1,19)
            PSPR        = V(1,20)
            PRHOKPR     = V(1,21:30)
            PP2PR       = V(1,31)
            PU2PR       = V(1,32)
            S           = V(1,33)
            
            CALL CKUMS ( T,UMS )                                                    ! Species internal energy in mass units
            CALL CKHMS ( T,HMS )                                                    ! Species enthropy in mass units
            CALL CKCVMS( T,CVMS )                                                   ! Species specific heats at constant volume in mass units
            CALL CKCPMS( T,CPMS )                                                   ! Species specific heats at constant pressure in mass units
            CALL CKRY2C( Rho, RHOK/Rho, C )                                         ! Given density, mass fractions, return molar concentrations
            CALL CKWC( T, C, Ct, Cdot, Ddot, Wdot )                                 ! Reaction sources in mole unit
            CALL AJ_dWdot_dZ( T, C, Ct, Cdot, Ddot, Wdot, dCdot, dDdot, dWdot )     ! Jacobians of Reaction Rates
            CALL AJ_dCVMS_dT( T, CVMS, dCVMS )                                      ! To Calculate  d(CVMS)/dT
            
            GAMA            = SUM(CPMS*RHOK)/SUM(CVMS*RHOK)                                 ! Specific heat
            
            X(1,1)          = -Rho*S                                                        ! DRho      Dt
            
            X(1,2)          = -PPPR/Rho                                                     ! DU        Dt
            
            X(1,3)          = -P*S/SUM(RHOK*CVMS)                                           ! DT        Dt1
            X(1,4)          = SUM(Wdot*WT*UMS)/SUM(RHOK*CVMS)                               ! DT        Dt2
            X(1,5)          = X(1,3)+X(1,4)                                                 ! DT        Dt
            
            X(1,6)          = -GAMA*P*S                                                     ! DP        Dt1
            X(1,7)          = SUM((HMS-GAMA*UMS)*Wdot*WT)                                   ! DP        Dt2
            X(1,8)          = X(1,6)+X(1,7)                                                 ! DP        Dt
            
            DO K            = 1,Ns             
                X(1,6+3*K)    = -RHOK(K)*S                                                  ! DRhok     Dt1  
                X(1,7+3*K)    = Wdot(K)*WT(K)                                               ! DRhok     Dt2
                X(1,8+3*K)    =X(1,6+3*K)+X(1,7+3*K)                                        ! DRhok     Dt
            END DO
            
            X(1,39)         = -PUPR**2                                                      ! DPUPR     Dt1
            X(1,40)         = PPPR*PRPR/Rho**2                                              ! DPUPR     Dt2
            X(1,41)         = -PP2PR/Rho                                                    ! DPUPR     Dt3
            X(1,42)         = SUM(X(1,39:41))                                               ! DPUPR     Dt
            
            X(1,43)         = -PRPR*(S+PUPR)                                                ! DPRPR     Dt1
            X(1,44)         = -Rho*PSPR                                                     ! DPRPR     Dt2
            X(1,45)         = X(1,43)+X(1,44)                                               ! DPRPR     Dt
            
            DO K            = 1,Ns
                X(1,42+4*K) = -PRHOKPR(K)*(S+PUPR)                                          ! DPRHOKPR   Dt1
                X(1,43+4*K) = -RHOK(K)*PSPR                                                 ! DPRHOKPR   Dt2
                X(1,44+4*K) = WT(K)*(dWdot(K,KK+2)*PTPR + dWdot(K,K)/WT(K)*PRHOKPR(K))      ! DPRHOKPR   Dt3
                X(1,45+4*K) = SUM(X(1,43+3*K:45+3*K))                                       ! DPRHOKPR   Dt
            END DO
            
            X(1,86)         = -PPPR*S/SUM(RHOK*CVMS)                                        ! DPTPR     Dt1
            X(1,87)         = -P*PSPR/SUM(RHOK*CVMS)                                        ! DPTPR     Dt2
            X(1,88)         = 0
            
            DO K = 1,Ns
                X(1,88)         = X(1,88)-WT(K)*WDot(K)*CVMS(K)*PTPR - WT(K)*UMS(K)*(dWdot(K,KK+2)*PTPR + dWdot(K,K)/WT(K)*PRHOKPR(K))
            ENDDO           
                                                                                            ! DPTPR     Dt3
            X(1,89)         = -X(1,5)*SUM(CVMS*PRHOKPR+RHOK*PTPR*dCVMS)/SUM(RHOK*CVMS)      ! DPTPR     Dt4
            X(1,90)         = -PUPR*PTPR                                                    ! DPTPR     Dt5
            X(1,91)         = SUM(X(1,86:90))                                               ! DPTPR     Dt
            
            X(1,92)         = SUM(RHOK*RU/WT)*X(1,81)                                       ! DPPPR     Dt1
            X(1,93)         = PTPR*SUM(RU/WT*X(1,9:38:3))                                   ! DPPPR     Dt2
            X(1,94)         = X(1,5)*SUM(RU/WT*PRHOKPR)                                     ! DPPPR     Dt3
            X(1,95)         = T*SUM(RU/WT*X(1,50:85:4))                                     ! DPPPR     Dt4
            X(1,96)         = SUM(X(1,92:95))                                               ! DPPPR     Dt
            
            X(1,97)         = P
            X(1,98)         = U
            X(1,99)         = T
            X(1,100)        = r

            
            WRITE(2,FMT2)   X
            
        ENDDO
        
        CLOSE(1)
        CLOSE(2)
        
    ENDDO

    END PROGRAM BUDGET

