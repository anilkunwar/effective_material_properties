    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Conductivity of bulk solid copper fitted as a function of temperature
    ! (kth_cu)solid = Ascu*(T-250)^2 + Bscu*(T-250) + Dscu, where Ascu = 0.00004 W/mK4 and Bscu = -0.10639 W/mK2 and Dscu= 405.75697 W/mK (250.0 < T < 1358.0 K)
    ! Reference for (kth_cu)solid: C.Y. Ho, R. W. Powell and P. E. Liley, Thermal Conductivity of the Elements: A Comprehensive Review,  Journal of Physical and Chemical Reference Data
    ! Volume 3 (1974) 1-756
    ! Conductivity of gaseous air fitted as a function of temperature
    ! Expression 1; (kth_air)gas = Agas*T^2 + Bgas*T    where Agas = 0.00012 W/mK3  and  Bgas=  -0.001557 W/mK2 (250.0 < T < 1358.0 K) 
    ! Reference for (kth_air)gas: C.Y. Ho, R. W. Powell and P. E. Liley, Thermal Conductivity of the Elements: A Comprehensive Review,  Journal of Physical and Chemical Reference Data
    ! Volume 3 (1974) 1-756
    ! The following function yields better result for thermal conductivity of air
    ! Expression 2; (kth)_air = Agas*T + Cgas , where Agas = 1.7082E-4 W/(m  K2) and Cgas = -7.488E-3 W/m K
    ! This function is adapted from the expression ((kth)_air = Agas*T + Bgas*T**2 + Dgas*T**3 + Egas*T**4 + Fgas*T**5+ Cgas ) provided in the following reference
    ! https://www.cambridge.org/core/books/abs/gas-turbines/equations-of-air-thermophysical-properties/9572106E068EFF1B7C0896124C17A196
    ! It means a value of 0 has been put in Bgas, Dgas, Egas and Fgas to enable thermal conducitivity to be well around the order of 10^{-2} for given temperature range
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Initial porosity (porosity) = 0.45
    ! Conductivity of bulk liquid copper fitted as a function of temperature 
    ! (kth_cu)liquid = Alcu*(T-1358)^2 + Blcu*(T-1358) + Dlcu, where Alcu = 2.51e-6 W/mK4 and Blcu = -1.6e-4 W/mK3 and Dlcu= 165.58 W/m  (1358.0 < T  K)
    ! Reference for (kth_cu)liquid: C.Y. Ho, R. W. Powell and P. E. Liley, Thermal Conductivity of the Elements: A Comprehensive Review,  Journal of Physical and Chemical Reference Data
    ! Volume 3 (1974) 1-756
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-08)
    !-----------------------------------------------------
    FUNCTION getThermalConductivity( model, n, temp ) RESULT(effcondt)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, effcondt, kthbulk, kthpowder, kthgas, tscaler

    ! variables needed inside function
    REAL(KIND=dp) ::  refMeltTemp, refStTemp, porosity,  &
    alphas, betas, deltas, alphal, betal, deltal, alphag, betag 
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getCuConductivity', 'No material found')
    END IF

    ! read in reference Conductivity at reference temperature
    !refSolidDenst = GetConstReal( material, 'Reference Conductivity Ccu Solid Cu',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getCuConductivity', 'Reference Conductivity Solid Cu not found')
    !END IF
    
    ! read in the coefficient of  Acu*(T-300)^2 term
    alphas = GetConstReal( material, 'Conductivity Coeff Acu Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Coefficientt of T3 term solid Cu not found')
    END IF
    
    ! read in the coefficient of  Bcu*(T-300)^2 term
    betas = GetConstReal( material, 'Conductivity Coeff Bcu Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Coefficientt of T2 term solid Cu not found')
    END IF
    
    ! read in the coefficient of Dcu*(T-300) term
    deltas = GetConstReal( material, 'Conductivity Coeff Dcu of Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'slope of Conductivity-temperature curve solid not found')
    END IF
    
    ! read in reference Conductivity of liquid Cu
    !refLiquidDenst = GetConstReal( material, 'Reference Conductivity Ccu Liquid Cu',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getCuConductivity', 'Reference Conductivity Liquid Cu not found')
    !END IF
    
    ! read in the coefficient of  Acu*(T-1330)^3 term
    alphal = GetConstReal( material, 'Conductivity Coeff Acu Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Coefficientt of T3 term solid Cu not found')
    END IF
    
    ! read in the coefficient of BCu*(T-1330)^2 term
    betal = GetConstReal( material, 'Conductivity Coeff Bcu Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Coefficient of T2 term Liquid Cu not found')
    END IF
    
    ! read in the coefficient of Dcu*(T-1330) term
    deltal = GetConstReal( material, 'Conductivity Coeff Dcu of Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'slope of Conductivity-temperature curve liquid not found')
    END IF
        
    ! read in pseudo reference Conductivity of gaseous air
    !refGasDenst = GetConstReal( material, 'Reference Conductivity Cgas of Air',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getCuConductivity', 'Reference Conductivity for Gaseous Air not found')
    !END IF
        
    ! read in the coefficient of  Agas*(T)^2 term
    alphag = GetConstReal( material, 'Conductivity Coeff Agas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Coefficient of T^2 term for air Conductivity not found')
    END IF
    
    ! read in the coefficient of Bgas*T term for air Conductivity
    betag = GetConstReal( material, 'Conductivity Coeff Bgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Slope of rhogas-T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    !deltag = GetConstReal( material, 'Conductivity Coeff Dgas of Air', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getCuConductivity', 'Coefficient of lnT term Gaseous Air not found')
    !END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of Cu powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Porosity of Cu powder not found')
    END IF
    
    ! read in reference sintering temperature of Ag powder
    refStTemp = GetConstReal( material, 'Sintering Temperature of Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Reference Sintering Temperature Cu not found')
    END IF

    ! read in reference melting temperature of Ag powder
    refMeltTemp = GetConstReal( material, 'Melting Point Temperature of Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Reference Melting Temperature Cu not found')
    END IF
    
     ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Scaling Factor for T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of Cu powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuConductivity', 'Porosity of Cu powder not found')
    END IF

     
    ! compute Conductivity conductivity
    ! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/else-if.html
    IF (refMeltTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getCuConductivity', 'The Cu material is in liquid state.')
            !CALL Warn('getCuConductivity', 'Using Conductivity reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    effcondt = alphal*(tscaler*temp)**3 + betal*(tscaler*temp)**2 +deltal*(tscaler*temp)
    ELSE IF (refMeltTemp > temp .AND. refStTemp < temp) THEN
       CALL Warn('getCuConductivity', 'The Cu material is being sintered.')
    kthbulk = alphas*(tscaler*temp)**3 + betas*(tscaler*temp)**2 + deltas*(tscaler*temp)
    !kthgas = alphag*(tscaler*temp)**2 + betag*(tscaler*temp) !Expression 1
    kthgas = alphag*(tscaler*temp)  + betag !Expression 2
    kthpowder = (1-porosity)*kthbulk + porosity*kthgas
    effcondt = ((kthbulk - kthpowder)*(temp-refStTemp))/(refMeltTemp-refStTemp)+kthpowder
    ELSE
    effcondt = (1-porosity)*(alphas*(tscaler*temp)**3 + betas*(tscaler*temp)**2 + deltas*(tscaler*temp))+ &
    porosity*(alphag*(tscaler*temp)**2 + betag*tscaler*temp)
    END IF

    END FUNCTION getThermalConductivity

