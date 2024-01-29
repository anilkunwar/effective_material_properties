    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Density of bulk solid copper fitted as a function of temperature
    ! (rho_cu)solid = Ascu*(T-300)^3 + Bscu*(T-300)^2 + Dscu*(T-300) + Cscu, where Ascu = -2.77E-11 kg/m3K3 and Bscu = -6.12E-8 kg/m3K2 and Dscu=-4.28E-04 kg/m3K and Ccu = 881.0 kg/m3 
    ! Reference for (rho_cu)solid: Demin et al. MATHEMATICA MONTISNIGRI (2020), Vol. XLVII:137-151
    ! Density of gaseous air fitted as a function of temperature
    ! (rho_air)gas = Agas*T^2 + Bgas*T  + Dgas*ln(Egas*T) +Cgas, where Agas = -7.65e-07 kg/m3K2  and  Bgas=  3.16e-03 kg/m3K and Dgas=-1.94 kg/m3 and Egas = 1.0 K-1 and Cgas = 1.14e+01  kg/m3
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Initial porosity (porosity) = 0.45
    ! Density of bulk liquid copper fitted as a function of temperature 
    ! (rho_cu)liquid = Alcu*(T-1330)^3 + Blcu*(T-1330)^2 + Dlcu*(T-1330) + Clcu, where Alcu = -2.383E-11 kg/m3K3 and Blcu = +8.69E-8 kg/m3K2 and Fcu=-7.96E-04 and Ccu = 7.89 kg/m3 
    ! Reference for (rho_cu)liquid: Demin et al. MATHEMATICA MONTISNIGRI (2020), Vol. XLVII:137-151
    ! https://www.montis.pmf.ac.me/allissues/47/Mathematica-Montisnigri-47-12.pdf
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-08)
    !-----------------------------------------------------
    FUNCTION getDensity( model, n, temp ) RESULT(effdenst)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, effdenst, rhobulk, rhopowder, rhogas, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolidDenst, refLiquidDenst, refGasDenst,  &
    alphas, betas, deltas, alphal, betal, deltal, alphag, betag, deltag,  & 
    porosity,  refStTemp, refMeltTemp
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getCuDensity', 'No material found')
    END IF

    ! read in reference density at reference temperature
    refSolidDenst = GetConstReal( material, 'Reference Density Ccu Solid Cu',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Reference Density Solid Cu not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-300)^2 term
    alphas = GetConstReal( material, 'Density Coeff Acu Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficientt of T3 term solid Cu not found')
    END IF
    
    ! read in the coefficient of  Bcu*(T-300)^2 term
    betas = GetConstReal( material, 'Density Coeff Bcu Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficientt of T2 term solid Cu not found')
    END IF
    
    ! read in the coefficient of Dcu*(T-300) term
    deltas = GetConstReal( material, 'Density Coeff Dcu of Solid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'slope of Density-temperature curve solid not found')
    END IF
    
    ! read in reference density of liquid Cu
    refLiquidDenst = GetConstReal( material, 'Reference Density Ccu Liquid Cu',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Reference Density Liquid Cu not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-1330)^3 term
    alphal = GetConstReal( material, 'Density Coeff Acu Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficientt of T3 term solid Cu not found')
    END IF
    
    ! read in the coefficient of BCu*(T-1330)^2 term
    betal = GetConstReal( material, 'Density Coeff Bcu Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficient of T2 term Liquid Cu not found')
    END IF
    
    ! read in the coefficient of Dcu*(T-1330) term
    deltal = GetConstReal( material, 'Density Coeff Dcu of Liquid Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'slope of Density-temperature curve liquid not found')
    END IF
        
    ! read in pseudo reference density of gaseous air
    refGasDenst = GetConstReal( material, 'Reference Density Cgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Reference Density for Gaseous Air not found')
    END IF
        
    ! read in the coefficient of  Agas*(T)^2 term
    alphag = GetConstReal( material, 'Density Coeff Agas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficient of T^2 term for air density not found')
    END IF
    
    ! read in the coefficient of Bgas*T term for air density
    betag = GetConstReal( material, 'Density Coeff Bgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Slope of rhogas-T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    deltag = GetConstReal( material, 'Density Coeff Dgas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Coefficient of lnT term Gaseous Air not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of Cu powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Porosity of Cu powder not found')
    END IF
    
    ! read in reference sintering temperature of Ag powder
    refStTemp = GetConstReal( material, 'Sintering Temperature of Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Reference Sintering Temperature Cu not found')
    END IF

    ! read in reference melting temperature of Ag powder
    refMeltTemp = GetConstReal( material, 'Melting Point Temperature of Cu', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Reference Melting Temperature Cu not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Scaling Factor for T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of Cu powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getCuDensity', 'Porosity of Cu powder not found')
    END IF

     
    ! compute density conductivity
    ! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/else-if.html
    IF (refMeltTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getCuDensity', 'The Cu material is in liquid state.')
            !CALL Warn('getCuDensity', 'Using density reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    effdenst = refLiquidDenst + alphal*(tscaler*temp-1330.00)**3 + betal*(tscaler*temp-1330.00)**2 +deltal*(tscaler*temp-1330.0)
    ELSE IF (refMeltTemp > temp .AND. refStTemp < temp) THEN
       CALL Warn('getCuDensity', 'The Cu material is being sintered.')
    rhobulk = refSolidDenst +  alphas*(tscaler*temp-300.0)**3 + betas*(tscaler*temp-300.0)**2 + deltas*(tscaler*temp-300.0)
    rhogas = refGasDenst +  alphag*(tscaler*temp)**2 + betag*tscaler*temp + deltag*LOG(tscaler*temp)
    rhopowder = (1-porosity)*rhobulk + porosity*rhogas
    effdenst = ((rhobulk - rhopowder)*(temp-refStTemp))/(refMeltTemp-refStTemp)+rhopowder
    ELSE
    effdenst = (1-porosity)*(refSolidDenst +  alphas*(tscaler*temp-300.00)**3 + betas*(tscaler*temp-300)**2 + &
    deltas*(tscaler*temp-300))+porosity*(refGasDenst +  alphag*(tscaler*temp)**2 + betag*tscaler*temp + deltag*LOG(tscaler*temp))
    END IF

    END FUNCTION getDensity

