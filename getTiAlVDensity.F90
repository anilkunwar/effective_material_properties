    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Density of solid Ti6Al4V fitted as a function of temperature
    ! (rho_tialv)solid = As*(T) +  Bs*T^2 +Ds*ln(T)+Cs, where As = 1.152 kg/m3K  and Bs = -4.1197e-04 kg/m3K2 and Ds = -4.254e+02 and Cs = 6.54263e+03 kg/m3
    ! 298 < T < Tm where Tm = 1923.0 K
    ! Reference: Lu et al., Additive Manufacturing, 26 (2019) 166-179.
    ! https://www.sciencedirect.com/science/article/pii/S2214860418305955
    ! Density of liquid Ti6Al4V fitted as a function of temperature
    !(rho_tialv)liquid  = Al*(T) + Bl, where Al = -0.452 kg/m3K  and Bl = 4955.0 kg/m3 (1923.0 K < T < 2500.0 K)
    ! Reference: Schmon et al, EPJ Web of Conferences, 151 (2017) 04003.
    ! https://www.epj-conferences.org/articles/epjconf/pdf/2017/20/epjconf_lam2017_04003.pdf
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-12)
    !-----------------------------------------------------
    FUNCTION getDensity( model, n, temp ) RESULT(denst)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, denst, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolDenst, refLiqDenst, refTemp,  &
    alphas, betas, deltas, alphal
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getDensity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolDenst = GetConstReal( material, 'Reference Density Cs Solid TiAlV',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Density Solid TiAlV not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Density Coeff As Solid TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'slope of Density-temperature curve solid not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    betas = GetConstReal( material, 'Density Coeff Bs Solid TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Coefficientt of Bs*T2 term solid TiAlV not found')
    END IF
    
    ! read in  Ds in Ds*ln(T) term
    deltas = GetConstReal( material, 'Density Coeff Ds Liquid TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Coefficient of logT term liquid TiAlV not found')
    END IF
    
    ! read in reference density at reference temperature
    refLiqDenst = GetConstReal( material, 'Reference Density Bl Liquid TiAlV',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Density Solid TiAlV not found')
    END IF
    
    ! read in pseudo reference conductivity at reference temperature of liquid
    alphal = GetConstReal( material, 'Density Coefficient Al Liquid TiAlV',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Density Coefficient Al of Liquid TiAlV not found')
    END IF

    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Melting Temperature of TiAlV not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    END IF

    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getDensity', 'The TiAlV material is in liquid state.')
            !CALL Warn('getDensity', 'Using density reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    denst = refLiqDenst + alphal*((tscaler)*temp)
    ELSE
    denst = refSolDenst + alphas*((tscaler)*temp) + betas*((tscaler)*temp)**2 + deltas*LOG((tscaler)*temp)
    END IF

    END FUNCTION getDensity

