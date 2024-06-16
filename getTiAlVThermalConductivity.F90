    !-----------------------------------------------------
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-12)
    ! material property user defined function for ELMER:
    ! Thermal conductivity of Ti6Al4V fitted as a function of temperature
    ! (kth_tialv)solid = As*(T)  +  Cs, where As = 1.46E-02 W /mK2 and    Cs = -0.32 W/mK (1400.0 K < T < 1800.0 K)
    ! Reference: Boivineau et al., International Journal of Thermophysics, 27 (2006) 507-529.
    ! https://link.springer.com/article/10.1007/PL00021868
    ! (kth_tialv)liquid = Al*(T) + Cl, where Al = 1.83E-02  W/mK2 and E = -6.66 W/m K (1950.0 K < T < 2700.0 K) 
    ! Reference: Boivineau et al., International Journal of Thermophysics, 27 (2006) 507-529.
    ! https://link.springer.com/article/10.1007/PL00021868
    !-----------------------------------------------------
    FUNCTION getThermalConductivity( model, n, temp ) RESULT(thcondt)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, thcondt, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolThCond, refLiqThCond,refTemp, &
    alphas, alphal 
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalConductivity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolThCond = GetConstReal( material, 'Reference Thermal Conductivity C Solid TiAlV',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Thermal Conductivity Solid TiAlV not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Cond Coeff A Solid TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'slope of thermal conductivity-temperature curve solid not found')
    END IF
    
    ! read in reference conductivity at reference temperature
    refLiqThCond = GetConstReal( material, 'Reference Thermal Conductivity C Liquid TiAlV',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Thermal Conductivity Liquid TiAlV not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    alphal = GetConstReal( material, 'Cond Coeff A Liquid TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Coefficientt of T  term Liquid TiAlV not found')
    END IF
    
    
    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature TiAlV', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Temperature not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    END IF


    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalConductivity', 'The Ti material is in liquid state.')
            !CALL Warn('getThermalConductivity', 'Using density reference value')
    !thcondt = 1.11*(refThCond + alpha*(temp))
    thcondt = refLiqThCond + alphal*((tscaler)*temp)
    ELSE
    thcondt = refSolThCond + alphas*((tscaler)*temp) 
    END IF

    END FUNCTION getThermalConductivity

