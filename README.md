[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://phasefieldcalculator.streamlit.app/)

# Ti6Al4V ALLOY

The Ti6Al4V material is assumed as a continuum material with zero porosity


# Cu MATERIAL
  The Cu material is assumed to be a part of a porous layer with initial porosity of 0.45, and this porosity value is temperature dependent.
  
 *a. Conductivity:*
  
  Liquid State:
  
If the reference melting temperature (refMeltTemp) is less than or equal to the current temperature (temp), it considers the material to be in a liquid state and calculates the effective thermal conductivity based on a model:

effcondt = α_l * (tscaler * temp)^3 + β_l * (tscaler * temp)^2 + δ_l * (tscaler * temp)

Sintering State:

If the reference melting temperature is greater than the current temperature (temp) and the reference sintering temperature (refStTemp) is less than the current temperature, it considers the material to be in a sintering state. It calculates the effective thermal conductivity based on a sintering model:

effcondt = ((α_bulk - kthpowder) * (temp - refStTemp)) / (refMeltTemp - refStTemp) + kthpowder

Other States:
If the above conditions are not met, it calculates the effective thermal conductivity based on a default model:

effcondt = (1 - porosity) * (α_s * (tscaler * temp)^3 + β_s * (tscaler * temp)^2 + δ_s * (tscaler * temp)) + porosity * (α_g * (tscaler * temp)^2 + β_g * tscaler * temp)

*b. Density:*

Liquid State:

If the reference melting temperature (refMeltTemp) is less than or equal to the current temperature (temp), it considers the material to be in a liquid state and calculates the effective density conductivity based on a model:

effdenst = refLiquidDenst + alphal * (tscaler * temp - 1330.00)**3 + betal * (tscaler * temp - 1330.00)**2 + deltal * (tscaler * temp - 1330.0)


Sintering State:
If the reference melting temperature is greater than the current temperature (temp) and the reference sintering temperature (refStTemp) is less than the current temperature, it considers the material to be in a sintering state. It calculates the effective density conductivity based on a model:



effdenst = ((rhobulk - rhopowder) * (temp - refStTemp)) / (refMeltTemp - refStTemp) + rhopowder


where rhobulk, rhogas, and rhopowder are densities of bulk Cu, gaseous air and powdered Cu respectively.

Other States:
If the above conditions are not met, it calculates the effective density conductivity based on a default model:

effdenst = (1 - porosity) * (refSolidDenst + alphas * (tscaler * temp - 300.00)**3 + betas * (tscaler * temp - 300)**2 + deltas * (tscaler * temp - 300)) + porosity * (refGasDenst + alphag * (tscaler * temp)**2 + betag * tscaler * temp + deltag * LOG(tscaler * temp))





