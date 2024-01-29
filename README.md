# The Ti6Al4V material is assumed as a continuum material with zero porosity


#Cu MATERIAL
  The Cu material is assumed to be a part of a porous layer with initial porosity of 0.45, and this porosity value is temperature dependent.
  Conductivity:
  Liquid State:
If the reference melting temperature (refMeltTemp) is less than or equal to the current temperature (temp), it considers the material to be in a liquid state and calculates the effective thermal conductivity based on a model:
effcondt = α_l * (tscaler * temp)^3 + β_l * (tscaler * temp)^2 + δ_l * (tscaler * temp)

