# Code based on:
# https://biocycle.atmos.colostate.edu/shiny/carbonate/
# http://www-naweb.iaea.org/napc/ih/documents/global_cycle/vol%20I/cht_i_09.pdf

# INPUT VARIABLES
#       T      #  Temperature (Celsius)
#       DIC    #  total dissolved inorganic carbon (mg/L)
#       pH     #  pH of the water
#       output #  units. Default to "mg". Option = "mol"

# CHEMICAL & THERMODYNAMIC COEFFICIENTS  
#       K0     #  Henry's Law constant
#       K1     #  first dissociation coefficient for H2CO3
#       K2     #  second dissociation coefficient for H2CO3
#       Kw     #  dissociation constant of water 


carbonate <- function(TEMP=16, DIC=2002, pH=8, output = "mg") {
  
  DIC = DIC / (12.01*1000) # Conversion from mg/L to mol
  TEMP = TEMP + 273.15 # temperature from Celsius to Kelvin
  
  # Most freshwaters can be considered as ideal solutions (extrapolated to zero ionic strength)

  # all values of the solubilities and dissociation constants are temperature dependent. 
  # Carbonate equilibrium constants as functions of temp (K)
  pK0 = -2622.38/TEMP - 0.0178471*TEMP + 15.5873 
  
  pK1 = 3404.71/TEMP + 0.032786*TEMP - 14.8435
  
  pK2 = 2902.39/TEMP + 0.02379*TEMP - 6.4980
  
  lnKw = 148.9802 - 13847.26/TEMP - 23.6521*log(TEMP)
  Kw = exp(lnKw)
  
  # Convert pH to H+
  Hplus = 10^(-pH)           
  OHneg = Kw/Hplus
    
  K0 = 10^(-pK0)
  K1 = 10^(-pK1)
  K2 = 10^(-pK2)
  
  # The total concentration of dissolved inorganic carbon is defined by
  # DIC = CO2aq + H2CO3 + HCO3 + CO3 
  
  H2CO3 = DIC * Hplus^2 / (Hplus^2 + (Hplus*K1) + (K1*K2))
  HCO3 = DIC * (Hplus*K1) / (Hplus^2 + (Hplus*K1) + (K1*K2))
  CO3 = DIC * (K1*K2) / (Hplus^2 + (Hplus*K1) + (K1*K2))
  
  # titration alkalinity
  ALK = HCO3 + 2*CO3 + OHneg - Hplus #eq/kg
  # To find the alkalinity in terms of mg/L of calcium carbonate, a commonly used measure of alkalinity, multiply by 50,000:
  #ALK eq/kg × 50,000 mg/eq = mg/kg as CaCO3
  
  if (output == "mol") {
    return(data.frame(co2_molkg = H2CO3, bicarbonate_molkg = HCO3, carbonate = CO3_molkg, alkalinity_eqkg = ALK)) # mol/kg
  } else if (output == "mg") {
  # To convert moles per kg to milligrams per kg (mg/kg), 
  # multiply the bicarbonate result by 61,017.1, the carbonate result by 60,009.2, 
  # and the hydroxide result by 17,007.3. 
    return(data.frame(co2_mgkg = H2CO3 * 62024.8, bicarbonate_mgkg = HCO3*61017.1, 
           carbonate_mgkg = CO3*60009.2, alkalinity_uEqkg = ALK * 1e6)) # mg/kg
  }
}
