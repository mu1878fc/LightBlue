# LightBlue
# Seagrass Carbon Model
#
# Project LightBlue - Shedding light on the European Blue Carbon Budget
# Units S4 and D2 (Ocean and Water)
# European Commission, Directorate-General Joint Research Centre (JRC), Ispra (Italy)
#
# Aim: simulate long-term carbon dynamics in seagrass meadows.
#   It tracks dissolved inorganic carbon (DIC), carbon burial in sediments,
#   and partial pressure of CO₂ of water (pCO₂)water in response to fluxes such as
#   air–sea gas exchange, burial, export, resuspension, and dilution.
#
# Dependencies:
#   - CO2SYS.m ((van Heuven (2011), with uncertainty propagation (Orr et al., 2018),
#     based on the original CO2SYS program by Lewis and Wallace (1998);
#     calculate carbonate chemistry from alkalinity and DIC)
#   - Temp_timeseries.mat  (daily temperature time series, °C)
#   - Sal_timeseries.mat   (daily salinity time series, )
#
# Outputs:
#   - simulation_output.xlsx: full time series of states and fluxes
#
# Authors:
#   Ana DE AZEVEDO E COSTA and Luca POLIMENE
#   last version: August 2025
Seagrass Carbon Model
