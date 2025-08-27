%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seagrass Carbon Model
%
% Project LightBlue - Shedding light on the European Blue Carbon Budget
% Units S4 and D2 (Ocean and Water)
% European Commission, Directorate-General Joint Research Centre (JRC), Ispra (Italy)
%
% Aim: simulate long-term carbon dynamics in seagrass meadows.
%   It tracks dissolved inorganic carbon (DIC), carbon burial in sediments,
%   and partial pressure of CO₂ of water (pCO₂)water in response to fluxes such as
%   air–sea gas exchange, burial, export, resuspension, and dilution.
%
% Dependencies:
%   - CO2SYS.m ((van Heuven (2011), with uncertainty propagation (Orr et al., 2018),
%     based on the original CO2SYS program by Lewis and Wallace (1998);
%     calculate carbonate chemistry from alkalinity and DIC)
%   - Temp_timeseries.mat  (daily temperature time series, °C)
%   - Sal_timeseries.mat   (daily salinity time series, )
%
% Outputs:
%   - simulation_output.xlsx: full time series of states and fluxes
%
% Authors:
%   Ana DE AZEVEDO E COSTA and Luca POLIMENE
%   last version: August 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
tic;

%% =======================================================================
% 1. MODEL SETUP
% ========================================================================
% Simulation time
days        = 365;        % Days per year
numberYears = 1500;       % Simulation length (years)
extent      = days * numberYears;  % Total number of timesteps
tstep       = 1.0;        % Model timestep (days)

% -------------------- STATE VARIABLES -----------------------------------
% Units: mmol C m⁻²
dic      = zeros(extent, 1); % Dissolved Inorganic Carbon in water column
cburied  = zeros(extent, 1); % Carbon buried in sediments

% Initial conditions
dic(1)     = 2050 * 10;   % mmol/m², assuming 10 m water depth
cburied(1) = 0.0;         % No initial sediment carbon

% -------------------- FLUX VARIABLES ------------------------------------
% Units: mmol C m⁻² d⁻¹
cuptake      = zeros(extent, 1); % Air–sea CO₂ exchange
burial       = zeros(extent, 1); % carbon burial
resuspension = zeros(extent, 1); % Carbon released from sediments
dilution     = zeros(extent, 1); % Mixing with reference DIC
export       = zeros(extent, 1); % Organic carbon exported from system

% -------------------- ADDITIONAL FLUXES ---------------------------------
% (µatm)
PCO2     = zeros(extent, 1); % Partial pressure of CO₂

%% =======================================================================
% 2. LOAD FORCING DATA
% ========================================================================
% External drivers: daily time series of temperature and salinity.
load 'Temp_timeseries.mat';   % [°C]
load 'Sal_timeseries.mat';    % [ ]

% Repeat input forcing to cover all simulation years
T = repmat(Temp_timeseries, [numberYears, 1]); % Temperature forcing
S = repmat(Sal_timeseries, [numberYears, 1]);  % Salinity forcing

%% =======================================================================
% 3. MODEL PARAMETERS
% ========================================================================
rho         = 1025;   % Seawater density (kg/m³)
CO2atm      = 280;    % Atmospheric CO₂ concentration (µatm) - pre-industrial levels
alkalinity  = 2300;   % Ocean alkalinity (µmol/kg)
wnd         = 4;      % Mean wind speed (m/s)

% Seagrass community carbon fluxes (Duarte, 2017), converted to mmol
ncp         = 460 / 12;   % Net community production (mmol C m⁻² d⁻¹)
burialrate  = 184 / 12;   % Burial rate (mmol C m⁻² d⁻¹)
exportrate  = 276 / 12;   % Export rate (mmol C m⁻² d⁻¹)

% Other rates
resrate     = 0.000012; % Sediment resuspension rate (d⁻¹)
dilrate     = 0.01;     % Dilution/mixing rate (d⁻¹)

% Reference DIC concentration (mmol C m⁻²)
dicr        = 2050 * 10;

%% =======================================================================
% 4. MAIN SIMULATION LOOP
% ========================================================================
for i = 2:extent
    % Convert to Kelvin
    TK = T(i) + 273.15;
    tk100 = TK / 100; % For Weiss solubility equation

    % -------------------------------------------------------------------
    % 4.1 Carbonate chemistry (via CO2SYS)
    % -------------------------------------------------------------------
    % Convert DIC from mmol/m² to µmol/kg (assuming 10 m depth & 1025 kg/m³)
    dic_mmol_kg = dic(i-1) / 10;

    % CO2SYS input parameters
    par1       = alkalinity; % Total alkalinity (µmol/kg)
    par2       = dic_mmol_kg; % DIC (µmol/kg)
    par1_type  = 1;          % Input type 1 = TA
    par2_type  = 2;          % Input type 2 = DIC
    sal        = S(i);       % Salinity
    tempin     = T(i);       % Temperature at input
    presin     = 0;          % Pressure input (dbar)
    tempout    = T(i);       % Output temperature (same as input)
    presout    = 0;          % Output pressure (dbar)
    sil        = 0;          % Silicate concentration (µmol/kg)
    po4        = 0;          % Phosphate concentration (µmol/kg)
    pHscale    = 1;          % pH scale (1 = Total scale)
    k1k2c      = 4;          % Constants K1/K2 (Mehrbach refit by Dickson & Millero)
    kso4c      = 1;          % HSO₄⁻ dissociation constant (Dickson)

    % Call CO2SYS
    result = CO2SYS(par1, par2, par1_type, par2_type, ...
        sal, tempin, tempout, presin, presout, ...
        sil, po4, pHscale, k1k2c, kso4c);

    % Extract pCO₂ from CO2SYS output
    pco2 = result(4);
    PCO2(i) = pco2;

    % -------------------------------------------------------------------
    % 4.2 Air–sea CO₂ solubility
    % -------------------------------------------------------------------
    k0 = exp(93.4517 / tk100 - 60.2409 + 23.3585 * log(tk100) ...
        + S(i) * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk100^2));

    % -------------------------------------------------------------------
    % 4.3 Gas exchange
    % -------------------------------------------------------------------
    sc = 2116.8 - 136.25*T(i) + 4.7353*T(i)^2 ...
         - 0.092307*T(i)^3 + 0.0007555*T(i)^4; % Schmidt number
    fwind = (0.222 * wnd^2 + 0.333 * wnd) * (sc / 600)^(-0.5);
    fwind = fwind * 24 / 100; % Convert to m/day

    % -------------------------------------------------------------------
    % 4.4 Flux calculations
    % -------------------------------------------------------------------
    UPTAKE      = fwind * k0 * (CO2atm - pco2) * rho / 1000; % mmol/m²/d
    DIL         = (dicr - dic(i-1)) * dilrate;               % Dilution flux
    RESUSP      = cburied(i-1) * resrate;                    % Sediment resusp.
    EXPORT      = exportrate;                                % Constant export
    BURIAL      = burialrate;                                % Constant burial

    % -------------------------------------------------------------------
    % 4.5 State variable calculations
    % -------------------------------------------------------------------
    SDIC   = UPTAKE + DIL + RESUSP - BURIAL - EXPORT; % Change in DIC
    SCBUR  = BURIAL - RESUSP;                         % Change in buried C

    % Update state variables
    dic(i)     = dic(i-1) + SDIC * tstep;
    cburied(i) = cburied(i-1) + SCBUR * tstep;

    % Save fluxes
    cuptake(i)      = UPTAKE;
    burial(i)       = BURIAL;
    resuspension(i) = RESUSP;
    dilution(i)     = DIL;
    export(i)       = EXPORT;
end

%% =======================================================================
% 5. MASS BALANCE CHECK
% ========================================================================
tot_C = dic(end) - dic(1) + ...
        cburied(end) - cburied(1) - ...
        sum(cuptake) - sum(dilution) + sum(export);

fprintf('Mass balance check (net C): %.2f mmol/m²\n', tot_C);

%% =======================================================================
% 6. PLOTTING
% ========================================================================
time = (1:extent)'; % Simulation time (days)

figure;
plot(time, dic, 'b');
xlabel('Time (days)'); ylabel('DIC (mmol/m²)');
title('Dissolved Inorganic Carbon (DIC)');

figure;
plot(time, cburied, 'g');
xlabel('Time (days)'); ylabel('Carbon Buried (mmol/m²)');
title('Carbon Burial');

figure;
plot(time, PCO2, 'r');
xlabel('Time (days)'); ylabel('pCO₂ (µatm)');
title('Partial Pressure of CO₂');

figure;
plot(time, tot_C, 'k');
xlabel('Time (days)'); ylabel('Total Carbon Balance (mmol/m²)');
title('Carbon Mass Balance');

%% =======================================================================
% 7. SAVE OUTPUTS
% ========================================================================

% --- Save full simulation output ---
outputTable = table(time, dic, cburied, PCO2, ...
    cuptake, burial, resuspension, dilution, export);
writetable(outputTable, 'simulation_output.xlsx');

toc;
return;

