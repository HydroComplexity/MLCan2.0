function [VARIABLES] = CN_deposition (VARIABLES,PARAMS,CONSTANTS,VERTSTRUC)

% Applying soil layer
ApLayer = 1;

% AMM & Nit are time step! Not daily!

dtime            = CONSTANTS.dtime; %[s]
dz_mm            = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
dz_m             = dz_mm./1000;
N_Adepo_amm      = PARAMS.CN.N_Adepo_amm; % [g/m2/yr]
N_Adepo_nit      = PARAMS.CN.N_Adepo_nit; % [g/m2/yr]
N_Adepo_amm_gm3s = N_Adepo_amm/dz_m(ApLayer)/(365/(24*60*60/dtime)); % [g/m3/dtime]
N_Adepo_nit_gm3s = N_Adepo_nit/dz_m(ApLayer)/(365/(24*60*60/dtime)); % [g/m3/dtime]

Nit=VARIABLES.Nit; % Nitrate concentration in soil [gN / m^3]
Amm=VARIABLES.Amm; % Ammonium concentration in soil [gN / m^3]
%Ammdd = VARIABLES.Ammdd;
%Nitdd = VARIABLES.Nitdd;

% N deposition is simply added into the "applied soil layer "
Nit(ApLayer)=Nit(ApLayer)+N_Adepo_nit_gm3s;
Amm(ApLayer)=Amm(ApLayer)+N_Adepo_amm_gm3s;
%Nitdd(ApLayer)=Nitdd(ApLayer)+N_Adepo_nit_gm3s;
%Ammdd(ApLayer)=Ammdd(ApLayer)+N_Adepo_amm_gm3s;

VARIABLES.Nit=Nit;
VARIABLES.Amm=Amm;
%VARIABLES.Nitdd=Nitdd;
%VARIABLES.Ammdd=Ammdd;