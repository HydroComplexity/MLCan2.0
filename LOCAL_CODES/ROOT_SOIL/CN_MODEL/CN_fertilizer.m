function [VARIABLES] = CN_fertilizer (VARIABLES,PARAMS,CONSTANTS,VERTSTRUC)

% Applying soil layer
ApLayer = 1:2;

dtime       = CONSTANTS.dtime; %[s]
dz_mm       = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
dz_m        = dz_mm./1000;

N_Fert_DOY  = PARAMS.CN.N_Fert_DOY; % [DOY]
N_Fert_amm  = PARAMS.CN.N_Fert_amm; % [gN/m2]
N_Fert_nit  = PARAMS.CN.N_Fert_nit; % [gN/m2]
N_Fert_urea = PARAMS.CN.N_Fert_urea; % [gN/m2]
VolAmm      = PARAMS.CN.VolAmm;
VolUrea     = PARAMS.CN.VolUrea;


Nit=VARIABLES.Nit; % Nitrate concentration in soil [gN / m^3]
Amm=VARIABLES.Amm; % Ammonium concentration in soil [gN / m^3]
%Ammdd = VARIABLES.Ammdd;
%Nitdd = VARIABLES.Nitdd;

% Compute fertilizer amount with volatilization factor
N_Fert_amm_gm3s(ApLayer)  = N_Fert_amm.*(1-VolAmm)./dz_m(ApLayer).*(dz_m(ApLayer)./sum(dz_m(ApLayer))); % [g/m3/dtime]
N_Fert_nit_gm3s(ApLayer)  = N_Fert_nit./dz_m(ApLayer).*(dz_m(ApLayer)./sum(dz_m(ApLayer))); % [g/m3/dtime]
N_Fert_urea_gm3s(ApLayer) = N_Fert_urea.*(1-VolUrea)./dz_m(ApLayer).*(dz_m(ApLayer)./sum(dz_m(ApLayer))); % [g/m3/dtime]
N_Fert_amm_gm3s=N_Fert_amm_gm3s(:);
N_Fert_nit_gm3s=N_Fert_nit_gm3s(:);
N_Fert_urea_gm3s=N_Fert_urea_gm3s(:);

% Check
%sum(N_Fert_amm_gm3s(ApLayer).*dz_m(ApLayer)); ==
%N_Fert_amm.*(1-VolAmm)

% N fertilizer is simply added into the "applied soil layer "
Nit(ApLayer)=Nit(ApLayer)+N_Fert_nit_gm3s;
Amm(ApLayer)=Amm(ApLayer)+N_Fert_amm_gm3s+N_Fert_urea_gm3s;
%Nitdd(ApLayer)=Nitdd(ApLayer)+N_Fert_nit_gm3s;
%Ammdd(ApLayer)=Ammdd(ApLayer)+N_Fert_amm_gm3s+N_Fert_urea_gm3s;

VARIABLES.Nit=Nit;
VARIABLES.Amm=Amm;
%VARIABLES.Nitdd=Nitdd;
%VARIABLES.Ammdd=Ammdd;
