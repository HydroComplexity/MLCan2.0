function [VARIABLES] = ...
        CN_Cycle_Porporato2(rootfr,PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC)

%   Calculate Dynamics of Soil Carbon and Nitrogen
%
%   INPUTS:
%       INITIAL CONDITIONS:
%       sm = soil moisture profile [m^3 / m^3]
%       porsl = porosity profile [m^3 / m^3]
%       Ts = soil temperature profile [C]
%       Ts_max = Maximum temperature [C]. Used to compute environmental factors
%       qq = water fluxes between layers [mm / s]
%       layeruptake = Uptake from each layer [mm/s]
%       layeruptake_All = Uptake from each layer by each species separate [mm/s]
%       rootfr  = Root distribution of fine roots.

%       ku_Amm = 
%       ku_Nit = 
%       ki_Amm = partitioning coefficient for ammonium immobilization [m^3 d / gC]
%       ki_Nit = partitioning coefficient for nitrate immobilization [m^3 d / gC]

% DE REFERENCE BLOCKS   
% water and temperature in the soil
nspecies = PARAMS.CanStruc.nspecies; % number of species
if SWITCHES.CN.Bioturbation;
    sm = [VARIABLES.SOIL.volliqli ; VARIABLES.SOIL.volliq];
    Ts = [VARIABLES.SOIL.Tli ; VARIABLES.SOIL.Ts];
    qq = [VARIABLES.SOIL.qinflL ; VARIABLES.SOIL.qlayer];
    layeruptake_all = [zeros(1,nspecies) ; VARIABLES.SOIL.layeruptake_all];
    smp = [VARIABLES.SOIL.psil ; VARIABLES.SOIL.smp];
    nl_soil=PARAMS.nl_soil + 1;%       nl_soil = # soil layers
    rootfr=[zeros(1,nspecies) ; rootfr];
    dz_mm = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
else
    sm = VARIABLES.SOIL.volliq;
    Ts = VARIABLES.SOIL.Ts;
    qq = VARIABLES.SOIL.qlayer;
    layeruptake_all = VARIABLES.SOIL.layeruptake_all;
    smp = VARIABLES.SOIL.smp;
    nl_soil=PARAMS.nl_soil ;%       nl_soil = # soil layers
    rootfr=rootfr;
    dz_mm = VERTSTRUC.dzsmm;%       dz_mm = grid depth [mm]    
end

a_Amm = PARAMS.CN.a_Amm;%       a_Amm = fraction of dissolved ammonium [-]
a_Nit = PARAMS.CN.a_Nit;%       a_Nit = fraction of dissolved nitrate [-]
ki_Amm = PARAMS.CN.ki_Amm; % constant that determines the partition between NO3 and NH4 for immobilization
ki_Nit = PARAMS.CN.ki_Nit; % constant that determines the partition between NO3 and NH4 for immobilization
rr = PARAMS.CN.rr;%       rr = fraction of decomposed organic carbon that goes to respiration [-]
Cl = VARIABLES.Cl;%       Cl = carbon concentration in litter pool [gC / m^3]
Ch = VARIABLES.Ch;%       Ch = carbon concentration in humus pool [gC / m^3]
Cb = VARIABLES.Cb;%       Cb = carbon concentration in biomass pool [gC / m^3]
Nl = VARIABLES.Nl;%       Nl = nitrogen concentration in litter pool [gN / m^3]
Amm = VARIABLES.Amm;%       Amm = ammonium concentration in soil [gN / m^3]
Nit = VARIABLES.Nit;%       Nit = nitrate concentration in soil [gN / m^3]
TR_can_store = VARIABLES.CANOPY.TR_can_all;%        TR total transpiration in each species [mm/s]      
dtime = CONSTANTS.dtime;%       dtime = Time step [1800 s]
%*************************************************************************
%                       PREALLOCATE VECTORS 
  PHI = zeros(nl_soil,1);
  phi = ones(nl_soil,1);  
%*************************************************************************
%                       CHANGE OF UNITS AND DE-REFERENCE BLOCK

qq=qq*86400/1000;  % Change units from [mm/s] to [m/d]

%layeruptake=layeruptake*86400/1000;  % Change units from [mm/s] to [m/d]
layeruptake_all=layeruptake_all*86400/1000;  % Change units from [mm/s] to [m/d]
TR_can_store = TR_can_store*86400/1000;

% check that evaporation fluxes are not included
 if SWITCHES.CN.Bioturbation
     qq(1) = max(qq(1),0);
     qq(2) = max(qq(2),0);
 end
%*************************************************************************
%     COMPUTE ENVIRONMENTAL FACTORS 
%*************************************************************************            
CN_envfactor ();
%*************************************************************************
%   BIOTURBATION (take out)
%*************************************************************************
if SWITCHES.CN.Bioturbation
    [Cl, VARIABLES] = CN_bioturbation (PARAMS, VARIABLES, CONSTANTS, FORCING, VERTSTRUC, SWITCHES, Cl,fTd);
end
%*************************************************************************
%     COMPUTE ADITTION OF LITTER. (LITTER, HUMUS, BIOMASS) 
%*************************************************************************            
[CNa, ADD, CNo, OUT, ADD_bio, ADD_ex, ADD_net] = CN_addlitter (FORCING, PARAMS,...
    VERTSTRUC, VARIABLES, CONSTANTS, SWITCHES, rootfr);
%*************************************************************************
%     COMPUTE RATES OF CHANGES IN POOLS CONCENTRATION. (LITTER, HUMUS, BIOMASS) 
%*************************************************************************            
[dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, VARIABLES] = ...
    CN_dynamics(VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo);
            
%*************************************************************************
%     COMPUTE NET MINERALIZATION OR NET IMMOBILIZATION
%*************************************************************************
% Preallocate vectors        
    [phi, PHI, MIN_net, IMM_net, MIN_gross, IMM_gross, Nreg, DECl] = ...
        CN_computephi (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, ADD, CNa); 

                        
    [dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, VARIABLES] =...
        CN_dynamics (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo);    
    
    % PARTITION IMMOBILIZATION BETWEEN NO3- AND NH4+
    indn = Amm ~= 0;
    IMM_Amm(indn) = (ki_Amm*Amm(indn)./ (ki_Amm*Amm(indn) + ki_Nit*Nit(indn))).*IMM_net(indn);
    IMM_Amm(~indn) = 0;
    IMM_Amm = IMM_Amm(:);
    indn = Nit ~= 0;
    IMM_Nit(indn) = (ki_Nit*Nit(indn)./ (ki_Amm*Amm(indn) + ki_Nit*Nit(indn))).*IMM_net(indn);
    IMM_Nit(~indn) = 0;
    IMM_Nit = IMM_Nit(:);
            
%*************************************************************************
%     COMPUTE NITRATE AND AMMONIUM LEACHING FLUX [gN / m^3 / d]
%*************************************************************************        
    [LCH_nit_m2, TLCH_nit_m2] = CN_leach(PARAMS, SWITCHES, Nit, a_Nit, qq, sm);    
    [LCH_amm_m2, TLCH_amm_m2] = CN_leach(PARAMS, SWITCHES, Amm, a_Amm, qq, sm);    

    % Compute leaching in units of [gr / m^3 / d]
    LCH_nit = LCH_nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    LCH_amm = LCH_amm_m2./(dz_mm/1000);   % [gr/m^3/d]
    LCH_N_m2 = TLCH_nit_m2 + TLCH_amm_m2;  %[gr/m2/d]
        
%*************************************************************************
%     COMPUTE NITRATE AND AMMONIUM WATER UPTAKE FLUX [gN / m^2 / d]
%*************************************************************************            
    [UP_amm_m2, UP_nit_m2, UP_amm_all_m2, UP_nit_all_m2, F_HR_N] = CN_nuptake (VARIABLES, PARAMS, SWITCHES, a_Nit,...
         a_Amm, sm, layeruptake_all, TR_can_store);
    
    % Compute uptake in units of [gr / m^3 / d]
    UP_nit = UP_nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    UP_amm = UP_amm_m2./(dz_mm/1000);   % [gr/m^3/d]
    UP_nit_all = UP_nit_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]        
    UP_amm_all = UP_amm_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]
        
    % Compute total Nitrogen Uptake [gr/m2/d].        
    UP_N_m2 = sum(UP_amm_m2) + sum(UP_nit_m2);   %[gr/m^2/d]
                 
%*************************************************************************
%       CHANGES IN MINERAL CONCENTRATIONS
%*************************************************************************
    [Nitrif, Denitrif, Volat] = CN_nfluxes (VARIABLES, PARAMS, VERTSTRUC, SWITCHES, fSn, fTn);            
                        
%*************************************************************************
%      STATE UPDATES
%*************************************************************************
    dAmm_dt = MIN_net - IMM_Amm - Nitrif - Volat - UP_amm;  %[gr N / m^3 / d]
    dNit_dt = Nitrif - IMM_Nit - Denitrif - LCH_nit - UP_nit;
   
    Cl_new = Cl + dCl_dt*dtime/86400;          %[gC/m3]
    if SWITCHES.CN_type
        Ch_new = Ch + dCh_dt*dtime/86400;      %[gC/m3]
    else
        Ch_new = NaN;                          %[gC/m3]
    end        
    Cb_new = Cb + dCb_dt*dtime/86400;          %[gC/m3]
    Nl_new = Nl + dNl_dt*dtime/86400;          %[gC/m3] 
    Amm_new = Amm + dAmm_dt*dtime/86400;       %[gC/m3] 
    indneg = Amm_new < 0;                       
    if sum(indneg)==1
        stop =4;
    end
    Amm_new(indneg) = 0;                       %[gC/m3] 
    Nit_new = Nit + dNit_dt*dtime/86400;       %[gC/m3] 
    indneg = Nit_new < 0;       
    if sum(indneg)==1
        stop =4;
    end
    Nit_new(indneg) = 0;
    CNl_new = Cl_new./Nl_new;
    
    indcb = Cb_new<1; 
    Cb_new(indcb) = 1;
%*************************************************************************
    % STORE IN STRUCTURE 
%*************************************************************************
    CN_massbalance ();
    
    
%*************************************************************************
    % STORE IN STRUCTURE 
%*************************************************************************
        
    VARIABLES.Cl=Cl_new; 
    VARIABLES.Ch=Ch_new; 
    VARIABLES.Cb=Cb_new;
    VARIABLES.Nl=Nl_new; 
    VARIABLES.Amm=Amm_new;
    VARIABLES.Nit=Nit_new; 
    VARIABLES.CNl=CNl_new;
    VARIABLES.dCl_dt = dCl_dt;
    VARIABLES.dCb_dt = dCb_dt;
    
    VARIABLES.UP_amm=UP_amm;
    VARIABLES.UP_amm_all=UP_amm_all;
    VARIABLES.UP_amm_all_m2=UP_amm_all_m2;    
    VARIABLES.UP_nit=UP_nit;
    VARIABLES.UP_nit_all=UP_nit_all;
    VARIABLES.UP_nit_all_m2=UP_nit_all_m2;    
    VARIABLES.UP_N_m2=UP_N_m2;
    VARIABLES.F_HR_N=F_HR_N;
    
    VARIABLES.LCH_amm=LCH_amm;
    VARIABLES.LCH_nit=LCH_nit;
    VARIABLES.LCH_N_m2=LCH_N_m2;    
    VARIABLES.TLCH_amm_m2=TLCH_amm_m2;
    VARIABLES.TLCH_nit_m2=TLCH_nit_m2;

    VARIABLES.MIN_net=MIN_net;    
    VARIABLES.MIN_gross=MIN_gross;
    VARIABLES.IMM_net=IMM_net; 
    VARIABLES.IMM_gross=IMM_gross;
    VARIABLES.Nreg=Nreg;
    VARIABLES.DECl=DECl;
    VARIABLES.PHI=PHI;
    VARIABLES.phi=phi; 
    VARIABLES.fSd=fSd;
    VARIABLES.fTd=fTd;
    VARIABLES.mberrorN=mberrorN;
    VARIABLES.mberrorC=mberrorC;
    
    VARIABLES.ADD = ADD;
    VARIABLES.ADD_bio = ADD_bio;
    VARIABLES.ADD_ex = ADD_ex;
    VARIABLES.OUT = OUT;
    VARIABLES.ADD_net = ADD_net;
        
    