function [dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, VARIABLES]...
    = CN_dynamics (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo)


Cl = VARIABLES.Cl;%       Cl = carbon concentration in litter pool [gC / m^3]
Ch = VARIABLES.Ch;%       Ch = carbon concentration in humus pool [gC / m^3]
Cb = VARIABLES.Cb;%       Cb = carbon concentration in biomass pool [gC / m^3]
kd = PARAMS.CN.kd;%       kd = slope of linear dependence of biomass death on biomass concentration [d^-1]
kl = PARAMS.CN.kl;%       kl = rate of decomposition of the litter pool [m^3 d / gC]
kh = PARAMS.CN.kh;%       kh = rate of decomposition of the humus pool [m^3 d / gC]
koae = PARAMS.CN.koae;% Organic Assimilation Efficiency parameter
rr = PARAMS.CN.rr;%       rr = fraction of decomposed organic carbon that goes to respiration [-]
rhmin = PARAMS.CN.rhmin;%       rhmin = minimum value of isohumic coefficient [-]
CNl = VARIABLES.CNl;%       CNl = carbon/nitrogen ratio of litter pool [gC / gN]
CNb = PARAMS.CN.CNb;%       CNb = carbon/nitrogen ratio of biomass pool [gC / gN]
CNh = PARAMS.CN.CNh;%       CNh = carbon/nitrogen ratio of humus pool [gC / gN]
if SWITCHES.CN.Bioturbation;
    nl_soil=PARAMS.nl_soil + 1;%       nl_soil = # soil layers
else
    nl_soil=PARAMS.nl_soil ;%       nl_soil = # soil layers
end

% LITTER POOL
% Rate of microbial death and return to litter pool
BD = kd .* Cb;             
% Decomposition rate of litter
ratel = phi .* fTd .* fSd .* kl .* Cb;
DECl = ratel .* Cl;

if SWITCHES.CN.Bioturbation 
    dCl_dt = ADD -OUT + BD - DECl;    %[gC/m3/d]
    dNl_dt = ADD./CNa - OUT./CNo + BD./CNb - DECl./CNl;
    dNl_dt = dNl_dt(:);          %[gN/m3/d]
else    
    dCl_dt = ADD + BD - DECl;    %[gC/m3/d]
    dNl_dt = ADD./CNa + BD./CNb - DECl./CNl;
    dNl_dt = dNl_dt(:);          %[gN/m3/d]
end

% HUMUS POOL
% Decomposition rate of humus        
rh = min(rhmin, CNh./CNl);
        
if SWITCHES.CN_type 
    % Decomposition rate of humus
    rateh = phi .* fTd .* fSd .* kh .* Cb;
    DECh = rateh .* Ch;

    dCh_dt = rh.*DECl - DECh;     %[gC/m3/d]

    Nladd = min(DECl./CNl, rh.*DECl./CNh);
    dNh_dt = Nladd - DECh./CNh;   %[gN/m3/d]
else
    dCh_dt = nan(nl_soil,1);
    dNh_dt = nan(nl_soil,1);
end            
            
% BIOMASS POOL        
if SWITCHES.CN_type == 1
    dCb_dt = (1-rh-rr).*DECl + (1-rr).*DECh - BD;                      %[gC/m3/d]
    dNb_dt = (1-rh.*CNl./CNh).*DECl./CNl + DECh./CNh - BD./CNb - PHI;  %[gN/m3/d]
elseif SWITCHES.CN_type == 2             
     dCb_dt = (1-rr-rr*rh).*DECl - BD;                      %[gC/m3/d]    
     dNb_dt = DECl./CNl - BD./CNb - PHI;  %[gN/m3/d]
else        
     dCb_dt = (1-rr).*DECl - BD;                      %[gC/m3/d]    
     dNb_dt = koae*DECl./CNl - BD./CNb - PHI;  %[gN/m3/d]        
end
  

VARIABLES.dCl_dT = dCl_dt;              % rate of change in gr/[m3]
VARIABLES.BD = BD;                      % rate of change in gr/[m3]

