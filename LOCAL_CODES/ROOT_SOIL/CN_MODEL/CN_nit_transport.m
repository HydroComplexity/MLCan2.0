function [VARIABLES] = CN_nit_transport(PARAMS, SWITCHES, VARIABLES, VERTSTRUC,...
    CONSTANTS, FORCING, rootarq, CC, CCdd, aa, qq, sm, layeruptake_all)


ioncase = 1;
% Dereference blocks
% *************************************************************************
if SWITCHES.CN.Bioturbation;
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
    dz = [VARIABLES.SOIL.litterthickness ; VERTSTRUC.dzsmm/1000];%       dz = grid depth [m]
    % computation of zn only for conductivities
    nz = [VARIABLES.SOIL.litterthickness/2 ; VERTSTRUC.znsmm/1000 + VARIABLES.SOIL.litterthickness];%       dz = grid node [m]    
    zh = [VARIABLES.SOIL.litterthickness ; VERTSTRUC.zhsmm/1000 + VARIABLES.SOIL.litterthickness];%       dz = grid node [m]    
else
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
    dz = VERTSTRUC.dzsmm/1000; % dz = grid depth [m]
    % computation of zn only for conductivities
    nz = [VERTSTRUC.znsmm/1000];%       dz = grid node [m]    
    zh = [VERTSTRUC.znsmm/1000];%       dz = grid node [m]        
end
dtime = CONSTANTS.dtime;%       dtime = Time step [1800 s]
% *************************************************************************
CCnew = CC;
CCddnew = CCdd;
% Compute advection in the soil
% *************************************************************************
 [ADV_m2, ADV_m3, TLCADV_m2, TLCADV_m3] = CN_advection(nl_soil, dz, CCnew, aa, qq, sm, dtime);
 CCnew = CCnew + ADV_m3*dtime/86400;  % [g/m3/dtime] 
CCnew(CCnew<0) = 0;                    % Test for negative values that may arise 
% *************************************************************************

% Compute diffusion in the soil
% *************************************************************************
[DIFF_m3, DIFF_m2, TLCDFF_m2, TLCDFF_m3, De] = CN_diffusion(nl_soil, dz, nz,...
    zh, CCnew, sm, dtime, ioncase);
%DIFF_m3 = 0;
%DIFF_m2 = 0;
CCnew = CCnew + DIFF_m3*dtime/86400;  % [g/m3/dtime] 
CCnew(CCnew<0) = 0;                    % Test for negative values that may arise 


% Compute diffusion by water uptake
% *************************************************************************
[DIFFr_m2, DIFFr_m3, DIFFrdd_m3] = CN_diffwu (PARAMS, rootarq, nl_soil, dz, nz, zh, CCnew, CCdd, sm,...
        dtime, De);  
%DIFFr_m2 = 0;    
%DIFFr_m3 = 0;
%DIFFrdd_m3 = 0;
CCnew = CCnew + DIFFr_m3*dtime/86400;  % [g/m3/dtime] 
CCnew(CCnew<0) = 0;                    % Test for negative values that may arise 
CCddnew = CCdd + DIFFrdd_m3*dtime/86400;  % [g/m3/dtime]
CCddnew(CCddnew<0) = 0;                    % Test for negative values that may arise 


% Compute advection by water uptake
% ************************************************************************* 
[ADVr_m2, ADVr_m3, ADVrdd_m3] = CN_ruptakeadvection...
            (CCnew, CCddnew, PARAMS, SWITCHES, aa,...
         sm, dz, rootarq.RVID, layeruptake_all, dtime);       
%ADVr_m2 = 0;
%ADVr_m3 = 0;
%ADVrdd_m3 = 0;
CCnew = CCnew + ADVr_m3*dtime/86400;  % [g/m3/dtime] 
CCnew(CCnew<0) = 0;                    % Test for negative values that may arise 
CCddnew = CCddnew + ADVrdd_m3*dtime/86400;  % [g/m3/dtime] 
CCddnew(CCddnew<0) = 0;                    % Test for negative values that may arise 
% *************************************************************************

    
% Compute nitrogen uptake
% *************************************************************************
[UP_nit_m3, UP_nit_m2, UP_nit_all_m2] = CN_nituptake (PARAMS, rootarq, CCddnew, nl_soil, dtime, ioncase);     
%UP_nit_m3 = 0;
%UP_nit_m2 = 0;
%UP_nit_all_m2 = 0;
CCddnew = CCddnew - UP_nit_m3*dtime/86400;  % [g/m3/dtime] 
CCddnew(CCddnew<0) = 0;                    % Test for negative values that may arise 


% CHECK OF MASS BALANCE 
% *************************************************************************
% check in the soil first
% I. layer by layer
%errorly = (((ADV_m2 + frradvection_m2 + DIFF_m2 + rrdiff_m2)./dz/48 - (CCnew-CC))./CC)*100;
% II. entire layer

% check in the ENTIRE VOLUME 
RVID = sum(rootarq.RVID,2);
errormbnit = sum(CCnew.*dz + CCddnew.*RVID )-sum(CC.*dz + CCdd.*RVID) + (sum(UP_nit_m2) + TLCADV_m2)/48;

% COMPUTE DIFFERENCES IN CC AND CCdD
% *************************************************************************
% 
deltanit = (CCnew - CC)*86400/dtime;
deltanitdd = (CCddnew - CCdd)*86400/dtime;
fluxzd_m2 = -(DIFFr_m2 + ADVr_m2);                  %[g/m^2/d]  % Net loss from soil to ZD
fluxzd_m3 = -(DIFFr_m2 + ADVr_m2)./dz;              %[g/m^3/d]  % Net loss from soil to ZD
% you can check that 
% 1. sum(fluxzd_m2) - sum(UP_nit_m2) - sum(deltanitdd.*RVID) = 0
% 2. sum(fluxzd_m2) + sum(TLCADV_m2 + TLCDFF_m2) + sum(deltanit.*dz) = 0  % the 3 main losses from the soil
% - sum(UP_nit_m2)- sum(TLCADV_m2 + TLCDFF_m2) -(sum(deltanit.*dz) + sum(deltanitdd.*RVID))
% *************************************************************************

% STORE NEW VARIABLES
VARIABLES.ADV_nit_m2 = ADV_m2;                      %[g/m^2/d]
VARIABLES.ADV_nit_m3 = ADV_m3;                      %[g/m^3/d]
VARIABLES.ADVr_nit_m2 = ADVr_m2;                    %[g/m^2/d]
VARIABLES.ADVr_nit_m3 = ADVr_m3;                    %[g/m^3/d]
VARIABLES.ADVrdd_nit_m3 = ADVrdd_m3;                %[g/m^3/d]
VARIABLES.DIFF_nit_m3 = DIFF_m3;                    %[g/m^3/d]
VARIABLES.DIFF_nit_m2 = DIFF_m2;                    %[g/m^2/d]
VARIABLES.LCH_nit_m3 = TLCADV_m3 + TLCDFF_m3;       %[g/m^3/d]
VARIABLES.LCH_nit_m2 = TLCADV_m2 + TLCDFF_m2;       %[g/m^2/d]
VARIABLES.DIFFr_nit_m2 = DIFFr_m2;                      %[g/m^2/d]
VARIABLES.DIFFr_nit_m3 = DIFFr_m3;                      %[g/m^3/d]
VARIABLES.DIFFrdd_nit_m3 = DIFFrdd_m3;                  %[g/m^3/d]
VARIABLES.UP_nit_m3 = UP_nit_m3;                    %[g/m^3/d]
VARIABLES.UP_nit_m2 = UP_nit_m2;                    %[g/m^2/d]
VARIABLES.UP_nit_all_m2 = UP_nit_all_m2;            %[g/m^2/d]
VARIABLES.errormbnit = errormbnit; 
VARIABLES.deltanit = deltanit;          %[g/m^3/d]
VARIABLES.Nitdd = CCddnew;                       %[g/m^3/d]
VARIABLES.deltanitdd = deltanitdd;          %[g/m^3/d]
VARIABLES.fluxzd_m2 = fluxzd_m2 ;                   %[g/m^2/d]
VARIABLES.fluxzd_m3 = fluxzd_m3 ;                   %[g/m^3/d]



