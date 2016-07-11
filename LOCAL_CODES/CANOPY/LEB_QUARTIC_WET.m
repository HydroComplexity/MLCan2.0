function [Tl, H, LE, gv, vinds, nvinds] = LEB_QUARTIC_WET (VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, sunlit, cntspecies)
                                  
%=========================================================================
% Solve exact energy balance equation as presented in Nikolov et al (1995)
% using an expansion of the saturation vapor pressure at leaf temperature
% using a fourth order polynomial. 
% 
% Latent Heat occurs only from the wet part of leaves. Water that is
% located over the leaves.
%
%   Rabs = LE + H + Me + LWout
%
%   where   Rabs is the bi-directional absorbed radiation
%           Me is the energy stored in biochemical reactions
%           LW is the emitted LW radiation by the leaf layer
%
%   all terms are in [W / m^2 leaf area]
%
% Written By: Darren Drewry, Modified by Juan Quijano
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES        % VARIABLES structure
%       PARAMS           % PARAMS structure
%       VERTSTRUC        % VERTSTRUC structure
%       CONSTANTS        % CONSTANTS structure
%       sunlit           % [1] Sunlit [0] Shade
%       cntspecies       % [] Number of species
%------------------------- Output Variables ------------------------------
%       Tl               % [C] Leaf temperature profile  
%       H                % [W / m^2 leaf] Sensible heat flux profile   
%       LE               % [W / m^2 leaf] Latent heat flux profile   
%       gv               % [W / m^2 leaf area] Total vapor conductance profile                
%========================================================================              
%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES FOR EACH SPECIES 
%*************************************************************************
    nl_can = PARAMS.CanStruc.nl_can;    
    if (sunlit)
        Rabs = VARIABLES.CANOPY.Rabs_sun;                                  % [W/m^2 leaf area] Rabs - Total absorbed radiation Sunlit Fraction
        An =  VARIABLES.CANOPY.An_sun(:,cntspecies);                       % [umol/m^2 la/s] An - net CO2 uptake  
        gsv =  VARIABLES.CANOPY.gsv_sun(:,cntspecies);                     % [mol/m^2 la/s] gsv - stomatal conductance to vapor transport 
        gbv = VARIABLES.CANOPY.gbv_sun(:,cntspecies);                      % [mol/m^2 la/s] gbv - boundary layer conductance to vapor transport
        gbh = VARIABLES.CANOPY.gbh_sun(:,cntspecies);                      % [mol/m^2 la/s] gbh - boundary layer conductance to heat transport 
        leaffrac = VARIABLES.CANOPY.fsun;                                  % [] Fraction of LAI in each species 
        LAIfrac = VARIABLES.CANOPY.LAIsun;                                 % [] Fraction of LAI in sunlit     
    else
        Rabs = VARIABLES.CANOPY.Rabs_shade;                                % [W/m^2 leaf area] Rabs - Total absorbed radiation shade Fraction
        An =  VARIABLES.CANOPY.An_shade(:,cntspecies);                     % [umol/m^2 la/s] An - net CO2 uptake  
        gsv =  VARIABLES.CANOPY.gsv_shade(:,cntspecies);                   % [mol/m^2 la/s] gsv - stomatal conductance to vapor transport            
        gbv =  VARIABLES.CANOPY.gbv_shade(:,cntspecies);                   % [mol/m^2 la/s] gbv - boundary layer conductance to vapor transport
        gbh =  VARIABLES.CANOPY.gbh_shade(:,cntspecies);                   % [mol/m^2 la/s] gbh - boundary layer conductance to heat transport 
        leaffrac = VARIABLES.CANOPY.fshade;                                % [] Fraction of LAI in each species     
        LAIfrac = VARIABLES.CANOPY.LAIshade;                               % [] Fraction of LAI in shade 
    end
    wetfrac = VARIABLES.CANOPY.wetfrac;                                    % [] Fraction of canopy that is wet  
    
    nvinds_all = VERTSTRUC.nvinds_all;                                     % [] Indicator of Layers where LAI is zero of lower for all species 
    nvinds1 = nvinds_all{cntspecies};                                      % [] Indicator 1 for species  
    nvinds2 = find(LAIfrac <= 0);                                          % [] Indicator 2 for LAI 
    nvinds3 = find(wetfrac == 0);                                          % [] Indicator 3 for LAI
    
    nvinds = sort(unique([nvinds1 ; nvinds2 ; nvinds3]));                  %[] Final Indicator 
    
    all = (1:nl_can)';                                                     %[] Vector of numbers 
    vinds = all(~ismember(all,nvinds));                                    %[] Fill vinds  
%    vinds_all = VERTSTRUC.vinds_all;
%    vinds = vinds_all{cntspecies}; 
    % check logic
    nvinds = nvinds(:);                                                    %[] Convert to column vector 
    vinds = vinds(:);                                                      %[] Convert to column vector 
    
    LWmethod = PARAMS.LWmethod;                                            %[] Method to compute in LW_ATTENUATION FUNCITON LW 0. Original Darren. 1. Corrected  
    retainLW = PARAMS.retainLW;                                            %[] If LWmethod = 1,  retainLW is an extra LW retain the canopy 
    taud = VARIABLES.CANOPY.taud;                                          %[] Transmitted fraction  
                                        
    eaz = VARIABLES.CANOPY.EAz;                                            %[kPa] Vapor Pressure Air, for all layers  
    Taz = VARIABLES.CANOPY.TAz;                                            %[C] Air Temperature, for all layers 
    Paz = VARIABLES.CANOPY.PAz;                                            %[umol/mol] Atmospheric Concentration of CO2, for all layers 

    epsv = PARAMS.Rad.epsv;                                                %[] Vegetation emissivity  
    LWfact = PARAMS.CanStruc.LWfact(cntspecies);                           %[] Leaf factor for emission of LW 
    Hfact = PARAMS.CanStruc.Hfact(cntspecies);                             %[] Leaf factor for emission of H 
    LEfact = PARAMS.CanStruc.LEfact(cntspecies);                           %[] Leaf factor for emission of LE 
    
    
    Lv = CONSTANTS.Lv;                                                     %[J / mol] latent heat of vaporization 
    boltz = CONSTANTS.boltz;                                               %[W / m^2 / K^4] Stefan-Boltzmann constant 
    cp = CONSTANTS.cp_mol;                                                 %[J/mol/K] specific heat of air at constant pressure
                                      
%************************************************************************* 
%*************************************************************************
    %          ALLOCATE MEMORY
    Tl = nan(nl_can,1);
%*************************************************************************


    gbh = gbh * 0;
    
    % Energy stored in biochemical reactions (Me)
        Me = 0.506 * An;                    % [W / m^2] - energy stored in biochemical reactions
        
    % Total vapor conductance (gv)
        gv = gbv;   % [mol/m^2/s]

    
    aa = 273.15;
    if LWmethod == 0
        bb = LWfact * epsv * boltz * (1-taud) .* leaffrac ./ LAIfrac;
    else
        bb = LWfact * epsv * boltz;
    end
%    bb(LAIfrac==0)=0;
    bb = bb(:);
    dd = Hfact * cp * gbh;
    ee = LEfact * gv * Lv ./ Paz;
    
    c1 = (5.82436*10^-4) / 1000;
    c2 = (1.5842*10^-2) / 1000;
    c3 = 1.55186 / 1000;
    c4 = 44.513596 / 1000;
    c5 = 607.919 / 1000;
    
    % Terms of Quartic
    T4 = ee.*c1 + bb;                                       % * Tl^4
    T3 = ee.*c2 + 4*aa.*bb;                                 % * Tl^3
    T2 = ee.*c3 + 6*bb.*aa.^2;                              % * Tl^2
    T1 = ee.*c4 + dd + 4*bb.*aa.^3;                         % * Tl
    T0 = ee.*(c5-eaz) - dd.*Taz + bb.*aa.^4 - Rabs + Me;    % constant term
            
    
for ii = vinds'%1:length(Taz)

    Tl_roots = roots([T4(ii), T3(ii), T2(ii), T1(ii), T0(ii)]);    
    rcnt = 0;
%   Initialize var vector with zeros, later we select  max(var);
     var=zeros(4,1);%*inf;
     for jj = 1:length(Tl_roots)
         if (isreal(Tl_roots(jj)));% && Tl_roots(jj)>0)
             rcnt = rcnt + 1;
             var(rcnt)=Tl_roots(jj);            
         end
     end
    Tl_all = max(var(1:rcnt)); 
        
    % Error Check
    if (rcnt == 0)
        Tlout = fzero(@(Tlin) LEB_residual(ee(ii),dd(ii),bb(ii),Tlin,Taz(ii),eaz(ii),Rabs(ii),Me(ii)), Taz(ii));
        if isnan(Tlout)
            disp(['NO REAL, POSITIVE ROOTS FOUND IN LEB WET!!!'])
            Tl(ii) = Taz(ii);                    
        else
        Tl(ii) = Tlout;
        [res LEout Hout LWout Tlout]=LEB_residual(ee(ii),dd(ii),bb(ii),Tlout,Taz(ii),eaz(ii),Rabs(ii), Me(ii));        
        end
    else    
        Tl(ii) = Tl_all;
    end    
end

%estarTl = 0.611 * exp( (17.502*Tl) ./ (Tl + 240.97) );
estarTl = 0.6107 * exp( (17.38*Tl) ./ (Tl + 239) );
LE = LEfact.*(Lv*gv./Paz).*(estarTl - eaz);
H = Hfact.*cp.*gbh.*(Tl - Taz);

H = H(:);
LE = LE(:);
Tl = Tl(:);
gv = gv(:);

% Set H and LE to zero when there is not LAI
H(nvinds) = 0;
LE(nvinds) = 0;



% Compute normalized energy balance terms for checking:
%     LE_check = ee.*(c1.*Tl.^4 + c2.*Tl.^3 + c3.*Tl.^2 + c4.*Tl + c5 - eaz);
%     H_check = dd.*(Tl - Taz);
%     LWout_check = bb.*(Tl.^4 + 4*aa.*Tl.^3 + (6*aa.^2).*Tl.^2 + (4*aa.^3).*Tl + aa.^4);




