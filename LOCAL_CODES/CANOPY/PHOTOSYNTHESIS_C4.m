function [Ph, An, Ph_limit, Jc, Jj, Js] = PHOTOSYNTHESIS_C4(VARIABLES, PARAMS, VERTSTRUC, sunlit, cntspecies) 
                                                   
 % Calculates leaf photosynthesis according to the Collatz et al (AJPP,
 % 1992) model for C4 photosynthesis
 %
 %  INPUTS:
 %      Qabs = absorbed PAR [umol/m^2 leaf area/s]
 %      Tl = leaf temperature [C]
 %      Ci = internal CO2 concentration [ppm]
 %      O = oxygen concentration
 %      Vmax = Maximum rubisco capacity [umol/m^2 leaf area/s]
 %      Q10 = proportional increase for parameter value for 10C increase in Tl [-]
 %      kk = initial slope of photosynthetic CO2 response [mol/m^2 leaf area/s]
 %      theta = curvature parameter [-]
 %      beta = curvature parameter [-]
 %      Vz = decay of photosynthetic capacity with canopy depth [-]
 %
 %  OUTPUTS:
 %      Ph = leaf photosynthesis            [umol/m^2 leaf area/s]
 %      An = net leaf CO2 uptake rate       [umol/m^2 leaf area/s]
 
 
 %*************************************************************************
 %                          DE-REFERENCE BLOCK
 %*************************************************************************
 % Dongkook Woo - Edit
%     if (sunlit)  
%         Qabs = VARIABLES.CANOPY.PARabs_sun;
%         Tl =  VARIABLES.CANOPY.Tl_sun;
%         Ci =  VARIABLES.CANOPY.Ci_sun;
%     else
%         Qabs = VARIABLES.CANOPY.PARabs_shade;
%         Tl =  VARIABLES.CANOPY.Tl_shade;
%         Ci =  VARIABLES.CANOPY.Ci_shade;
%     end
    if (sunlit)  
        Qabs = VARIABLES.CANOPY.PARabs_sun;
        Tl =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);
        Ci =  VARIABLES.CANOPY.Ci_sun(:,cntspecies); 
    else
        Qabs = VARIABLES.CANOPY.PARabs_shade;
        Tl =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);
        Ci =  VARIABLES.CANOPY.Ci_sun(:,cntspecies); 
    end
 
%     Vmax = PARAMS.Photosyn.Vmax_C4(1);        
%     kk = PARAMS.Photosyn.kk_C4(1);
%     Q10 = PARAMS.Photosyn.Q10_C4(1);
%     theta = PARAMS.Photosyn.theta_C4(1);
%     beta = PARAMS.Photosyn.beta_C4(1);
%     Rd = PARAMS.Photosyn.Rd_C4(1);
%     al = PARAMS.Photosyn.al_C4(1);    
%     nl_can = PARAMS.CanStruc.nl_can;
    Vmax = PARAMS.Photosyn.Vmax_C4(cntspecies);        
    kk = PARAMS.Photosyn.kk_C4(cntspecies);
    Q10 = PARAMS.Photosyn.Q10_C4(cntspecies);
    theta = PARAMS.Photosyn.theta_C4(cntspecies);
    beta = PARAMS.Photosyn.beta_C4(cntspecies);
    Rd = PARAMS.Photosyn.Rd_C4(cntspecies);
    al = PARAMS.Photosyn.al_C4(cntspecies);    
    nl_can = PARAMS.CanStruc.nl_can;
    
    Vz = VARIABLES.CANOPY.Vz;
    
    nvinds_all = VERTSTRUC.nvinds_all;
    nvinds1 = nvinds_all{cntspecies}; 
    
    nvinds_all = VERTSTRUC.nvinds_all;
    nvinds1 = nvinds_all{cntspecies}; 
    nvinds2 = find(isnan(Tl));
    nvinds = sort(unique([nvinds1 ; nvinds2]));
        
    all = (1:nl_can)';
    vinds = all(~ismember(all,nvinds));
        

 %*************************************************************************
 % ALLOCATE MEMORY
 M = nan(nl_can,1);
 Ph = nan(nl_can,1);
  
 Vmaxs = Vz * Vmax;
 
 % Eqns 5B
 Q10s = Q10.^((Tl-25)./10);
 VT = (Vmaxs .* Q10s)./( (1+exp(0.3*(13-Tl))).*(1+exp(0.3*(Tl-36))) );
 RT = Rd.*Q10s./(1 + exp(1.3*(Tl-55)));
 kT = kk*Q10s;
 
 % Eqn 2B
 aa = theta * ones(size(Qabs));
 bb = -(VT + al*Qabs);
 cc = VT*al.*Qabs;
 for ii = vinds';%1:length(aa)
     Mroots = roots([aa(ii), bb(ii), cc(ii)]);
     M(ii) = min(Mroots);
 end
 M = M(:);
  
 % Eqn 3B
 aa = beta * ones(size(Qabs));
 bb = -(M+kT.*Ci);
 cc = M.*kT.*Ci;
  for ii = vinds';%1:length(aa)
     Aroots = roots([aa(ii), bb(ii), cc(ii)]);
     Ph(ii) = min(Aroots);
  end
  
Ph = Ph(:);
An = Ph - RT;
Ph(nvinds) = 0;
An(nvinds) = 0;   


% Calculate Limiting Rates
    % Light-Limited (Eqn 3)
    Jj = al*Qabs;  
    % CO2 Limited (Eqn 4)
    Jc = kT .* Ci;
    % (Eqn 5)
    Js = VT;
        
    Ph_limit = NaN(size(Jc));
    Jcinds = find(Jc<Jj & Jc<Js);
    Ph_limit(Jcinds) = 1;
    Jjinds = find(Jj<=Jc & Jj<=Js);
    Ph_limit(Jjinds) = 2;
    Jsinds = find(isnan(Ph_limit));
    Ph_limit(Jsinds) = 3;
    
    Ph(isnan(Tl)) = 0;
    An(isnan(Tl)) = 0;

    
% Dongkook Woo - Edit End
 