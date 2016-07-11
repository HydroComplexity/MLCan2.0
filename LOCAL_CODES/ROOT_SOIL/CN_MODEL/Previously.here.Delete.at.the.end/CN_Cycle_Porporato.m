function [Cl_new, Ch_new, Cb_new, Nl_new, Amm_new, Nit_new, CNl_new, ...
          UP_amm, UP_amm_all, UP_nit, UP_nit_all, UP_N, LCH_amm, LCH_nit, TLCH_amm, TLCH_nit, MIN, IMM] = ...
        CN_Cycle_Porporato(Cl, Ch, Cb, Nl, Amm, Nit, ...
                           CNa, CNl, CNb, CNh, ADD,   ...
                           sm, porsl, sm_min, Ts, Ts_max, Ts_reduc_on, ...
                           dz_mm, dtime, nl_soil, ...
                           a_Amm, a_Nit, kd, kl, kh, kn, ...
                           ku_Amm, ku_Nit, ki_Amm, ki_Nit, rhmin, rr, qq, leach_on,...
                           layeruptake,layeruptake_all,rootfr, VARIABLES, nspecies)
                                    
                       
                        
%
%   Calculate Dynamics of Soil Carbon and Nitrogen
%
%   INPUTS:
%       INITIAL CONDITIONS:
%       Cl = carbon concentration in litter pool [gC / m^3]
%       Ch = carbon concentration in humus pool [gC / m^3]
%       Cb = carbon concentration in biomass pool [gC / m^3]
%       Nl = nitrogen concentration in litter pool [gN / m^3]
%       Amm = ammonium concentration in soil [gN / m^3]
%       Nit = nitrate concentration in soil [gN / m^3]
%
%       CNl = carbon/nitrogen ratio of litter pool [gC / gN]
%       CNb = carbon/nitrogen ratio of biomass pool [gC / gN]
%       CNh = carbon/nitrogen ratio of humus pool [gC / gN]
%
%       sm = soil moisture profile [m^3 / m^3]
%       porsl = porosity profile [m^3 / m^3]
%       sm_min = minimum soil moisture content [m^3 / m^3]
%       Ts = soil temperature profile [C]
%       Ts_max = Maximum temperature [C]. Used to compute environmental factors
%       Ts_reduc_on = soil temperature reduction function switch [-]
%
%       dz_mm = grid depth [mm]
%       dtime = Time step [1800 s]
%       nl_soil = # soil layers
%       
%       a_Amm = fraction of dissolved ammonium [-]
%       a_Nit = fraction of dissolved nitrate [-]
%       kd = slope of linear dependence of biomass death on biomass
%           concentration [d^-1]
%       kl = rate of decomposition of the litter pool [m^3 d / gC]
%       kh = rate of decomposition of the humus pool [m^3 d / gC]
%       kn = rate of nitrification [m^3 d / gC]
%       ku_Amm = 
%       ku_Nit = 
%       ki_Amm = partitioning coefficient for ammonium immobilization [m^3 d / gC]
%       ki_Nit = partitioning coefficient for nitrate immobilization [m^3 d / gC]
%       rhmin = minimum value of isohumic coefficient [-]
%       rr = fraction of decomposed organic carbon that goes to respiration [-]
%       qq = water fluxes between layers [mm / s]
%       leach_on = Leaching ON or OF [1 0]
%       layeruptake = Uptake from each layer [mm/s]
%       rootfr  = Root distribution of fine roots.
%       VARIABLES = Structure with information of some variables

%       DEM_Amm = plant demand for ammonium under active uptake [gN / m^3 / d]
%       DEM_Nit = plant demand for nitrate under active uptake [gN / m^3 / d]
%
%   OUTPUTS:
%       Cl_new = Updated Carbon concentration in litter pool [gC / m^3]
%       Ch_new = Updated Carbon concentration in humus pool [gC / m^3]
%       Cb_new = Updated Carbon concentration in biomass pool [gC / m^3]
%       Nl_new = Updated nitrogen concentration in litter pool [gN / m^3]
%       Nh_new = Updated Nitrogen concentration in humus pool [gC / m^3]
%       Nb_new = Updated Nitrogen concentration in biomass pool [gC / m^3]
%       Amm_new = Updated ammonium concentration in soil [gN / m^3]
%       Nit_new = Updated nitrate concentration in soil [gN / m^3]
%       CNl_new = Updated of the C/N ratio concentration in litter [gC / gN]
%       UP_amm = Uptake N in the of Ammonium [gN / m3 / d]
%       UP_nit = Uptake of N in the form of Nitrate [gN / m3 / d]
%       UP_N = Total Nitrogen Uptake [gr/m^2/d]
%       LCH_amm = Leaching of ammoniun in each layer[gN/m^3/d]
%       LCH_nit = Leaching of nitrige in each layer [gN/m^3/d]
%       TLCH_amm = Leaching of ammoniun out of the root zone [gr/m^2/d]
%       TLCH_nit = Leaching of nitrate out of the root zone [gr/m^2/d]
%       MIN = Net N mineralized
%       IMM = Net N immobilized

%*************************************************************************
% %                       De reference block  
% % Above ground
% CNaa = VARIABLES.NITROALLO.CNad ;                  % C/N ratio  above ground at end of time step [gr/m^2]
% ADDa = VARIABLES.NITROALLO.TBMla;                 % Total biomass Drop in this time step as litter from above [gr / m^2]     
% ADDa = ADDa * 86400/dtime;                        % Convert from [gr/m^2/dtime s] to [gr/m^2/d] 
% % Below ground
% CNab = 50;         
% ADDb = 1.04;                                      % Total biomass Drop in this time step as drop from below [gr / m^2]  [gC/m^2/d]
%*************************************************************************
%                       PREALLOCATE VECTORS 

IMM_Amm = zeros(nl_soil,1);
IMM_Nit = zeros(nl_soil,1);
IMM = zeros(nl_soil,1);
MIN = zeros(nl_soil,1);
kind = zeros(nl_soil,1);
PHI = zeros(nl_soil,1);
LCH_nit = zeros (nl_soil,1);
UP_amm=zeros(nl_soil,1);
UP_nit=zeros(nl_soil,1);
LCH_amm=zeros(nl_soil,1);
%*************************************************************************
%                       CHANGE OF UNITS AND DE-REFERENCE BLOCK

qq=qq*86400/1000;  % Change units from [mm/s] to [m/d]

layeruptake=layeruptake*86400/1000;  % Change units from [mm/s] to [m/d]

%*************************************************************************
%                   COMPUTE THE LITTER TO ADD AGAIN
%*************************************************************************
% % Litter from plant below ground is distributed based on root distribution
%     ADD2=ADDb.*rootfr./(dz_mm/1000); 
% % Litter from plant above ground is distributed in uniformly in first 2 layers     
%     ADD1=zeros(nl_soil,1); 
%     ADD1(1)= ADDa/2/(dz_mm(1)/1000);    
%     ADD1(2)= ADDa/2/(dz_mm(2)/1000);    
% % Supersposition of litter from both
%     ADD=ADD1+ADD2;
% % Compute CNadd vector based on weights 
%     CNa=(ADD1*CNaa+ADD2*CNab)./(ADD1+ADD2); %

%*************************************************************************

%*************************************************************************
%                   COMPUTE ENVIRONMENTAL FACTORS
%*************************************************************************
    % Relative Water Content
    rwc = (sm-sm_min)./(porsl-sm_min);
    fc= 0.39;     %Saxton et al (1986)
    
    % for Decomposition
        % Temp
          %  fTd = exp(3.36 * (Ts - Ts_max) ./ (Ts + 31.79));
            fTd = ones(nl_soil,1);
        % Moisture
          %  fSd = 1./(1+4*exp(-6*rwc));
            % using equations from porporato et al 2003
             findn=rwc<fc;
             find=rwc>=fc;
             fSd(findn)=rwc(findn)./fc;
             fSd(find)=fc./rwc(find);
             fSd=fSd(:);
            
    % for Nitrification
        fTn = fTd;
        fSn = fSd;
    % for Volatilization
        fTv = fTd;
    % for Denitrification
        % Temp
            fTdn = fTd;
        % Moisture
            finds = find( sm./porsl <= 0.8 );
            fSdn(finds) = 0;
            finds = find( (sm./porsl>0.8) & (sm./porsl<=0.9) );
            fSdn(finds) = -1.6 + 2*sm(finds)./porsl(finds);
            finds = find( sm./porsl > 0.9 );
            fSdn(finds) = -7 + 8*sm(finds)./porsl(finds);

%*************************************************************************
%     COMPUTE RATES OF CHANGES IN POOLS CONCENTRATION. (LITTER, HUMUS, BIOMASS) 
%*************************************************************************            

    phi = 1;    % ""
    
    % LITTER POOL
        % Rate of microbial death and return to litter pool
        BD = kd .* Cb;         
    
        % Decomposition rate of litter
        ratel = phi .* fTd .* fSd .* kl .* Cb;
        DECl = ratel .* Cl;
        
        dCl_dt = ADD + BD - DECl;    %[gC/m3/d]
        cninds = find(CNa>0);
        dNl_dt(cninds) = ADD(cninds)./CNa(cninds) + BD(cninds)./CNb - DECl(cninds)./CNl(cninds);
        cninds = find(CNa<=0);
        dNl_dt(cninds) = BD(cninds)./CNb - DECl(cninds)./CNl(cninds);
        dNl_dt = dNl_dt(:);          %[gN/m3/d]

    % HUMUS POOL
        % Decomposition rate of humus
        rateh = phi .* fTd .* fSd .* kh .* Cb;
        DECh = rateh .* Ch;
        
        rh = min(rhmin, CNh./CNl);
        dCh_dt = rh.*DECl - DECh;     %[gC/m3/d]
                
        Nladd = min(DECl./CNl, rh.*DECl./CNh);
        dNh_dt = Nladd - DECh./CNh;   %[gN/m3/d]


    % BIOMASS POOL
        dCb_dt = (1-rh-rr).*DECl + (1-rr).*DECh - BD;                      %[gC/m3/d]
        dNb_dt = (1-rh.*CNl./CNh).*DECl./CNl + DECh./CNh - BD./CNb - PHI;  %[gN/m3/d]
        [y,in]=max(dNb_dt);
        if in~=1
            stop=3;
        end
%*************************************************************************
%     COMPUTE NET MINERALIZATION OR NET IMMOBILIZATION
%*************************************************************************
% Preallocate vectors        
    for ii = 1:nl_soil    
        
        % Eqn. 18. Compute PHI
        phi = 1;
        PHI(ii) = phi*fTd(ii)*fSd(ii)*Cb(ii)*(kh(ii)*Ch(ii)*(1./CNh - (1-rr)./CNb) + kl(ii)*Cl(ii)*(1./CNl(ii) - rh(ii)./CNh - (1-rh(ii)-rr)./CNb));
        if (PHI(ii)>0)   % NET MINERALIZATION
                        
            IMM_Amm(ii) = 0;
            IMM_Nit(ii) = 0;
            IMM(ii) = 0;
            
            MIN(ii) = PHI(ii);
            
            kind(ii) = 1;
            
        else        % NET IMMOBILIZATION        
            
            MIN(ii) = 0;
            
            IMMmax = fTd(ii)*fSd(ii)*(ku_Amm*Amm(ii) + ku_Nit*Nit(ii));
            
            if (abs(PHI(ii))<IMMmax) % UNRESTRICTED IMMOBILIZATION
                
                phi = 1;
                IMM(ii) = abs(PHI(ii));
                
                kind(ii) = 2;   
                
            else            % RESTRICTED IMMOBILIZATION
                
                num = IMMmax;
                denom = fTd(ii)*fSd(ii)*Cb(ii)*(kh(ii)*Ch(ii)*(1/CNh - (1-rr)/CNb) + kl(ii)*Cl(ii)*(1/CNl(ii) - rh(ii)/CNh - (1-rh(ii)-rr)/CNb));
                phi = -num / denom;
                
                IMM(ii) = IMMmax;
                
                kind(ii) = 3;
                
                % Re-compute Pool Changes With New 'phi' Value
                    % Decomposition rate of litter
                        ratel(ii) = phi * fTd(ii) * fSd(ii) * kl(ii) * Cb(ii);
                        DECl(ii) = ratel(ii) * Cl(ii);

                        dCl_dt(ii) = ADD(ii) + BD(ii) - DECl(ii);
                        if (CNa(ii)>0)
                            dNl_dt(ii) = ADD(ii)/CNa(ii) + BD(ii)/CNb - DECl(ii)/CNl(ii);
                        else
                            dNl_dt(ii) =                   BD(ii)/CNb - DECl(ii)/CNl(ii);
                        end


                    % HUMUS POOL
                        % Decomposition rate of humus
                        rateh(ii) = phi * fTd(ii) * fSd(ii) * kh(ii) * Cb(ii);
                        DECh(ii) = rateh(ii) * Ch(ii);

                        rh(ii) = min(rhmin, CNh/CNl(ii));
                        dCh_dt(ii) = rh(ii)*DECl(ii) - DECh(ii);

                        Nladd(ii) = min(DECl(ii)/CNl(ii), rh(ii)*DECl(ii)/CNh);
                        dNh_dt(ii) = Nladd(ii) - DECh(ii)/CNh;


                    % BIOMASS POOL
                        dCb_dt(ii) = (1-rh(ii)-rr)*DECl(ii) + (1-rr)*DECh(ii) - BD(ii);        
                        dNb_dt(ii) = (1-rh(ii)*CNl(ii)/CNh)*DECl(ii)/CNl(ii) + DECh(ii)/CNh - BD(ii)/CNb - PHI(ii);
        
            end
            
            % PARTITION IMMOBILIZATION BETWEEN NO3- AND NH4+
            IMM_Amm(ii) = (ki_Amm*Amm(ii) / (ki_Amm*Amm(ii) + ki_Nit*Nit(ii)) ) * IMM(ii);
            IMM_Nit(ii) = (ki_Nit*Nit(ii) / (ki_Amm*Amm(ii) + ki_Nit*Nit(ii)) ) * IMM(ii);
            
        end
        PHI(ii) = MIN(ii) - IMM(ii);

%*************************************************************************
%     COMPUTE NITRATE AND AMMONIUM LEACHING FLUX [gN / m^3 / d]
%*************************************************************************        
        
        %   positive value = loss (subtracted from budget)
        %   NOTE: first H2O flux value is from infiltration into soil column
            LCH_nit(ii) = 0;
            if (leach_on)
                if (ii>1)   % transport through top layer interface
                    if (qq(ii)>0)  % downward transport into cell
                        LCH_nit(ii) = LCH_nit(ii) - (a_Nit/sm(ii-1)*Nit(ii-1)*qq(ii));% [gr/m^2/d]
                    else           % upward transport out of cell
                        LCH_nit(ii) = LCH_nit(ii) - (a_Nit/sm(ii)*Nit(ii)*qq(ii));% [gr/m^2/d]
                    end
                end
                % transport through bottom layer interface
                if (qq(ii+1)>0)  % downward transport out of cell
                    LCH_nit(ii) = LCH_nit(ii) + (a_Nit/sm(ii)*Nit(ii)*qq(ii+1));% [gr/m^2/d]
                else           % upward transport into cell
                    LCH_nit(ii) = LCH_nit(ii) + (a_Nit/sm(ii+1)*Nit(ii+1)*qq(ii+1));% [gr/m^2/d]
                end
                if (ii==nl_soil)
                    TLCH_nit = LCH_nit(ii) + (a_Nit/sm(ii)*Nit(ii)*qq(ii+1));% [gr/m^2/d]
                    TLCH_amm = 0;
                end   
            end            
            
            % Assume ammonium does not leach due to bonding
            LCH_amm(ii) = 0;
            
    end
    
    % Compute leaching in units of [gr / m^3 / d]
    LCH_nit = LCH_nit./(dz_mm/1000);   % [gr/m^3/d]
    LCH_amm = LCH_amm./(dz_mm/1000);   % [gr/m^3/d]
           
    IMM_Amm = IMM_Amm(:);
    IMM_Nit = IMM_Nit(:);
    MIN = MIN(:);
    PHI = PHI(:);
    LCH_nit = LCH_nit(:);
    LCH_amm = LCH_amm(:);
    
%*************************************************************************
%     COMPUTE NITRATE AND AMMONIUM WATER UPTAKE FLUX [gN / m^2 / d]
%*************************************************************************            
    
    % PLANT UPTAKE
        % AMMONIUM
        for i=1:1:nspecies
            layeruptn=layeruptake(:,i);
            wtind=layeruptn>0;
            wtnind=layeruptn<=0;
            % First for the case in which the flow is from the soil to the root
            UP_amm_n(wtind) = (a_Amm./sm(wtind).*Amm(wtind)).*layeruptn(wtind);% [gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system. 
            conAmmroots=(sum(a_Amm./sm(wtind).*Amm(wtind).*layeruptn(wtind)))/(sum(layeruptn(wtind))); 
            UP_amm_n(wtnind) = (conAmmroots).*layeruptn(wtnind);% [gr/m^2/d]           
            UP_amm_n = max(0, UP_amm_n);
            UP_amm_n=UP_amm_n(:);
            UP_amm_all(:,i) = UP_amm_n;

            % NITRATE
            % First for the case in which the flow is from the soil to the root        
            UP_nit_n(wtind) = (a_Nit./sm(wtind).*Nit(wtind)).*layeruptn(wtind);%[gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system. 
            conNitroots=(sum(a_Nit./sm(wtind).*Nit(wtind).*layeruptn(wtind)))/(sum(layeruptn(wtind)));
            UP_nit_n(wtnind) = (conNitroots).*layeruptn(wtnind);% [gr/m^2/d];           
            UP_nit_n = max(0, UP_nit_n);
            UP_nit_n=UP_nit_n(:);
            UP_nit_all(:,i) = UP_nit_n;
        end
        UP_amm = sum(UP_amm_all,2);
        UP_nit = sum(UP_nit_all,2);
        
        % Compute total Nitrogen Uptake.
        UP_N = sum(UP_amm) + sum(UP_nit);   %[gr/m^2/d]
        
        % Compute uptake in units of [gr / m^3 / d]
        UP_nit = UP_nit./(dz_mm/1000);   % [gr/m^3/d]
        UP_amm = UP_amm./(dz_mm/1000);   % [gr/m^3/d]

         
%*************************************************************************
%       CHANGES IN MINERAL CONCENTRATIONS
%*************************************************************************
    
        % Nitrification Rate
            Nitrif = (fSn.*fTn*kn.*Cb) .* Amm;
            
        % Rate of CO2 Production (Eqn. 23)
%            Prod_CO2 = (1-fe)*(1-fh)*kl*fTd*fSd*Cl + (1-fe)*kh*fTd*fSd*Ch;

        % Denitrification Rate
%            Denitrif = min(alph*fTdn*fSdn*Prod_CO2, bet*NO3);
            Denitrif = zeros(nl_soil,1);

        % Volatilization Rate
            %Volat = kv*fTv*Amm;
            Volat = zeros(nl_soil,1);
            
            
        if min(size(MIN)==[nl_soil,1])~=1    
            st=1;
        end    
        if min(size(IMM_Amm)==[nl_soil,1])~=1    
            st=1;
        end    
        if min(size(Nitrif)==[nl_soil,1])~=1    
            st=1;
        end    
        if min(size(Volat)==[nl_soil,1])~=1    
            st=1;
        end    
        if min(size(UP_amm)==[nl_soil,1])~=1    
            st=1;
        end    
        dAmm_dt = MIN - IMM_Amm - Nitrif - Volat - UP_amm;  %[gr Amm / m^3 / d]
        dNit_dt = Nitrif - IMM_Nit - Denitrif - LCH_nit - UP_nit;
            
%*************************************************************************
%      STATE UPDATES
%*************************************************************************


    Cl_new = Cl + dCl_dt*dtime/86400;
    Ch_new = Ch + dCh_dt*dtime/86400;
    Cb_new = Cb + dCb_dt*dtime/86400;
    Nl_new = Nl + dNl_dt*dtime/86400;
    Amm_new = Amm + dAmm_dt*dtime/86400;
    Nit_new = Nit + dNit_dt*dtime/86400;
    CNl_new = Cl_new./Nl_new;
