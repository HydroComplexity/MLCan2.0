function [Nitrif, Denitrif, Volat, Total_N20, N2O_nit] = CN_nfluxes (VARIABLES, PARAMS, VERTSTRUC, SWITCHES, fSn, fTn, DECl, DECh)
    
        Cb = VARIABLES.Cb;%       Cb = carbon concentration in biomass pool [gC / m^3]
        Ch = VARIABLES.Ch;%       Cb = carbon concentration in biomass pool [gC / m^3]
        Amm = VARIABLES.Amm;%       Cb = carbon concentration in biomass pool [gC / m^3]
        if SWITCHES.CN.Bioturbation;
            nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
            dz_mm = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
        else
            nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
        end
        
        kd = PARAMS.CN.kd;%       kd = slope of linear dependence of biomass death on biomass concentration [d^-1]
        kl = PARAMS.CN.kl;%       kl = rate of decomposition of the litter pool [m^3 d / gC]
        kh = PARAMS.CN.kh;%       kh = rate of decomposition of the humus pool [m^3 d / gC]
        kn = PARAMS.CN.kn;%       kh = rate of nitrification of the humus pool [m^3 d / gC]

        % Dongkook Woo - Edit
        N_denit=SWITCHES.CN.Denitrification;
        Nit = VARIABLES.Nit;
        rr=PARAMS.CN.rr;
        K2=PARAMS.CN.K2;
        
        porosity=[VERTSTRUC.porsl(1);VERTSTRUC.porsl]; 
        real_sm = [VARIABLES.SOIL.volliqli ; VARIABLES.SOIL.volliq];
        BulkDensity=VERTSTRUC.rhod; % [kg/m3]
        BulkDensity=[BulkDensity(1).*1000;BulkDensity.*1000]; % [g/m3]
        
        % Prepare to compute denitrification
        WFPS=(real_sm./porosity).*100; % Water-Filled Pore Space
        indOver100=WFPS>100;
        WFPS(indOver100)=100;
        WFPS=WFPS(:);
        
        if SWITCHES.CN_type == 1
            PCO2=(DECl.*rr+DECh.*rr);
        elseif SWITCHES.CN_type == 0
            PCO2=(DECl.*rr);
        end
        
        ClayContent=[VERTSTRUC.clayColum(1);VERTSTRUC.clayColum]./100; % [%]
        SoilAir=porosity-real_sm;
        indSoilAir = SoilAir < 0;
        SoilAir(indSoilAir)=0;
        Dfc=porosity.^2.*(SoilAir./porosity).^(2+3./(13.6.*ClayContent+3.5));
        Dfc=Dfc(:);
        % Dongkook Woo - Edit End
        
        % Nitrification Rate
            Nitrif = (fSn.*fTn*kn.*Cb).* Amm;   %[gr N/m3/d]
            
            indnit = Nitrif>Amm;
            Nitrif(indnit) = Amm(indnit);
        % Rate of CO2 Production (Eqn. 23)
%            Prod_CO2 = (1-fe)*(1-fh)*kl*fTd*fSd*Cl + (1-fe)*kh*fTd*fSd*Ch;

        % Volatilization Rate
            %Volat = kv*fTv*Amm;
            Volat = zeros(nl_soil,1);

        % Dongkook Woo - Edit            
        % Denitrification Rate
%            Denitrif = min(alph*fTdn*fSdn*Prod_CO2, bet*NO3);
        if N_denit == 0
            Denitrif = zeros(nl_soil,1);
            Total_N20 = zeros(nl_soil,1);
            N2O_nit = zeros(nl_soil,1);
        elseif N_denit == 1

        % N2O flux from nitrification
            N2O_nit=Nitrif.*K2; % gN/m3/d
            
        % Nit [g/m3] ./ Bulkdensity [g/m3] .* g to ug [1g = 1000000ug]
            FdNO3=1.15.*(((Nit./BulkDensity*(1000000/1))).^0.57); % [ugN/g/d]
        % PCO2 [g/m3] ./ Bulkdensity [g/m3] .* g to ug [1g = 1000000ug]
            FdCO2=0.1.*(((PCO2./BulkDensity*(1000000/1))).^1.3); % [ugN/g/d]
            
                %Mmin=ones(size(dz,1),1).*0.113;
                Mmin=ones(nl_soil,1).*0.113;
                
            M=min(Mmin,Dfc).*(-3.05)+0.36;
            a=0.90-M.*(PCO2./BulkDensity*(1000000/1));
            % Paper mistake correction based on new one
            %%% FdWFPS=0.5+(atan(0.6.*pi.*0.1.*WFPS-a))./(pi);
            %FdWFPS=0.5.*(atan(0.6.*pi.*0.1.*WFPS-a))./(pi);
            % Due to FdWFPS's error from the paper! they arbitrariy change
            % this eqn. why not me?
            FdWFPS=(0.5).*(atan(0.6.*pi.*(0.1.*WFPS-a)))./(pi);
            isneg=FdWFPS<0;
            FdWFPS(isneg)=0;
%            FdWFPS=(0.5.*PARAMS.dnitrification_factor)+(atan(0.6.*pi.*(0.1.*WFPS-a)))./(pi);
            
            Dt=min(FdNO3,FdCO2).*FdWFPS; % [ugN/g/d]
            
                %kllmax=ones(size(dz,1),1).*1.7;
                kllmax=ones(nl_soil,1).*1.7;
            kll=max(kllmax,38.4-350.*Dfc);
            FrNO3CO2=max(0.16.*kll,kll.*exp(-0.8.*(Nit./PCO2)));
            FrWFPS=max(0.1,0.015.*(WFPS)-0.32);
            % Test
            %FrWFPS=max(0.1,0.5.*(WFPS)-0.32);
            RN2N2O=FrNO3CO2.*FrWFPS;
            
            DN2O=Dt./(1+RN2N2O); % [ugN/g/d]
           
            
        % Dt [ugN/g/d] *  Bulkdensity [g/m3]  * ug to g [1g = 1000000ug]
            Denitrif=Dt.*BulkDensity.*(1/1000000); % gN/m3/d
            N20_denitrif=DN2O.*BulkDensity.*(1/1000000); % gN/m3/d
            Total_N20=N2O_nit+N20_denitrif; % gN/m3/d
        end
        % Dongkook Woo - Edit End