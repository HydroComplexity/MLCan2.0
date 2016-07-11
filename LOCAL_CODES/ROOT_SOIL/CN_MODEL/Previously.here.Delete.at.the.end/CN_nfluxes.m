function [Nitrif, Denitrif, Volat] = CN_nfluxes (VARIABLES, PARAMS, VERTSTRUC, SWITCHES, fSn, fTn)
    
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


        % Nitrification Rate
            Nitrif = (fSn.*fTn*kn.*Cb).* Amm;   %[gr N/m3/d]
            
            indnit = Nitrif>Amm;
            Nitrif(indnit) = Amm(indnit);
        % Rate of CO2 Production (Eqn. 23)
%            Prod_CO2 = (1-fe)*(1-fh)*kl*fTd*fSd*Cl + (1-fe)*kh*fTd*fSd*Ch;

        % Denitrification Rate
%            Denitrif = min(alph*fTdn*fSdn*Prod_CO2, bet*NO3);
            Denitrif = zeros(nl_soil,1);

        % Volatilization Rate
            %Volat = kv*fTv*Amm;
            Volat = zeros(nl_soil,1);


