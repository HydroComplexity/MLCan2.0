function [UP_amm, UP_nit, UP_amm_all, UP_nit_all, F_HR_N] = CN_nuptake (VARIABLES, PARAMS, SWITCHES, a_Nit,...
         a_Amm, sm, layeruptake_all, TR_can_store)        

        if SWITCHES.CN.Bioturbation;
            nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
        else
            nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
        end
        nspecies = PARAMS.CanStruc.nspecies; % number of species
        Nit = VARIABLES.Nit;
        Amm = VARIABLES.Amm;        
        UP_amm_n=zeros(nl_soil,1);
        UP_nit_n=zeros(nl_soil,1);
        UP_amm_all=zeros(nl_soil,nspecies);
        UP_nit_all=zeros(nl_soil,nspecies);
        factorrec = PARAMS.CN.factorrec; % fraction of Nitrate that is recycled
        Recycling = SWITCHES.Recycling; % Switch that decides if recycling of nitrate is on or off
        
        
        % PLANT UPTAKE
        % AMMONIUM
        for i=1:1:nspecies
            layeruptn=layeruptake_all(:,i);
            TR = TR_can_store(i);
            wtind=layeruptn>0;
            wtnind=layeruptn<=0;
            wtzind=layeruptn==0;
%            F_HR_N = TR./sum(layeruptn(wtind));
            F_HR_N = 1;
            % First for the case in which the flow is from the soil to the root
            UP_amm_n(wtind) = (a_Amm./sm(wtind).*Amm(wtind)).*layeruptn(wtind)*F_HR_N;% [gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system. 
            if Recycling
                waterrec = sum(UP_amm_n(wtind))*factorrec;
                UP_amm_n(wtnind) = - waterrec.*layeruptn(wtnind)/(sum(layeruptn(wtnind)));% [gr/m^2/d]           
            else
                UP_amm_n(wtnind) = 0;
            end
            UP_amm_n(wtzind) = 0;            
             
            UP_amm_n=UP_amm_n(:);
            UP_amm_all(:,i) = UP_amm_n;

            % NITRATE
            % First for the case in which the flow is from the soil to the root        
            UP_nit_n(wtind) = (a_Nit./sm(wtind).*Nit(wtind)).*layeruptn(wtind)*F_HR_N;%[gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system. 
            if Recycling
                waterrec = sum(UP_nit_n(wtind))*factorrec;
                UP_nit_n(wtnind) = - waterrec.*layeruptn(wtnind)/(sum(layeruptn(wtnind)));% [gr/m^2/d]           
            else
                UP_nit_n(wtnind) = 0;
            end
            UP_nit_n(wtzind) = 0;            
 
            UP_nit_n=UP_nit_n(:);
            UP_nit_all(:,i) = UP_nit_n;
            clear UP_nit_n;
            clear UP_amm_n;            
        end
        UP_amm = sum(UP_amm_all,2);         % [gr/m^2/d]
        UP_nit = sum(UP_nit_all,2);         % [gr/m^2/d]
        