function    [VARIABLES, SUP_amm_pas_all_m2, SUP_nit_pas_all_m2, SUP_amm_act_all_m2, SUP_nit_act_all_m2, ...
    FUP_amm_all_m2, FUP_nit_all_m2, RUP_amm_all_m2, RUP_nit_all_m2] = ...
    CN_nuptake (VARIABLES, PARAMS, SWITCHES, VERTSTRUC, FORCING, STORAGE, CONSTANTS, a_Nit, a_Amm, sm, layeruptake_all, rootfr)

%%%%% LOAD PARAMETERS %%%%%
%nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers
if SWITCHES.CN.Bioturbation;   % include litter
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
    dz_mm = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm];%       dz_mm = grid depth [mm]
else
    nl_soil = PARAMS.nl_soil;%       nl_soil = # soil layers        
    dz_mm = VERTSTRUC.dzsmm;%       dz_mm = grid depth [mm]    
end

nspecies = PARAMS.CanStruc.nspecies; % number of species
%n=PARAMS.porosity;
n=VERTSTRUC.eff_poros;
n=[n(1);n];
%F=PARAMS.F;
%d=PARAMS.d;
F=PARAMS.CN.fRd;
d=PARAMS.CN.fTortuo;
% Nremobilization_rate=PARAMS.Nremobilization_rate;
Nremobilization_rate=PARAMS.CN.NrI;

%nfix_rate=PARAMS.nfix_rate;
%NFr=nfix_rate(day);
NFr=PARAMS.CN.NfI;

fDEMamm=PARAMS.CN.fDEMamm;
fDEMnit=PARAMS.CN.fDEMnit;

dtime=CONSTANTS.dtime;
DOY_start=PARAMS.CN.DOY_start;
DOY_end=PARAMS.CN.DOY_end;
%%%%%% Switchs %%%%%
%factorrec = PARAMS.CN.factorrec; % fraction of Nitrate that is recycled
Recycling = SWITCHES.Recycling; % Switch that decides if recycling of nitrate is on or off
%dz_in_m=VERTSTRUC.dzsmm./1000;

Active_Nuptake=SWITCHES.Active_Nuptake;

% Nremobil=SWITCHES.Nremobil;
Nremobil=SWITCHES.CN.N_Remo;

%%%%% LOAD VARIABLES %%%%%
Nit = VARIABLES.Nit;
Amm = VARIABLES.Amm;
%Nit=Nit_Update; 
%Amm=Amm_Update;
timestep = VARIABLES.timestep;          % timestep = Current time step

%%%%% LOAD Forcings %%%%%
dtime = CONSTANTS.dtime;
DEM=FORCING.CNdem(:,timestep).*86400/dtime; % Making /d

%%%%% Allocation %%%%%
UP_amm_n=zeros(nl_soil,1); %0
UP_nit_n=zeros(nl_soil,1); %0
UP_amm_all=zeros(nl_soil,nspecies); %0
UP_nit_all=zeros(nl_soil,nspecies); %0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixation swiches is not working. If one do not want to consider the
% fixation, just put "0" for fixation rate
if SWITCHES.CN.N_Fix == 1;
    NFr=NFr;
elseif SWITCHES.CN.N_Fix == 0;
    NFr=zeros(size(NFr));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*************************************************************************
%       Nitrogen Uptake from soil (Passive)                              *
%*************************************************************************

for i=1:1:nspecies
    % AMMONIUM
    layeruptn=layeruptake_all(:,i); % plant water uptake
    wtind=layeruptn>0;
    wtnind=layeruptn<0;
    wtzind=layeruptn==0;
    % First for the case in which the flow is from the soil to the root
    UP_amm_n(wtind) = (a_Amm./sm(wtind).*Amm(wtind)).*layeruptn(wtind);% [gr/m^2/d];
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
    UP_nit_n(wtind) = (a_Nit./sm(wtind).*Nit(wtind)).*layeruptn(wtind);%[gr/m^2/d];
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
    
    UP_nit_active=nan;
    UP_amm_active=nan;
%    Amm_DEM=zeros(nl_soil,nspecies);
%    Nit_DEM=zeros(nl_soil,nspecies);
end
SUP_amm = sum(UP_amm_all,2);         % [gr/m^2/d]
SUP_nit = sum(UP_nit_all,2);         % [gr/m^2/d]
% Organize the output
SUP_amm_pas_all_m2=UP_amm_all;
SUP_nit_pas_all_m2=UP_nit_all;

%*************************************************************************
%       Nitrogen Uptake from soil (Active)                               *
%*************************************************************************
if Active_Nuptake
%     % Nit, and Amm in the unit of [g/m3]
%     % carbon content's units from PALMS [kg/m2]
%     grain_carbon=CARBONS.grain_carbon_day_point(day)*1000;     % [g/m2]
%     leaf_carbon=CARBONS.leaf_carbon_day_point(day)*1000;       % [g/m2]
%     rhizome_carbon=CARBONS.rhizome_carbon_day_point(day)*1000; % [g/m2]
%     root_carbon=CARBONS.root_carbon_day_point(day)*1000;       % [g/m2]
%     stem_carbon=CARBONS.stem_carbon_day_point(day)*1000;       % [g/m2]
%     
%     %CNratio=ALL_CNveg(day);
%     %CNratio_root=ALL_CNveg_root(day);
%     CNleaf=CNratio_leaf(day);
%     CNstem=CNratio_stem(day);
%     CNgrain=CNratio_grain(day);
%     CNroot=CNratio_root(day);
%     
%     if day>1
%         pre_grain_carbon=CARBONS.grain_carbon_day_point(day-1)*1000;
%         pre_leaf_carbon=CARBONS.leaf_carbon_day_point(day-1)*1000;
%         pre_rhizome_carbon=CARBONS.rhizome_carbon_day_point(day-1)*1000;
%         pre_root_carbon=CARBONS.root_carbon_day_point(day-1)*1000;
%         pre_stem_carbon=CARBONS.stem_carbon_day_point(day-1)*1000;
%         
%         %pre_CNratio=ALL_CNveg(day-1);
%         %pre_CNratio_root=ALL_CNveg_root(day-1);
%         pre_CNleaf=CNratio_leaf(day-1);
%         pre_CNstem=CNratio_stem(day-1);
%         pre_CNgrain=CNratio_grain(day-1);
%         pre_CNroot=CNratio_root(day-1);
%     else
%         pre_grain_carbon=CARBONS.grain_carbon_day_point(day)*1000;
%         pre_leaf_carbon=CARBONS.leaf_carbon_day_point(day)*1000;
%         pre_rhizome_carbon=CARBONS.rhizome_carbon_day_point(day)*1000;
%         pre_root_carbon=CARBONS.root_carbon_day_point(day)*1000;
%         pre_stem_carbon=CARBONS.stem_carbon_day_point(day)*1000;
%         
%         %pre_CNratio=ALL_CNveg(day);
%         %pre_CNratio_root=ALL_CNveg_root(day);
%         pre_CNleaf=CNratio_leaf(day);
%         pre_CNstem=CNratio_stem(day);
%         pre_CNgrain=CNratio_grain(day);
%         pre_CNroot=CNratio_root(day);
%     end
%     
%     F=PARAMS.F;
%     d=PARAMS.d;
%     n=PARAMS.porosity;
%     
    % Calculate ku: dependence of the diffusion process on
    % the soil moisture level [-/d]
    for i=1:1:nspecies
        Nit_ku(:,i)=(a_Nit./(sm)).*F(i).*((sm./n).^d(i)); % [m/d]
        Amm_ku(:,i)=(a_Amm./(sm)).*F(i).*((sm./n).^d(i)); % [m/d]
    end
%     
%     % Calculate DEM based on carbon content change.
%     %Above_Current_total_C=grain_carbon+leaf_carbon+stem_carbon; % [g/m2]
%     Below_Current_total_C=rhizome_carbon+root_carbon; % [g/m2]
%     
%     %Above_Previous_total_C=pre_grain_carbon+pre_leaf_carbon+pre_stem_carbon; % [g/m2]
%     Below_Previous_total_C=pre_rhizome_carbon+pre_root_carbon; % [g/m2]
%     
%     %Above_C_change=Above_Current_total_C-Above_Previous_total_C;
%     leaf_C_change=leaf_carbon-pre_leaf_carbon;
%     stem_C_change=stem_carbon-pre_stem_carbon;
%     grain_C_change=grain_carbon-pre_grain_carbon;
%     Below_C_change=Below_Current_total_C-Below_Previous_total_C;
%     
%     % Calibaration
%     %if Above_C_change > 0
%     %    Above_C_change=(Above_Current_total_C-Above_Previous_total_C)*Add_Carbon; % [g/m2]
%     %end
%     %if Below_C_change > 0
%     %    Below_C_change=(Below_Current_total_C-Below_Previous_total_C)*Add_Carbon; % [g/m2]
%     %end
%     
%     % in order to have 0 carbon change for every first day of year
%     for i=1:how_many_year
%         if day == detect_first_day_of_year(i)
%             %Above_C_change=0;
%             leaf_C_change=0;
%             stem_C_change=0;
%             grain_C_change=0;
%             Below_C_change=0;
%         end
%     end
%     
%     % average CN ratio
%     %ave_CNratio=(pre_CNratio+CNratio)/2;
%     ave_leaf_CNratio=(pre_CNleaf+CNleaf)/2;
%     ave_stem_CNratio=(pre_CNstem+CNstem)/2;
%     ave_grain_CNratio=(pre_CNgrain+CNgrain)/2;
%     ave_CNratio_root=(pre_CNroot+CNroot)/2;
%     
%     % Calculate Demand of N
%     %if Above_C_change > 0
%     %    Above_DEM=(Above_C_change/ave_CNratio).*rootfr; % [g/m2]
%     %else
%     %    Above_DEM=zeros(nl_soil,1);
%     %end
%     if leaf_C_change > 0
%         leaf_DEM=(leaf_C_change/ave_leaf_CNratio).*rootfr; % [g/m2]
%     else
%         leaf_DEM=zeros(nl_soil,1);
%     end
%     if stem_C_change > 0
%         stem_DEM=(stem_C_change/(ave_stem_CNratio)).*rootfr; % [g/m2]
%     else
%         stem_DEM=zeros(nl_soil,1);
%     end
%     if grain_C_change > 0
%         grain_DEM=(grain_C_change/ave_grain_CNratio).*rootfr; % [g/m2]
%     else
%         grain_DEM=zeros(nl_soil,1);
%     end
%     
%     if Below_C_change > 0
%         Below_DEM=(Below_C_change/ave_CNratio_root).*rootfr; % [g/m2]
%     else
%         Below_DEM=zeros(nl_soil,1);
%     end
%     %DEM=Above_DEM+Below_DEM;
%     DEM=leaf_DEM+stem_DEM+grain_DEM+Below_DEM;
%     
%     % For GUI
%     VARIABLES.DEM=sum(DEM); % [g/m2]
%     
%     % DEM fraction
%     if sum(DEM) == 0
%         DEM_fraction_l_s_g_b = [0 0 0 0];
%     else
%         DEM_fraction_l_s_g_b=[sum(leaf_DEM) sum(stem_DEM) sum(grain_DEM) sum(Below_DEM)]./sum(DEM);
%     end
    
    % Due to low DEM
    %             if choose_crop ==4
    %                 if day < 365*3
    %                     DEM=DEM*3;
    %                 end
    %                 if day < 365*2
    %                     DEM=DEM/3*10;
    %                 end
    %                 if day < 365
    %                     DEM=DEM/3/10*30;
    %                 end
    %             end
    
    % Corn and soy
%     if choose_crop ==2
%         if iscorn == 1 % Corn
%             Nit_DEM=DEM.*0.00;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
%             Amm_DEM=DEM.*1.00;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
%         else % Soy
%             Nit_DEM=DEM.*0.00;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
%             Amm_DEM=DEM.*1.00;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
%         end
%     end
%     % SG: 0.85, 0.15
%     if choose_crop ==3
%         Nit_DEM=DEM.*0.9;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
%         Amm_DEM=DEM.*0.1;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
%     end
%     % MG
%     if choose_crop ==4
%         Nit_DEM=DEM.*0.45;%0.45; % The values (0.2 and 0.5) came from D'Odorico et al., 2003
%         Amm_DEM=DEM.*0.55;%0.55; % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
%     end
    for i=1:1:nspecies
        Amm_DEM(:,i)=DEM(i).*fDEMamm(i).*rootfr(:,i);
        Nit_DEM(:,i)=DEM(i).*fDEMnit(i).*rootfr(:,i);
    end

%     % for NP, I choose Nit_DEM and Amm_DEM based on the paper by Porporato et al. (2003).
%     if choose_crop ==1
%         %if day>1
%         %    pre_totlail=totlail_layer_day_point(day-1);
%         %else
%         %    pre_totlail=totlail_layer_day_point(day);
%         %end
%         %totlail=totlail_layer_day_point(day);
%         %lai_change=totlail-pre_totlail;
%         
%         % 0.8m -> 0.5, 0.2 when rootfr = 0.9760
%         % 0.5*(1+1-0.9760)=0.512
%         % 0.2*(1+1-0.9760)=0.2048
%         
%         %if lai_change > 0
%         Nit_DEM=(rootfr.*0.075).*dz_in_m; % 0.5 [g/m3/d] * soil depth or 0.06
%         Amm_DEM=(rootfr.*0.000).*dz_in_m; % 0.2 [g/m3/d] * soil depth or 0.01
%         %else
%         %    Nit_DEM=zeros(nl_soil,1);
%         %    Amm_DEM=zeros(nl_soil,1);
%         %end
%     end
    
    % Porporato et al. (2003) 's three cases in eqn. 31.
    % UP_amm, UP_nit [gr/m^2/d]
    % Nit, and Amm [g/m3]
    % DEM [g/m2]
    % Nit_ku, Amm_ku [m/d]
    
    
    [UP_nit_active, Fix_Nit] = CN_passive_fixation_nuptake_Nit (nl_soil, nspecies, Nit_DEM, SUP_nit_pas_all_m2, Nit_ku, Nit, NFr);
    [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, nspecies, Amm_DEM, SUP_amm_pas_all_m2, Amm_ku, Amm, NFr);
    
    %SUP_nit_active=UP_nit_active(:); % [g/m2/d]
    %SUP_amm_active=UP_amm_active(:); % [g/m2/d]
    SUP_nit_active=UP_nit_active;
    SUP_amm_active=UP_amm_active;
    
    %SUP_nit_passive=SUP_nit(:); % [g/m2/d]
    %SUP_amm_passive=SUP_amm(:); % [g/m2/d]
    SUP_nit_passive=SUP_nit_pas_all_m2; % [g/m2/d]
    SUP_amm_passive=SUP_amm_pas_all_m2; % [g/m2/d]
    
%     if size(SUP_nit_passive) ~= size(SUP_nit_passive)
%         stop=1;
%         
%     end
%     if size(SUP_amm_passive) ~= size(SUP_amm_active)
%         stop=1; 
%     end
%     if sum(isnan(SUP_amm_passive))
%        stop=1; 
%     end
%     if sum(isnan(SUP_amm_active))
%        stop=1; 
%     end
%     if size(SUP_amm_passive,1) ~=13
%         stop=1;
%     end
%     if size(SUP_amm_active,1) ~=13
%         stop=1;
%     end
%     if sum(isinf(SUP_amm_passive))
%        stop=1;
%     end
%     if sum(isinf(SUP_amm_active))
%        stop=1;
%     end
    
    % N Uptake from soil = passive upatke + active uptake
    SUP_nit=sum(SUP_nit_passive+SUP_nit_passive,2); % [g/m2/d]
    SUP_amm=sum(SUP_amm_passive+SUP_amm_active,2); % [g/m2/d]
    
    %SUP_nit_all=SUP_nit; % [g/m2/d]
    %SUP_amm_all=SUP_amm; % [g/m2/d]
    SUP_nit_all=SUP_nit_passive+SUP_nit_active; % [g/m2/d]
    SUP_amm_all=SUP_amm_passive+SUP_amm_active; % [g/m2/d]
    
    % N Uptake from N fixation
    %Fix_Nit=Fix_Nit(:);
    %Fix_Amm=Fix_Amm(:);
    Fix_Nit;
    Fix_Amm;
end
% Organize the output
SUP_nit_act_all_m2=SUP_nit_active;
SUP_amm_act_all_m2=SUP_amm_active;
FUP_nit_all_m2=Fix_Nit;
FUP_amm_all_m2=Fix_Amm;



RUP_amm_all_m2=zeros(size(FUP_amm_all_m2));
RUP_nit_all_m2=zeros(size(FUP_amm_all_m2));
UPr_amm=zeros(nl_soil,nspecies);
UPr_nit=zeros(nl_soil,nspecies);

if Nremobil
    if DOY_start == 1 && (DOY_end == 365 || DOY_end == 366)
%     for i=1:how_many_year
%         if day == detect_first_day_of_year(i)
%             if Check_year == 1
%                 TUPr_yrAmm_pool=0;
%                 TUPr_yrNit_pool=0;
%                 
%                 VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
%                 VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
%             else
%                 DEM_fraction_Abov_p_below=(1-DEM_fraction_l_s_g_b_rate(:,4))';
%                 TUPr_yrAmm_pool=sum(sum(Tot_UP_Amm_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                     (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                     (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1))); % g/m2
%                 TUPr_yrNit_pool=sum(sum(Tot_UP_Nit_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                     (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                     (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)));
%                 
%                 VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
%                 VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
%             end
%         end
%     end

        % This only works for 100 years
        HOW_MANY_YEAR_DO_YOU_WANT=100;
        timestep_inayear=365*24*60*60/dtime;
        detect_first_day_of_year=[1,(1:HOW_MANY_YEAR_DO_YOU_WANT).*timestep_inayear+1];
                
        for i=1:HOW_MANY_YEAR_DO_YOU_WANT
            if timestep == detect_first_day_of_year(i)
                if timestep == 1
                    TUPr_yrAmm_pool=zeros(1,nspecies);
                    TUPr_yrNit_pool=zeros(1,nspecies);
                    
                    VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
                    VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
                else
                    UP_nit_m2_store=STORAGE.UP_nit_m2_store;
                    UP_amm_m2_store=STORAGE.UP_amm_m2_store;
%                     DEM_fraction_Abov_p_below=(1-DEM_fraction_l_s_g_b_rate(:,4))';
%                     TUPr_yrAmm_pool=sum(sum(Tot_UP_Amm_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                         (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                         (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1))); % g/m2
%                     TUPr_yrNit_pool=sum(sum(Tot_UP_Nit_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                         (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
%                         (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)));
                    TUPr_yrAmm_pool=sum((sum(UP_amm_m2_store(:,:,timestep-timestep_inayear:timestep-1))),3); % g/m2
                    TUPr_yrNit_pool=sum((sum(UP_nit_m2_store(:,:,timestep-timestep_inayear:timestep-1))),3);
                    VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
                    VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
                end
            end
        end

    TUPr_yrAmm_pool=VARIABLES.TUPr_yrAmm_pool;
    TUPr_yrNit_pool=VARIABLES.TUPr_yrNit_pool;
    UPr_amm=zeros(nl_soil,nspecies);
    UPr_nit=zeros(nl_soil,nspecies);
    
%     % Ammonium
%     if TUPr_yrAmm_pool > 0
%         if TUPr_yrAmm_pool >= sum(Amm_DEM)
%             UPr_amm=Amm_DEM;
%             TUPr_yrAmm_pool=TUPr_yrAmm_pool-sum(Amm_DEM); % g/m2
%             
%             % Update active, passive, fixation uptakes
%             SUP_amm_active=zeros(nl_soil,1); % [g/m2/d]
%             SUP_amm_passive=zeros(nl_soil,1); % [g/m2/d]
%             
%             % N Uptake from soil = passive upatke + active uptake
%             SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
%             SUP_amm_all=SUP_amm; % [g/m2/d]
%             
%             % N Uptake from N fixation
%             Fix_Amm=zeros(nl_soil,1);
%         else
%             UPr_amm=TUPr_yrAmm_pool.*rootfr;
%             TUPr_yrAmm_pool=0; % g/m2
%             
%             % Update active, passive, fixation uptakes
%             Amm_DEM_updated=(sum(Amm_DEM)- sum(UPr_amm)).*rootfr;
%             [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, Amm_DEM_updated, SUP_amm, Amm_ku, Amm, NFr);
%             
%             SUP_amm_active=UP_amm_active(:); % [g/m2/d]
%             SUP_amm_passive=SUP_amm_passive; % [g/m2/d]
%             
%             % N Uptake from soil = passive upatke + active uptake
%             SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
%             SUP_amm_all=SUP_amm; % [g/m2/d]
%             
%             % N Uptake from N fixation
%             Fix_Amm=Fix_Amm(:);
%         end
%     end
%     VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
 % Ammonium   
    
        for j=1:nspecies
            if TUPr_yrAmm_pool(j) > 0
                if TUPr_yrAmm_pool(j) >= sum(Amm_DEM(:,j))
                    UPr_amm(:,j)=Amm_DEM(:,j);
                    TUPr_yrAmm_pool(j)=TUPr_yrAmm_pool(j)-sum(Amm_DEM(:,j)); % g/m2

                    % Update active, passive, fixation uptakes
                    SUP_amm_active(:,j)=zeros(nl_soil,1); % [g/m2/d]
                    SUP_amm_passive(:,j)=zeros(nl_soil,1); % [g/m2/d]

                    % N Uptake from soil = passive upatke + active uptake
                    %SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
                    %SUP_amm_all=SUP_amm; % [g/m2/d]
                    SUP_amm_all(:,j)=(SUP_amm_passive(:,j)+SUP_amm_active(:,j)); % [g/m2/d]

                    % N Uptake from N fixation
                    Fix_Amm(:,j)=zeros(nl_soil,1);
                else
                    UPr_amm(:,j)=TUPr_yrAmm_pool(j).*rootfr(:,j);
                    TUPr_yrAmm_pool(j)=0; % g/m2

                    % Update active, passive, fixation uptakes
                    Amm_DEM_updated(:,j)=(sum(Amm_DEM(:,j))- sum(UPr_amm(:,j))).*rootfr(:,j);
                    [UP_amm_active(:,j), Fix_Amm(:,j)] = CN_passive_fixation_nuptake_Amm (nl_soil, 1, Amm_DEM_updated(:,j), SUP_amm_pas_all_m2(:,j), Amm_ku(:,j), Amm, NFr(j));

                    %SUP_amm_active=UP_amm_active(:); % [g/m2/d]
                    %SUP_amm_passive=SUP_amm_passive; % [g/m2/d]
                    SUP_amm_active(:,j)=UP_amm_active(:,j); % [g/m2/d]
                    SUP_amm_passive(:,j)=SUP_amm_passive(:,j); % [g/m2/d]

                    % N Uptake from soil = passive upatke + active uptake
                    %SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
                    %SUP_amm_all=SUP_amm; % [g/m2/d]
                    SUP_amm_all(:,j)=(SUP_amm_passive(:,j)+SUP_amm_active(:,j)); % [g/m2/d]

                    % N Uptake from N fixation
                    %Fix_Amm=Fix_Amm(:);
                    Fix_Amm(:,j)=Fix_Amm(:,j);
                end
            end

            VARIABLES.TUPr_yrAmm_pool(j)=TUPr_yrAmm_pool(j);

            % Nitrate
            if TUPr_yrNit_pool(j) > 0
                if TUPr_yrNit_pool(j) >= sum(Nit_DEM(:,j))
                    UPr_nit(:,j)=Nit_DEM(:,j);
                    TUPr_yrNit_pool(j)=TUPr_yrNit_pool(j)-sum(Nit_DEM(:,j)); % g/m2

                    % Update active, passive, fixation uptakes
                    SUP_nit_active(:,j)=zeros(nl_soil,1); % [g/m2/d]
                    SUP_nit_passive(:,j)=zeros(nl_soil,1); % [g/m2/d]

                    % N Uptake from soil = passive upatke + active uptake
                    %SUP_nit=(SUP_nit_passive+SUP_nit_active); % [g/m2/d]
                    %SUP_nit_all=SUP_nit; % [g/m2/d]
                    SUP_nit_all=(SUP_nit_passive(:,j)+SUP_nit_active(:,j)); % [g/m2/d]

                    % N Uptake from N fixation
                    Fix_Nit(:,j)=zeros(nl_soil,1);
                else
                    UPr_nit(:,j)=TUPr_yrNit_pool(j).*rootfr(:,j);
                    TUPr_yrNit_pool(j)=0; % g/m2

                    % Update active, passive, fixation uptakes
                    Nit_DEM_updated(:,j)=(sum(Nit_DEM(:,j))- sum(UPr_nit(:,j))).*rootfr(:,j);
                    [UP_nit_active(:,j), Fix_Nit(:,j)] = CN_passive_fixation_nuptake_Nit (nl_soil, 1, Nit_DEM_updated(:,j), SUP_nit_pas_all_m2(:,j), Nit_ku(:,j), Nit, NFr(j));

                    %SUP_nit_active=UP_nit_active(:); % [g/m2/d]
                    %SUP_nit_passive=SUP_nit_passive; % [g/m2/d]
                    SUP_nit_active(:,j)=UP_nit_active(:,j); % [g/m2/d]
                    SUP_nit_passive(:,j)=SUP_nit_passive(:,j); % [g/m2/d]

                    % N Uptake from soil = passive upatke + active uptake
                    %SUP_nit=(SUP_nit_passive+SUP_nit_active); % [g/m2/d]
                    %SUP_nit_all=SUP_nit; % [g/m2/d]
                    SUP_nit_all(:,j)=(SUP_nit_passive(:,j)+SUP_nit_active(:,j)); % [g/m2/d]

                    % N Uptake from N fixation
                    %Fix_Nit=Fix_Nit(:);
                    Fix_Nit=Fix_Nit;
                end
            end
            VARIABLES.TUPr_yrNit_pool(j)=TUPr_yrNit_pool(j);
        end
    
    end
end
% Organize the output
SUP_nit_pas_all_m2=SUP_nit_passive;
SUP_amm_pas_all_m2=SUP_amm_passive;
SUP_nit_act_all_m2=SUP_nit_active;
SUP_amm_act_all_m2=SUP_amm_active;
FUP_nit_all_m2=Fix_Nit;
FUP_amm_all_m2=Fix_Amm;
RUP_nit_all_m2=UPr_nit;
RUP_amm_all_m2=UPr_amm;


