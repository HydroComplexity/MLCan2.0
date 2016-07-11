% Dongkook Woo - Comment
% function []= core_N(VARIN)
% % Generate char for root cut type
% % Decode info for root cutting
% if VARIN(4) == 1
%    strcut = 'NN';
% elseif VARIN(4) == 2
%    strcut = 'OH';
% elseif VARIN(4) == 3
%    strcut = 'AH1';
% elseif VARIN(4) == 4
%    strcut = 'AH2';
% end
% 
% 
% % Generate VAR
% infile=[num2str(VARIN(1))  '_CUT' strcut '_HR' num2str(VARIN(2)) '_nsp' num2str(VARIN(3))];
% outfile1 = ['CN_' infile];
% outfile2 = ['CNSS_' infile];
% 
% %clc;
% filename = ['./MLCAN_SA/' infile '.mat'];
% % Code Library Paths
%     addpath('./CN_MODEL/');
% 
% if isempty(dir(filename))
%      disp(sprintf('file missing: %s\n',filename));
%      return;
% else    
%     load (filename,'volliq_store','Ts_store','qlayer_store','wuptake_all_store',...
%     'wuptake_store','TR_can_all_store','smp_store','SWC1_in','SWC2_in','SWC3_in','Ts1_in',...
%     'Ts2_in','Ts3_in','Ts4_in','Ts5_in','Ts6_in','Ts7_in','Ts8_in','Ts9_in',...
%     'Ts10_in','dzsmm','zhsmm','znsmm','nl_soil','theta_dry',...
%     'porsl', 'rootfr','CONSTANTS','PARAMS','VERTSTRUC',...
%     'volliqli_store','Tli_store','qinflL_store','psili_store','Ts_store');
%     load ('tsoil_corr.mat','fTdl');
% end
% Dongkook Woo - Comment end

% Dongkook Woo insert
function  [VARIABLES, SWITCHES, PARAMS, FORCING]  = ...
                core_N(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE); 
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%               Loading Parameters and Forcings for CN model            %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Juan Quijano                                           %
%   Modified by  : Dongkook Woo                                           %
%   Date         : January 11, 2014                                       %
%-------------------------------------------------------------------------%
%                                                                         %    
%   Loading parameters and forcings Soil Carbon and Nitrogen Model        %
%                                                                         %
%-------------------------------------------------------------------------%
%%
if VARIABLES.timestep == 1
    % Dongkook Woo insert end
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'Soil_C_pool1'
    load (fullpath_forcings)
    
    nl_soil=PARAMS.nl_soil;
    N=length(CNratio_species1_crop);
    nspecies=PARAMS.CanStruc.nspecies;
    
    %VARIABLES.SOIL.litterthickness = 0.03;
    VARIABLES.SOIL.litterthickness = VARIABLES.SOIL.dzlit_m;
    %N = 17520;
    
    %*************************************************************************%
    %                             SWITCHES                                   *%
    %*************************************************************************%
    % switches ();
    % SWITCHES.CN_type = 0;             % 1 = 3 pools Include humus
    % % 2 = 3 pools steady stated in humus
    % % 0 = 2 pools humus and litter are coupled
    if Soil_C_pool1 == 1
        SWITCHES.CN_type = 0;   % 1 = 3 pools Include humus
        % 2 = 3 pools steady stated in humus
        % 0 = 2 pools humus and litter are coupled
    elseif Soil_C_pool1 == 0
        SWITCHES.CN_type = 1;
    end
    
    %SWITCHES.initialcond = 3;         % 1 = load data with file 2. fill with desire numbers 3. same as Cl_bloi
    SWITCHES.CN.Bioturbation = 1;     % 0 No Bioturbation, 1 Including Bioturbation (include litter)
    SWITCHES.Recycling = 0;           % 1 = yes, there is recycling of NO3 in the soil by HR
    % 0 = No, there is not recycling of NO3 in the soil by HR
    %SWITCHES.CN.consmethod = 2;       % Method to calulate bioturbation in the computation of constants kl,kd
    % 1. using bioturbation fluxes
    % 2. using an exponential function
    %SWITCHES.comconstants = 0;        % 1. Compute Constants. 0. Does not compute Constants (load from file)
    
    % SWITCHES.save = 0;
    % SWITCHES.savetype = 0;
    % SWITCHES.fitlitter = 0;
    load './Temps/temp_variable.mat'...
        'N_denit' 'N_Adepo' 'N_Fert' 'N_Uptake_RB' 'N_Fix' 'N_Remo'
    SWITCHES.CN.Denitrification = N_denit; % 1 = Including denitrification 0= No denitrification
    SWITCHES.CN.AtmosDeposition = N_Adepo; % 1 = Including Atmospheric Nitrogen Deposition 0= No Atmospheric Nitrogen Deposition
    SWITCHES.CN.Fertilizer      = N_Fert; % 1 = Including Fertilizer application 0= No Fertilizer application              = 0;
    SWITCHES.CN.NupRootBiomass  = N_Uptake_RB; % 1 = Computing N uptake with root biomass 
                                               % 0 = Computing N uptake with plant C dynamics
    SWITCHES.CN.leach_on        = 1;%       leach_on = Leaching ON or OF [1 0]     
    SWITCHES.CN.N_Fix           = N_Fix; % 1= including fixation 0 no
    SWITCHES.CN.N_Remo          = N_Remo; % 1=including fixation 0 no
    SWITCHES.Active_Nuptake     = 1;
    SWITCHES.Recycling          = 0; % Switch that decides if recycling of nitrate is on or off
    %*************************************************************************%
    %                            PARAMETERS                                  *%
    %*************************************************************************%
    
    % FOR CARBON NITROGEN SOIL MODEL
    %PARAMS.CN.factorrec = 0.3;     %
    %Fraction of C and N in leaves and fine roots
    %PARAMS.CN.fCleaves=0.55;
    %PARAMS.CN.fCroots=0.5;
    %PARAMS.CN.Nrecyc=0.0;      % Amount of Nitrogen that is recycled before leaves drop
    %PARAMS.CN.CR_ratio=0.36;   % Portion of C in roots compared to leaves
    
    load './Temps/temp_variable.mat'...
        'DOY_start' 'DOY_end'    
    PARAMS.CN.DOY_start = DOY_start;
    PARAMS.CN.DOY_end = DOY_end;
    
    %  MODIFICATIONS TO INCLUDE IN BIG CODE
    load './Temps/temp_variable.mat'...
        'para_soilCN'
    
    %    PARAMS.CN.a_Amm = 0.05;       % [-]
    %    PARAMS.CN.a_Nit = 1;          % [-]
    if isempty(cell2mat(para_soilCN(8,2)))
        PARAMS.CN.a_Amm = cell2mat(para_soilCN(8,3)); % [cm2/year]
    else
        PARAMS.CN.a_Amm = cell2mat(para_soilCN(8,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(9,2)))
        PARAMS.CN.a_Nit = cell2mat(para_soilCN(9,3)); % [cm2/year]
    else
        PARAMS.CN.a_Nit = cell2mat(para_soilCN(9,2)); % [cm2/year]
    end
    %PARAMS.CN.DEMp = 0.2;         % [gN / m^3 / d]
    %PARAMS.CN.DEMm = 0.5;         % [gN / m^3 / d]
    %    PARAMS.CN.kd = ones(12,1)*8.5*10^-3;     % [d^-1]
    %    PARAMS.CN.kl = ones(12,1)*6.5*10^-6;     % [m^3 / d / gC]
    %    PARAMS.CN.kh = ones(12,1)*2.5*10^-7;     % [m^3 / d / g?]
    %PARAMS.CN.kn = 0.6/100;       % [m^3 / d / gC]
    
    %    PARAMS.CN.ki_Amm = 1;         % [m^3 / d / gC]
    %    PARAMS.CN.ki_Nit = 1;         % [m^3 / d / gC]
    %    PARAMS.CN.ku_Amm = 0.1;       % Maximum Immobilization [m^3 / d / gC]
    %    PARAMS.CN.ku_Nit = 0.1;       % Maximum Immobilization [m^3 / d / gC]
    %    PARAMS.CN.rhmin = 0.2;
    if isempty(cell2mat(para_soilCN(4,2)))
        PARAMS.CN.ki_Amm = cell2mat(para_soilCN(4,3)); % [cm2/year]
    else
        PARAMS.CN.ki_Amm = cell2mat(para_soilCN(4,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(5,2)))
        PARAMS.CN.ki_Nit = cell2mat(para_soilCN(5,3)); % [cm2/year]
    else
        PARAMS.CN.ki_Nit = cell2mat(para_soilCN(5,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(6,2)))
        PARAMS.CN.ku_Amm = cell2mat(para_soilCN(6,3)); % [cm2/year]
    else
        PARAMS.CN.ku_Amm = cell2mat(para_soilCN(6,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(7,2)))
        PARAMS.CN.ku_Nit = cell2mat(para_soilCN(7,3)); % [cm2/year]
    else
        PARAMS.CN.ku_Nit = cell2mat(para_soilCN(7,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(14,2)))
        PARAMS.CN.rhmin = cell2mat(para_soilCN(14,3)); % [cm2/year]
    else
        PARAMS.CN.rhmin = cell2mat(para_soilCN(14,2)); % [cm2/year]
    end
    
    %PARAMS.CN.rr = 0.7;
    PARAMS.CN.koae = 1 ; % Organic Assimilation Efficiency parameter
    %PARAMS.CN.Ts_reduc_on = 0;
    %PARAMS.CN.Ts_max = 10;
    %PARAMS.CN.leach_on = 1;
    %     if (SWITCHES.soilCN_on && SWITCHES.CN.Bioturbation)
    %         CNb(1:PARAMS.nl_soil+1) = 11.8;
    %         load ('tsoil_corr.mat','fTdl');
    %         PARAMS.CN.fTdl = fTdl;
    %     else
    %         CNb(1:PARAMS.nl_soil) = 11.8;
    %     end
    %     PARAMS.CN.CNb= CNb(:);
    %PARAMS.CN.CNh = 22;
    
    
    %bioturbation
    %PARAMS.CN.bio.biofactor = 0.5; % Factor used to compute Bioturbation rate[gr/m2]
    % Is a factor of the Total Biomass that comes
    % into the soil
    %PARAMS.kbio = 0.0000215;
    %PARAMS.CN.bio.dsoil = 2000*1000;%         % Bulk density Mineral Soil [gr/m3]
    %PARAMS.CN.bio.lini = 2;%         % layer where difussion starts
    %PARAMS.CN.bio.lfin = 12;%         % layer where difussion finish
    %PARAMS.CN.Clitter = 60000;        % [gr/m3]
    %PARAMS.CN.layerbio = 7;   % Maximum layer for bioturbation
    %PARAMS.CN.D = 4;         % [cm2/year]
    
    
    %PARAMS.CN.leafspan = 3;               % leaf span shrubs [years]
    
    
    
    
    %PARAMS.K2=0.02;
    if isempty(cell2mat(para_soilCN(12,2)))
        PARAMS.CN.K2 = cell2mat(para_soilCN(12,3)); % [cm2/year]
    else
        PARAMS.CN.K2 = cell2mat(para_soilCN(12,2)); % [cm2/year]
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--Needed to change--%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    PARAMS.CN.fTdl = fTdl;
    
    % Try to chage! -> in CN_bioturbation
    %PARAMS.CN.fTdl = ones(17520*5,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Dongkook Woo - Comment
    %    permat = permute(wuptake_all_store,[2,3,1]);
    % Dongkook Woo - Comment end
    %********************************************************************
    
    %  MODIFICATIONS TO INCLUDE IN BIG CODE
    load './Temps/temp_variable.mat'...
        'Soil_C_pool1' 'opt_root_para' 'N_Uptake_RB' 'para_soilCN'
    
    % PARAMS.CN.Clitter = 60000;
    % PARAMS.CN.D = 4;         % [cm2/year]
    % PARAMS.Soil.rootarq.diam = [0.25 0.3 0.35 0.45]; % Root Diameter [mm]
    % PARAMS.Soil.rootarq.dist = [0.5 0.24 0.13 0.13]; % Distribution
    % PARAMS.Soil.rootarq.SRL = [22 28 38 20];        % Specific Root Length [m/gr]
    % PARAMS.Soil.rootarq.deltarZD = 0.1; % Definition zone of disturbance [mm]
    
    if isempty(cell2mat(para_soilCN(10,2)))
        PARAMS.CN.D = cell2mat(para_soilCN(10,3)); % [cm2/year]
    else
        PARAMS.CN.D = cell2mat(para_soilCN(10,2)); % [cm2/year]
    end
    if N_Uptake_RB == 1
        PARAMS.Soil.rootarq.diam = cell2mat(opt_root_para(1,2:end)); % Root Diameter [mm]
        PARAMS.Soil.rootarq.dist = cell2mat(opt_root_para(2,2:end)); % Distribution
        PARAMS.Soil.rootarq.SRL  = cell2mat(opt_root_para(3,2:end));        % Specific Root Length [m/gr]
        if isempty(cell2mat(para_soilCN(19,2)))
            PARAMS.Soil.rootarq.deltarZD = cell2mat(para_soilCN(17,3)); % Definition zone of disturbance [mm]
        else
            PARAMS.Soil.rootarq.deltarZD = cell2mat(para_soilCN(17,3)); % Definition zone of disturbance [mm]
        end
    elseif N_Uptake_RB == 0
        PARAMS.Soil.rootarq.diam = nan; % Root Diameter [mm]
        PARAMS.Soil.rootarq.dist = nan; % Distribution
        PARAMS.Soil.rootarq.SRL = nan;        % Specific Root Length [m/gr]
        PARAMS.Soil.rootarq.deltarZD = nan; % Definition zone of disturbance [mm]
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--Needed to change--%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load ('.\LOCAL_CODES\ROOT_SOIL\CN_MODEL\jaeger.mat');
    %PARAMS.CN.jaeger = jaeger_1942;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PARAMS.vmaxnit = 0.5;%0.037; %[mumolNO3/g_rr/hr]
    % PARAMS.kmnit = 72;%15; %[mmolNO3/m^3]
    % PARAMS.vmaxamm = 0.053; %[mumolNO3/g_rr/hr]
    % PARAMS.kmamm = 15; %[mmolNO3/m^3]
    
    if isempty(cell2mat(para_soilCN(18,2)))
        PARAMS.vmaxnit = cell2mat(para_soilCN(18,3)); %[mumolNO3/g_rr/hr]
    else
        PARAMS.vmaxnit = cell2mat(para_soilCN(18,2)); %[mumolNO3/g_rr/hr]
    end
    if isempty(cell2mat(para_soilCN(19,2)))
        PARAMS.kmnit = cell2mat(para_soilCN(19,3)); %[mmolNO3/m^3]
    else
        PARAMS.kmnit = cell2mat(para_soilCN(19,2)); %[mmolNO3/m^3]
    end
    if isempty(cell2mat(para_soilCN(20,2)))
        PARAMS.vmaxamm = cell2mat(para_soilCN(20,3)); %[mumolNO3/g_rr/hr]
    else
        PARAMS.vmaxamm = cell2mat(para_soilCN(20,2)); %[mumolNO3/g_rr/hr]
    end
    if isempty(cell2mat(para_soilCN(21,2)))
        PARAMS.kmamm = cell2mat(para_soilCN(21,3)); %[mmolNO3/m^3]
    else
        PARAMS.kmamm = cell2mat(para_soilCN(21,2)); %[mmolNO3/m^3]
    end
    
    
    
    
    
    
    
    
    
    
    
    load './Temps/temp_variable.mat'...
        'para_soilCN'
    % PARAMS.CN.kn = 0.0001;
    % PARAMS.CN.CNb = ones(13,1)*11.8;
    if isempty(cell2mat(para_soilCN(2,2)))
        PARAMS.CN.kn = cell2mat(para_soilCN(2,3));
    else
        PARAMS.CN.kn = cell2mat(para_soilCN(2,2));
    end
    if isempty(cell2mat(para_soilCN(1,2)))
        PARAMS.CN.CNb = ones(PARAMS.nl_soil+1,1)*cell2mat(para_soilCN(1,3));
    else
        PARAMS.CN.CNb = ones(PARAMS.nl_soil+1,1)*cell2mat(para_soilCN(1,2));
    end
    
    nspecies = PARAMS.CanStruc.nspecies;
    
    
    % PARAMS.kbio = 0.0000215;
    % PARAMS.CN.layerbio = 7;   % Maximum layer for bioturbation
    if isempty(cell2mat(para_soilCN(11,2)))
        PARAMS.kbio = cell2mat(para_soilCN(11,3));
    else
        PARAMS.kbio = cell2mat(para_soilCN(11,2));
    end
    
    
    %CN_nutrients_parameters;
    %TBMla = FORCING.TBMla;
    %CNh = PARAMS.CN.CNh;
    if isempty(cell2mat(para_soilCN(13,2)))
        PARAMS.CN.CNh = cell2mat(para_soilCN(13,3));
    else
        PARAMS.CN.CNh = cell2mat(para_soilCN(15,2));
    end
    if isempty(cell2mat(para_soilCN(3,2)))
        PARAMS.CN.rr = cell2mat(para_soilCN(3,3));
    else
        PARAMS.CN.rr = cell2mat(para_soilCN(3,2));
    end
    
    
    
    
    load './Temps/temp_variable.mat'...
        'para_soilCN' 'num_root1' 'para_soilCN_crop1' 'para_soilCN_crop2' 'para_soilCN_crop3' 'para_soilCN_crop4' 'Sim_species' ...
        'para_soilCN_crop11' 'para_soilCN_crop22' 'para_soilCN_crop33' 'para_soilCN_crop44' 'Sim_species_con'
    % PARAMS.CN.CR_ratio = 0.8;
    if N_Uptake_RB == 1
        if Sim_species == 1
            if isempty(cell2mat(para_soilCN_crop1(1,2)))
                PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop1(1,3));
            else
                PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop1(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop2(1,2)))
                PARAMS.CN.CR_ratio(2) = cell2mat(para_soilCN_crop2(1,3));
            else
                PARAMS.CN.CR_ratio(2) = cell2mat(para_soilCN_crop2(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop3(1,2)))
                PARAMS.CN.CR_ratio(3) = cell2mat(para_soilCN_crop3(1,3));
            else
                PARAMS.CN.CR_ratio(3) = cell2mat(para_soilCN_crop3(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop4(1,2)))
                PARAMS.CN.CR_ratio(4) = cell2mat(para_soilCN_crop4(1,3));
            else
                PARAMS.CN.CR_ratio(4) = cell2mat(para_soilCN_crop4(1,2));
            end
        elseif Sim_species == 0
            if Sim_species_con == 1
                if isempty(cell2mat(para_soilCN_crop1(1,2)))
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop1(1,3));
                else
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop1(1,2));
                end
            elseif Sim_species_con == 2
                if isempty(cell2mat(para_soilCN_crop2(1,2)))
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop2(1,3));
                else
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop2(1,2));
                end
            elseif Sim_species_con == 3
                if isempty(cell2mat(para_soilCN_crop3(1,2)))
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop3(1,3));
                else
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop3(1,2));
                end
            elseif Sim_species_con == 4
                if isempty(cell2mat(para_soilCN_crop4(1,2)))
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop4(1,3));
                else
                    PARAMS.CN.CR_ratio(1) = cell2mat(para_soilCN_crop4(1,2));
                end
            end
        end
    elseif N_Uptake_RB == 0
        if Sim_species == 1
            if isempty(cell2mat(para_soilCN_crop11(1,2)))
                PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop11(1,3));
            else
                PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop11(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(1,2)))
                PARAMS.CN.fDEMamm(2) = cell2mat(para_soilCN_crop22(1,3));
            else
                PARAMS.CN.fDEMamm(2) = cell2mat(para_soilCN_crop22(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(1,2)))
                PARAMS.CN.fDEMamm(3) = cell2mat(para_soilCN_crop33(1,3));
            else
                PARAMS.CN.fDEMamm(3) = cell2mat(para_soilCN_crop33(1,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(1,2)))
                PARAMS.CN.fDEMamm(4) = cell2mat(para_soilCN_crop44(1,3));
            else
                PARAMS.CN.fDEMamm(4) = cell2mat(para_soilCN_crop44(1,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop11(2,2)))
                PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop11(2,3));
            else
                PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop11(2,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(2,2)))
                PARAMS.CN.fDEMnit(2) = cell2mat(para_soilCN_crop22(2,3));
            else
                PARAMS.CN.fDEMnit(2) = cell2mat(para_soilCN_crop22(2,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(2,2)))
                PARAMS.CN.fDEMnit(3) = cell2mat(para_soilCN_crop33(2,3));
            else
                PARAMS.CN.fDEMnit(3) = cell2mat(para_soilCN_crop33(2,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(2,2)))
                PARAMS.CN.fDEMnit(4) = cell2mat(para_soilCN_crop44(2,3));
            else
                PARAMS.CN.fDEMnit(4) = cell2mat(para_soilCN_crop44(2,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop11(3,2)))
                PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop11(3,3));
            else
                PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop11(3,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(3,2)))
                PARAMS.CN.fTortuo(2) = cell2mat(para_soilCN_crop22(3,3));
            else
                PARAMS.CN.fTortuo(2) = cell2mat(para_soilCN_crop22(3,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(3,2)))
                PARAMS.CN.fTortuo(3) = cell2mat(para_soilCN_crop33(3,3));
            else
                PARAMS.CN.fTortuo(3) = cell2mat(para_soilCN_crop33(3,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(3,2)))
                PARAMS.CN.fTortuo(4) = cell2mat(para_soilCN_crop44(3,3));
            else
                PARAMS.CN.fTortuo(4) = cell2mat(para_soilCN_crop44(3,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop11(4,2)))
                PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop11(4,3));
            else
                PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop11(4,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(4,2)))
                PARAMS.CN.fRd(2) = cell2mat(para_soilCN_crop22(4,3));
            else
                PARAMS.CN.fRd(2) = cell2mat(para_soilCN_crop22(4,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(4,2)))
                PARAMS.CN.fRd(3) = cell2mat(para_soilCN_crop33(4,3));
            else
                PARAMS.CN.fRd(3) = cell2mat(para_soilCN_crop33(4,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(4,2)))
                PARAMS.CN.fRd(4) = cell2mat(para_soilCN_crop44(4,3));
            else
                PARAMS.CN.fRd(4) = cell2mat(para_soilCN_crop44(4,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop11(5,2)))
                PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop11(5,3));
            else
                PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop11(5,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(5,2)))
                PARAMS.CN.NfI(2) = cell2mat(para_soilCN_crop22(5,3));
            else
                PARAMS.CN.NfI(2) = cell2mat(para_soilCN_crop22(5,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(5,2)))
                PARAMS.CN.NfI(3) = cell2mat(para_soilCN_crop33(5,3));
            else
                PARAMS.CN.NfI(3) = cell2mat(para_soilCN_crop33(5,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(5,2)))
                PARAMS.CN.NfI(4) = cell2mat(para_soilCN_crop44(5,3));
            else
                PARAMS.CN.NfI(4) = cell2mat(para_soilCN_crop44(5,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop11(6,2)))
                PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop11(6,3));
            else
                PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop11(6,2));
            end
            if isempty(cell2mat(para_soilCN_crop22(6,2)))
                PARAMS.CN.NrI(2) = cell2mat(para_soilCN_crop22(6,3));
            else
                PARAMS.CN.NrI(2) = cell2mat(para_soilCN_crop22(6,2));
            end
            if isempty(cell2mat(para_soilCN_crop33(6,2)))
                PARAMS.CN.NrI(3) = cell2mat(para_soilCN_crop33(6,3));
            else
                PARAMS.CN.NrI(3) = cell2mat(para_soilCN_crop33(6,2));
            end
            if isempty(cell2mat(para_soilCN_crop44(6,2)))
                PARAMS.CN.NrI(4) = cell2mat(para_soilCN_crop44(6,3));
            else
                PARAMS.CN.NrI(4) = cell2mat(para_soilCN_crop44(6,2));
            end
        elseif Sim_species == 0
            if Sim_species_con == 1
                if isempty(cell2mat(para_soilCN_crop11(1,2)))
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop11(1,3));
                else
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop11(1,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop11(2,2)))
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop11(2,3));
                else
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop11(2,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop11(3,2)))
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop11(3,3));
                else
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop11(3,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop11(4,2)))
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop11(4,3));
                else
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop11(4,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop11(5,2)))
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop11(5,3));
                else
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop11(5,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop11(6,2)))
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop11(6,3));
                else
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop11(6,2));
                end
                
            elseif Sim_species_con == 2                
                if isempty(cell2mat(para_soilCN_crop22(1,2)))
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop22(1,3));
                else
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop22(1,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop22(2,2)))
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop22(2,3));
                else
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop22(2,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop22(3,2)))
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop22(3,3));
                else
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop22(3,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop22(4,2)))
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop22(4,3));
                else
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop22(4,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop22(5,2)))
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop22(5,3));
                else
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop22(5,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop22(6,2)))
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop22(6,3));
                else
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop22(6,2));
                end
            elseif Sim_species_con == 3                
                if isempty(cell2mat(para_soilCN_crop33(1,2)))
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop33(1,3));
                else
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop33(1,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop33(2,2)))
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop33(2,3));
                else
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop33(2,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop33(3,2)))
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop33(3,3));
                else
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop33(3,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop33(4,2)))
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop33(4,3));
                else
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop33(4,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop33(5,2)))
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop33(5,3));
                else
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop33(5,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop33(6,2)))
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop33(6,3));
                else
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop33(6,2));
                end
            elseif Sim_species_con == 4                                
                if isempty(cell2mat(para_soilCN_crop44(1,2)))
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop44(1,3));
                else
                    PARAMS.CN.fDEMamm(1) = cell2mat(para_soilCN_crop44(1,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop44(2,2)))
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop44(2,3));
                else
                    PARAMS.CN.fDEMnit(1) = cell2mat(para_soilCN_crop44(2,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop44(3,2)))
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop44(3,3));
                else
                    PARAMS.CN.fTortuo(1) = cell2mat(para_soilCN_crop44(3,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop44(4,2)))
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop44(4,3));
                else
                    PARAMS.CN.fRd(1) = cell2mat(para_soilCN_crop44(4,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop44(5,2)))
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop44(5,3));
                else
                    PARAMS.CN.NfI(1) = cell2mat(para_soilCN_crop44(5,2));
                end
                
                if isempty(cell2mat(para_soilCN_crop44(6,2)))
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop44(6,3));
                else
                    PARAMS.CN.NrI(1) = cell2mat(para_soilCN_crop44(6,2));
                end
            end
        end
    end
    % Reduce according to the number of species
    if Sim_species == 1
        if N_Uptake_RB == 1
            if nspecies == 1
                PARAMS.CN.CR_ratio(2:4) = [];
            elseif nspecies == 2
                PARAMS.CN.CR_ratio(3:4) = [];
            elseif nspecies == 3
                PARAMS.CN.CR_ratio(4)   = [];
            elseif nspecies == 4
                PARAMS.CN.CR_ratio      = PARAMS.CN.CR_ratio;
            end
        elseif N_Uptake_RB == 0
            if nspecies == 1
                PARAMS.CN.fDEMamm(2:4) = [];
                PARAMS.CN.fDEMnit(2:4) = [];
                PARAMS.CN.fTortuo(2:4) = [];
                PARAMS.CN.fRd(2:4) = [];
                PARAMS.CN.NfI(2:4) = [];
                PARAMS.CN.NrI(2:4) = [];
            elseif nspecies == 2
                PARAMS.CN.fDEMamm(3:4) = [];
                PARAMS.CN.fDEMnit(3:4) = [];
                PARAMS.CN.fTortuo(3:4) = [];
                PARAMS.CN.fRd(3:4) = [];
                PARAMS.CN.NfI(3:4) = [];
                PARAMS.CN.NrI(3:4) = [];
            elseif nspecies == 3
                PARAMS.CN.fDEMamm(4)   = [];
                PARAMS.CN.fDEMnit(4)   = [];
                PARAMS.CN.fTortuo(4)   = [];
                PARAMS.CN.fRd(4)   = [];
                PARAMS.CN.NfI(4)   = [];
                PARAMS.CN.NrI(4)   = [];
            elseif nspecies == 4
                PARAMS.CN.fDEMamm = PARAMS.CN.fDEMamm;
                PARAMS.CN.fDEMnit = PARAMS.CN.fDEMnit;
                PARAMS.CN.fTortuo = PARAMS.CN.fTortuo;
                PARAMS.CN.fRd = PARAMS.CN.fRd;
                PARAMS.CN.NfI = PARAMS.CN.NfI;
                PARAMS.CN.NrI = PARAMS.CN.NrI;
            end
        end
    end
    
    
    
    
    
    
    
    %fc= 0.39;     %Saxton et al (1986)
    %fHorO = 0.;
    %load constants.mat;
    PARAMS.CN.fHorO =0;
    load './Temps/temp_variable.mat'...
        'dat_decom' 'dat_decom_litter'
    % kd=[kd(1); kd(:)];
    % PARAMS.CN.kd = kd(:);
    % kl=[kl(1); kl(:)];
    % PARAMS.CN.kl = kl(:);
    % kh=[kh(1); kh(:)];
    % PARAMS.CN.kh = kh(:);
    kd=[dat_decom_litter(1,2); dat_decom(:,2)];
    PARAMS.CN.kd = kd(:);
    kl=[dat_decom_litter(1,3); dat_decom(:,3)];
    PARAMS.CN.kl = kl(:);
    kh=[dat_decom_litter(1,4); dat_decom(:,4)];
    PARAMS.CN.kh = kh(:);
    
    load './Temps/temp_variable.mat'...
        'N_Adepo' 'N_Adepo_amm' 'N_Adepo_nit'
    if N_Adepo == 1
        PARAMS.CN.N_Adepo_amm = N_Adepo_amm; % [g/m2/yr]
        PARAMS.CN.N_Adepo_nit = N_Adepo_nit; % [g/m2/yr]
    elseif N_Adepo == 0
        PARAMS.CN.N_Adepo_amm = nan;
        PARAMS.CN.N_Adepo_nit = nan;
    end
    
    load './Temps/temp_variable.mat'...
        'N_Fert' 'N_Fert_DOY' 'N_Fert_amm' 'N_Fert_nit' 'N_Fert_urea'
    if N_Fert == 1
        PARAMS.CN.N_Fert_DOY  = N_Fert_DOY; % [DOY]
        PARAMS.CN.N_Fert_amm  = N_Fert_amm; % [gN/m2]
        PARAMS.CN.N_Fert_nit  = N_Fert_nit; % [gN/m2]
        PARAMS.CN.N_Fert_urea = N_Fert_urea; % [gN/m2]
    elseif N_Fert == 0
        PARAMS.CN.N_Fert_DOY  = nan; % [DOY]
        PARAMS.CN.N_Fert_amm  = nan; % [gN/m2]
        PARAMS.CN.N_Fert_nit  = nan; % [gN/m2]
        PARAMS.CN.N_Fert_urea = nan; % [gN/m2]
    end

         
    if isempty(cell2mat(para_soilCN(15,2)))
        PARAMS.CN.VolUrea = cell2mat(para_soilCN(15,3));
    else
        PARAMS.CN.VolUrea = cell2mat(para_soilCN(15,2));
    end
    if isempty(cell2mat(para_soilCN(16,2)))
        PARAMS.CN.VolAmm = cell2mat(para_soilCN(16,3));
    else
        PARAMS.CN.VolAmm = cell2mat(para_soilCN(16,2));
    end
    
    %*************************************************************************%
    %                             Forcings                                   *%
    %*************************************************************************%
    
    %********************************************************************
    %  JUST FOR THIS CODE
    load './Temps/temp_variable.mat'...
        'fullpath_forcings' 'num_species' 'Sim_species' 'Sim_species_con'
    load (fullpath_forcings)
    
    %load ('CNinfo.mat','CNveg');
    %load ('FORCING2','FORCING');
    %PARAMS.CN.CNveg{1} = CNveg{1};
    %PARAMS.CN.CNveg{2} = CNveg{2};
    
    % Organize the data for simulation
    load './Temps/temp_variable.mat'...
        'working_forcings' 'DOY_start' 'DOY_end'
    load(working_forcings,...
        'year_crop','doy_crop');
    
    years = unique(year_crop)';
    doys = [DOY_start:DOY_end];
    inds = find(ismember(year_crop, years) & ismember(floor(round((doy_crop+1).*10^10)./10^10), doys));
    FORCING.doy_all=doy_crop(inds);
    
    if N_Uptake_RB == 1
        if Sim_species == 1
            if nspecies == 1
                PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
                FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species1_crop(inds);
            elseif nspecies == 2
                PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
                PARAMS.CN.CNveg{2} = CNratio_species2_crop(inds);
                FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
                FORCING.TBMla(2,:) = BMla_species2_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species1_crop(inds);
                FORCING.LMI(2,:)   = LMI_species2_crop(inds);
            elseif nspecies == 3
                PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
                PARAMS.CN.CNveg{2} = CNratio_species2_crop(inds);
                PARAMS.CN.CNveg{3} = CNratio_species3_crop(inds);
                FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
                FORCING.TBMla(2,:) = BMla_species2_crop(inds)';
                FORCING.TBMla(3,:) = BMla_species2_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species1_crop(inds);
                FORCING.LMI(2,:)   = LMI_species2_crop(inds);
                FORCING.LMI(3,:)   = LMI_species3_crop(inds);
            elseif nspecies == 4
                PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
                PARAMS.CN.CNveg{2} = CNratio_species2_crop(inds);
                PARAMS.CN.CNveg{3} = CNratio_species3_crop(inds);
                PARAMS.CN.CNveg{4} = CNratio_species4_crop(inds);
                FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
                FORCING.TBMla(2,:) = BMla_species2_crop(inds)';
                FORCING.TBMla(3,:) = BMla_species2_crop(inds)';
                FORCING.TBMla(4,:) = BMla_species2_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species1_crop(inds);
                FORCING.LMI(2,:)   = LMI_species2_crop(inds);
                FORCING.LMI(3,:)   = LMI_species3_crop(inds);
                FORCING.LMI(4,:)   = LMI_species4_crop(inds);
            end
            
        elseif Sim_species == 0
            
            if Sim_species_con == 1
                PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
                FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species1_crop(inds);
            elseif Sim_species_con == 2
                PARAMS.CN.CNveg{1} = CNratio_species2_crop(inds);
                FORCING.TBMla(1,:) = BMla_species2_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species2_crop(inds);
            elseif Sim_species_con == 3
                PARAMS.CN.CNveg{1} = CNratio_species3_crop(inds);
                FORCING.TBMla(1,:) = BMla_species3_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species3_crop(inds);
            elseif Sim_species_con == 4
                PARAMS.CN.CNveg{1} = CNratio_species4_crop(inds);
                FORCING.TBMla(1,:) = BMla_species4_crop(inds)';
                FORCING.LMI(1,:)   = LMI_species4_crop(inds);
            end
            
        end
        

    elseif N_Uptake_RB == 0
        
        if Sim_species == 1
            if nspecies == 1       
                FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
                
            elseif nspecies == 2
                FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
                
                FORCING.BMla(2,:)    = BMlaDEM_species2_crop(inds);
                FORCING.BMrb(2,:)    = BMrbDEM_species2_crop(inds);
                FORCING.CNabove(2,:) = CNratioA_species2_crop(inds);
                FORCING.CNbelow(2,:) = CNratioB_species2_crop(inds);
                FORCING.CNdem(2,:)   = NDEM_species2_crop(inds);
            elseif nspecies == 3
                FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
                
                FORCING.BMla(2,:)    = BMlaDEM_species2_crop(inds);
                FORCING.BMrb(2,:)    = BMrbDEM_species2_crop(inds);
                FORCING.CNabove(2,:) = CNratioA_species2_crop(inds);
                FORCING.CNbelow(2,:) = CNratioB_species2_crop(inds);
                FORCING.CNdem(2,:)   = NDEM_species2_crop(inds);
                
                FORCING.BMla(3,:)    = BMlaDEM_species3_crop(inds);
                FORCING.BMrb(3,:)    = BMrbDEM_species3_crop(inds);
                FORCING.CNabove(3,:) = CNratioA_species3_crop(inds);
                FORCING.CNbelow(3,:) = CNratioB_species3_crop(inds);
                FORCING.CNdem(3,:)   = NDEM_species3_crop(inds);
            elseif nspecies == 4
                FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
                
                FORCING.BMla(2,:)    = BMlaDEM_species2_crop(inds);
                FORCING.BMrb(2,:)    = BMrbDEM_species2_crop(inds);
                FORCING.CNabove(2,:) = CNratioA_species2_crop(inds);
                FORCING.CNbelow(2,:) = CNratioB_species2_crop(inds);
                FORCING.CNdem(2,:)   = NDEM_species2_crop(inds);
                
                FORCING.BMla(3,:)    = BMlaDEM_species3_crop(inds);
                FORCING.BMrb(3,:)    = BMrbDEM_species3_crop(inds);
                FORCING.CNabove(3,:) = CNratioA_species3_crop(inds);
                FORCING.CNbelow(3,:) = CNratioB_species3_crop(inds);
                FORCING.CNdem(3,:)   = NDEM_species3_crop(inds);
                
                FORCING.BMla(4,:)    = BMlaDEM_species4_crop(inds);
                FORCING.BMrb(4,:)    = BMrbDEM_species4_crop(inds);
                FORCING.CNabove(4,:) = CNratioA_species4_crop(inds);
                FORCING.CNbelow(4,:) = CNratioB_species4_crop(inds);
                FORCING.CNdem(4,:)   = NDEM_species4_crop(inds);
            end
            
        elseif Sim_species == 0
            
            if Sim_species_con == 1
                FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
            elseif Sim_species_con == 2
                FORCING.BMla(1,:)    = BMlaDEM_species2_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species2_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species2_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species2_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species2_crop(inds);
            elseif Sim_species_con == 3
                FORCING.BMla(1,:)    = BMlaDEM_species3_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species3_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species3_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species3_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species3_crop(inds);
            elseif Sim_species_con == 4
                FORCING.BMla(1,:)    = BMlaDEM_species4_crop(inds);
                FORCING.BMrb(1,:)    = BMrbDEM_species4_crop(inds);
                FORCING.CNabove(1,:) = CNratioA_species4_crop(inds);
                FORCING.CNbelow(1,:) = CNratioB_species4_crop(inds);
                FORCING.CNdem(1,:)   = NDEM_species4_crop(inds);
            end
        end
        
    end

    
    
    %*************************************************************************%
    %                        Initial Condition                               *%
    %*************************************************************************%
    
    % INITIALIZE NUTRIENTS. POOLS SIZE AND REQUIRED VARIABLES.
    VARIABLES.dCl_dT = 0;
    VARIABLES.DECl = 0;
    VARIABLES.BD = 0;
    
    
    
    load './Temps/temp_variable.mat'...
        'Soil_C_pool1' 'nutrient_int' 'nutrient_int_litter'
    %CN_IC();
    % Cl =[PARAMS.CN.Clitter ; Cl_bloi];
    % VARIABLES.Cl = Cl(:);
    % Cb =[Cb_bloi(1) ; Cb_bloi];
    % VARIABLES.Cb = Cb(:);
    % Ch =[30000 30000 5000 5000 5000 5000 5000 1000 1000 1000 100 100 100];
    % VARIABLES.Ch =Ch(:);
    % VARIABLES.Amm = 10*ones(nl_soil+1,1);%0.5*ones(nl_soil+1,1);
    % VARIABLES.Ammdd = 10*ones(nl_soil+1,1);%0.5*ones(nl_soil+1,1);
    % VARIABLES.Nit = 3*ones(nl_soil+1,1);%2*ones(nl_soil+1,1);
    % VARIABLES.Nitdd = 3*ones(nl_soil+1,1);%2*ones(nl_soil+1,1);
    % VARIABLES.CNl(1:nl_soil+1,1) = PARAMS.CN.CNb./(1-PARAMS.CN.rr)-15;  % Initialize with CN min
    % VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
    if Soil_C_pool1 == 1
        Cl =[nutrient_int_litter(2) ; nutrient_int(:,2)];
        VARIABLES.Cl = Cl(:);
        PARAMS.CN.Clitter=nutrient_int_litter(2);
        Cb =[nutrient_int_litter(3) ; nutrient_int(:,3)];
        VARIABLES.Cb = Cb(:);
        Ch =ones(PARAMS.nl_soil+1,1)*nan;
        VARIABLES.Ch =Ch(:);
        VARIABLES.Amm = [nutrient_int_litter(4) ; nutrient_int(:,4)];
        VARIABLES.Ammdd = [nutrient_int_litter(4) ; nutrient_int(:,4)];
        VARIABLES.Nit = [nutrient_int_litter(5) ; nutrient_int(:,5)];
        VARIABLES.Nitdd = [nutrient_int_litter(5) ; nutrient_int(:,5)];
        VARIABLES.CNl = [nutrient_int_litter(6) ; nutrient_int(:,6)];  % Initialize with CN min
        VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
    elseif Soil_C_pool1 ==0
        Cl =[nutrient_int_litter(2) ; nutrient_int(:,2)];
        VARIABLES.Cl = Cl(:);
        PARAMS.CN.Clitter=nutrient_int_litter(2);
        Cb =[nutrient_int_litter(4) ; nutrient_int(:,4)];
        VARIABLES.Cb = Cb(:);
        Ch =[nutrient_int_litter(3) ; nutrient_int(:,3)];
        VARIABLES.Ch =Ch(:);
        VARIABLES.Amm = [nutrient_int_litter(5) ; nutrient_int(:,5)];
        VARIABLES.Ammdd = [nutrient_int_litter(5) ; nutrient_int(:,5)];
        VARIABLES.Nit = [nutrient_int_litter(6) ; nutrient_int(:,6)];
        VARIABLES.Nitdd = [nutrient_int_litter(6) ; nutrient_int(:,6)];
        VARIABLES.CNl = [nutrient_int_litter(7) ; nutrient_int(:,7)];  % Initialize with CN min
        VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--Needed to change--%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to fix it! -> ok
    %VARIABLES.SOIL.psil=VARIABLES.SOIL.smp(1);
    %VARIABLES.SOIL.psil;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %*************************************************************************%
    %                               Allocaition                              *%
    %*************************************************************************%
    addpath('./LOCAL_CODES/ROOT_SOIL/CN_MODEL');
    ALLOCATE ();
    
end
% if SWITCHES.fitlitter;
% % dzlitter = fzero(@(kbio) fitlitter_dz (TBMla, kl, kd, Cb(1), PARAMS.CN.Clitter,...
% %            smp_store, Ts_store, PARAMS, VARIABLES, CONSTANTS, kbio), 0.000004); 
% load tsoil_corr;
% fitlitter_dz (TBMla, kl, kd, VARIABLES.Cb(2), PARAMS.CN.Clitter,...
%             psili_store, Tli_store, tsoil_corr, PARAMS, VARIABLES, CONSTANTS, PARAMS.kbio);
% 
% end

% for yy =1:1:500
%     if mod(yy, 50) == 0
%         yy
%     end
%     if (yy == 5 || yy == 10 || yy==20)
%         stop = 23;
%     end
%    tic
%     for tt=1:1:17520
%         if tt == 17500
%             toc
%             stop = 17500;
%         end
        % obtain info for this time step
%                 VARIABLES.SOIL.volliq = volliq_store(:,tt);
%                 VARIABLES.SOIL.smp = smp_store(:,tt);
%                 VARIABLES.SOIL.qinfl = qlayer_store(1,tt);
%                 VARIABLES.SOIL.qlayer = qlayer_store(:,tt);
%                 VARIABLES.SOIL.layeruptake = wuptake_store(:,tt);
%                 VARIABLES.CANOPY.TR_can_all = TR_can_all_store(:,tt);    
% %                vec = wuptake_all_store(100,:,1:2);
% %                VARIABLES.SOIL.layeruptake_all = squeeze(vec); 
% %                vec = permat(:,1:2,tt); 
%                 vec = permat(:,1:nspecies,tt); 
%                 VARIABLES.SOIL.layeruptake_all = vec;
%                 VARIABLES.SOIL.volliqli = volliqli_store(:,tt);
%                 VARIABLES.SOIL.qinflL = qinflL_store(tt);
%                 VARIABLES.SOIL.Tli = Tli_store(tt);
%                 VARIABLES.SOIL.psil = psili_store(tt);
%                 VARIABLES.SOIL.Ts = Ts_store(:,tt);  
%                 
%                 VARIABLES.timestep = tt;
%                 VARIABLES.CANOPY.TR_can_all=VARIABLES.CANOPY.TR_can_all';
                
                



%*************************************************************************%
%             Calculate Dynamics of Soil Carbon and Nitrogen             *%
%*************************************************************************%
[VARIABLES] = ...
    CN_Cycle_Porporato2(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE);
  
%                CN_STORE_DATA ();
%     end
%     toc
%     store_means ();    
        
%     if SWITCHES.save
%         if yy<=10
%             SAVE_DATA;
%             varsave =['save (''./SS/' num2str(yy)  outfile2 ''', ''strdata'')'];
%             eval(varsave);
%         else
%             if mod(yy,20) == 0
%                 SAVE_DATA;
%                 varsave =['save (''./SS/' num2str(yy)  outfile2 ''', ''strdata'')'];
%                 eval(varsave);
%             end
%         end
%     end
% end
