nl_soil=PARAMS.nl_soil;

load './Temps/temp_variable.mat'...
    'dat_root1' 'dat_root2' 'dat_root3' 'dat_root4'
zns = dat_root1(:,1);

% Soil layer thicknesses
dzs(1)  = 0.5*(zns(1)+zns(2));
dzs(nl_soil)= zns(nl_soil)-zns(nl_soil-1);
for j = 2:nl_soil-1
    dzs(j)= 0.5*(zns(j+1)-zns(j-1));
end
dzs=dzs';

% Soil layer interface depths from the surface [m]
zhs(nl_soil) = zns(nl_soil) + 0.5*dzs(nl_soil);
for j = 1:nl_soil-1
    zhs(j)= 0.5*(zns(j)+zns(j+1));
end
zhs=zhs';

znsmm = zns(:)*1000;      % [mm]
dzsmm = dzs(:)*1000;      % [mm]
zhsmm = zhs(:)*1000;      % [mm]

if (kspecies == 5 && Sim_species_con == 1 )
    rootfr = dat_root1(:,2);
    roottr = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root = sum(count_nl_root1 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 2 )
    rootfr = dat_root2(:,2);
    roottr = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root = sum(count_nl_root2 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 3 )
    rootfr = dat_root3(:,2);
    roottr = dat_root3(:,2);
    count_nl_root3=dat_root3(:,2);
    count_nl_root3(1)=1;
    nl_root = sum(count_nl_root3 ~= 0);
elseif (kspecies == 5 && Sim_species_con == 4 )
    rootfr = dat_root4(:,2);
    roottr = dat_root4(:,2);
    count_nl_root4=dat_root4(:,2);
    count_nl_root4(1)=1;
    nl_root = sum(count_nl_root4 ~= 0);
elseif (kspecies == 1)
    rootfr = dat_root1(:,2);
    roottr = rootfr;
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root = sum(count_nl_root1 ~= 0);
elseif (kspecies == 2)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root(2) = sum(count_nl_root2 ~= 0);
elseif (kspecies == 3)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
elseif (kspecies == 4)
    rootfr(:,1) = dat_root1(:,2);
    roottr(:,1) = dat_root1(:,2);
    count_nl_root1=dat_root1(:,2);
    count_nl_root1(1)=1;
    nl_root(1) = sum(count_nl_root1 ~= 0);
    rootfr(:,2) = dat_root2(:,2);
    roottr(:,2) = dat_root2(:,2);
    count_nl_root2=dat_root2(:,2);
    count_nl_root2(1)=1;
    nl_root(2) = sum(count_nl_root2 ~= 0);
    rootfr(:,3) = dat_root3(:,2);
    roottr(:,3) = dat_root3(:,2);
    count_nl_root3=dat_root3(:,2);
    count_nl_root3(1)=1;
    nl_root(3) = sum(count_nl_root3 ~= 0);
    rootfr(:,4) = dat_root4(:,2);
    roottr(:,4) = dat_root4(:,2);
    count_nl_root4=dat_root4(:,2);
    count_nl_root4(1)=1;
    nl_root(4) = sum(count_nl_root4 ~= 0);
end

% ASSIGN
VERTSTRUC.zns = zns;
VERTSTRUC.dzs = dzs;
VERTSTRUC.zhs = zhs;
VERTSTRUC.znsmm = znsmm;
VERTSTRUC.dzsmm = dzsmm;
VERTSTRUC.zhsmm = zhsmm;
VERTSTRUC.rootfr = rootfr;
VERTSTRUC.roottr = roottr;
VERTSTRUC.nl_root = nl_root;
VERTSTRUC.nl_soil = nl_soil;
PARAMS.Soil.nl_soil = nl_soil;

[VERTSTRUC] = SOIL_PROPERTIES(PARAMS, VERTSTRUC);
theta_dry = VERTSTRUC.theta_dry;
porsl = VERTSTRUC.porsl;
TK_dry = VERTSTRUC.TK_dry;
TK_sol = VERTSTRUC.TK_sol;
HC_sol = VERTSTRUC.HC_sol;

% ALLOCATE STORAGE FOR MODELLED VARIABLES
ALLOCATE_STORAGE;

% INITIALIZE CANOPY STATES
VARIABLES.CANOPY.gsv_sun = 0.01*ones(nl_can,nspecies);
VARIABLES.CANOPY.gsv_shade = 0.01*ones(nl_can,nspecies);
VARIABLES.CANOPY.TR = zeros(length(znc),nspecies);
VARIABLES.CANOPY.Sh2o_prof = zeros(length(znc),1);
VARIABLES.CANOPY.Tl_prev_dt = Ta_in(1) * ones(nl_can,1);

% INITIALIZE SOIL STATES
% initialize soil moisture
VARIABLES.SOIL.volliq = volliqinit;
VARIABLES.SOIL.smp = VERTSTRUC.psi0 .* (VARIABLES.SOIL.volliq ./ VERTSTRUC.porsl).^(-VERTSTRUC.bsw);

% Message box for soil moisture
if sum(volliqinit>=porsl) > 0
    msgbox({'Initial soil moisture is higher than saturated soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
    error('Initial soil moisture is higher than saturated soil moisture!')
end
if sum(volliqinit<=theta_dry) > 0
    msgbox({'Initial soil moisture is lower than residual soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
    error('Initial soil moisture is lower than residual soil moisture!')
end

% Initialize snow moisture variables
VARIABLES.SOIL.voltotsn = 1;
VARIABLES.SOIL.voltotli = 1;
VARIABLES.SOIL.volliqli = volliqliinit;
VARIABLES.SOIL.voliceli = 0;
VARIABLES.SOIL.volliqsn = 0;
VARIABLES.SOIL.volicesn = 0;
VARIABLES.SOIL.zliqsl = (VARIABLES.SOIL.dzlit_m*1000)*volliqliinit;
VARIABLES.SOIL.zicesl = 0;
VARIABLES.SOIL.zliqsl_prev = VARIABLES.SOIL.zliqsl;
VARIABLES.SOIL.zicesl_prev = VARIABLES.SOIL.zicesl;
VARIABLES.SOIL.zsn = VARIABLES.SOIL.zliqsl;
VARIABLES.SOIL.wliqsl = (VARIABLES.SOIL.zliqsl/1000)*PARAMS.Soil.rho_liq;
VARIABLES.SOIL.wicesl = 0;
VARIABLES.SOIL.wsn = VARIABLES.SOIL.wliqsl;
VARIABLES.SOIL.rhosn = 1000;

% Fixing the constant soil layer problem
VARIABLES.SOIL.volice = zeros(nl_soil,1);
VARIABLES.SOIL.snow_tcount = 0;
VARIABLES.SOIL.volliqli = volliqliinit ;   % Initial value of litter soil moisture



VARIABLES.CANOPY.Sh2o_can_prev = 0;

% initialize ice content & soil temperature
Ts = Tsinit;
VARIABLES.SOIL.Tli = Ta_in(1);
VARIABLES.SOIL.Tsl = Tslint;
VARIABLES.SOIL.Tlprev = Ta_in(1);
VARIABLES.SOIL.TKsoil=VERTSTRUC.TK_sol;

% INITIALIZE ROOT POTENTIAL
for ii=1:nspecies
    VARIABLES.ROOT.rpp_wgt(:,ii) =  VARIABLES.SOIL.smp(1);
    VARIABLES.ROOT.rpp(:,ii)= VARIABLES.SOIL.smp;
end

% PEDOTRANSFER FUNCTIONS
if SWITCHES.Pedofunctions
    [VERTSTRUC] = PEDOSOIL_PROPERTIES(PARAMS, VERTSTRUC, VARIABLES);
    porsl = VERTSTRUC.porsl;                                         % POROSITY
    psi0 = VERTSTRUC.psi0;                                           % MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
    bsw = VERTSTRUC.bsw;                                             % B PARAMETER BROKS AND COREY SHAPE PARAMETER
    Ksat = VERTSTRUC.HKsat;                                          % HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s]
    eff_poros = VERTSTRUC.eff_poros;
end


% RUN CANOPY-ROOT-SOIL MODEL
% LOOP OVER EACH YEAR TO RE-INITIALIZE CANOPY/SOIL STATES FOR EACH YEAR

% TimeBar 1/3
tff=yendinds(end);
t00=ybeginds(1);
hh = timebar('Progress','MLCan Simulation');

tic
for yy = 1:length(Each_year)
    % compute the range of time steps in current year
    yy;
    ybind = ybeginds(yy);
    yeind = yendinds(yy);
    
    % LOOP OVER EACH TIME PERIOD IN YEAR yy
    %for tt = ybind:yeind
    for tt = ybind:1:yeind
        
        tt;
        % TimeBar 2/3
        timebar(hh,(tt-t00)/(tff-t00))

        timestep = tt-ybind + 1;
        VARIABLES.timestep = timestep;
        [VERTSTRUC VARIABLES rootfr] = ROOT_RESPONSE_DRY(VARIABLES,...
            SWITCHES, VERTSTRUC, CONSTANTS, PARAMS, doy, smp_store);
        
        % FORCING CONDITIONS
        FORCING.doy = doy(tt);
        FORCING.Rg = Rg_in(tt);
        FORCING.Pa = Pa_in;
        if PARAMS.LWcom == 1
            FORCING.LWdn = LWdn_in(tt);
        end
        FORCING.zen = ZEN_in(tt);
        FORCING.U = U_in(tt);
        FORCING.ppt = PPT_in(tt);    % [mm]
        FORCING.Ta = Ta_in(tt);
        FORCING.ea = ea_in(tt);
        FORCING.Ca = CO2base;
        FORCING.ELEV=ELEV;
        if (~SWITCHES.soilheat_on)
            VARIABLES.SOIL.Ts = (Ta_in(tt)-5)*ones(size(zns));
        else
            VARIABLES.SOIL.Ts = Ts;
        end
        VARIABLES.SOIL.Ts = Ts;
        VARIABLES.SOIL.Tsurf=Ts(1);
        
        % CANOPY STRUCTURE
        if SWITCHES.plants
            for kk=1:1:nspecies
                if SWITCHES.LT == 1
                    LAILT = LAI_in(timestep,kk);
                    VERTSTRUC.LAIzall(:,kk) = LAILT*LADnorm_all(:,kk);
                else
                    VERTSTRUC.LAIzall(:,kk) = LAI_in(tt,kk)*LADnorm_all(:,kk);
                end
            end
        else
            VERTSTRUC.LAIzall(:,1:nspecies) = zeros(nl_can,nspecies);
        end
        LADnorm = sum(VERTSTRUC.LAIzall,2)/sum(LAI_in(tt,:));
        VERTSTRUC.LAIz = sum(VERTSTRUC.LAIzall,2);
        VERTSTRUC.LADz = VERTSTRUC.LAIz ./ dzc; % Total LAD distribution
        fLAIz =VERTSTRUC.LAIzall./(repmat(sum(VERTSTRUC.LAIzall,2),1,nspecies));
        fLAIz(isnan(fLAIz)) = 0; % Set to zero whenever there is not LAI at a given layer
        VERTSTRUC.fLAIz = fLAIz; % Fraction of LAI in each species at each relative height level

        % create vinds
        % 1. For the total canopy
        LADmax = (max(VERTSTRUC.LAIzall,[],2)); % Maximum LAD
        nvinds = find(LADmax<=0);
        vinds = find(LADmax>0);
        VERTSTRUC.vinds = vinds;
        VERTSTRUC.nvinds = nvinds;
        
        % 2. For All the species
        for kk=1:nspecies
            nvinds_all{kk} = find(VERTSTRUC.LAIzall(:,kk) <= 0);
            vinds_all{kk} = find(VERTSTRUC.LAIzall(:,kk) > 0);
        end
        VERTSTRUC.nvinds_all = nvinds_all;
        VERTSTRUC.vinds_all = vinds_all;
        
        % INITIALIZE CANOPY ENVIRONMENT
        VARIABLES.CANOPY.TAz = Ta_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.CAz = CO2base * ones(nl_can,1);
        VARIABLES.CANOPY.EAz = ea_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.PAz = Pa_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.Uz = U_in(tt) * ones(nl_can,1);
        
        VARIABLES.CANOPY.TR_sun = zeros(nl_can,nspecies);
        VARIABLES.CANOPY.TR_shade = zeros(nl_can,nspecies);
        
        % INITIALIZE CANOPY STATES
        VARIABLES.CANOPY.Tl_can_sun = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_can_shade = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_sun = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Tl_shade = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Ci_sun = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        VARIABLES.CANOPY.Ci_shade = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        
        
        % CANOPY MODEL SOLUTION      
        [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
            Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil,remaincan,remaineco,...
            Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
            An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
            Tl_sun, Tl_shade, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun, fsvm_sun,...
            fsvg_shade, fsvm_shade, Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
            LAI_sun, LAI_shade, fsun, fshade, ...
            Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
            Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
            PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
            LWabs_can, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
            Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
            dryfrac, wetfrac, Vz, VARIABLES, FORCING,...
            SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out,...
            SWsoildif_in, SWsoildif_out, LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON] = ...
            CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS);
        
        % SOLUTION OF SNOW-LITTER PACK DYNAMICS
        if SWITCHES.litter
            [VARIABLES] = FLUXES_WATER_SOIL_LITTER (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES);
        else
            [VARIABLES] = FLUXES_WATER_SOIL (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES);
        end
        
        % assign
        qinfl = VARIABLES.SOIL.qinfl;
        qinflL = VARIABLES.SOIL.qinfl;
        net_qinflL = VARIABLES.SOIL.net_qinflL;
        drainlitter = VARIABLES.SOIL.drainlitter;
        volliqli = VARIABLES.SOIL.volliqli;
        
        
        % Implicit Solution
        if (SWITCHES.ns)
            if tt == 321
                stop = 1;
            end
            [rpp,rpp_wgt,krad,kax,dwat,smp,bk,hk, ...
                qlayer,layeruptake,layeruptake_all,mberrormm,type, hor_drainage, hor_drainage_lay,flux_Ss]=ROOTSOIL(SWITCHES, VERTSTRUC,...
                PARAMS, VARIABLES, CONSTANTS, nspecies);
            VARIABLES.SOIL.flux_Ss =flux_Ss;
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            VARIABLES.SOIL.type = type;
        else
            if (SWITCHES.HR_on)
                [rpp, rpp_wgt, krad, kax] = ROOTS_HR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES);
            else
                [rpp, rpp_wgt, krad, kax] = ROOTS_NOHR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES);
            end
            
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            
            % Soil Moisture Solution
            [dwat, smp, hk, smp_wgt, thsatfrac_wgt, qlayer] = ...
                SOILMOISTURE(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, CONSTANTS);
            
            hor_drainage = nan;
            layeruptake = (smp - rpp).*krad;
            layeruptake_all = layeruptake;
            hor_drainage = 0;
            hor_drainage_lay = zeros(length(smp),1);
            mberrormm = nan;
        end
        
        % Update Volumetric Liquid Content
        volliq = VARIABLES.SOIL.volliq;
        volliq = volliq(:) + dwat(:);
        volliq = max(theta_dry, volliq);
        volliq = min(VERTSTRUC.eff_poros, volliq);
        
        % ASSIGN
        VARIABLES.SOIL.dwat = dwat;
        VARIABLES.SOIL.volliq = volliq;
        VARIABLES.SOIL.smp = smp;
        VARIABLES.SOIL.qlayer = qlayer;
        VARIABLES.SOIL.hor_drainage = hor_drainage;
        VARIABLES.SOIL.hor_drainage_lay = hor_drainage_lay;
        VARIABLES.SOIL.layeruptake = layeruptake;
        VARIABLES.SOIL.layeruptake_all = layeruptake_all;
        
        % RECOMPUTE MASS BALANCE INCLUDING THE FLUX BACK FROM
        % INFILTRATION SOLUTION
        if SWITCHES.litter
            [VARIABLES] = FLUXES_WATER_SOIL_LITTER_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);
        else
            [VARIABLES] = FLUXES_WATER_SOIL_BACK (VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);
        end
        % assign
        qinfl = VARIABLES.SOIL.qinfl;
        qinflL = VARIABLES.SOIL.qinfl;
        net_qinflL = VARIABLES.SOIL.net_qinflL;
        drainlitter = VARIABLES.SOIL.drainlitter;
        volliqli = VARIABLES.SOIL.volliqli;
        
        if (SWITCHES.soilheat_on)
            Ginto = VARIABLES.SOIL.Gs;
            
            % Soil Temperature Solution
            volice = 0;
            [VARIABLES] = SOILHEAT (Ginto, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS);
            Ts = VARIABLES.SOIL.Ts;
            cpv = VARIABLES.SOIL.cpv;
        end
        
        % ************************************************************************
        %     [PARAMS, VARIABLES] = Nitrogen_Plant(PARAMS, FORCING, VARIABLES, CONSTANTS,nspecies);
        % ************************************************************************
        if SWITCHES.soilCN_on
            % For nitrogen remobilization
            if tt == 1
                STORAGE.UP_nit_m2_store = zeros(nl_soil,nl_soil,1);
                STORAGE.UP_amm_m2_store = zeros(nl_soil,nl_soil,1);
            end

            [VARIABLES, SWITCHES, PARAMS, FORCING] = ...
                core_N(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE);
            
            CN_STORE_DATA ();
        end

        if (SWITCHES.entropy_on)            
            [SSresults] = ...
                COMPUENTROPY (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
                SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                SWout, fdiff,LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                LWout, Tsurf, FORCING,SWITCHES, CONSTANTS, PARAMS,...
                VARIABLES, VERTSTRUC);
        end
        
        [VARIABLES] = MASS_BALANCE (VARIABLES, CONSTANTS, PARAMS, FORCING, SWITCHES, VERTSTRUC, tt);
        
        STORE_DATA;
        
    end
    
    toc
end

% Timebar 3/3
close(hh);
