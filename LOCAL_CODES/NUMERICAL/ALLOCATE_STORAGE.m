% Script to allocate storage vectors / matrices

 N = length(doy);

% CANOPY VARIABLES
    
    % Radiation Absorption
        PARabs_sun_prof = NaN(nl_can,N);         %#ok<AGROW>
        PARabs_shade_prof = NaN(nl_can,N);       %#ok<AGROW>
        PARabs_canopy_prof = NaN(nl_can,N);      %#ok<AGROW>
        
        PARabs_sun_norm_prof = NaN(nl_can,N);    %#ok<AGROW>
        PARabs_shade_norm_prof = NaN(nl_can,N);  %#ok<AGROW>
        PARabs_canopy_norm_prof = NaN(nl_can,N); %#ok<AGROW>
                
        NIRabs_sun_prof = NaN(nl_can,N);         %#ok<AGROW>
        NIRabs_shade_prof = NaN(nl_can,N);       %#ok<AGROW>
        NIRabs_canopy_prof = NaN(nl_can,N);      %#ok<AGROW>
        
        NIRabs_sun_norm_prof = NaN(nl_can,N);    %#ok<AGROW>
        NIRabs_shade_norm_prof = NaN(nl_can,N);  %#ok<AGROW>
        NIRabs_canopy_norm_prof = NaN(nl_can,N); %#ok<AGROW>
        
        LWabs_can_prof = NaN(nl_can,N);    %#ok<AGROW>
        LWemit_can_prof = NaN(nl_can,N);   %#ok<AGROW> 
        
        LWemit_can_store = NaN(1,N);        
        LWemit_soil_store = NaN(1,N);
                
        SWout_store = NaN(1,N);             %#ok<AGROW>
        LWout_store = NaN(1,N);             %#ok<AGROW>
        fdiff_store = NaN(1,N);             %#ok<AGROW>
        LWdn_in = NaN(1,N);
%        LWdn2_in = NaN(1,N);
        refl_soil_store = NaN(1,N);
        
% FLUXES: 
    % Ecosystem Fluxes (Canopy + Soil)    
        Fc_eco_store = NaN(1,N);            %#ok<AGROW>
        LE_eco_store = NaN(1,N);            %#ok<AGROW> 
        H_eco_store = NaN(1,N);             %#ok<AGROW> 
        Rnrad_eco_store = NaN(1,N);         %#ok<AGROW>
        
    % Canopy Fluxes
		Ph_can_store = NaN(1,N);
		Ph_can_all_store = NaN(1, N, nspecies);
        An_can_store = NaN(1,N);            %#ok<AGROW>
        An_can_all_store = NaN(1, N, nspecies);            %#ok<AGROW        
        LE_can_store = NaN(1,N);            %#ok<AGROW> 
        LE_can_all_store = NaN(1, N, nspecies);            %#ok<AGROW        
        H_can_store = NaN(1,N);             %#ok<AGROW>
        H_can_all_store = NaN(1, N, nspecies);            %#ok<AGROW        
        TR_can_store = NaN(1,N);            %#ok<AGROW>
        TR_can_all_store = NaN(1, N, nspecies);            %#ok<AGROW                
        Rnrad_can_store = NaN(1,N);         %#ok<AGROW>
        dH_can_store = NaN(1,N);
        dH_can_all_store = NaN(1, N, nspecies);
        
    % Soil Fluxes
        Fc_soil_store = NaN(1,N);           %#ok<AGROW>
        H_soil_store = NaN(1,N);            %#ok<AGROW>
        LE_soil_store = NaN(1,N);           %#ok<AGROW>
        G_store = NaN(1,N);                 %#ok<AGROW>
        Tsurf_store = NaN(1,N);
        Rnrad_soil_store = NaN(1,N);        %#ok<AGROW>
        RH_soil_store = NaN(1,N);           %#ok<AGROW>
        E_soil_store = NaN(1,N);            %#ok<AGROW>
        Remain_soil_store = NaN(1,N);
        Remain_can_store = NaN(1,N);
		Remain_eco_store = NaN(1,N);

   % Mass Balance 
        MBerror_soil = NaN(1,N);                            % [mm/s]
        MBerror_littersoil = NaN(1,N);                      % [mm/s]
        MBerror_mbcan = NaN(1,N);                       % [mm/s]
        MBerror_mbcanlittersoil = NaN(1,N);                        % [mm/s]
        
        
   % Flux Profiles
        An_sun_prof = NaN(nl_can,N);        %#ok<AGROW>
        LE_sun_prof = NaN(nl_can,N);        %#ok<AGROW>
        H_sun_prof = NaN(nl_can,N);         %#ok<AGROW>
        Rnrad_sun_prof = NaN(nl_can,N);     %#ok<AGROW>
        TR_sun_prof = NaN(nl_can,N);        %#ok<AGROW>
        An_shade_prof = NaN(nl_can,N);      %#ok<AGROW>
        LE_shade_prof = NaN(nl_can,N);      %#ok<AGROW>
        H_shade_prof = NaN(nl_can,N);       %#ok<AGROW>
        Rnrad_shade_prof = NaN(nl_can,N);   %#ok<AGROW>
        TR_shade_prof = NaN(nl_can,N);      %#ok<AGROW>        
        % In some variables only I save prof info for all species
        Ph_sun_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        An_sun_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        Cs_sun_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        Ph_shade_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        An_shade_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        Cs_shade_prof_all = NaN(nl_can,N,nspecies);   %#ok<AGROW>
  
        
    % Mean Flux Profiles
        An_canopy_prof = NaN(nl_can,N);     %#ok<AGROW>
        LE_canopy_prof = NaN(nl_can,N);     %#ok<AGROW>
        H_canopy_prof = NaN(nl_can,N);      %#ok<AGROW>
        Rnrad_canopy_prof = NaN(nl_can,N);  %#ok<AGROW> 
        TR_canopy_prof = NaN(nl_can,N);     %#ok<AGROW>
        
    % Normalized Flux Profiles (ie. per unit LAI)    
    %   Sunlit
        An_sun_norm_prof = NaN(nl_can,N, nspecies);   %#ok<AGROW>
        LE_sun_norm_prof = NaN(nl_can,N, nspecies);   %#ok<AGROW>
        H_sun_norm_prof = NaN(nl_can,N, nspecies);    %#ok<AGROW>
        Rnrad_sun_norm_prof = NaN(nl_can,N);%#ok<AGROW>
        
    %   Shaded
        An_shade_norm_prof = NaN(nl_can,N,nspecies);  %#ok<AGROW>
        LE_shade_norm_prof = NaN(nl_can,N,nspecies);  %#ok<AGROW> 
        H_shade_norm_prof = NaN(nl_can,N,nspecies);   %#ok<AGROW>
        Rnrad_shade_norm_prof = NaN(nl_can,N); %#ok<AGROW>
        
    %   Canopy
        An_canopy_norm_prof = NaN(nl_can,N);  %#ok<AGROW>
        LE_canopy_norm_prof = NaN(nl_can,N);  %#ok<AGROW>
        H_canopy_norm_prof = NaN(nl_can,N);   %#ok<AGROW>
        Rnrad_canopy_norm_prof = NaN(nl_can,N); %#ok<AGROW>
         
        
    % Leaf States
        Tl_sun_prof = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        Tl_shade_prof = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Tl_sun_Ta_Diff = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Tl_shade_Ta_Diff = NaN(nl_can, N, nspecies);   %#ok<AGROW>
        
        psil_sun_prof = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        psil_shade_prof = NaN(nl_can, N, nspecies);   %#ok<AGROW>
        
        fsvg_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        fsvm_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        fsvg_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        fsvm_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        
        gsv_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        gsv_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        
        Ci_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        Ci_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        
        gbv_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        gbh_sun_prof = NaN(nl_can, N, nspecies);      %#ok<AGROW>
        gbv_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        gbh_shade_prof = NaN(nl_can, N, nspecies);    %#ok<AGROW>

        LAI_sun_prof = NaN(nl_can,N);      %#ok<AGROW>
        LAI_shade_prof = NaN(nl_can,N);    %#ok<AGROW>
        
        fsun_prof = NaN(nl_can,N);         %#ok<AGROW>
        fshade_prof = NaN(nl_can,N);       %#ok<AGROW>
        LSsunCON_store=NaN(1, N, nspecies);
        LSshaCON_store=NaN(1, N, nspecies);

     % Leaf convergence
     	cntleafsolution_store=NaN(1, N, nspecies);
        
    % Photosynthetic Biochemistry     
        Ph_limit_sun_store = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        Jc_C3_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        Jj_C3_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        Js_C3_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        Jc_C4_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        Jj_C4_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>  
        Js_C4_sun_store = NaN(nl_can, N, nspecies);       %#ok<AGROW>
        
        Ph_limit_shade_store = NaN(nl_can, N, nspecies);    %#ok<AGROW>
        Jc_C3_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>  
        Jj_C3_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Js_C3_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Jc_C4_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Jj_C4_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        Js_C4_shade_store = NaN(nl_can, N, nspecies);     %#ok<AGROW>
        
        dryfrac_prof = NaN(nl_can,N);           %#ok<AGROW>
        wetfrac_prof = NaN(nl_can,N);           %#ok<AGROW>
                
        ppt_ground_store = NaN(1,N);         %#ok<AGROW>
        qinfl_store = NaN(1,N);              %#ok<AGROW>
        
    % Mean Canopy Profiles
        Ci_canopy_prof = NaN(nl_can,N);     %#ok<AGROW>
        Tl_canopy_prof = NaN(nl_can,N);     %#ok<AGROW>
        gsv_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
        psil_canopy_prof = NaN(nl_can,N);   %#ok<AGROW>   
        fsvg_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
        fsvm_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
        gbv_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
        gbh_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
       
        
    % Mean Canopy States
        Ci_mean = NaN(1,N);                 %#ok<AGROW>
        Tl_mean = NaN(1,N);                 %#ok<AGROW>
        gsv_mean = NaN(1,N);                %#ok<AGROW>
        psil_mean = NaN(1,N);               %#ok<AGROW>
        fsv_mean = NaN(1,N);                %#ok<AGROW>
        fsvg_mean = NaN(1,N);
        fsvm_mean = NaN(1,N);
    % Canopy Microenvironment
        CAz_prof = NaN(nl_can,N);          %#ok<AGROW>
        TAz_prof = NaN(nl_can,N);          %#ok<AGROW>
        EAz_prof = NaN(nl_can,N);          %#ok<AGROW>
        Uz_prof = NaN(nl_can,N);           %#ok<AGROW>
 
    % canopy h2o storage
        Sh2o_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>
        Sh2o_canopy = NaN(1,N);              %#ok<AGROW>
        
        % Condensation
        Ch2o_canopy_prof = NaN(nl_can,N);    %#ok<AGROW>                                
        Ch2o_canopy = NaN(1,N);              %#ok<AGROW>
        
        % evaporation
        Evap_canopy_prof = NaN(nl_can,N);    %#ok<AGROW> 
        Evap_canopy = NaN(1,N);              %#ok<AGROW>
             
    
% SOIL VARIABLES
        volliq_store = NaN(nl_soil,N);    
        krad_store = NaN(nl_soil, N, nspecies);
        kax_store = NaN(nl_soil, N, nspecies);        
        hk_store = NaN(nl_soil,N);
        Ts_store = NaN(nl_soil,N);
        smp_store = NaN(nl_soil,N);
        rpp_store = NaN(nl_soil, N, nspecies);
		mberrormm_store = NaN(nl_soil,N);

        smpMPA_store = NaN(nl_soil,N);
        rppMPA_store = NaN(nl_soil, N, nspecies);

        smp_weight_store = NaN(1,N); 
        rpp_weight_store = NaN(1, N, nspecies); 
        volliq_weight_store = NaN(1,N);

        qlayer_store = NaN(nl_soil+1,N);
        bot_drainage_store = NaN(1,N);
        hor_drainage_store = NaN(1,N);
        %hor_drainage_lay_store = NaN(12,N);
        % Fixing the constant soil layer problem
        hor_drainage_lay_store = NaN(nl_soil,N);
        
        if (SWITCHES.soilheat_on)
            Ts_store = NaN(nl_soil,N);
			Gnew_store = NaN(nl_soil, N); 
			cpv_store = NaN(nl_soil, N); 
			TKsoil_store = NaN(nl_soil, N); 
			TKsol_store = NaN(nl_soil, N); 
        end

        if (SWITCHES.ns)
            wuptake_store = NaN(nl_soil,N);
			wuptake_all_store = NaN(nl_soil,N,nspecies);            
        end
        
        runoff_store = NaN(1,N);
        
% All of the CN model switch and parameters are in core_N 
%         if (SWITCHES.soilCN_on && SWITCHES.CN.Bioturbation)
%             Cl_store = NaN(nl_soil+1,N);
%             Ch_store = NaN(nl_soil+1,N);
%             Cb_store = NaN(nl_soil+1,N);
%             Nl_store = NaN(nl_soil+1,N);
%             Amm_store = NaN(nl_soil+1,N);
%             Nit_store = NaN(nl_soil+1,N);
%             CNl_store = NaN(nl_soil+1,N);
%             dCb_dt_store = NaN(nl_soil+1,N);
%             dCl_dt_store = NaN(nl_soil+1,N);
% 
%             UP_amm_store = NaN(nl_soil+1,N);
%             UP_amm_all_store = NaN(nl_soil+1,N, nspecies);
%             UP_amm_all_m2_store = NaN(nl_soil+1,N, nspecies);    
%             UP_nit_store = NaN(nl_soil+1,N);
%             UP_nit_all_store = NaN(nl_soil+1,N, nspecies);
%             UP_nit_all_m2_store = NaN(nl_soil+1,N, nspecies);
%             UP_N_m2_store = NaN(nl_soil+1,N);
%             F_HR_N_store = NaN(N);   
% 
% 
%             LCH_amm_store = NaN(nl_soil+1,N);
%             LCH_nit_store = NaN(nl_soil+1,N);
%             LCH_N_m2_store = NaN(nl_soil+1,N);
%             TLCH_amm_m2_store = NaN(nl_soil+1,N);
%             TLCH_nit_m2_store = NaN(nl_soil+1,N);
% 
%             MIN_net_store = NaN(nl_soil+1,N);
%             MIN_gross_store = NaN(nl_soil+1,N);
%             IMM_net_store= NaN(nl_soil+1,N);
%             IMM_gross_store = NaN(nl_soil+1,N);
%             Nreg_store = NaN(nl_soil+1,N);
%             DECl_store = NaN(nl_soil+1,N);
%             PHI_store = NaN(nl_soil+1,N);
%             phi_store = NaN(nl_soil+1,N);
%             fSd_store = NaN(nl_soil+1,N);
%             fTd_store = NaN(nl_soil+1,N);
%             mberrorN_store = NaN(1,N);
%             mberrorC_store = NaN(1,N); 
% 
%             %new
% 
%             litterthickness_store = NaN(1,N);
%             Cl_change_store = NaN(nl_soil+1,N);
%             Cflux_store=NaN(1,N);
%             Cl_bio_change_store=NaN(nl_soil+1,N);
%             Cl_bio_in_store=NaN(nl_soil+1,N);
%             Cl_bio_out_store=NaN(nl_soil+1,N);
%             bioNerror_store=NaN(1,N);
%             bioCerror_store=NaN(1,N);
%             dCb_dt_store = NaN(nl_soil+1,N);
%             dCl_dt_store = NaN(nl_soil+1,N);
%             ADD_store = NaN(nl_soil+1,N);
%             ADD_bio_store = NaN(nl_soil+1,N);
%             ADD_ex_store = NaN(nl_soil+1,N);
%             OUT_store = NaN(nl_soil+1,N);
%             ADD_net_store = NaN(nl_soil+1,N);
% 
%         elseif SWITCHES.soilCN_on
%             Cl_store = NaN(nl_soil,N);
%             Ch_store = NaN(nl_soil,N);
%             Cb_store = NaN(nl_soil,N);
%             Nl_store = NaN(nl_soil,N);
%             Amm_store = NaN(nl_soil,N);
%             Nit_store = NaN(nl_soil,N);
%             CNl_store = NaN(nl_soil,N);
%             dCb_dt_store = NaN(nl_soil,N);
%             dCl_dt_store = NaN(nl_soil,N);
% 
%             UP_amm_store = NaN(nl_soil,N);
%             UP_amm_all_store = NaN(nl_soil, N, nspecies);
%             UP_amm_all_m2_store = NaN(nl_soil, N, nspecies);    
%             UP_nit_store = NaN(nl_soil,N);
%             UP_nit_all_store = NaN(nl_soil, N, nspecies);
%             UP_nit_all_m2_store = NaN(nl_soil, N, nspecies);
%             UP_N_m2_store_store = NaN(nl_soil,N);
% 
%             LCH_amm_store = NaN(nl_soil,N);
%             LCH_nit_store = NaN(nl_soil,N);
%             LCH_N_m2_store = NaN(nl_soil,N);
%             TLCH_amm_m2_store = NaN(nl_soil,N);
%             TLCH_nit_m2_store = NaN(nl_soil,N);
% 
%             MIN_net_store = NaN(nl_soil,N);
%             MIN_gross_store = NaN(nl_soil,N);
%             IMM_net_store= NaN(nl_soil,N);
%             IMM_gross_store = NaN(nl_soil,N);
%             Nreg_store = NaN(nl_soil,N);
%             DECl_store = NaN(nl_soil,N);
%             PHI_store = NaN(nl_soil,N);
%             phi_store = NaN(nl_soil,N);
%             fSd_store = NaN(nl_soil,N);
%             fTd_store = NaN(nl_soil,N);
%             mberrorN_store = NaN(nl_soil,N);
%             mberrorC_store = NaN(nl_soil,N);
% 
%             %new
%             Cflux_store = NaN(1,N);
%             litterthickness_store = NaN(1,N);        
%             ADD_store = NaN(nl_soil+1,N);
%         end        

            Gs_store = NaN(1,N);
            Gsl_store = NaN(1,N);
            LE_soli_store = NaN(1,N);
            LE_liat_store = NaN(1,N);
            volliqli_store = NaN(1,N);
            qback_store = NaN(1,N);
            psili_MPa_store = NaN(1,N);
            psili_store = NaN(1,N);
            qinflL_store = NaN(1,N);
            net_qinflL_store = NaN(1,N);
            drainlitter_store = NaN(1,N);
            Tli_store = NaN(1,N);
            Tsl_store = NaN(1,N);
                       
            qadflux_store = NaN(1,N);
            qback_lit_store = NaN(1,N);
        
            
            % rhosn
            rhosn_store = NaN(1,N);                        % [kg/m3]
            % w
            wliqsl_store = NaN(1,N); 
            wicesl_store = NaN(1,N);
            deltawice_store = NaN(1,N);
            dH_store = NaN(1,N);            
            dS_store = NaN(1,N);            
            wsn_store = NaN(1,N);
            % z
            zicesl_store = NaN(1,N); 
            zliqsl_store = NaN(1,N);
            zliqsl_prev_store = NaN(1,N);            
            zicesl_prev_store = NaN(1,N);                        
            zsn_store = NaN(1,N);
            % vols
            volliqli_store = NaN(1,N);     
            voliceli_store = NaN(1,N);     
            volliqsn_store = NaN(1,N);     
            volicesn_store = NaN(1,N);         
            voltotsn_store = NaN(1,N);
            voltotli_store = NaN(1,N); 
            %Other variables 
            CRm_store = NaN(1,N);
            CRo_store = NaN(1,N);
            drhosn_store = NaN(1,N);            
            qinflL_store = NaN(1,N);
            net_qinflL_store = NaN(1,N);
            qinfl_store = NaN(1,N);
            drainlitter_store = NaN(1,N);
            pptrate_ground_store = NaN(1,N);
            Esoil_store = NaN(1,N);
            Esl_store = NaN(1,N);                
            snow_tcount_store = NaN(1,N);

        % Mass Balance
        
MB_store.PPTf = nan(1,N);                  % Total Rainfall  [mm/s]
MB_store.PPTgr = nan(1,N);                % Total Rainfall Reaching Ground  [mm/s]
MB_store.INFsoil = nan(1,N);            % Infiltration into Soil, Drainage From Litter [mm/s]
MB_store.EVli = nan(1,N);                  % Evaporation from litter [mm/s]
MB_store.EVsoil = nan(1,N);              % Evaporation from soil [mm/s]
MB_store.EVcanf = nan(1,N);              % Evaporation from canopy [mm/s]
MB_store.drainage = nan(1,N);          % Bottom Drainage from soil [mm/s]
MB_store.TR = nan(1,N,nspecies);                      % Transpiration all species [mm/s]
MB_store.runoff = nan(1,N);              % Transpiration all species [mm/s]
MB_store.COcanf = nan(1,N);              % Total canopy condensation [mm/s]

MB_store.dws_f = nan(1,N);                % Rate of change in soil water content [mm/s]
MB_store.dwc_f = nan(1,N);                % Rate of change in canopy water content [mm/s]
MB_store.dwl_f = nan(1,N);                % Rate of change in litter water content [mm/s]

MB_store.mbsoil = nan(1,N);  % rate of change in soil mass balance content [mm/day]
MB_store.mblittersoil = nan(1,N);  % rate of change in soil - litter mass balance content [mm/day]
MB_store.mbcan = nan(1,N);  % rate of change in canopy mass balance content [mm/day]
MB_store.mbcanlittersoil = nan(1,N);  % rate of change in soil - litter - canopy mass balance content [mm/day]



if (SWITCHES.entropy_on)
        % CANOPY
        % in
        SScandir_in_all_lay_store = nan(nl_can,N,nspecies);
        SScandir_in_all_store = nan(1,N,nspecies);
        SScandir_in_tot_store = nan(1,N);

        SScandif_in_all_lay_store = nan(nl_can,N,nspecies) ;
        SScandif_in_all_store = nan(1,N,nspecies);
        SScandif_in_tot_store = nan(1,N);

        SScanLW_in_all_lay_store = nan(nl_can,N,nspecies);
        SScanLW_in_all_store = nan(1,N,nspecies);
        SScanLW_in_tot_store = nan(1,N);

        SScan_in_all_lay_store = nan(nl_can,N,nspecies);
        SScan_in_all_store = nan(1,N,nspecies);
        SScan_in_store = nan(1,N);

        % out
        SScandir_out_all_lay_store = nan(nl_can,N,nspecies);
        SScandir_out_all_store = nan(1,N,nspecies);
        SScandir_out_tot_store = nan(1,N);

        SScandif_out_all_lay_store = nan(nl_can,N,nspecies);
        SScandif_out_all_store = nan(1,N,nspecies);
        SScandif_out_tot_store = nan(1,N);

        SScanLW_out_all_lay_store = nan(nl_can,N,nspecies);
        SScanLW_out_all_store = nan(1,N,nspecies);
        SScanLW_out_tot_store = nan(1,N);

        SScanLE_out_all_lay_store = nan(nl_can,N,nspecies);
        SScanLE_out_all_store = nan(1,N,nspecies);
        SScanLE_out_tot_store = nan(1,N);

        SScanH_out_all_lay_store = nan(nl_can,N,nspecies);
        SScanH_out_all_store = nan(1,N,nspecies);
        SScanH_out_tot_store = nan(1,N);

        % Total OUT
        SScan_out_all_lay_store = nan(nl_can,N,nspecies);
        SScan_out_all_store = nan(1,N,nspecies);
        SScan_out_store = nan(1,N);

        % TOTAL CANOPY
        SScan_all_lay_store = nan(nl_can,N,nspecies);
        SScan_all_store = nan(1,N,nspecies);
        SScan_tot_store = nan(1,N);
        
        % PHOTOSYNTHESIS        
        Epho_dif_store = nan(1,N);
        Epho_dir_store = nan(1,N);
        Epho_tot_store = nan(1,N);
        SSphodir_in_store = nan(1,N);
        SSphodif_in_store = nan(1,N);        


        % SOIL
        % in
        SSsoildir_in_tot_store = nan(1,N);
        SSsoildif_in_tot_store = nan(1,N);
        SSsoilLW_in_tot_store = nan(1,N);
        SSsoilG_in_tot_store = nan(1,N);
        SSsoil_in_store = nan(1,N);

        % out 
        SSsoildir_out_tot_store = nan(1,N);
        SSsoildif_out_tot_store = nan(1,N);
        SSsoilLW_out_tot_store = nan(1,N);
        SSsoilLE_out_tot_store = nan(1,N);
        SSsoilH_out_tot_store = nan(1,N);
        SSsoil_out_store = nan(1,N);

        % G, dS, dH SOIL
        SSsoilG_store = nan(1,N);
        SSsoildH_store = nan(1,N);
        SSsoildS_store = nan(1,N);        

        % dH CANOPY
        SSdHcan_store = nan(1,N);
        
        % TOTAL SOIL
        SSsoil_tot_store = nan(1,N);

        
        % TOTALS BY TYPE OF ENERGY        
        SSeco_totSWdir_store = nan(1,N);
        SSeco_totSWdif_store = nan(1,N);
        SSeco_totLW_store = nan(1,N);
        SSeco_totH_store = nan(1,N);
        SSeco_totLE_store = nan(1,N);
        SSeco_tot_store_test = nan(1,N);
        
                  
        % TOTAL ECOSYSTEM
        SSeco_tot_store = nan(1,N);
        
        % NET CANOPY
        
        SSdir_net_in_store = nan(1,N); 
        SSdir_net_out_store = nan(1,N);
        SSdif_net_in_store = nan(1,N);
        SSdif_net_out_store = nan(1,N);

        SSLW_net_in_store = nan(1,N);
        SSLW_net_inX_store = nan(1,N);        
        SSLW_net_out_store = nan(1,N);
        SSLW_net_outX_store = nan(1,N);        
        
        SSSW_net_inX_store = nan(1,N);        
        SSSW_net_outX_store = nan(1,N);        
        
        
        SSLE_net_store = nan(1,N);
        SSH_net_store = nan(1,N);
        SSG_net_store = nan(1,N);
        SSdHcan_net_store = nan(1,N);
        
        SSeco_net_in_store = nan(1,N);
        SSeco_net_out_store = nan(1,N); 
        SSeco_net_tot_store = nan(1,N);
        SSeco_net_ratio_store = nan(1,N);
        
        SSTl_net_store = nan(1,N);
        SSTl_net2_store = nan(1,N);
        
        % ENERGY AND TEFF
        
        SSeco_tot_store = nan(1,N);
        SWdir_in_store = nan(1,N);
        SWdir_out_store = nan(1,N);
        SWdif_in_store = nan(1,N);
        SWdif_out_store = nan(1,N);
        LWin_net_store = nan(1,N);
        LWout_net_store = nan(1,N);
        LWemi_net_store = nan(1,N);
        LE_net_store = nan(1,N);
        H_net_store = nan(1,N);
        G_net_store = nan(1,N);
        dHcan_net_store = nan(1,N);
        SSTeffent_store = nan(1,N);  
        SSXeffent_store = nan(1,N);
        
        SH2O_store = nan(nl_soil-1,N);
        EH2O_store = nan(nl_soil-1,N);
        
end


        