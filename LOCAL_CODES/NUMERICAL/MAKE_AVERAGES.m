% CALCULATE DIURNAL AVERAGES
    hr = unique(hour);
    hour=hour';
    % Modeled Canopy States
    [Ci_diurnal_mod, Ci_std_mod] = DIURNAL_AVERAGE (hr, hour, Ci_mean);
    [Tl_diurnal_mod, Tl_std_mod] = DIURNAL_AVERAGE (hr, hour, Tl_mean);
    [gsv_diurnal_mod, gsv_std_mod] = DIURNAL_AVERAGE (hr, hour, gsv_mean);
    [psil_diurnal_mod, psil_std_mod] = DIURNAL_AVERAGE (hr, hour, psil_mean);
    [Sh2o_diurnal_mod, Sh2o_std_mod] = DIURNAL_AVERAGE (hr, hour, Sh2o_canopy);
    [Ch2o_diurnal_mod, Ch2o_std_mod] = DIURNAL_AVERAGE (hr, hour, Ch2o_canopy);
    
    % Modeled Ecosystem Fluxes
    [Fc_diurnal_eco, Fc_std_eco] = DIURNAL_AVERAGE (hr, hour, Fc_eco_store);
    [LE_diurnal_eco, LE_std_eco] = DIURNAL_AVERAGE (hr, hour, LE_eco_store);
    [H_diurnal_eco, H_std_eco] = DIURNAL_AVERAGE (hr, hour, H_eco_store);
    [Rn_diurnal_eco, Rn_std_eco] = DIURNAL_AVERAGE (hr, hour, Rnrad_eco_store);
    
    % Modeled Canopy Fluxes
    [An_diurnal_can, An_std_can] = DIURNAL_AVERAGE (hr, hour, An_can_store);
    [LE_diurnal_can, LE_std_can] = DIURNAL_AVERAGE (hr, hour, LE_can_store);
    [H_diurnal_can, H_std_can] = DIURNAL_AVERAGE (hr, hour, H_can_store);
    [Rn_diurnal_can, Rn_std_can] = DIURNAL_AVERAGE (hr, hour, Rnrad_can_store);
    [Evap_diurnal_mod, Evap_std_mod] = DIURNAL_AVERAGE (hr, hour, Evap_canopy);
    
    % Modeled Soil Fluxes
    [Fc_diurnal_soil, Fc_std_soil] = DIURNAL_AVERAGE (hr, hour, Fc_soil_store);
    [LE_diurnal_soil, LE_std_soil] = DIURNAL_AVERAGE (hr, hour, LE_soil_store);
    [H_diurnal_soil, H_std_soil] = DIURNAL_AVERAGE (hr, hour, H_soil_store);
    [Rn_diurnal_soil, Rn_std_soil] = DIURNAL_AVERAGE (hr, hour, Rnrad_soil_store);
    [G_diurnal_mod, G_std_mod] = DIURNAL_AVERAGE (hr, hour, G_store);    
    
    % Observed Fluxes
    %[Fc_diurnal_obs, Fc_std_obs] = DIURNAL_AVERAGE (hr, hour, Fc_in);
    % Dongkook Woo - Edit
    %[LE_diurnal_obs, LE_std_obs] = DIURNAL_AVERAGE (hr, hour, LE_in);
    %[H_diurnal_obs, H_std_obs] = DIURNAL_AVERAGE (hr, hour, H_in);
    %[Rn_diurnal_obs, Rn_std_obs] = DIURNAL_AVERAGE (hr, hour, Rn_in);
    %[G_diurnal_obs, G_std_obs] = DIURNAL_AVERAGE (hr, hour, G_in);    
    % Dongkook Woo - Edit End
    
    % Observed Met
    [Rg_diurnal_obs, Rg_std_obs] = DIURNAL_AVERAGE (hr, hour, Rg_in');
    [LW_diurnal_obs, LW_std_obs] = DIURNAL_AVERAGE (hr, hour, LWdn_in);
    [VPD_diurnal_obs, VPD_std_obs] = DIURNAL_AVERAGE (hr, hour, VPD_in');
    [ea_diurnal_obs, ea_std_obs] = DIURNAL_AVERAGE (hr, hour, ea_in');
    [Ta_diurnal_obs, Ta_std_obs] = DIURNAL_AVERAGE (hr, hour, Ta_in');
    [U_diurnal_obs, U_std_obs] = DIURNAL_AVERAGE (hr, hour, U_in');
%    [Tskin_diurnal_obs, Tskin_std_obs] = DIURNAL_AVERAGE (hr, hour, Tskin_in);
    

% Dongkook - Edit

    if size(An_sun_norm_prof,3) == 1
        An_sun_norm_prof1=An_sun_norm_prof(:,:,1);
        An_sun_norm_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_sun_norm_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_sun_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        An_shade_norm_prof1=An_shade_norm_prof(:,:,1);
        An_shade_norm_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_shade_norm_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_shade_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_sun_prof1=Tl_sun_prof(:,:,1);
        Tl_sun_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_sun_prof1=gsv_sun_prof(:,:,1);
        gsv_sun_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_sun_prof1=gbv_sun_prof(:,:,1);
        gbv_sun_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_sun_prof1=Ci_sun_prof(:,:,1);
        Ci_sun_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_shade_prof1=Tl_shade_prof(:,:,1);
        Tl_shade_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_shade_prof1=gsv_shade_prof(:,:,1);
        gsv_shade_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_shade_prof1=gbv_shade_prof(:,:,1);
        gbv_shade_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_shade_prof1=Ci_shade_prof(:,:,1);
        Ci_shade_prof2=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
    elseif size(An_sun_norm_prof,3) == 2
        An_sun_norm_prof1=An_sun_norm_prof(:,:,1);
        An_sun_norm_prof2=An_sun_norm_prof(:,:,2);
        An_sun_norm_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_sun_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        An_shade_norm_prof1=An_shade_norm_prof(:,:,1);
        An_shade_norm_prof2=An_shade_norm_prof(:,:,2);
        An_shade_norm_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        An_shade_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_sun_prof1=Tl_sun_prof(:,:,1);
        Tl_sun_prof2=Tl_sun_prof(:,:,2);
        Tl_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_sun_prof1=gsv_sun_prof(:,:,1);
        gsv_sun_prof2=gsv_sun_prof(:,:,2);
        gsv_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_sun_prof1=gbv_sun_prof(:,:,1);
        gbv_sun_prof2=gbv_sun_prof(:,:,2);
        gbv_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_sun_prof1=Ci_sun_prof(:,:,1);
        Ci_sun_prof2=Ci_sun_prof(:,:,2);
        Ci_sun_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_shade_prof1=Tl_shade_prof(:,:,1);
        Tl_shade_prof2=Tl_shade_prof(:,:,2);
        Tl_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Tl_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_shade_prof1=gsv_shade_prof(:,:,1);
        gsv_shade_prof2=gsv_shade_prof(:,:,2);
        gsv_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gsv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_shade_prof1=gbv_shade_prof(:,:,1);
        gbv_shade_prof2=gbv_shade_prof(:,:,2);
        gbv_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        gbv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_shade_prof1=Ci_shade_prof(:,:,1);
        Ci_shade_prof2=Ci_shade_prof(:,:,2);
        Ci_shade_prof3=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
        Ci_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
    elseif size(An_sun_norm_prof,3) == 3
        An_sun_norm_prof1=An_sun_norm_prof(:,:,1);
        An_sun_norm_prof2=An_sun_norm_prof(:,:,2);
        An_sun_norm_prof3=An_sun_norm_prof(:,:,3);
        An_sun_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        An_shade_norm_prof1=An_shade_norm_prof(:,:,1);
        An_shade_norm_prof2=An_shade_norm_prof(:,:,2);
        An_shade_norm_prof3=An_shade_norm_prof(:,:,3);
        An_shade_norm_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_sun_prof1=Tl_sun_prof(:,:,1);
        Tl_sun_prof2=Tl_sun_prof(:,:,2);
        Tl_sun_prof3=Tl_sun_prof(:,:,3);
        Tl_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_sun_prof1=gsv_sun_prof(:,:,1);
        gsv_sun_prof2=gsv_sun_prof(:,:,2);
        gsv_sun_prof3=gsv_sun_prof(:,:,3);
        gsv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_sun_prof1=gbv_sun_prof(:,:,1);
        gbv_sun_prof2=gbv_sun_prof(:,:,2);
        gbv_sun_prof3=gbv_sun_prof(:,:,3);
        gbv_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_sun_prof1=Ci_sun_prof(:,:,1);
        Ci_sun_prof2=Ci_sun_prof(:,:,2);
        Ci_sun_prof3=Ci_sun_prof(:,:,3);
        Ci_sun_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Tl_shade_prof1=Tl_shade_prof(:,:,1);
        Tl_shade_prof2=Tl_shade_prof(:,:,2);
        Tl_shade_prof3=Tl_shade_prof(:,:,3);
        Tl_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gsv_shade_prof1=gsv_shade_prof(:,:,1);
        gsv_shade_prof2=gsv_shade_prof(:,:,2);
        gsv_shade_prof3=gsv_shade_prof(:,:,3);
        gsv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        gbv_shade_prof1=gbv_shade_prof(:,:,1);
        gbv_shade_prof2=gbv_shade_prof(:,:,2);
        gbv_shade_prof3=gbv_shade_prof(:,:,3);
        gbv_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));

        Ci_shade_prof1=Ci_shade_prof(:,:,1);
        Ci_shade_prof2=Ci_shade_prof(:,:,2);
        Ci_shade_prof3=Ci_shade_prof(:,:,3);
        Ci_shade_prof4=nan.*ones(size(An_sun_norm_prof,1),size(An_sun_norm_prof,2));
    elseif size(An_sun_norm_prof,3) == 4
        An_sun_norm_prof1=An_sun_norm_prof(:,:,1);
        An_sun_norm_prof2=An_sun_norm_prof(:,:,2);
        An_sun_norm_prof3=An_sun_norm_prof(:,:,3);
        An_sun_norm_prof4=An_sun_norm_prof(:,:,4);

        An_shade_norm_prof1=An_shade_norm_prof(:,:,1);
        An_shade_norm_prof2=An_shade_norm_prof(:,:,2);
        An_shade_norm_prof3=An_shade_norm_prof(:,:,3);
        An_shade_norm_prof4=An_shade_norm_prof(:,:,4);

        Tl_sun_prof1=Tl_sun_prof(:,:,1);
        Tl_sun_prof2=Tl_sun_prof(:,:,2);
        Tl_sun_prof3=Tl_sun_prof(:,:,3);
        Tl_sun_prof4=Tl_sun_prof(:,:,4);

        gsv_sun_prof1=gsv_sun_prof(:,:,1);
        gsv_sun_prof2=gsv_sun_prof(:,:,2);
        gsv_sun_prof3=gsv_sun_prof(:,:,3);
        gsv_sun_prof4=gsv_sun_prof(:,:,4);

        gbv_sun_prof1=gbv_sun_prof(:,:,1);
        gbv_sun_prof2=gbv_sun_prof(:,:,2);
        gbv_sun_prof3=gbv_sun_prof(:,:,3);
        gbv_sun_prof4=gbv_sun_prof(:,:,4);

        Ci_sun_prof1=Ci_sun_prof(:,:,1);
        Ci_sun_prof2=Ci_sun_prof(:,:,2);
        Ci_sun_prof3=Ci_sun_prof(:,:,3);
        Ci_sun_prof4=Ci_sun_prof(:,:,4);

        Tl_shade_prof1=Tl_shade_prof(:,:,1);
        Tl_shade_prof2=Tl_shade_prof(:,:,2);
        Tl_shade_prof3=Tl_shade_prof(:,:,3);
        Tl_shade_prof4=Tl_shade_prof(:,:,4);

        gsv_shade_prof1=gsv_shade_prof(:,:,1);
        gsv_shade_prof2=gsv_shade_prof(:,:,2);
        gsv_shade_prof3=gsv_shade_prof(:,:,3);
        gsv_shade_prof4=gsv_shade_prof(:,:,4);

        gbv_shade_prof1=gbv_shade_prof(:,:,1);
        gbv_shade_prof2=gbv_shade_prof(:,:,2);
        gbv_shade_prof3=gbv_shade_prof(:,:,3);
        gbv_shade_prof4=gbv_shade_prof(:,:,4);

        Ci_shade_prof1=Ci_shade_prof(:,:,1);
        Ci_shade_prof2=Ci_shade_prof(:,:,2);
        Ci_shade_prof3=Ci_shade_prof(:,:,3);
        Ci_shade_prof4=Ci_shade_prof(:,:,4);
    end
% Dongkook - Edit End




% Mean Profiles
    for ii = 1:nl_can
        
    % Radiation
        [PARabs_sun_diurnal_prof(ii,:), PARabs_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_sun_prof(ii,:));
        [PARabs_shade_diurnal_prof(ii,:), PARabs_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_shade_prof(ii,:));
        [PARabs_canopy_diurnal_prof(ii,:), PARabs_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_canopy_prof(ii,:));
        
        [PARabs_sun_norm_diurnal_prof(ii,:), PARabs_sun_norm_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_sun_norm_prof(ii,:));
        [PARabs_shade_norm_diurnal_prof(ii,:), PARabs_shade_norm_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_shade_norm_prof(ii,:));
        [PARabs_canopy_norm_diurnal_prof(ii,:), PARabs_canopy_norm_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, PARabs_canopy_norm_prof(ii,:));
        
        [NIRabs_sun_diurnal_prof(ii,:), NIRabs_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, NIRabs_sun_prof(ii,:));
        [NIRabs_shade_diurnal_prof(ii,:), NIRabs_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, NIRabs_shade_prof(ii,:));
        [NIRabs_canopy_diurnal_prof(ii,:), NIRabs_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, NIRabs_canopy_prof(ii,:));
        
        [LWabs_can_diurnal_prof(ii,:), LWabs_can_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, LWabs_can_prof(ii,:));
        [LWemit_can_diurnal_prof(ii,:), LWemit_can_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, LWemit_can_prof(ii,:));
        
    % Fluxes
        [An_sun_diurnal_prof(ii,:), An_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_prof(ii,:));
        [LE_sun_diurnal_prof(ii,:), LE_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, LE_sun_prof(ii,:));
        [H_sun_diurnal_prof(ii,:), H_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, H_sun_prof(ii,:));
        [Rnrad_sun_diurnal_prof(ii,:), Rnrad_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Rnrad_sun_prof(ii,:));
        
        [An_shade_diurnal_prof(ii,:), An_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_prof(ii,:));
        [LE_shade_diurnal_prof(ii,:), LE_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, LE_shade_prof(ii,:));
        [H_shade_diurnal_prof(ii,:), H_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, H_shade_prof(ii,:));
        [Rnrad_shade_diurnal_prof(ii,:), Rnrad_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Rnrad_shade_prof(ii,:));
        
        [An_canopy_diurnal_prof(ii,:), An_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_canopy_prof(ii,:));
        [LE_canopy_diurnal_prof(ii,:), LE_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, LE_canopy_prof(ii,:));
        [H_canopy_diurnal_prof(ii,:), H_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, H_canopy_prof(ii,:));
        [Rnrad_canopy_diurnal_prof(ii,:), Rnrad_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Rnrad_canopy_prof(ii,:));
        
        % Dongkook Edit 
        %[An_sun_norm_diurnal_prof(ii,:), An_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_norm_prof(ii,:)');
        %[An_shade_norm_diurnal_prof(ii,:), An_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_norm_prof(ii,:)');
        [An_sun_norm_diurnal_prof1(ii,:), An_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_norm_prof1(ii,:));
        [An_sun_norm_diurnal_prof2(ii,:), An_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_norm_prof2(ii,:));
        [An_sun_norm_diurnal_prof3(ii,:), An_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_norm_prof3(ii,:));
        [An_sun_norm_diurnal_prof4(ii,:), An_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_sun_norm_prof4(ii,:));
        
        [An_shade_norm_diurnal_prof1(ii,:), An_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_norm_prof1(ii,:));
        [An_shade_norm_diurnal_prof2(ii,:), An_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_norm_prof2(ii,:));
        [An_shade_norm_diurnal_prof3(ii,:), An_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_norm_prof3(ii,:));
        [An_shade_norm_diurnal_prof4(ii,:), An_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, An_shade_norm_prof4(ii,:));
     % Canopy States
        %[Tl_sun_diurnal_prof(ii,:), Tl_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_sun_prof(ii,:)');
        %[gsv_sun_diurnal_prof(ii,:), gsv_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_sun_prof(ii,:)');
        %[gbv_sun_diurnal_prof(ii,:), gbv_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_sun_prof(ii,:)');
        %[Ci_sun_diurnal_prof(ii,:), Ci_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_sun_prof(ii,:)');
        
        %[Tl_shade_diurnal_prof(ii,:), Tl_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_shade_prof(ii,:)');
        %[gsv_shade_diurnal_prof(ii,:), gsv_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_shade_prof(ii,:)');
        %[gbv_shade_diurnal_prof(ii,:), gbv_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_shade_prof(ii,:)');
        %[Ci_shade_diurnal_prof(ii,:), Ci_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_shade_prof(ii,:)');
        
        [Tl_sun_diurnal_prof1(ii,:), Tl_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_sun_prof1(ii,:));
        [Tl_sun_diurnal_prof2(ii,:), Tl_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_sun_prof2(ii,:));
        [Tl_sun_diurnal_prof3(ii,:), Tl_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_sun_prof3(ii,:));
        [Tl_sun_diurnal_prof4(ii,:), Tl_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_sun_prof4(ii,:));
        
        [gsv_sun_diurnal_prof1(ii,:), gsv_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_sun_prof1(ii,:));
        [gsv_sun_diurnal_prof2(ii,:), gsv_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_sun_prof2(ii,:));
        [gsv_sun_diurnal_prof3(ii,:), gsv_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_sun_prof3(ii,:));
        [gsv_sun_diurnal_prof4(ii,:), gsv_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_sun_prof4(ii,:));
        
        [gbv_sun_diurnal_prof1(ii,:), gbv_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_sun_prof1(ii,:));
        [gbv_sun_diurnal_prof2(ii,:), gbv_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_sun_prof2(ii,:));
        [gbv_sun_diurnal_prof3(ii,:), gbv_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_sun_prof3(ii,:));
        [gbv_sun_diurnal_prof4(ii,:), gbv_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_sun_prof4(ii,:));
        
        [Ci_sun_diurnal_prof1(ii,:), Ci_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_sun_prof1(ii,:));
        [Ci_sun_diurnal_prof2(ii,:), Ci_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_sun_prof2(ii,:));
        [Ci_sun_diurnal_prof3(ii,:), Ci_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_sun_prof3(ii,:));
        [Ci_sun_diurnal_prof4(ii,:), Ci_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_sun_prof4(ii,:));

        
        [Tl_shade_diurnal_prof1(ii,:), Tl_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_shade_prof1(ii,:));
        [Tl_shade_diurnal_prof2(ii,:), Tl_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_shade_prof2(ii,:));
        [Tl_shade_diurnal_prof3(ii,:), Tl_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_shade_prof3(ii,:));
        [Tl_shade_diurnal_prof4(ii,:), Tl_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_shade_prof4(ii,:));
        
        [gsv_shade_diurnal_prof1(ii,:), gsv_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_shade_prof1(ii,:));
        [gsv_shade_diurnal_prof2(ii,:), gsv_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_shade_prof2(ii,:));
        [gsv_shade_diurnal_prof3(ii,:), gsv_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_shade_prof3(ii,:));
        [gsv_shade_diurnal_prof4(ii,:), gsv_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_shade_prof4(ii,:));
        
        [gbv_shade_diurnal_prof1(ii,:), gbv_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_shade_prof1(ii,:));
        [gbv_shade_diurnal_prof2(ii,:), gbv_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_shade_prof2(ii,:));
        [gbv_shade_diurnal_prof3(ii,:), gbv_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_shade_prof3(ii,:));
        [gbv_shade_diurnal_prof4(ii,:), gbv_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_shade_prof4(ii,:));
        
        [Ci_shade_diurnal_prof1(ii,:), Ci_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_shade_prof1(ii,:));
        [Ci_shade_diurnal_prof2(ii,:), Ci_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_shade_prof2(ii,:));
        [Ci_shade_diurnal_prof3(ii,:), Ci_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_shade_prof3(ii,:));
        [Ci_shade_diurnal_prof4(ii,:), Ci_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_shade_prof4(ii,:));
        % Dongkook Edit - End
        
        [Tl_canopy_diurnal_prof(ii,:), Tl_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Tl_canopy_prof(ii,:));
        [gsv_canopy_diurnal_prof(ii,:), gsv_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gsv_canopy_prof(ii,:));
        [gbv_canopy_diurnal_prof(ii,:), gbv_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, gbv_canopy_prof(ii,:));
        [Ci_canopy_diurnal_prof(ii,:), Ci_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ci_canopy_prof(ii,:));
        
        [Sh2o_canopy_diurnal_prof(ii,:), Sh2o_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Sh2o_canopy_prof(ii,:));
        [Ch2o_canopy_diurnal_prof(ii,:), Ch2o_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Ch2o_canopy_prof(ii,:));
        [dryfrac_diurnal_prof(ii,:), dryfrac_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, dryfrac_prof(ii,:));
        [wetfrac_diurnal_prof(ii,:), wetfrac_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, wetfrac_prof(ii,:));
        [Evap_canopy_diurnal_prof(ii,:), Evap_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Evap_canopy_prof(ii,:));
        
                
     % Canopy Microenvironment
        [CAz_diurnal_prof(ii,:), CAz_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, CAz_prof(ii,:));
        [TAz_diurnal_prof(ii,:), TAz_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, TAz_prof(ii,:));
        [EAz_diurnal_prof(ii,:), EAz_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, EAz_prof(ii,:));
        [Uz_diurnal_prof(ii,:), Uz_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, Uz_prof(ii,:));
        
     % Tl-Ta   
     % Dongkook Woo - Edit
        %[TlmTa_sun_diurnal_prof(ii,:), TlmTa_sun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_sun_prof(ii,:)-TAz_prof(ii,:)));
        %[TlmTa_shade_diurnal_prof(ii,:), TlmTa_shade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_shade_prof(ii,:)-TAz_prof(ii,:)));
        [TlmTa_sun_diurnal_prof1(ii,:), TlmTa_sun_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_sun_prof1(ii,:)-TAz_prof(ii,:)));
        [TlmTa_sun_diurnal_prof2(ii,:), TlmTa_sun_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_sun_prof2(ii,:)-TAz_prof(ii,:)));
        [TlmTa_sun_diurnal_prof3(ii,:), TlmTa_sun_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_sun_prof3(ii,:)-TAz_prof(ii,:)));
        [TlmTa_sun_diurnal_prof4(ii,:), TlmTa_sun_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_sun_prof4(ii,:)-TAz_prof(ii,:)));
        
        [TlmTa_shade_diurnal_prof1(ii,:), TlmTa_shade_std_prof1(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_shade_prof1(ii,:)-TAz_prof(ii,:)));
        [TlmTa_shade_diurnal_prof2(ii,:), TlmTa_shade_std_prof2(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_shade_prof2(ii,:)-TAz_prof(ii,:)));
        [TlmTa_shade_diurnal_prof3(ii,:), TlmTa_shade_std_prof3(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_shade_prof3(ii,:)-TAz_prof(ii,:)));
        [TlmTa_shade_diurnal_prof4(ii,:), TlmTa_shade_std_prof4(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_shade_prof4(ii,:)-TAz_prof(ii,:)));
     % Dongkook Woo - Edit End           
        [TlmTa_canopy_diurnal_prof(ii,:), TlmTa_canopy_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, (Tl_canopy_prof(ii,:)-TAz_prof(ii,:)));
        
        
     % Sunlit / Shaded Fractions
        [fsun_diurnal_prof(ii,:), fsun_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, fsun_prof(ii,:));
        [fshade_diurnal_prof(ii,:), fshade_std_prof(ii,:)] = DIURNAL_AVERAGE (hr, hour, fshade_prof(ii,:));
    end
    
    

    