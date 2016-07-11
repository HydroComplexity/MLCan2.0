% Plot Environmental Variable Profiles against meteorological conditions


ninds = find(Rg_in < 10);
dinds = find(Rg_in >= 10);


% PLOT LIMITS
    Fcmin = -55;
    Fcmax = 20;
    LEmin = -70;
    LEmax = 630;
    Hmin = -100;
    Hmax = 230;

    Rgmin = 0;
    Rgmax = 980;
    Tamin = 8;
    Tamax = 33;
    VPDmin = 0;
    VPDmax = 2.2;
    Umin = 0;
    Umax = 4.5;
    

% BINNED AVERAGES
    Rg_binwidth = 100;
    VPD_binwidth = 0.15;
    Ta_binwidth = 2;
    U_binwidth = 0.4;
    
    Rg_binmin = 50;
    Rg_binmax = 850;
    
    Ta_binmin = 12;
    Ta_binmax = 30;
    
    VPD_binmin = VPD_binwidth;
    VPD_binmax = 1.8;
    
    U_binmin = U_binwidth;
    U_binmax = 4;
    
    Rg_bins = [Rg_binmin:Rg_binwidth:Rg_binmax];
    Ta_bins = [Ta_binmin:Ta_binwidth:Ta_binmax];
    VPD_bins = [VPD_binmin:VPD_binwidth:VPD_binmax];
    U_bins = [U_binmin:U_binwidth:U_binmax];
    
    dinds = find(Rg_in >= 10);
    ninds = find(Rg_in < 10);
       
        
  % BINNED MODELED PROFILE FLUXES - DAY
  inds = [1:length(An_canopy_prof)];
  for ii = 1:length(znc)
    % Rg Binning
        [An_canopy_Rgbins_avg(ii,:), An_canopy_Rgbins_std(ii,:), Rgbin_avg] = Bin_Average (Rg_in, An_canopy_prof(ii,inds), Rg_bins);
        [An_sun_Rgbins_avg(ii,:), An_sun_Rgbins_std(ii,:)] = Bin_Average (Rg_in, An_sun_prof(ii,inds), Rg_bins);
        [An_shade_Rgbins_avg(ii,:), An_shade_Rgbins_std(ii,:)] = Bin_Average (Rg_in, An_shade_prof(ii,inds), Rg_bins);

        [LE_canopy_Rgbins_avg(ii,:), LE_canopy_Rgbins_std(ii,:)] = Bin_Average (Rg_in, LE_canopy_prof(ii,inds), Rg_bins);
        [LE_sun_Rgbins_avg(ii,:), LE_sun_Rgbins_std(ii,:)] = Bin_Average (Rg_in, LE_sun_prof(ii,inds), Rg_bins);
        [LE_shade_Rgbins_avg(ii,:), LE_shade_Rgbins_std(ii,:)] = Bin_Average (Rg_in, LE_shade_prof(ii,inds), Rg_bins);

        [H_canopy_Rgbins_avg(ii,:), H_canopy_Rgbins_std(ii,:)] = Bin_Average (Rg_in, H_canopy_prof(ii,inds), Rg_bins);
        [H_sun_Rgbins_avg(ii,:), H_sun_Rgbins_std(ii,:)] = Bin_Average (Rg_in, H_sun_prof(ii,inds), Rg_bins);
        [H_shade_Rgbins_avg(ii,:), H_shade_Rgbins_std(ii,:)] = Bin_Average (Rg_in, H_shade_prof(ii,inds), Rg_bins);

        % Ta Binning
        [An_canopy_Tabins_avg(ii,:), An_canopy_Tabins_std(ii,:), Tabin_avg] = Bin_Average (Ta_in, An_canopy_prof(ii,inds), Ta_bins);
        [An_sun_Tabins_avg(ii,:), An_sun_Tabins_std(ii,:)] = Bin_Average (Ta_in, An_sun_prof(ii,inds), Ta_bins);
        [An_shade_Tabins_avg(ii,:), An_shade_Tabins_std(ii,:)] = Bin_Average (Ta_in, An_shade_prof(ii,inds), Ta_bins);

        [LE_canopy_Tabins_avg(ii,:), LE_canopy_Tabins_std(ii,:)] = Bin_Average (Ta_in, LE_canopy_prof(ii,inds), Ta_bins);
        [LE_sun_Tabins_avg(ii,:), LE_sun_Tabins_std(ii,:)] = Bin_Average (Ta_in, LE_sun_prof(ii,inds), Ta_bins);
        [LE_shade_Tabins_avg(ii,:), LE_shade_Tabins_std(ii,:)] = Bin_Average (Ta_in, LE_shade_prof(ii,inds), Ta_bins);

        [H_canopy_Tabins_avg(ii,:), H_canopy_Tabins_std(ii,:)] = Bin_Average (Ta_in, H_canopy_prof(ii,inds), Ta_bins);
        [H_sun_Tabins_avg(ii,:), H_sun_Tabins_std(ii,:)] = Bin_Average (Ta_in, H_sun_prof(ii,inds), Ta_bins);
        [H_shade_Tabins_avg(ii,:), H_shade_Tabins_std(ii,:)] = Bin_Average (Ta_in, H_shade_prof(ii,inds), Ta_bins);

        % VPD Binning
        [An_canopy_VPDbins_avg(ii,:), An_canopy_VPDbins_std(ii,:), VPDbin_avg] = Bin_Average (VPD_in, An_canopy_prof(ii,inds), VPD_bins);
        [An_sun_VPDbins_avg(ii,:), An_sun_VPDbins_std(ii,:)] = Bin_Average (VPD_in, An_sun_prof(ii,inds), VPD_bins);
        [An_shade_VPDbins_avg(ii,:), An_shade_VPDbins_std(ii,:)] = Bin_Average (VPD_in, An_shade_prof(ii,inds), VPD_bins);

        [LE_canopy_VPDbins_avg(ii,:), LE_canopy_VPDbins_std(ii,:)] = Bin_Average (VPD_in, LE_canopy_prof(ii,inds), VPD_bins);
        [LE_sun_VPDbins_avg(ii,:), LE_sun_VPDbins_std(ii,:)] = Bin_Average (VPD_in, LE_sun_prof(ii,inds), VPD_bins);
        [LE_shade_VPDbins_avg(ii,:), LE_shade_VPDbins_std(ii,:)] = Bin_Average (VPD_in, LE_shade_prof(ii,inds), VPD_bins);

        [H_canopy_VPDbins_avg(ii,:), H_canopy_VPDbins_std(ii,:)] = Bin_Average (VPD_in, H_canopy_prof(ii,inds), VPD_bins);
        [H_sun_VPDbins_avg(ii,:), H_sun_VPDbins_std(ii,:)] = Bin_Average (VPD_in, H_sun_prof(ii,inds), VPD_bins);
        [H_shade_VPDbins_avg(ii,:), H_shade_VPDbins_std(ii,:)] = Bin_Average (VPD_in, H_shade_prof(ii,inds), VPD_bins);

        % U Binning
        [An_canopy_Ubins_avg(ii,:), An_canopy_Ubins_std(ii,:), Ubin_avg] = Bin_Average (U_in, An_canopy_prof(ii,inds), U_bins);
        [An_sun_Ubins_avg(ii,:), An_sun_Ubins_std(ii,:)] = Bin_Average (U_in, An_sun_prof(ii,inds), U_bins);
        [An_shade_Ubins_avg(ii,:), An_shade_Ubins_std(ii,:)] = Bin_Average (U_in, An_shade_prof(ii,inds), U_bins);

        [LE_canopy_Ubins_avg(ii,:), LE_canopy_Ubins_std(ii,:)] = Bin_Average (U_in, LE_canopy_prof(ii,inds), U_bins);
        [LE_sun_Ubins_avg(ii,:), LE_sun_Ubins_std(ii,:)] = Bin_Average (U_in, LE_sun_prof(ii,inds), U_bins);
        [LE_shade_Ubins_avg(ii,:), LE_shade_Ubins_std(ii,:)] = Bin_Average (U_in, LE_shade_prof(ii,inds), U_bins);

        [H_canopy_Ubins_avg(ii,:), H_canopy_Ubins_std(ii,:)] = Bin_Average (U_in, H_canopy_prof(ii,inds), U_bins);
        [H_sun_Ubins_avg(ii,:), H_sun_Ubins_std(ii,:)] = Bin_Average (U_in, H_sun_prof(ii,inds), U_bins);
        [H_shade_Ubins_avg(ii,:), H_shade_Ubins_std(ii,:)] = Bin_Average (U_in, H_shade_prof(ii,inds), U_bins);

  end         
      
minAn = min([ min(min(An_canopy_Rgbins_avg)), min(min(An_canopy_Tabins_avg)), min(min(An_canopy_VPDbins_avg)), min(min(An_canopy_Ubins_avg)) ]);  
maxAn = max([ max(max(An_canopy_Rgbins_avg)), max(max(An_canopy_Tabins_avg)), max(max(An_canopy_VPDbins_avg)), max(max(An_canopy_Ubins_avg)) ]);  

minLE = min([ min(min(LE_canopy_Rgbins_avg)), min(min(LE_canopy_Tabins_avg)), min(min(LE_canopy_VPDbins_avg)), min(min(LE_canopy_Ubins_avg)) ]);  
maxLE = max([ max(max(LE_canopy_Rgbins_avg)), max(max(LE_canopy_Tabins_avg)), max(max(LE_canopy_VPDbins_avg)), max(max(LE_canopy_Ubins_avg)) ]);  

minH = min([ min(min(H_canopy_Rgbins_avg)), min(min(H_canopy_Tabins_avg)), min(min(H_canopy_VPDbins_avg)), min(min(H_canopy_Ubins_avg)) ]);  
maxH = max([ max(max(H_canopy_Rgbins_avg)), max(max(H_canopy_Tabins_avg)), max(max(H_canopy_VPDbins_avg)), max(max(H_canopy_Ubins_avg)) ]);  

      
figure(fignum); clf
% VARIATION WITH Rg
    subplot(3,4,1)
        hold on
        pcolor(Rgbin_avg, znc./znc(end), An_canopy_Rgbins_avg)
        shading interp
        colorbar
        caxis([minAn maxAn])
        title('\itA_n \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
        box on
        %axis([Rgmin Rgmax Fcmin Fcmax])
        %set(gca, 'XTickLabel', [])
        axis tight
        
    subplot(3,4,5)
        hold on
        pcolor(Rgbin_avg, znc./znc(end), LE_canopy_Rgbins_avg)
        shading interp
        colorbar
        caxis([minLE maxLE])
        title('\itLE \rm[W m^{-2}]', 'FontSize', 12)
        box on
        %axis([Rgmin Rgmax LEmin LEmax])
        %set(gca, 'XTickLabel', [])
        axis tight 
        
    subplot(3,4,9)
        hold on
        pcolor(Rgbin_avg, znc./znc(end), H_canopy_Rgbins_avg)
        shading interp
        colorbar
        caxis([minH maxH])
        title('\itH \rm[W m^{-2}]', 'FontSize', 12)
        xlabel('\itR_g \rm[W m^{-2}]', 'FontSize', 12)
        box on
        %axis([Rgmin Rgmax Hmin Hmax])
        axis tight

        
% VARIATION WITH Ta
    subplot(3,4,2)
        hold on
        pcolor(Tabin_avg, znc./znc(end), An_canopy_Tabins_avg)
        shading interp
        colorbar
        caxis([minAn maxAn])
        title('\itA_n \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
        box on
        %axis([Tamin Tamax Fcmin Fcmax])
        %set(gca, 'XTickLabel', [])
        axis tight
        
    subplot(3,4,6)
        hold on
        pcolor(Tabin_avg, znc./znc(end), LE_canopy_Tabins_avg)
        shading interp
        colorbar
        caxis([minLE maxLE])
        title('\itLE \rm[W m^{-2}]', 'FontSize', 12)
        box on
        %axis([Tamin Tamax LEmin LEmax])
        %set(gca, 'XTickLabel', [])
        axis tight 
        
    subplot(3,4,10)
        hold on
        pcolor(Tabin_avg, znc./znc(end), H_canopy_Tabins_avg)
        shading interp
        colorbar
        caxis([minH maxH])
        title('\itH \rm[W m^{-2}]', 'FontSize', 12)
        xlabel('\itT_a \rm[\circC]', 'FontSize', 12)
        box on
        %axis([Tamin Tamax Hmin Hmax])
        axis tight
        
                
% VARIATION WITH VPD
    subplot(3,4,3)
        hold on
        pcolor(VPDbin_avg, znc./znc(end), An_canopy_VPDbins_avg)
        shading interp
        colorbar
        caxis([minAn maxAn])
        title('\itA_n \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
        box on
        %axis([VPDmin VPDmax Fcmin Fcmax])
        %set(gca, 'XTickLabel', [])
        axis tight
        
    subplot(3,4,7)
        hold on
        pcolor(VPDbin_avg, znc./znc(end), LE_canopy_VPDbins_avg)
        shading interp
        colorbar
        caxis([minLE maxLE])
        title('\itLE \rm[W m^{-2}]', 'FontSize', 12)
        box on
        %axis([VPDmin VPDmax LEmin LEmax])
        %set(gca, 'XTickLabel', [])
        axis tight 
        
    subplot(3,4,11)
        hold on
        pcolor(VPDbin_avg, znc./znc(end), H_canopy_VPDbins_avg)
        shading interp
        colorbar
        caxis([minH maxH])
        title('\itH \rm[W m^{-2}]', 'FontSize', 12)
        xlabel('\itVPD \rm[kPa]', 'FontSize', 12)
        box on
        %axis([VPDmin VPDmax Hmin Hmax])
        axis tight

        
% VARIATION WITH U
    subplot(3,4,4)
        hold on
        pcolor(Ubin_avg, znc./znc(end), An_canopy_Ubins_avg)
        shading interp
        colorbar
        caxis([minAn maxAn])
        title('\itA_n \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
        box on
        %axis([Umin Umax Fcmin Fcmax])
        %set(gca, 'XTickLabel', [])
        axis tight
        
    subplot(3,4,8)
        hold on
        pcolor(Ubin_avg, znc./znc(end), LE_canopy_Ubins_avg)
        shading interp
        colorbar
        caxis([minLE maxLE])
        title('\itLE \rm[W m^{-2}]', 'FontSize', 12)
        box on
        %axis([Umin Umax LEmin LEmax])
        %set(gca, 'XTickLabel', [])
        axis tight 
        
    subplot(3,4,12)
        hold on
        pcolor(Ubin_avg, znc./znc(end), H_canopy_Ubins_avg)
        shading interp
        colorbar
        caxis([minH maxH])
        title('\itH \rm[W m^{-2}]', 'FontSize', 12)
        xlabel('\itU \rm[m s^{-1}]', 'FontSize', 12)
        box on
        %axis([Umin Umax Hmin Hmax])
        axis tight
        
         
%RePosition_3x4;        

