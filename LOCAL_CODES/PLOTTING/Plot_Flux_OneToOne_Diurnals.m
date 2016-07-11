% CANOPY-TOP FLUXES

    Rnvar = Rnrad_eco_store;
    Rg_min = 10;

    % CALCULATE REGRESSION COEFFICIENTS
  %  Fcqc_good_inds = find(Fc_qc==1);
    LEqc_good_inds = find(LE_qc==1);
    Hqc_good_inds = find(H_qc==1);
    Rn_good_inds = find(~isnan(Rn_in));
    
    % Total Sum of Squares (Data)
  %      TSSE_Fc = sum( (Fc_in(Fcqc_good_inds) - mean(Fc_in(Fcqc_good_inds))).^2 );
        TSSE_LE = sum( (LE_in(LEqc_good_inds) - mean(LE_in(LEqc_good_inds))).^2 );
        TSSE_H = sum( (H_in(Hqc_good_inds) - mean(H_in(Hqc_good_inds))).^2 );
        TSSE_Rn = sum( (Rn_in(Rn_good_inds) - mean(Rn_in(Rn_good_inds))).^2 );

    % Sum of Squared Errors
  %      SSE_Fc = sum( (Fc_in(Fcqc_good_inds) - Fc_eco_store(Fcqc_good_inds)).^2 );
        SSE_LE = sum( (LE_in(LEqc_good_inds) - LE_eco_store(LEqc_good_inds)).^2 );
        SSE_H = sum( (H_in(Hqc_good_inds) - H_eco_store(Hqc_good_inds)).^2 );
        SSE_Rn = sum( (Rn_in(Rn_good_inds) - Rnvar(Rn_good_inds)).^2 );
        
    % R-Squared Values
%        RS_Fc = 1 - SSE_Fc/TSSE_Fc;
        RS_LE = 1 - SSE_LE/TSSE_LE;
        RS_H = 1 - SSE_H/TSSE_H;
        RS_Rn = 1 - SSE_Rn/TSSE_Rn;
    %******************************

    
%    Fcqc_best_inds = find(Fc_qc==0);
    LEqc_best_inds = find(LE_qc==0);
    Hqc_best_inds = find(H_qc==0);
    Rn_good_inds = find(~isnan(Rn_in));
    
    pfFc = polyfit(Fc_in(Fcqc_best_inds), Fc_eco_store(Fcqc_best_inds), 1);
    pfLE = polyfit(LE_in(LEqc_best_inds), LE_eco_store(LEqc_best_inds), 1);
    pfH = polyfit(H_in(Hqc_best_inds), H_eco_store(Hqc_best_inds), 1);
    pfRn = polyfit(Rn_in(Rn_good_inds), Rnvar(Rn_good_inds), 1);
    
    Fcqc_good_inds = find(Fc_qc==1);
    LEqc_good_inds = find(LE_qc==1);
    Hqc_good_inds = find(H_qc==1);
    
    Fcqc_bad_inds = find(Fc_qc>1);
    LEqc_bad_inds = find(LE_qc>1);
    Hqc_bad_inds = find(H_qc>1); 
        
    
    minFc = min(min(Fc_in),min(Fc_eco_store)) - 5;
    maxFc = max(max(Fc_in),max(Fc_eco_store)) + 5;
    minLE = min(min(LE_in),min(LE_eco_store)) - 20;
    maxLE = max(max(LE_in),max(LE_eco_store)) + 20;
    minH = min(min(H_in),min(H_eco_store));
    maxH = max(max(H_in),max(H_eco_store));
    minRn = min(min(H_in(Rn_good_inds)),min(Rnvar));
    maxRn = max(max(H_in(Rn_good_inds)),max(Rnvar));
    
    minFc = -65; %min(minFc, -70);
    maxFc = 23; %max(maxFc, 20);
    minLE = -200; %min(minLE, -200);
    maxLE = 620; %max(maxLE, 700);
    minH = -200; %min(minH, -200);
    maxH = 620; %max(maxH, 700);
    minRn = -200;
    maxRn = 850;
    
    
    % TRUNCATE slopes, intercepts, R2
        pfFc = round(pfFc*100) / 100;
        pfLE = round(pfLE*100) / 100;
        pfH = round(pfH*100) / 100;
        pfRn = round(pfRn*100) / 100;
        
        RS_Fc = round(RS_Fc*100) / 100;
        RS_LE = round(RS_LE*100) / 100;
        RS_H = round(RS_H*100) / 100;
        RS_Rn = round(RS_Rn*100) / 100;
    
    figure(fignum); clf
            
    % FLUX ONE-TO-ONES
        subplot(4,2,1)
            hold on
            plot([minFc maxFc], [minFc maxFc], '-k')
            plot(Fc_in(Fcqc_best_inds), Fc_eco_store(Fcqc_best_inds), '.k', 'MarkerSize', 8)
            plot([minFc, maxFc], pfFc(1)*[minFc, maxFc]+pfFc(2), '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
            %plot([minFc, maxFc], pfFc_day(1)*[minFc, maxFc]+pfFc_day(2), '--r', 'LineWidth', 2)
            axis([minFc maxFc minFc maxFc])
            xlabel('\itF_c \rm(Obs) [\mumol m^{-2} s^{-1}]', 'FontSize', 12)
            ylabel('\itF_c \rm(Mod) [\mumol m^{-2} s^{-1}]', 'FontSize', 12)
            text(-60, 15, ['slope = ', num2str(pfFc(1))], 'FontSize', 12)
            text(-60, 5, ['int = ', num2str(pfFc(2))], 'FontSize', 12)
            text(-60, -5, ['R^2 = ', num2str(RS_Fc)], 'FontSize', 12)
            box on
            
        subplot(4,2,3)
            hold on
            plot([minLE maxLE], [minLE maxLE], '-k')
            plot(LE_in(LEqc_best_inds), LE_eco_store(LEqc_best_inds), '.k', 'MarkerSize', 8)
            plot([minLE, maxLE], pfLE(1)*[minLE, maxLE]+pfLE(2), '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
            %plot([minLE, maxLE], pfLE_day(1)*[minLE, maxLE]+pfLE_day(2), '--r', 'LineWidth', 2)            
            axis([minLE maxLE minLE maxLE])
            xlabel('\itLE \rm(Obs) [W m^{-2}]', 'FontSize', 12)
            ylabel('\itLE \rm(Mod) [W m^{-2}]', 'FontSize', 12)
            text(-150, 500, ['slope = ', num2str(pfLE(1))], 'FontSize', 12)
            text(-150, 400, ['int = ', num2str(pfLE(2))], 'FontSize', 12)
            text(-150, 300, ['R^2 = ', num2str(RS_LE)], 'FontSize', 12)
            box on
            
        subplot(4,2,5)
            hold on       
            plot([minH maxH], [minH maxH], '-k')
            plot(H_in(Hqc_best_inds), H_eco_store(Hqc_best_inds), '.k', 'MarkerSize', 8)
            plot([minH, maxH], pfH(1)*[minH, maxH]+pfH(2), '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
            %plot([minH, maxH], pfH_day(1)*[minH, maxH]+pfH_day(2), '--r', 'LineWidth', 2)
            axis([minH maxLE minH maxLE])
            xlabel('\itH \rm(Obs) [W m^{-2}]', 'FontSize', 12)
            ylabel('\itH \rm(Mod) [W m^{-2}]', 'FontSize', 12)
            text(-150, 500, ['slope = ', num2str(pfH(1))], 'FontSize', 12)
            text(-150, 400, ['int = ', num2str(pfH(2))], 'FontSize', 12)
            text(-150, 300, ['R^2 = ', num2str(RS_H)], 'FontSize', 12)
            box on

        subplot(4,2,7)
            hold on       
            plot([minRn maxRn], [minRn maxRn], '-k')
            plot(Rn_in(Rn_good_inds), Rnvar(Rn_good_inds), '.k', 'MarkerSize', 8)
            plot([minRn, maxRn], pfRn(1)*[minRn, maxRn]+pfRn(2), '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5])
            %plot([minRn, maxRn], pfRn_day(1)*[minRn, maxRn]+pfRn_day(2), '--r', 'LineWidth', 2)
            axis([minRn maxRn minRn maxRn])
            xlabel('\itR_n \rm(Obs) [W m^{-2}]', 'FontSize', 12)
            ylabel('\itR_n \rm(Mod) [W m^{-2}]', 'FontSize', 12)
            text(-150, 700, ['slope = ', num2str(pfRn(1))], 'FontSize', 12)
            text(-150, 570, ['int = ', num2str(pfRn(2))], 'FontSize', 12)
            text(-150, 440, ['R^2 = ', num2str(RS_Rn)], 'FontSize', 12)
            box on
        
            
    % FLUX DIURNALS        
        subplot(4,2,2)
            hold on
            plot(hr, zeros(size(hr)), ':k')
            errorbar(hr, Fc_diurnal_obs, Fc_std_obs, '-ok', 'MarkerFaceColor', 'k');
            errorbar(hr, Fc_diurnal_eco, Fc_std_eco, '-sr');
            axis([0 23.5 -50 15])
            box on
            xlabel('Hour', 'FontSize',12)
            ylabel('\itF_c \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
            
        subplot(4,2,4)
            hold on
            plot(hr, zeros(size(hr)), ':k')
            errorbar(hr, LE_diurnal_obs, LE_std_obs, '-ok', 'MarkerFaceColor', 'k')
            errorbar(hr, LE_diurnal_eco, LE_std_eco, '-sr')   
            axis([0 23.5 -100 450])
            box on
            xlabel('Hour', 'FontSize',12)
            ylabel('\itLE \rm[W m^{-2}]', 'FontSize', 12)
            
        subplot(4,2,6)
            hold on
            plot(hr, zeros(size(hr)), ':k')
            p1 = errorbar(hr, H_diurnal_obs, H_std_obs, '-ok', 'MarkerFaceColor', 'k');
            p2 = errorbar(hr, H_diurnal_eco, H_std_eco, '-sr');
            axis([0 23.5 -100 450])
            box on
            legend([p1,p2], 'Observed', 'Modeled')
            xlabel('Hour', 'FontSize',12)
            ylabel('\itH \rm[W m^{-2}]', 'FontSize', 12)
            
        subplot(4,2,8)
            hold on
            plot(hr, zeros(size(hr)), ':k')
            errorbar(hr, Rn_diurnal_obs, Rn_std_obs, '-ok', 'MarkerFaceColor', 'k')
            errorbar(hr, Rn_diurnal_eco, Rn_std_eco, '-sr')
            axis([0 23.5 -100 800])
            box on
            xlabel('Hour', 'FontSize',12)
            ylabel('\itR_n \rm[W m^{-2}]', 'FontSize', 12)
            
            