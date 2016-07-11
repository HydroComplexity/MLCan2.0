

% LEAF POTENTIAL AND HYDRAULIC STOMATAL CONTROL

canind = 13;

    minpsi = min( min(min(psil_sun_prof)), min(min(psil_shade_prof)) );  
    maxpsi = max( max(max(psil_sun_prof)), max(max(psil_shade_prof)) );  
    figure(fignum); clf
        subplot(6,1,1)
            pcolor(timevect, zhc, psil_canopy_prof)
            colormap(flipdim(jet,1))
            shading interp
            colorbar
            caxis([minpsi 0])
            ylabel('\bf z / h', 'FontSize', 12)
            title('\bf \Psi_l Canopy', 'FontSize', 12)
        subplot(6,1,2)
            pcolor(timevect, zhc, fsv_canopy_prof)
            colormap(flipdim(jet,1))
            caxis([0 1])
            shading interp
            colorbar
            title('\bf fsv Canopy', 'FontSize', 12)
            ylabel('\bf z / h', 'FontSize', 12)
        subplot(6,1,3)    
            pcolor(timevect, zhc, TR_canopy_prof)
            colormap(flipdim(jet,1))
            shading interp
            colorbar
            title('\bf Canopy Transpiration (TR)', 'FontSize', 12)
            ylabel('\bf z / h', 'FontSize', 12)
        subplot(6,1,4)    
            pcolor(timevect, zhc, TR_canopy_prof.*PARAMS.StomCond.Rp)
            colormap(flipdim(jet,1))
            shading interp
            colorbar
            title('\bf Transpiration Factor in \Psi_l Calculation', 'FontSize', 12)
            ylabel('\bf z / h', 'FontSize', 12)    
        subplot(6,1,5)
            hold on
            plot(timevect, smp_weight_store, '-k', 'LineWidth', 2)
            plot(timevect, rpp_weight_store, '-g', 'LineWidth', 2)
            plot(timevect, psil_sun_prof(canind,:), '-r', 'LineWidth', 1)
            plot(timevect, PARAMS.StomCond.psif.*ones(size(timevect)), '-', 'Color', [0.5 0.5 0.5])
            ylabel('\Psi [MPa]')
            box on
            axis([timevect(1) timevect(end) -1 0])
            h = legend('\Psi_s', '\Psi_r', '\Psi_l (top)');
            set(h, 'FontSize', 14)
        subplot(6,1,6)
            hold on
            plot(timevect, LE_in, '.-k', 'LineWidth', 2)
            plot(timevect, LE_can_store, '-b')
            plot(timevect, H_in, '.-r', 'LineWidth', 2)
            plot(timevect, H_can_store, '-g')
            xlabel('\bf DOY', 'FontSize', 12)
            box on 
            axis tight