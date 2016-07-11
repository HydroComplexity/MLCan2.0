

% fsv - reduction function for stomatal conductance
figure(fignum); clf
    subplot(3,1,1)
        pcolor(timevect, zhc, fsv_sun_prof)
        colormap(flipdim(jet,1))
        caxis([0 1])
        shading interp
        colorbar  
        title('\bf fsv Sun', 'FontSize', 12)
        ylabel('\bf z / h', 'FontSize', 12)
    subplot(3,1,2)
        pcolor(timevect, zhc, fsv_shade_prof)
        colormap(flipdim(jet,1))
        caxis([0 1])
        shading interp
        colorbar
        title('\bf fsv Shade', 'FontSize', 12)
        ylabel('\bf z / h', 'FontSize', 12)
    subplot(3,1,3)
        pcolor(timevect, zhc, fsv_canopy_prof)
        colormap(flipdim(jet,1))
        caxis([0 1])
        shading interp
        colorbar
        title('\bf fsv Canopy (Sunlit,Shaded Fraction Weighted)', 'FontSize', 12)
        ylabel('\bf z / h', 'FontSize', 12)
        xlabel('\bf DOY', 'FontSize', 12)