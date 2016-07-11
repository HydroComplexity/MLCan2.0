

figure(fignum); clf        
    hold on
    plot(timevect, smp_weight_store, '-k', 'LineWidth', 2)
    plot(timevect, rpp_weight_store, '-g', 'LineWidth', 2)
    plot(timevect, psil_sun_prof(nl_can,:), '-r', 'LineWidth', 1)
    plot(timevect, psil_shade_prof(nl_can,:), '-b', 'LineWidth', 1)
%    plot(timevect, psil_sun_prof(1,:), ':r', 'LineWidth', 1)
%    plot(timevect, psil_shade_prof(1,:), ':b', 'LineWidth', 1)
    ylabel('\Psi [MPa]')
    box on
    legend('soil (weighted)', 'root (weighted)', 'sun leaf (top)', 'shade leaf (top)')%, 'sun leaf (bottom)', 'shade leaf (bottom)')
