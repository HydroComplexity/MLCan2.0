
    [Rgupin_diurnal, Rgupin_std] = DIURNAL_AVERAGE (hr, hour, Rgup_in);
    [SWout_diurnal, SWout_std] = DIURNAL_AVERAGE (hr, hour, SWout_store);
    
    [LWupin_diurnal, LWupin_std] = DIURNAL_AVERAGE (hr, hour, LWup_in);
    [LWout_diurnal, LWout_std] = DIURNAL_AVERAGE (hr, hour, LWout_store);

% Difference Diurnals
    SWrefldiff = SWout_store - Rgup_in;
    [SWrefldiff_diurnal, SWrefldiff_std] = DIURNAL_AVERAGE (hr, hour, SWrefldiff);

    LWoutdiff = LWout_store - LWup_in;
    [LWoutdiff_diurnal, LWoutdiff_std] = DIURNAL_AVERAGE (hr, hour, LWoutdiff);


figure(fignum); clf
    subplot(4,2,1)
        hold on
        plot(timevect, Rgup_in, '.k')
        plot(timevect, SWout_store, 'r')
        box on 
        grid on
        ylabel('Rg Out [W m^{-2}]')
        xlabel('')
        title('SHORTWAVE')
    subplot(4,2,3)
        hold on
        plot(Rgup_in, SWout_store, '.k', 'MarkerSize', 1)
        plot([0 300], [0 300], 'Color', [0.5 0.5 0.5])
        box on 
        ylabel('Rg Out (MOD) [W m^{-2}]')
        xlabel('Rg Out (OBS) [W m^{-2}]')
    subplot(4,2,5)
        hold on
        errorbar(hr, Rgupin_diurnal, Rgupin_std, '-ob')
        errorbar(hr, SWout_diurnal, SWout_std, '.-r')
        axis([0 24 -Inf Inf])
        grid on
        legend('Obs', 'Mod')
        ylabel('SW_{out}(mod)-SW_{out}(obs)')
        xlabel('Hour')
        box on    
    subplot(4,2,7)
        hold on
        errorbar(hr, SWrefldiff_diurnal, SWrefldiff_std, '-ok')
        axis([0 24 -Inf Inf])
        grid on
        ylabel('SW_{out}(mod)-SW_{out}(obs)')
        xlabel('Hour')
        box on    
        
    % LONGWAVE    
    subplot(4,2,2)
        hold on
        plot(timevect, LWup_in, '.k')
        plot(timevect, LWout_store, 'r')
        legend('Obs', 'Mod')
        box on 
        grid on
        ylabel('LW out [W m^{-2}]')
        title('LONGWAVE')
    subplot(4,2,4)
        hold on
        plot(LWup_in, LWout_store, '.k', 'MarkerSize', 1)
        plot([300 600], [300 600], 'Color', [0.5 0.5 0.5])
        box on 
        ylabel('LW Out (MOD) [W m^{-2}]')
        xlabel('LW Out (OBS) [W m^{-2}]')
    subplot(4,2,6)
        hold on
        errorbar(hr, LWupin_diurnal, LWupin_std, '-ob')
        errorbar(hr, LWout_diurnal, LWout_std, '.-r')
        axis([0 24 -Inf Inf])
        grid on
        legend('Obs', 'Mod')
        ylabel('LW_{out}(mod)-LW_{out}(obs)')
        xlabel('Hour')
        box on            
    subplot(4,2,8)
        hold on
        errorbar(hr, LWoutdiff_diurnal, LWoutdiff_std, '-ob')
        axis([0 24 -Inf Inf])
        grid on
        ylabel('LW_{out}(mod)-LW_{out}(obs)')
        xlabel('Hour')
        box on    
    
        
        
        