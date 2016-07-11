
% Quality Control
    Fcqc_ginds = find(Fc_qc==1);
    LEqc_ginds = find(LE_qc==1);
    Hqc_ginds = find(H_qc==1);


dotsize = 6;    % Size of observed obs (dots)
msize = 8;      % Size of binned symbols

dclr1=1; dclr2=0.6; dclr3=0.6;      % Daytime observation dot color
nclr1=0.3; nclr2=0.3; nclr3=0.3;    % Nocturnal observation dot color

mclrbw1 = 0.4; mclrbw2 = 0.4; mclrbw3 = 0.4; % Model bin color for B&W plot

dclrbw1=0.8; dclrbw2=.8; dclrbw3=0.8;   % Daytime observation dot color for B&W plot


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
    
    % BINNED OBSERVATIONS - DAY
        [Fc_in_Rgbins_avg, Fc_in_Rgbins_std, Rgbin_avg] = Bin_Average (Rg_in(Fcqc_ginds), Fc_in(Fcqc_ginds), Rg_bins);
        [Fc_in_Tabins_avg, Fc_in_Tabins_std, Tabin_avg] = Bin_Average (Ta_in(Fcqc_ginds), Fc_in(Fcqc_ginds), Ta_bins);
        [Fc_in_VPDbins_avg, Fc_in_VPDbins_std, VPDbin_avg] = Bin_Average (VPD_in(Fcqc_ginds), Fc_in(Fcqc_ginds), VPD_bins);
        [Fc_in_Ubins_avg, Fc_in_Ubins_std, Ubin_avg] = Bin_Average (U_in(Fcqc_ginds), Fc_in(Fcqc_ginds), U_bins);

        [LE_in_Rgbins_avg, LE_in_Rgbins_std, Rgbin_avg] = Bin_Average (Rg_in(LEqc_ginds), LE_in(LEqc_ginds), Rg_bins);
        [LE_in_Tabins_avg, LE_in_Tabins_std, Tabin_avg] = Bin_Average (Ta_in(LEqc_ginds), LE_in(LEqc_ginds), Ta_bins);
        [LE_in_VPDbins_avg, LE_in_VPDbins_std, VPDbin_avg] = Bin_Average (VPD_in(LEqc_ginds), LE_in(LEqc_ginds), VPD_bins);
        [LE_in_Ubins_avg, LE_in_Ubins_std, Ubin_avg] = Bin_Average (U_in(LEqc_ginds), LE_in(LEqc_ginds), U_bins);

        [H_in_Rgbins_avg, H_in_Rgbins_std, Rgbin_avg] = Bin_Average (Rg_in(Hqc_ginds), H_in(Hqc_ginds), Rg_bins);
        [H_in_Tabins_avg, H_in_Tabins_std, Tabin_avg] = Bin_Average (Ta_in(Hqc_ginds), H_in(Hqc_ginds), Ta_bins);
        [H_in_VPDbins_avg, H_in_VPDbins_std, VPDbin_avg] = Bin_Average (VPD_in(Hqc_ginds), H_in(Hqc_ginds), VPD_bins);
        [H_in_Ubins_avg, H_in_Ubins_std, Ubin_avg] = Bin_Average (U_in(Hqc_ginds), H_in(Hqc_ginds), U_bins);

    % BINNED MODELED FLUXES - DAY
        [Fc_mod_Rgbins_avg, Fc_mod_Rgbins_std] = Bin_Average (Rg_in(Fcqc_ginds), Fc_eco_store(Fcqc_ginds), Rg_bins);
        [Fc_mod_Tabins_avg, Fc_mod_Tabins_std] = Bin_Average (Ta_in(Fcqc_ginds), Fc_eco_store(Fcqc_ginds), Ta_bins);
        [Fc_mod_VPDbins_avg, Fc_mod_VPDbins_std] = Bin_Average (VPD_in(Fcqc_ginds), Fc_eco_store(Fcqc_ginds), VPD_bins);
        [Fc_mod_Ubins_avg, Fc_mod_Ubins_std] = Bin_Average (U_in(Fcqc_ginds), Fc_eco_store(Fcqc_ginds), U_bins);

        [LE_mod_Rgbins_avg, LE_mod_Rgbins_std] = Bin_Average (Rg_in(LEqc_ginds), LE_eco_store(LEqc_ginds), Rg_bins);
        [LE_mod_Tabins_avg, LE_mod_Tabins_std] = Bin_Average (Ta_in(LEqc_ginds), LE_eco_store(LEqc_ginds), Ta_bins);
        [LE_mod_VPDbins_avg, LE_mod_VPDbins_std] = Bin_Average (VPD_in(LEqc_ginds), LE_eco_store(LEqc_ginds), VPD_bins);
        [LE_mod_Ubins_avg, LE_mod_Ubins_std] = Bin_Average (U_in(LEqc_ginds), LE_eco_store(LEqc_ginds), U_bins);

        [H_mod_Rgbins_avg, H_mod_Rgbins_std] = Bin_Average (Rg_in(Hqc_ginds), H_eco_store(Hqc_ginds), Rg_bins);
        [H_mod_Tabins_avg, H_mod_Tabins_std] = Bin_Average (Ta_in(Hqc_ginds), H_eco_store(Hqc_ginds), Ta_bins);
        [H_mod_VPDbins_avg, H_mod_VPDbins_std] = Bin_Average (VPD_in(Hqc_ginds), H_eco_store(Hqc_ginds), VPD_bins);
        [H_mod_Ubins_avg, H_mod_Ubins_std] = Bin_Average (U_in(Hqc_ginds), H_eco_store(Hqc_ginds), U_bins);
        
        

%*************
%   COLOR
%*************
figure(fignum); clf
    subplot(3,4,1)
        hold on
        plot(Rg_in(dinds), Fc_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Rg_in(ninds), Fc_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Rg_in, Fc_eco_store, '.r', 'MarkerSize', dotsize)
        p1 = errorbar(Rgbin_avg, Fc_in_Rgbins_avg, Fc_in_Rgbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2);
        p2 = errorbar(Rgbin_avg, Fc_mod_Rgbins_avg, Fc_mod_Rgbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize);
        ylabel('\itF_c \rm[\mumol m^{-2} s^{-1}]', 'FontSize', 12)
        box on
        axis([Rgmin Rgmax Fcmin Fcmax])
        set(gca, 'XTickLabel', [])
        legend([p1, p2], 'Observed', 'Modeled')
    subplot(3,4,5)
        hold on
        plot(Rg_in(dinds), LE_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Rg_in(ninds), LE_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Rg_in, LE_eco_store, '.r', 'MarkerSize', dotsize)        
        errorbar(Rgbin_avg, LE_in_Rgbins_avg, LE_in_Rgbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Rgbin_avg, LE_mod_Rgbins_avg, LE_mod_Rgbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)     
        ylabel('\itLE \rm[W m^{-2}]', 'FontSize', 12)
        box on
        axis([Rgmin Rgmax LEmin LEmax])
        set(gca, 'XTickLabel', [])
        
    subplot(3,4,9)
        hold on
        plot(Rg_in(dinds), H_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Rg_in(ninds), H_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Rg_in, H_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(Rgbin_avg, H_in_Rgbins_avg, H_in_Rgbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Rgbin_avg, H_mod_Rgbins_avg, H_mod_Rgbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        ylabel('\itH \rm[W m^{-2}]', 'FontSize', 12)
        xlabel('\itR_g \rm[W m^{-2}]', 'FontSize', 12)
        box on
        axis([Rgmin Rgmax Hmin Hmax])
        
% VARIATION WITH Ta
    subplot(3,4,2)
        hold on
        plot(Ta_in(dinds), Fc_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Ta_in(ninds), Fc_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Ta_in, Fc_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(Tabin_avg, Fc_in_Tabins_avg, Fc_in_Tabins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Tabin_avg, Fc_mod_Tabins_avg, Fc_mod_Tabins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        box on
        axis([Tamin Tamax Fcmin Fcmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,6)
        hold on
        plot(Ta_in(dinds), LE_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Ta_in(ninds), LE_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Ta_in, LE_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(Tabin_avg, LE_in_Tabins_avg, LE_in_Tabins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Tabin_avg, LE_mod_Tabins_avg, LE_mod_Tabins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        box on
        axis([Tamin Tamax LEmin LEmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,10)
        hold on
        plot(Ta_in(dinds), H_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(Ta_in(ninds), H_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(Ta_in, H_eco_store, '.r', 'MarkerSize', dotsize)    
        errorbar(Tabin_avg, H_in_Tabins_avg, H_in_Tabins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Tabin_avg, H_mod_Tabins_avg, H_mod_Tabins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        xlabel('\itT_a \rm[\circC]', 'FontSize', 12)
        box on
        axis([Tamin Tamax Hmin Hmax])
        set(gca, 'YTickLabel', [])
        
% VARIATION WITH VPD
    subplot(3,4,3)
        hold on
        plot(VPD_in(dinds), Fc_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(VPD_in(ninds), Fc_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(VPD_in, Fc_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(VPDbin_avg, Fc_in_VPDbins_avg, Fc_in_VPDbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(VPDbin_avg, Fc_mod_VPDbins_avg, Fc_mod_VPDbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)
        box on
        axis([VPDmin VPDmax Fcmin Fcmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,7)
        hold on
        plot(VPD_in(dinds), LE_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(VPD_in(ninds), LE_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(VPD_in, LE_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(VPDbin_avg, LE_in_VPDbins_avg, LE_in_VPDbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(VPDbin_avg, LE_mod_VPDbins_avg, LE_mod_VPDbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        box on
        axis([VPDmin VPDmax LEmin LEmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,11)
        hold on
        plot(VPD_in(dinds), H_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(VPD_in(ninds), H_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(VPD_in, H_eco_store, '.r', 'MarkerSize', dotsize)   
        errorbar(VPDbin_avg, H_in_VPDbins_avg, H_in_VPDbins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(VPDbin_avg, H_mod_VPDbins_avg, H_mod_VPDbins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)
        xlabel('\itVPD \rm[kPa]', 'FontSize', 12)
        box on
        axis([VPDmin VPDmax Hmin Hmax])
        set(gca, 'YTickLabel', [])
        
% VARIATION WITH U
    subplot(3,4,4)
        hold on
        plot(U_in(dinds), Fc_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(U_in(ninds), Fc_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(VPD_in, Fc_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(Ubin_avg, Fc_in_Ubins_avg, Fc_in_Ubins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Ubin_avg, Fc_mod_Ubins_avg, Fc_mod_Ubins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)
        box on
        axis([Umin Umax Fcmin Fcmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,8)
        hold on
        plot(U_in(dinds), LE_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(U_in(ninds), LE_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(U_in, LE_eco_store, '.r', 'MarkerSize', dotsize)
        errorbar(Ubin_avg, LE_in_Ubins_avg, LE_in_Ubins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Ubin_avg, LE_mod_Ubins_avg, LE_mod_Ubins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)  
        box on
        axis([Umin Umax LEmin LEmax])
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        
    subplot(3,4,12)
        hold on
        plot(U_in(dinds), H_in(dinds), '.', 'MarkerSize', dotsize, 'Color', [dclr1 dclr2 dclr3])
        plot(U_in(ninds), H_in(ninds), '.', 'MarkerSize', dotsize, 'Color', [nclr1 nclr2 nclr3])
        
        %plot(U_in, H_eco_store, '.r', 'MarkerSize', dotsize)   
        errorbar(Ubin_avg, H_in_Ubins_avg, H_in_Ubins_std, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        errorbar(Ubin_avg, H_mod_Ubins_avg, H_mod_Ubins_std, '-db', 'LineWidth', 2, 'MarkerSize', msize)
        xlabel('\itU \rm[m s^{-1}]', 'FontSize', 12)
        box on
        axis([Umin Umax Hmin Hmax])
        set(gca, 'YTickLabel', [])
 
RePosition_3x4;        