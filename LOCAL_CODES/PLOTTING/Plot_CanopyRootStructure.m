
% Compute conductance reduction function
    sf = PARAMS.StomCond.sf;
    psif = PARAMS.StomCond.psif;
    psil_all=[]; fsv_all=[];
    psil_all = [-3:0.1:0];  % [MPa]
    for ii = 1:length(psil_all)
        psil_MPa = psil_all(ii);

        fsv_all(ii) = (1 + exp(sf*psif)) ./ (1 + exp(sf*(psif - psil_MPa)));
    end

% Load Leaf Water Potential Data from Bunce
%    psil_bunce=[]; fgs_bunce=[];
%    if (soy_simulation)
%        X = load(['../../Data_Files/LeafWaterPotential_Bunce/Soy_Ambient_Points.rtf']);
%    else % corn
%        X = load(['../../Data_Files/LeafWaterPotential_Bunce/Maize_Ambient_Points.rtf']);
%    end
%    psil_bunce = X(:,1); % [MPa]
%    fgs_bunce = X(:,2);  % [fraction]
    

% CANOPY AND ROOT STRUCTURE
    rootfr_norm = rootfr./dzsmm;
    rootfr_norm = rootfr_norm./(sum(rootfr_norm));
    figure(fignum); clf
        subplot(3,1,1)
            plot(LADnorm, znc./PARAMS.CanStruc.hcan, '.-')
            ylabel('z / h [-]')
            xlabel('Normalized LAD Profile')
            axis([0 Inf 0 1])
            box on
        subplot(3,1,2)
            plot(rootfr_norm, -zns, '.-')
            ylabel('z [m]')
            xlabel('Root Fraction')
            box on
        subplot(3,1,3)
            hold on
%            plot(psil_bunce, fgs_bunce, '.r')
            plot(psil_all, fsv_all, '-k', 'LineWidth', 2)
            xlabel('\Psi_l [MPa]', 'FontSize', 12)
            ylabel('Fraction Conductivity', 'FontSize', 12)
            box on
            axis([-2 0 0 1])