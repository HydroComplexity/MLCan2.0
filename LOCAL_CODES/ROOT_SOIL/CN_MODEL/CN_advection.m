function [ADV_m2, ADV_m3, TLCADV_m2, TLCADV_m3] = CN_advection(nl_soil, dz, CC, aa, qq, sm, dtime)

    % CC Concentration variable [gr/m3]
    % aa solubility of element in study
        
    % Declare matrix
    Fo_adv = zeros (nl_soil,2);
    ADV_m2 = zeros (nl_soil,1);
     
    % Define Fo_adv which is a [nlayer X 2] matrix with out fluxes from each
    % layer. 1st column is upward flow while 2nd column is downward
    for ii = 1:nl_soil-1
        if (qq(ii+1)<0)             % upward transport out of cell
            Fo_adv(ii+1,1) = -(aa/sm(ii)*CC(ii+1)*qq(ii+1));% [gr/m^2/d]             
        end
        if (qq(ii+1)>0)           % downward transport out of cell
            Fo_adv(ii,2) = (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
        end
    end
    % Compute for the last layer
    if (qq(nl_soil+1)>0)
        Fo_adv(nl_soil,2) = (aa/sm(nl_soil)*CC(nl_soil)*qq(nl_soil+1));% [gr/m^2/d]
    end
    % It assumes that for those cases where there is an upward flow of water 
    % at the bottom there is not flux of ions coming into the control volume
    
    % Compute total fluxex and check there is enough mass in each layer to
    % satisfy that flux
    out_adv = sum(Fo_adv,2)./dz * dtime/86400;    % Total out due to advection in [gr/m^3/dtime]
    ind_mb = CC < out_adv;        % Check mass availability in that layer in [gr/m^3]
    F_adv = sum(Fo_adv,2);       % [gr/m^2/d]   
    F_adv(ind_mb) = 0;           % Those layers with no available are set to 0 flux  
    Fo_adv(ind_mb,:) = 0;     
    
    % Using the output fluxes matrix Fo_adv now we compute the flux into each layer 
    ADV_m2 = ADV_m2 - F_adv;% Substract the fluxes out from each layer first   
    % Add the fluxes from each layer from the 2 to the nlsoil -1
    for ii = 2:nl_soil-1                
        ADV_m2(ii-1) = ADV_m2(ii-1) + Fo_adv(ii,1);  %[g/m^2/d]
        ADV_m2(ii+1) = ADV_m2(ii+1) + Fo_adv(ii,2);  %[g/m^2/d] 
    end
    % Now compute for the top  
    ADV_m2(2) = ADV_m2(2) + Fo_adv(1,2);             %[g/m^2/d]
    % and botom
    ADV_m2(nl_soil-1) = ADV_m2(nl_soil-1) + Fo_adv(nl_soil,1); %[gr/m^2/d] 
    TLCADV_m2 = Fo_adv(nl_soil,2);           % Total leaching [gr/m^2/d]
    %
    %********************** WE CAN CHECK MASS BALANCE. COMPUTING THIS *****
    %sum(ADV_m2) + TLCADV_m2   % =0 ?????????????????????
    %*********************************************************************
    
    ADV_m3 = ADV_m2./dz;           % Compute the differences in concentration due to advection [g/m^3/d]
    TLCADV_m3 = TLCADV_m2/dz(nl_soil);           % Total leaching in [g/m^3/d]
    