function [LCH_CC, TLCH_CC] = CN_leach(PARAMS, SWITCHES, CC, aa, qq, sm)

    % CC Concentration variable [gr/m3]
    % aa solubility of element in study
    if SWITCHES.CN.Bioturbation;
        nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
    else
        nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
    end
    
    leach_on = PARAMS.CN.leach_on;%       leach_on = Leaching ON or OF [1 0]
    
    % Declare matrix
    LCH_CC = zeros (nl_soil,1);
    
    for ii = 1:nl_soil    
        %   positive value = loss (subtracted from budget)
        %   NOTE: first H2O flux value is from infiltration into soil column
            LCH_CC(ii) = 0;
            if (leach_on)
                if (ii>1)   % transport through top layer interface
                    if (qq(ii)>0)  % downward transport into cell
                        LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii-1)*CC(ii-1)*qq(ii));% [gr/m^2/d]
                    else           % upward transport out of cell
                        LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii)*CC(ii)*qq(ii));% [gr/m^2/d]
                    end
                end
                % transport through bottom layer interface
                if (qq(ii+1)>0)  % downward transport out of cell
                    LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
                else           % upward transport into cell
                    LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii+1)*CC(ii+1)*qq(ii+1));% [gr/m^2/d]
                end
                if (ii==nl_soil)
                    TLCH_CC = (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
                end   
            end                        
    end
    
    TLCH_CC = TLCH_CC(:);
    LCH_CC = LCH_CC(:);
    