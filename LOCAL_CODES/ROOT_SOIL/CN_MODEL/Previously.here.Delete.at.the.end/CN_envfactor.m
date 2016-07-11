%*************************************************************************

%*************************************************************************
%                   COMPUTE ENVIRONMENTAL FACTORS
%*************************************************************************
    
    % for Decomposition
        % Temp
          %   fTd = 0.9*(Ts./(Ts+exp(9.93-0.312*Ts)))+0.05;           
          %   fTd = exp(3.36 * (Ts - 40) ./ (Ts + 31.79));
          %   Use Q10 formulation
              indT = Ts > 25;
              fTd(indT) = 1;              
              fTd(~indT) = 1.9.^((Ts(~indT)-25)/10); 
              fTd = fTd(:);
          %  fTd = ones(nl_soil,1);
        % Moisture
          %  fSd = 1./(1+4*exp(-6*rwc));
            % using equations from porporato et al 2003
%              findn=rwc<fc;
%              findy=rwc>=fc;
%              fSd(findn)=rwc(findn)./fc;
%              fSd(findy)=fc./rwc(findy);
%              fSd=fSd(:);
            %  Use water potential in MPa 
            smp = smp * 9.8066e-06; % First convert to MPa
            fSd = log(7.58./abs(smp))/log(7.58/0.01);
            indWP = fSd < 0;
            fSd(indWP) = 0; 
            indWP = fSd > 1;
            fSd(indWP) = 1; 
            
    % for Nitrification
        fTn = fTd;
        fSn = fSd;
    % for Volatilization
        %fTv = fTd;
    % for Denitrification
        % Temp
            %fTdn = fTd;
        % Moisture
%             finds = sm./porsl <= 0.8;
%             fSdn(finds) = 0;
%             finds = (sm./porsl>0.8) & (sm./porsl<=0.9);
%             fSdn(finds) = -1.6 + 2*sm(finds)./porsl(finds);
%             finds = sm./porsl > 0.9;
%             fSdn(finds) = -7 + 8*sm(finds)./porsl(finds);
