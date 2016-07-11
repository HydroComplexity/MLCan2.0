function [VERTSTRUC VARIABLES rootfr] = ROOT_RESPONSE_DRY(VARIABLES,...
    SWITCHES, VERTSTRUC, CONSTANTS, PARAMS, doy, smp_store) 

% SAVE 
tt = VARIABLES.timestep;
nspecies = PARAMS.CanStruc.nspecies;
nl_soil = VERTSTRUC.nl_soil;
rootfr = VERTSTRUC.rootfr;
rpp = VARIABLES.ROOT.rpp;
nlc = CONSTANTS.nlc; % Number of layer to cut

if (tt> 1 && SWITCHES.cutroots ~=0 )
      if  SWITCHES.cutroots == 1
          rootfr(1:nlc,:) = 0;
      elseif and(doy(tt) ~= doy(tt-1),SWITCHES.cutroots == 2)                
          inddoy = doy==doy(tt-1);
          vecsmp1 = mean(abs(smp_store(:,inddoy))*9.8066e-06,2);
          indsmp = vecsmp1 >= abs(CONSTANTS.wilpoint);   
          for i=1:1:nspecies
              VERTSTRUC.rootfr(indsmp) = 0;
          end
      elseif SWITCHES.cutroots == 3
          VARIABLES.ROOT.PLC = 1./(1+exp(CONSTANTS.PLCa*((rpp*9.8066e-06)-CONSTANTS.PLCb)));  
          indplc = rpp*9.8066e-06 >= -0.05;
          VARIABLES.ROOT.PLC(indplc) = 0;                
      elseif SWITCHES.cutroots == 4
          % First cut roots in first layer
          rootfr(1:nlc,:) = 0; 
          % Second cut the roots based on water potential 
          vecsmp1 = abs(smp_store(:,tt-1))*9.8066e-06;
          indsmp = vecsmp1 >= abs(CONSTANTS.wilpoint);   
          for i=1:1:nspecies
              rootfr(indsmp) = 0;
          end                
          % Third the PLC
          VARIABLES.ROOT.PLC = 1./(1+exp(CONSTANTS.PLCa*((rpp*9.8066e-06)-CONSTANTS.PLCb)));  
          indplc = rpp*9.8066e-06 >= -0.05;
          VARIABLES.ROOT.PLC(indplc) = 0;
      end
      if ~VARIABLES.comroot     
          % Fixing the constant soil layer problem
%          rootfr = rootfr./repmat(sum(rootfr(nlc+1:12,:)),12,1);
          rootfr = rootfr./repmat(sum(rootfr(nlc+1:nl_soil,:)),nl_soil,1);
          VERTSTRUC.rootfr = rootfr; 
          VARIABLES.comroot = 1;              
      end    
else
    VARIABLES.ROOT.PLC = zeros(nl_soil,nspecies); 
    rootfr = VERTSTRUC.rootfr;
end        

