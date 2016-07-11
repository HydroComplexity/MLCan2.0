function [CCnew] = CN_updatecon(VERTSTRUC, SWITCHES, VARIABLES, PARAMS, RIVDf, CC, CCdd)

% Derefence blocks
if SWITCHES.CN.Bioturbation;
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
    dz = [VARIABLES.SOIL.litterthickness ; VERTSTRUC.dzsmm/1000];%       dz = grid depth [m]
    deltaCC = zeros(nl_soil,1);
else
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
    dz = VERTSTRUC.dzsmm/1000; % dz = grid depth [m]
    deltaCC = zeros(nl_soil,1);
end
RIVDi = VARIABLES.RVID;

% Compute delta volume in zone of disturbance
deltaRIVD = sum(RIVDf,2) - sum(RIVDi,2);

inddie = deltaRIVD < 0;
indgro = deltaRIVD > 0;
% Roots die, ammount in zone of disturbance returns to soil
deltaCC(inddie) = - (deltaRIVD(inddie).*CCdd(inddie)./dz(inddie));  % [gr/m3] is positive
% Roots 
deltaCC(indgro) = - (dz(indgro).*CCdd(indgro)./deltaRIVD(indgro));  % [gr/m3] is negative

CCnew = CC + deltaCC;

% check mass balance 

%CC.*dz + CCdd.*sum(RIVDi,2) - (CCnew.*dz + CCdd.*sum(RIVDf,2)); 