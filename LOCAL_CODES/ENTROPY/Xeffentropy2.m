% % This function is used to compute the Effective Temperature
% % from the Canopy Using the Entropy Calculated at each 
% % Level
% 
% function [remain] = Teffentropy (SSeco_tot, SSsoil_in, SScan_in, SSsoildif_out_tot, SScandif_out_tot, ...
%                 LWemi_net, LE_net, H_net, SSnewXout, Teffe)
% 
%     remain = - SSeco_tot + SSsoil_in + SScan_in ...
%         - (SSsoildif_out_tot + SScandif_out_tot + (4/3*SSnewXout)*LWemi_net/Teffe + ...
%                  LE_net/Teffe + H_net/Teffe);
%              
             
             
function [remain] = Xeffentropy2 (SSdirin, SSdifin, SSdirout, SSdifout, SSLWin,... 
SSecotot, LWout, LEnet, Hnet, X, Tefe)

% This function is used to compute the Effective Temperature
% from the Canopy Using the Entropy Calculated at each 
% Level


    remain = - SSecotot + SSdifin + SSdirin + SSLWin...
        - (SSdirout + SSdifout + (4/3)*X*LWout/Tefe + LEnet/Tefe + Hnet/Tefe);
                 