% This function is used to compute the Effective X
% from the Canopy Using the Entropy Calculated at each 
% Level

function [remain] = Xeffentropy(SScanLW_out_tot, SSsoilLW_out_tot, ...
                    SSnetLW_in, SScanLW_in_tot, SSsoilLW_in_tot,...            
                    LWout_net, Teffent, Xnet)

         
 remain = (SScanLW_in_tot + SSsoilLW_in_tot - SScanLW_out_tot - SSsoilLW_out_tot) -...
     (SSnetLW_in - (4/3)*Xnet*LWout_net/Teffent);
 
 
 
 