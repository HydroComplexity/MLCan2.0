
function [SSresults] = ENTROPY_net_results(SSresults)


% Compute the output of entropy and save in different variables


% SHORTWAVE
SScandir_net_in = SSresults.SSnetdir_in;
SScandir_net_out = SSresults.SSnetdir_out;
SScandif_net_in = SSresults.SSnetdif_in;
SScandif_net_out = SSresults.SSnetdif_out;


% LONGWAVE
SScanLW_net_in = SSresults.SSnetLW_in;
SScanLW_net_out = SSresults.SSnetLW_out;

% HEAT
SSLE_can_net = SSresults.SSLE_netLEcan;
SSH_can_net = SSresults.SSH_netLEcan;

