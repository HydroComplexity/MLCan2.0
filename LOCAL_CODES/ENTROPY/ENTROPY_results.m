
function [SSresults] = ENTROPY_results(SSresults, VERTSTRUC, PARAMS, VARIABLES)


% Compute the output of entropy and save in different variables


fLAIz = VERTSTRUC.fLAIz;
nspecies = PARAMS.CanStruc.nspecies;

% SHORTWAVE
SScandir_in = SSresults.SScandir_in;
SScandir_in_tot = SSresults.SScandir_in_tot;
SScandir_out = SSresults.SScandir_out;
SScandir_out_tot = SSresults.SScandir_out_tot;
SScandif_in = SSresults.SScandif_in;
SScandif_in_tot = SSresults.SScandif_in_tot;
SScandif_out = SSresults.SScandif_out;
SScandif_out_tot = SSresults.SScandif_out_tot;


SSsoildir_in = SSresults.SSsoildir_in;
SSsoildir_out = SSresults.SSsoildir_out;
SSsoildif_in = SSresults.SSsoildif_in;
SSsoildif_out = SSresults.SSsoildif_out;


% LONGWAVE
SScanLW_in = SSresults.SScanLW_in;
SScanLW_in_tot = SSresults.SScanLW_in_tot;
SScanLW_out = SSresults.SScanLW_out;
SScanLW_out_tot = SSresults.SScanLW_out_tot;

SSsoilLW_in = SSresults.SSsoilLW_in;
SSsoilLW_out = SSresults.SSsoilLW_out;


% HEAT
SSLE_can_all_lay = SSresults.SSLE_can_all_lay;
SSLE_can_all = SSresults.SSLE_can_all;
SSLE_can_tot = SSresults.SSLE_can_tot;

SSH_can_all_lay = SSresults.SSH_can_all_lay;
SSH_can_all = SSresults.SSH_can_all;
SSH_can_tot = SSresults.SSH_can_tot;

SSsoilLE_out = SSresults.SSsoilLE_out;
SSsoilH_out = SSresults.SSsoilH_out;
SSsoilG = SSresults.SSsoilG;

% NET ENERGY
LWout_net = SSresults.LWout_net;
LE_net = SSresults.LE_net;
H_net = SSresults.H_net;
G_net = SSresults.G_net;
LWemi_net = SSresults.LWemi_net;


% X factor
Xeffent = SSresults.SSnetLW_Xout;


% NET ENTROPY
SSnetdir_in = SSresults.SSnetdir_in;
SSnetdif_in = SSresults.SSnetdif_in;
SSnetdir_out = SSresults.SSnetdir_out;
SSnetdif_out = SSresults.SSnetdif_out;

SSnetLW_in = SSresults.SSnetLW_in;

zicesl = VARIABLES.SOIL.zicesl;
% CANOPY
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
% IN
% shortwave
% 1. direct
SScandir_in_all_lay = repmat(SScandir_in,1,nspecies).*fLAIz;
SScandir_in_all = sum(SScandir_in_all_lay,1);
SScandir_in_tot = SScandir_in_tot;
% 2. difusse
SScandif_in_all_lay = repmat(SScandif_in,1,nspecies).*fLAIz;
SScandif_in_all = sum(SScandif_in_all_lay,1);
SScandif_in_tot = SScandif_in_tot;
% Longwave
% 3.  
SScanLW_in_all_lay = repmat(SScanLW_in,1,nspecies).*fLAIz;
SScanLW_in_all = sum(SScanLW_in_all_lay,1);
SScanLW_in_tot = SScanLW_in_tot;
% Total IN
SScan_in_all_lay =  SScandir_in_all_lay + SScandif_in_all_lay + SScanLW_in_all_lay;
SScan_in_all =  SScandir_in_all + SScandif_in_all + SScanLW_in_all;
%SScan_in = sum(SScandir_in_tot + SScandif_in_tot + SScanLW_in_tot);
SScan_in = sum(SScan_in_all);
% *************************************************************************

% OUT
% shortwave
% 1. direct
SScandir_out_all_lay = repmat(SScandir_out,1,nspecies).*fLAIz;
SScandir_out_all = sum(SScandir_out_all_lay,1);
SScandir_out_tot = SScandir_out_tot;
% 2. difusse
SScandif_out_all_lay = repmat(SScandif_out,1,nspecies).*fLAIz;
SScandif_out_all = sum(SScandif_out_all_lay,1);
SScandif_out_tot = SScandif_out_tot;
% Longwave
% 3.  
SScanLW_out_all_lay = repmat(SScanLW_out,1,nspecies).*fLAIz;
SScanLW_out_all = sum(SScanLW_out_all_lay,1);
SScanLW_out_tot = SScanLW_out_tot;
% Heat out
% 4 H (latent heat)
SScanLE_out_all_lay = SSLE_can_all_lay;
SScanLE_out_all = SSLE_can_all;
SScanLE_out_tot = SSLE_can_tot;
% 5 LE (sensible heat)
SScanH_out_all_lay = SSH_can_all_lay;
SScanH_out_all = SSH_can_all;
SScanH_out_tot = SSH_can_tot;
% Total OUT
SScan_out_all_lay =  SScandir_out_all_lay + SScandif_out_all_lay + SScanLW_out_all_lay + ...
                     SScanLE_out_all_lay + SScanH_out_all_lay;
SScan_out_all =  SScandir_out_all + SScandif_out_all + SScanLW_out_all+...
                 SScanLE_out_all + SScanH_out_all;
SScan_out = sum(SScan_out_all);
% *************************************************************************

% TOTAL CANOPY
SScan_all_lay = SScan_in_all_lay - SScan_out_all_lay;
SScan_all = SScan_in_all - SScan_out_all;
SScan_tot = SScan_in - SScan_out;




% SOIL
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
% IN
% shortwave
% 1. direct
SSsoildir_in_tot = SSsoildir_in;
% 2. difusse
SSsoildif_in_tot = SSsoildif_in;
% Longwave
% 3.  
SSsoilLW_in_tot = SSsoilLW_in;
% Total IN
SSsoil_in = sum(SSsoildir_in_tot + SSsoildif_in_tot + SSsoilLW_in_tot);
% *************************************************************************

% OUT
% shortwave
% 1. direct
SSsoildir_out_tot = SSsoildir_out;
% 2. difusse
SSsoildif_out_tot = SSsoildif_out;
% Longwave
% 3.  
SSsoilLW_out_tot = SSsoilLW_out;
% Heat out
% 4 H (latent heat)
SSsoilLE_out_tot = SSsoilLE_out;
% 5 LE (sensible heat)
SSsoilH_out_tot = SSsoilH_out;
% 6 G (ground heat flux)
%SSsoilG = SSsoilG;
% Total OUT
SSsoil_out = sum(SSsoildir_out_tot + SSsoildif_out_tot + SSsoilLW_out_tot...
    + SSsoilLE_out_tot + SSsoilH_out_tot);
% *************************************************************************

% TOTAL SOIL
SSsoil_tot = SSsoil_in - SSsoil_out;


% SOIL AND CANOPY
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
SSeco_tot = SScan_tot + SSsoil_tot;
%Xeffent = SSnewXout; 
%dif = 1000;
%count = 0;
%Xeffent = 1.0073;

    %while dif > 0.01
        
    %    Teffent = fzero(@(Te) Teffentropy (SSnetdir_in, SSnetdif_in, SSnetdir_out, SSnetdif_out, SSnetLW_in,... 
    %        SSeco_tot, LWout_net, LE_net, H_net, G_net, Xeffent, Te), 273);
        
        Teffent = fzero(@(Te)  Teffentropy2(SScanLW_out_tot, SSsoilLW_out_tot, ...
                   SSnetLW_in, SScanLW_in_tot, SSsoilLW_in_tot,...            
                  LWout_net, Te, Xeffent), 273);                                                

    
        
    %    Xeffent = fzero(@(Xnet)  Xeffentropy(SScanLW_out_tot, SSsoilLW_out_tot, ...
    %               SSnetLW_in, SScanLW_in_tot, SSsoilLW_in_tot,...            
    %              LWout_net, Teffent, Xnet), Xeffent);                                                
    %    if count > 0
    %        dif = abs(Teffent - Teffentprev)/Teffent; 
   %      %   dif = abs(Xeffent - Xeffentprev)/Xeffent;             
    %    end
    %    Xeffentprev = Xeffent;
    %    Teffentprev = Teffent;        
    %    count = count + 1;
   % end
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
% *************************************************************************
% STORE INOFORMATION IN SSresults

        % CANOPY
        % in
SSresults.Total.SScandir_in_all_lay = SScandir_in_all_lay;
SSresults.Total.SScandir_in_all = SScandir_in_all;
SSresults.Total.SScandir_in_tot = SScandir_in_tot;

SSresults.Total.SScandif_in_all_lay = SScandif_in_all_lay;
SSresults.Total.SScandif_in_all = SScandif_in_all ;
SSresults.Total.SScandif_in_tot = SScandif_in_tot;

SSresults.Total.SScanLW_in_all_lay = SScanLW_in_all_lay;
SSresults.Total.SScanLW_in_all = SScanLW_in_all;
SSresults.Total.SScanLW_in_tot = SScanLW_in_tot;

SSresults.Total.SScan_in_all_lay = SScan_in_all_lay;
SSresults.Total.SScan_in_all = SScan_in_all;
SSresults.Total.SScan_in = SScan_in;

        % out
SSresults.Total.SScandir_out_all_lay = SScandir_out_all_lay;
SSresults.Total.SScandir_out_all = SScandir_out_all;
SSresults.Total.SScandir_out_tot = SScandir_out_tot;

SSresults.Total.SScandif_out_all_lay = SScandif_out_all_lay;
SSresults.Total.SScandif_out_all = SScandif_out_all;
SSresults.Total.SScandif_out_tot = SScandif_out_tot;

SSresults.Total.SScanLW_out_all_lay = SScanLW_out_all_lay;
SSresults.Total.SScanLW_out_all = SScanLW_out_all;
SSresults.Total.SScanLW_out_tot = SScanLW_out_tot;

SSresults.Total.SScanLE_out_all_lay = SScanLE_out_all_lay;
SSresults.Total.SScanLE_out_all = SScanLE_out_all;
SSresults.Total.SScanLE_out_tot = SScanLE_out_tot;

SSresults.Total.SScanH_out_all_lay = SScanH_out_all_lay;
SSresults.Total.SScanH_out_all = SScanH_out_all;
SSresults.Total.SScanH_out_tot = SScanH_out_tot;

        % Total OUT
SSresults.Total.SScan_out_all_lay = SScan_out_all_lay;
SSresults.Total.SScan_out_all = SScan_out_all;
SSresults.Total.SScan_out = SScan_out;

        % TOTAL CANOPY
SSresults.Total.SScan_all_lay = SScan_all_lay;
SSresults.Total.SScan_all = SScan_all;
SSresults.Total.SScan_tot = SScan_tot;


        % SOIL
        % in
SSresults.Total.SSsoildir_in_tot = SSsoildir_in_tot;
SSresults.Total.SSsoildif_in_tot = SSsoildif_in_tot;
SSresults.Total.SSsoilLW_in_tot = SSsoilLW_in_tot;
SSresults.Total.SSsoil_in = SSsoil_in;

        % out 
SSresults.Total.SSsoildir_out_tot = SSsoildir_out_tot;
SSresults.Total.SSsoildif_out_tot = SSsoildif_out_tot;
SSresults.Total.SSsoilLW_out_tot = SSsoilLW_out_tot;
SSresults.Total.SSsoilLE_out_tot = SSsoilLE_out_tot;
SSresults.Total.SSsoilH_out_tot = SSsoilH_out_tot;
SSresults.Total.SSsoil_out = SSsoil_out;
        
        
        % TOTAL SOIL
SSresults.Total.SSsoil_tot = SSsoil_tot;

        % TOTAL ECOSYSTEM
SSresults.Total.SSeco_tot = SSeco_tot;

        % EFF TEMPERATURE BASED ON ENTROPY
        
SSresults.Teffent = Teffent;
SSresults.Xeffent = Xeffent;
