function [LADnorm] = LAD_Weibull_pine(PARAMS, VERTSTRUC)

% Computes a normalized LAD profile using the beta distribution described
% in [Wang et al, Forest Ecology and Management, 1990, 32, 217-237]

    % De-reference Structure Values
%    hcan = PARAMS.CanStruc.hcan;
%    znc = VERTSTRUC.znc;     
    
%    Zopt = PARAMS.CanStruc.Zopt; 
%    sigsqd = PARAMS.CanStruc.sigsqd; 

%    zz = znc./hcan;
    
%    LAD = exp(-((zz-Zopt).^2) ./ (sigsqd));    
%    LAD = LAD(:);
%    LAD = flipdim(LAD,1);
%    LADnorm = LAD./sum(LAD);
    
    % De-reference Structure Values
    beta = PARAMS.CanStruc.beta;
    alpha = PARAMS.CanStruc.alpha;
     
    hcan = PARAMS.CanStruc.hcan;
    znc = VERTSTRUC.znc;     
     
    zz = znc./hcan;
    
    LAD = (beta\alpha)*(zz./alpha).^(beta-1).*exp(-(zz./alpha).^beta);
    LAD = LAD(:);
    LAD = flipdim(LAD,1);
    LADnorm = LAD./sum(LAD);

    
    
    