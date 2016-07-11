function [LADnorm] = LAD_Weibull_shrubs(PARAMS, VERTSTRUC)

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
    beta = PARAMS.CanStruc.beta(2);
    alpha = PARAMS.CanStruc.alpha(2);
     
    hcan = PARAMS.CanStruc.hcan;
    znc = VERTSTRUC.znc;     
     
    zz = znc./hcan;
    
    ind = zz >= 0.5 ;    
    nind = zz < 0.5 ;
    LAD(ind) = (beta\alpha)*(zz(ind)./alpha).^(beta-1).*exp(-(zz(ind)./alpha).^beta);
    
    LAD(nind) = 0;
        
    LAD = LAD(:);
    LAD = flipdim(LAD,1);
    LADnorm = LAD./sum(LAD);

    
    
    