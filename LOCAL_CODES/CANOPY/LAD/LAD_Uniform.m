function [cstruct] = LAD_Uniform (znc)

% Computes a normalized LAD profile using the beta distribution described
% in [Wang et al, Forest Ecology and Management, 1990, 32, 217-237]
    
    cstruct = ones(size(znc))./length(znc);
    cstruct = cstruct(:);
    