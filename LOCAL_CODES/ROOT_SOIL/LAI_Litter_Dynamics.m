function [FORCING] = LAI_Litter_Dynamics(PARAMS, CONSTANTS, LAI_in, doy, yeind) 

%**************************************************************************
nspecies = PARAMS.CanStruc.nspecies;
%    PREALLOCATE MATRICES TO IMPLEMENT
[dim1, dim2] = size(LAI_in);    
BMTt = nan(1,dim1);
TBMLa = nan(nspecies,dim1);
%**************************************************************************
% Time
dtime=CONSTANTS.dtime;
% ORGANIZE STORE OF SLAC
SLACstor = PARAMS.CN.SLAC;
% THIS FUNCTION IS MODIFIED TO ALLOW THE INCLUSION OF OTHER SHRUBS SPECIES
%************************   MODIFICATION  *********************************

for ii = 1:1:nspecies
    LAIvec = LAI_in(:,ii);
    SLACvec = SLACstor{ii};
    leafspan = PARAMS.CN.leafspan(ii);                    % leaf span [years]   
    for jj=1:1:yeind-1
        % LAI
        LAIi = LAIvec(jj);                                % LAI custom time step
        LAIf = LAIvec(jj+1);                                % LAI future time step   
        %SLAC Computed based on the day
        index = doy(jj)*48;
        SLAC = SLACvec(index);                              % Specific Leaf Area [m2 leaf/gr C dry leaf]   
        % compute normal turnover and phenologic changes in LAI units
        LAInt = LAIi/(leafspan*31536000)*dtime;
        LAIph = LAIf - LAIi + LAInt;
        % Compute LAI turnover biomass 
        if LAIph >= 0
            LAITt = LAInt;
        elseif LAIph < 0
            LAITt = abs(LAIph) + LAInt;
        end
        BMTt(jj) = LAITt/SLAC;
    end
    BMTt(length(LAIvec)) = BMTt(length(LAIvec)-1);
    TBMLa(ii,:) = BMTt;
end        

FORCING.TBMla = TBMLa;               % Total biomass Drop in this time step as litter from above [gr / m^2]     
