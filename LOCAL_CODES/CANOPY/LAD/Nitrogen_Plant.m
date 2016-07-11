function [FORCING] = Nitrogen_Plant(PARAMS, CONSTANTS, datafile, LAI1_in) 

%**************************************************************************
%    PREALLOCATE MATRICES TO IMPLEMENT

TBMLavec = nan(1,length(LAI1_in));
TBMLa = nan(3,length(LAI1_in));
%**************************************************************************

% ORGANIZE STORE LAI
LAIstor(1,:) = LAI1_in;
load (datafile,'shrubslai_2003');
LAIstor(2,:) = shrubslai_2003.manzanita;
LAIstor(3,:) = shrubslai_2003.ceanothus;

% ORGANIZE STORE OF SLA
SLAstor = PARAMS.CN.SLA;

% ORGANIZE STORE OF NL
Nalstor = PARAMS.CN.Nal;

%************************   MODIFICATION  *********************************
% THIS FUNCTION IS MODIFIED TO ALLOW THE INCLUSION OF OTHER SHRUBS SPECIES
nspecies = 3;
%************************   MODIFICATION  *********************************


for ii = 1:1:nspecies
    LAIivec = LAIstor(ii,1:length(LAI1_in)-1);
    LAIfvec = LAIstor(ii,2:length(LAI1_in));
    Nalvec = Nalstor{ii};
    Nalivec = Nalvec(1:length(LAI1_in)-1);
    Nalfvec = Nalvec(2:length(LAI1_in));
        
    SLAvec = SLAstor{ii};
    leafspan = PARAMS.CN.leafspan(ii);                    % leaf span [years]   
    for jj=1:1:length(LAI1_in)-1
        %******************************************************************
        %*******
        %                          DE-REFERENCE BLOCK
        %*************************************************************************
        SLA = SLAvec(jj);                              % Specific Leaf Area [m2 leaf/gr dry leaf]   

        % % LAI
        LAIi = LAIivec(jj);                                % LAI custom time step
        LAIf = LAIfvec(jj);                                % LAI future time step   
        % Time
        dtime=CONSTANTS.dtime;
        % Nitrogen Dynamics

        %*************************************************************************
        % COMPUTE LEAVES BIOMASS TO BUILD, LITTER

        % Compute LAI turnover biomass 
        BMLAIi = LAIi/SLA;                                      %[gr/m2]      
        BMLAIf = LAIf/SLA;                                      %[gr/m2]
        rate_NT = (1/(365*leafspan))*BMLAIi;                    % rate normal litter [gr/m2/d ]
        BMLAInl = rate_NT * (1/86400) * dtime;                  % Biomass from Normal Litter gr/m2  
        LAInl = BMLAInl * SLA;                                  % [m2 leaf / m2] 
        % compute deltaLAI and deltaLAI biomass
        BMLAId = abs(BMLAIf - BMLAIi + BMLAInl);                     % [gr/m2]
        LAId = BMLAId*SLA;                                      % Delta LAI [m2 leaf / m2]    
        % compute the total amount of biomass in the leaves remainning without taking
        % into account the the new leaves that are going to be built
        BMLAIrem = BMLAIi - BMLAInl;                                 % [gr/m2]
        LAIrem = BMLAIrem*SLA; 

        if (BMLAIf - BMLAIi + BMLAInl) < 0                                    % There is not construction of leaves in this time step
          if LAInl ~= LAId
              stop=3;
          end

        % Compute the total amount of biomass that is going as litter from above ground.    
           TBMLavec(jj) = BMLAInl + BMLAId; 
        elseif BMLAIf - BMLAIi + BMLAInl >= 0 
          if LAInl ~= LAId
              stop=3;
          end    
           TBMLavec(jj) = BMLAInl ; 
        end
        
        %*************************************************************************
        % COMPUTE NITROGEN DYNAMICS
        %delatN = 
        
        % NAL
        Nali = Nalivec(jj);                                % Nal custom time step
        Nalf = Nalfvec(jj);                                % Naf future time step
        deltaNal = Nalf - Nali; 
        deltaLAI = LAIf -LAIi;
        if (BMLAIf - BMLAIi + BMLAInl) < 0
            AddNdatavec(jj) = LAIi*deltaNal + deltaLAI*Nali +  deltaLAI*deltaLAI + Nali*(LAId + LAInl)*(1-PARAMS.CN.Nrecyc);
        else
            AddNdatavec(jj) = LAIi*deltaNal + deltaLAI*Nali +  deltaLAI*deltaLAI + Nali*LAInl*(1-PARAMS.CN.Nrecyc);            
        end
        AddNvec(jj) = AddNdatavec(jj) + LAInl*Nali;
        nlitter(jj) = TBMLavec(jj)*Nali;
    end
    LAIdm(ii,:) = LAId;
    LAIdnit(ii,:) = LAId.*Nalfvec;
    TBMLavec(length(LAI1_in)) = TBMLavec(length(LAI1_in)-1); 
    TBMLa(ii,:) =  TBMLavec;
    AddNdata(ii,:) = AddNdatavec; 
    AddN(ii,:) = AddNvec;    
    addl(ii,:) = nlitter; 
end

FORCING.TBMla = TBMLa;               % Total biomass Drop in this time step as litter from above [gr / m^2]     
