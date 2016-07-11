function [ADVr_m2, ADVr_m3, ADVrdd_m3] =...
    CN_ruptakeadvection (CC, CCdd, PARAMS, SWITCHES, aa,...
         sm, dz, RVID, layeruptake_all, dtime)        

        if SWITCHES.CN.Bioturbation;
            nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers    
        else
            nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers        
        end
        nspecies = PARAMS.CanStruc.nspecies; % number of species
        RVID = sum(RVID,2);
        
        rradvection = zeros(nl_soil,2);
        out_rradvection = zeros(nl_soil,1);
        in_rradvection = zeros(nl_soil,1);
        
        % PLANT UPTAKE
        for ii=1:1:nspecies
            layeruptn=layeruptake_all(:,ii);
            wtind=layeruptn>0;
            wtnind=layeruptn<=0;
            wtzind=layeruptn==0;
            % First for the case in which the flow is from the soil to the
            % root. By notation is positive toward the root and negative toward the soil  
            rradvection(wtind,ii) = (aa./sm(wtind).*CC(wtind)).*layeruptn(wtind);% [gr/m^2/d];                        
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed in the zone of disturbance
            rradvection(wtnind,ii) = (aa./sm(wtnind).*CCdd(wtnind)).*layeruptn(wtnind);% [gr/m^2/d];                                    
        end   
        % recompute indexes for net CC uptake
        ADVr_m2 = - sum(rradvection,2); % negative out of soil positive into
        wtind = ADVr_m2<0;
        wtnind = ADVr_m2>0;

       
        % check soil availability
        out_rradvection(wtind) = ADVr_m2(wtind)./dz(wtind)*dtime/86400;
        out_rradvection(~wtind) = 0;       
        ind_mb = CC < abs(out_rradvection);        % Check mass availability in that layer in [gr/m^3]
        ADVr_m2(ind_mb) = -CC(ind_mb).*dz(ind_mb)*86400/dtime;    % [gr/m^2/d]
        
        % check zone of disturbanace availability
        in_rradvection(wtnind) = ADVr_m2(wtnind)./(RVID(wtnind))*dtime/86400;
        in_rradvection(~wtnind) = 0;        
        ind_mb = CCdd < abs(in_rradvection);        % Check mass availability in that layer in [gr/m^3]
        ADVr_m2(ind_mb) = CCdd(ind_mb).*RVID(ind_mb)*86400/dtime;     % [gr/m^2/d]

        ADVr_m3 = ADVr_m2./dz;            % [gr/m^3/d]        
        ADVrdd_m3 = -ADVr_m2./RVID;            % [gr/m^3/d]

        
        %**** correct nan in frradvectionCCdd by zeros in RVID ************
        indnan = isnan(ADVrdd_m3);
        ADVrdd_m3(indnan) = 0; 
        
        