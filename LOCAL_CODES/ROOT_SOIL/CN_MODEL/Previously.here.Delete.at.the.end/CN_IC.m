

    if  SWITCHES.initialcond == 1  % load from file
            N = 365;
            load ('./RUNSNHR/year500.mat');
            VARIABLES.Cl = strdata.Cl(:,N); 
            VARIABLES.Ch = nan(13,1);
%            VARIABLES.Ch = strdata.Ch(:,N);
            VARIABLES.Cb = strdata.Cb(:,N);
            VARIABLES.Nl = strdata.Nl(:,N);    
            VARIABLES.Amm = strdata.Amm(:,N);
            VARIABLES.Nit = strdata.Nit(:,N); 
            VARIABLES.CNl =strdata.CNl(:,N); 
            VARIABLES.UP_amm = strdata.UP_amm(:,N);
            VARIABLES.Up_nit = strdata.UP_nit(:,N);
            VARIABLES.LCH_amm = strdata.LCH_amm(:,N);
            VARIABLES.LCH_nit = strdata.LCH_nit(:,N);
            VARIABLES.MIN_net = strdata.MIN_net(:,N);
            VARIABLES.IMM_gross = strdata.IMM_gross(:,N); 
            VARIABLES.MIN_net = strdata.MIN_net(:,N);
            VARIABLES.IMM_gross = strdata.IMM_gross(:,N); 
            VARIABLES.PHI = strdata.PHI(:,N);
            VARIABLES.phi = strdata.phi(:,N);    
            VARIABLES.fSd = strdata.fSd(:,N);
            VARIABLES.fTd = strdata.fTd(:,N);
    elseif SWITCHES.initialcond == 2        % fill wit numbers  
            if SWITCHES.litter;
            Cl =[PARAMS.CN.Clitter 1000 1000 1000 1000 1000 1000 10 10 10 10 10 10];
            VARIABLES.Cl = Cl(:);
            Cb =[500 500 500 500 100 100 100 40 40 40 10 10 10];
            VARIABLES.Cb = Cb(:);
            Ch =[30000 30000 5000 5000 5000 5000 5000 1000 1000 1000 100 100 100];
            VARIABLES.Ch =Ch(:);
            VARIABLES.Amm = 0.5*ones(nl_soil+1,1);
            VARIABLES.Nit = 2*ones(nl_soil+1,1);
            VARIABLES.CNl(1:nl_soil+1,1) = PARAMS.CN.CNb./(1-PARAMS.CN.rr);  % Initialize with CN min
            VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;             
            else
            Cl =[1000 1000 1000 1000 1000 1000 10 10 10 10 10 10];
            VARIABLES.Cl = Cl(:);
            Cb =[500 500 500 100 100 100 40 40 40 10 10 10];
            VARIABLES.Cb = Cb(:);
            Ch =[30000 5000 5000 5000 5000 5000 1000 1000 1000 100 100 100];
            VARIABLES.Ch =Ch(:);
            VARIABLES.Amm = 0.5*ones(nl_soil,1);
            VARIABLES.Nit = 2*ones(nl_soil,1);
            VARIABLES.CNl(1:nl_soil+1,1) = PARAMS.CN.CNb./(1-PARAMS.CN.rr);  % Initialize with CN min
            VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;             
            end       
    elseif SWITCHES.initialcond == 3        % use the same as Cl_bloi            
            if SWITCHES.CN.Bioturbation;
            Cl =[PARAMS.CN.Clitter ; Cl_bloi];
            VARIABLES.Cl = Cl(:);            
            Cb =[Cb_bloi(1) ; Cb_bloi];
            VARIABLES.Cb = Cb(:);
            Ch =[30000 30000 5000 5000 5000 5000 5000 1000 1000 1000 100 100 100];
            VARIABLES.Ch =Ch(:);
            VARIABLES.Amm = zeros(nl_soil+1,1);%0.5*ones(nl_soil+1,1);
            VARIABLES.Nit = zeros(nl_soil+1,1);%2*ones(nl_soil+1,1);
            VARIABLES.CNl(1:nl_soil+1,1) = PARAMS.CN.CNb./(1-PARAMS.CN.rr);  % Initialize with CN min
            VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;             
            else
            Cl =Cl_bloi;
            VARIABLES.Cl = Cl(:);
            Cb =Cb_bloi;
            VARIABLES.Cb = Cb(:);
            Ch =[30000 5000 5000 5000 5000 5000 1000 1000 1000 100 100 100];
            VARIABLES.Ch =Ch(:);
            VARIABLES.Amm = 0.5*ones(nl_soil,1);
            VARIABLES.Nit = 2*ones(nl_soil,1);
            PARAMS.CN.CNb = ones(12,1)*11.8;            
            VARIABLES.CNl(1:nl_soil,1) = PARAMS.CN.CNb./(1-PARAMS.CN.rr);  % Initialize with CN min
            VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;             
            end       
     end
                          