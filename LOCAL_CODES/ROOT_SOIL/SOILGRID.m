function [znode_all, dz_all, zlayer_all, nl_soil] = SOILGRID (PARAMS)
%
%   Written By: Darren Drewry 
%   Date:       1/22/08
%
%   This function constructs the soil grid, with 0.1m grid spacing down to 1m
%   depth, and 0.5 meter spacing below that, down to the specified maxdepth
%
%   Outputs: 
%       znode  - node depths (center of layers)     [m]
%       dznode - layer thicknesses                  [m]
%       zlayer - depths of soil layer interfaces    [m]
%       nl_soil - number of soil layers             [-]
%
%   Inputs:
%       dzs = layer spacings [m]
%       depths = depth of soil column to grid for each layer spacing (dzs) [m]
%                --> sum(depths) = total depth of soil column
%
% *************************************************************************
%  The code was changed. It uses the same distribution  as Amenu (Common land model)


% % De-reference Structure Values
%     dzs = PARAMS.Soil.dzs;
%     depths = PARAMS.Soil.depths;
     nl_soil = PARAMS.nl_soil;
     scaleza = PARAMS.Soil.scaleza;
     scalezb = PARAMS.Soil.scalezb;

% ALLOCATE MEMORY
    zsoi = nan(nl_soil,1);
    dzsoi = nan(nl_soil,1);
    zsoih = nan(nl_soil,1);

     % 

    % Soil layer node depths, i.e., depth of layer center from surface [m]
for j = 1:nl_soil
    zsoi(j) = scaleza*(exp(scalezb*(j-0.5))-1);  %node depths
end

% Soil layer thicknesses 
dzsoi(1)  = 0.5*(zsoi(1)+zsoi(2));        
dzsoi(nl_soil)= zsoi(nl_soil)-zsoi(nl_soil-1);
for j = 2:nl_soil-1
    dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1));   
end

% Soil layer interface depths from the surface [m]
zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil);
for j = 1:nl_soil-1
    zsoih(j)= 0.5*(zsoi(j)+zsoi(j+1));      
end
    % Make Column Vectors
        znode_all = zsoi(:);
        dz_all = dzsoi(:);
        zlayer_all = zsoih(:);    
    
      
%     if (length(dzs) ~= length(depths))
%         disp('*** ERROR: ''dzs'' and ''depths'' must have same length!');
%         return;
%     end
% 
%     
%     % Construct Grid
%     zl = 0;
%     znode_all=[]; dz_all=[]; zlayer_all=[];
%     for ii = 1:length(dzs)
%         znode=[]; dz=[]; zlayer=[];
%         
%         dzii = dzs(ii);
%         depthii = depths(ii);
%         
%         znode = linspace(dzii/2, depthii, (depthii-dzii)/dzii+1);
%         dz = ones(1,length(znode))*dzii;
%         zlayer = znode + dz/2;
%         
%         if (ii==1)
%             znode_all = [znode_all, znode];
%             zlayer_all = [zlayer_all, zlayer];
%         else
%             znode_all = [znode_all, znode+zlayer_all(end)];
%             zlayer_all = [zlayer_all, zlayer+zlayer_all(end)];
%         end
%         dz_all = [dz_all, dz];
%         
%     end
%     
%     nl_soil = length(zlayer_all);
%     
%     % Make Column Vectors
%         znode_all = znode_all(:);
%         dz_all = dz_all(:);
%         zlayer_all = zlayer_all(:);
    
          