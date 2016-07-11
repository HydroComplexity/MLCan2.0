function [Cl_final, dz_final, zh_final, zn_final, deltaz_final] = CN_layersdif (dz, Cl, layerbio, type, ind_join)

% Input variables

%  dz  Depth of each layer in the soil
%  zh  Bottom boundary of each layer in the soil
%  Cl  Carbon concentration
%  layerbio Maximum layer at which bioturbation happens

% Calculate Vector indicating the layers to use in diffusion



% CUT MATRIX
if type == 1
    ind_join = ind_join(1:layerbio);
    Cl = Cl(1:layerbio);
    dz = dz(1:layerbio);
    
    Cldz_vec = Cl.*dz;

    Cldz_final = accumarray(ind_join, Cldz_vec);
    dz_final = accumarray(ind_join, dz);
    Cl_final = Cldz_final./dz_final;
    zh_final = cumsum(dz_final);
    % compute nodes for these new layers
    zn_final = (zh_final + [0 ; zh_final(1:length(zh_final)-1)])/2;  
    deltaz_final = zn_final(2:length(zn_final))-zn_final(1:length(zn_final)-1);

elseif type == 2
    % compute vector with delta in carbon
    % compute vector of weights 
    dz_final = dz(:,1);
    dz_diff = dz(:,2);
    dweight = (dz_final./dz_diff(ind_join));
    Cl_add = Cl(ind_join).*dz_diff(ind_join);     %[gr/m2]
    Cl_final_m2 = Cl_add.*dweight;  %[gr/m2]
    Cl_final = Cl_final_m2./dz_final;
    zh_final=[];
    zn_final=[];    
    deltaz_final = [];
end
    