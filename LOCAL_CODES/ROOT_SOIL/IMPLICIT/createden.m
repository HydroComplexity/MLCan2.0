function [den]=createden(nl_soil,zsoi)

% Compute the distance between the nodes.
%-------------------------------------------------------------------------
% This function computes the distance between the different nodes in all 
% layers.
% Written by Juan Camilo Quijano, UIUC, 2008
%
%------------------------- Input Variables -------------------------------
%       nl_soil     % Number of layers 
%       zsoi         % Vector with depth of midpoint of each layer [mm]
%
%------------------------- Output Variables ------------------------------
%       den        % Distance the nodes of the layers.
%
%-------------------------------------------------------------------------

den=zeros(nl_soil-1,1);
for i=2:1:nl_soil;
    den(i-1)=zsoi(i)-zsoi(i-1);
end    