function [A, B]  = CN_biomatrices  (dz, deltaz, D, dtime, topflux, ndim)

% This function generates the matrices A and B to compute the linear system
% AxC(n-1) = (C+B) then
% C(n-1) = A^(-1)x(C+B) 
% middle layers

for ii=2:1:ndim-1
    A(ii,ii-1) = -(dtime*D)/(dz(ii)*deltaz(ii-1));
    A(ii,ii) = (1+(dtime*D)/(dz(ii)*deltaz(ii))+(dtime*D)/(dz(ii)*deltaz(ii-1)));
    A(ii,ii+1) = -(dtime*D)/(dz(ii)*deltaz(ii)); 
end

%top layer
A(1,1) = 1+(dtime*D)/(dz(1)*deltaz(1));
A(1,2) = -(dtime*D)/(dz(1)*deltaz(1));

%bottom layer
A(ndim,ndim) = 1+(dtime*D)/(dz(ndim)*deltaz(ndim-1));
A(ndim,ndim-1) = -(dtime*D)/(dz(ndim)*deltaz(ndim-1));

% VECTOR OF BOUNDARY CONDITIONS
B = topflux*dtime./dz;
%B = zeros(ndim,1);
%B(1) = (topflux(1)*dtime/dz(1));
%B(2) = (topflux(2)*dtime/dz(2));
%B(ndim)=0;