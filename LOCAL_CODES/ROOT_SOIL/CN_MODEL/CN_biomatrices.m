function [A, B]  = CN_biomatrices  (dz, deltaz, D, dtime, BC)

% This function generates the matrices A and B to compute the linear system
% AxC(n-1) = (C+B) then
% C(n-1) = A^(-1)x(C+B) 
% middle layers
ndim = length(dz);
% allocate matrix A
A = zeros(ndim,ndim);
for ii=2:1:ndim-1
    D1 = D(ii-1);
    D2 = D(ii);
    delz1 = deltaz(ii-1);
    delz2 = deltaz(ii);    
    A(ii,ii-1) = -(dtime*D1)/(dz(ii)*delz1);
    A(ii,ii) = (1+(dtime*D1)/(dz(ii)*delz1)+(dtime*D2)/(dz(ii)*delz2));
    A(ii,ii+1) = -(dtime*D2)/(dz(ii)*delz2); 
end

%top layer
D2 = D(1);
delz2 = deltaz(1);
A(1,1) = 1+(dtime*D2)/(dz(1)*delz2);
A(1,2) = -(dtime*D2)/(dz(1)*delz2);

%bottom layer
D1 = D(ndim-1);
delz1 = dz(ndim-1);
A(ndim,ndim) = 1+(dtime*D1)/(dz(ndim)*delz1);
A(ndim,ndim-1) = -(dtime*D1)/(dz(ndim)*delz1);

% VECTOR OF BOUNDARY CONDITIONS
B = BC*dtime./dz;
%B = zeros(ndim,1);
%B(1) = (topflux(1)*dtime/dz(1));
%B(2) = (topflux(2)*dtime/dz(2));
%B(ndim)=0;