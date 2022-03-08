function d = ddlm(z,rho,nnb,D)
% ddlm(z,rho,nnb,D) computes the density of local maxima using the Analytical DLM method.
%--------------------------------------------------------------------------
% ARGUMENTS
% z             a vector of local maxima value.
% rho           spatial correlation in each direction lattice. Input should be a n*d matrix with n voxels and d dimensions.
% nnb           number of neighbor in each direction lattice, can be 0,1,2. Input should be a n*d matrix with n voxels and d dimensions. 
% D             the dimension of the image.
%--------------------------------------------------------------------------
% OUTPUT
% d             the density of the discrete local maxima at value z.
%--------------------------------------------------------------------------
% EXAMPLES
% z = -2:0.001:2; 
% nz = length(z);
% D = 3;
% rho = ones(nz,D)*0.8;
% nnb = ones(nz,D)*2;
% dendlm = ddlm(z,rho,nnb,D);
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
switch D
    case 1
        ld = 1;
        PQd = asin(sqrt((1-rho(z==0,1).^2)./2))./pi; % correct for z is zero.
    case 2
        ld = 2;
        PQd = asin(sqrt((1-rho(z==0,1).^2)./2))./pi.*asin(sqrt((1-rho(z==0,2).^2)./2))./pi;
    case 3
        ld = 3;
        PQd = asin(sqrt((1-rho(z==0,1).^2)./2))./pi.*asin(sqrt((1-rho(z==0,2).^2)./2))./pi.*asin(sqrt((1-rho(z==0,3).^2)./2))./pi;
end
phi=exp(-z.^2./2)./sqrt(2*pi);
nz=length(z);
PQ=zeros(1,nz);
for i=1:nz
    Q=1;
    zi=z(i);
    for l=1:ld
        % loop in all dimensions
        alpha=asin(sqrt((1-rho(i,l)^2)/2));
        da=alpha/1000;
        h=sqrt((1-rho(i,l))/(1+rho(i,l)));
        zplus=max(zi,0);
        zall=0:da:alpha;
        f=exp(-1/2*h^2*zi^2./(sin(zall)).^2); % third part in equation (2) in manuscript
        if nnb(i,l)==1
            % Whether the neighbor in one dimension includes 1 or 2 voxels
            Q=Q*(1-erfc(h*zplus/sqrt(2))/2);
        end
        if nnb(i,l)==2
            Q=Q*(1-erfc(h*zplus/sqrt(2))+sum(f)*da/pi);
        end
    end
    PQ(i)=Q;
end
if any(z == 0)
    PQ(z==0) = PQd;
end
d=PQ.*phi;
return