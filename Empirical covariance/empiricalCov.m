function covv = empiricalCov(D, field)
% empiricalCov(D, field) calculate the empirical covariance function for
% both stationary and non-stationary field.
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the field.
% field         the filed that the covariance function calculated for,
%               dimension should be d1*d2*n.
%--------------------------------------------------------------------------
% OUTPUT
% the 3^D by 3^D covariance matrix

dim = size(field);
nsim = dim(end);

switch D
    case 1
        cov0 = sum(sum(field(1:(dim(1)),:).*field(1:dim(1),:)))/(dim(1)*nsim);
        cov1 = sum(sum(field(1:(dim(1)-1),:).*field(2:dim(1),:)))/((dim(1)-1)*nsim);
        cov2 = sum(sum(field(1:(dim(1)-2),:).*field(3:dim(1),:)))/((dim(1)-2)*nsim);
        c = [cov0,cov1,cov2];
        covv = toeplitz(c);
    case 2
        cov0 = sum(sum(sum(field(1:dim(1),1:dim(2),:).*field(1:dim(1),1:dim(2),:))))/...
            (dim(1)*dim(2)*nsim);
        cov1v = sum(sum(sum(field(1:(dim(1)-1),1:dim(2),:).*field(2:dim(1),1:dim(2),:))))/...
            ((dim(1)-1)*dim(2)*nsim);
        cov1h = sum(sum(sum(field(1:dim(1),1:(dim(2)-1),:).*field(1:dim(1),2:dim(2),:))))/...
            (dim(1)*(dim(2)-1)*nsim);
        cov2v = sum(sum(sum(field(1:(dim(1)-2),1:dim(2),:).*field(3:dim(1),1:dim(2),:))))/...
            ((dim(1)-2)*dim(2)*nsim);
        cov2h = sum(sum(sum(field(1:dim(1),1:(dim(2)-2),:).*field(1:dim(1),3:dim(2),:))))/...
            (dim(1)*(dim(2)-2)*nsim);
        cov3 = sum(sum(sum(field(1:(dim(1)-1),1:(dim(2)-1),:).*field(2:dim(1),2:dim(2),:))))/...
            ((dim(1)-1)*(dim(2)-1)*nsim);
        cov4 = sum(sum(sum(field(1:(dim(1)-2),1:(dim(2)-1),:).*field(3:dim(1),2:dim(2),:))))/...
            ((dim(1)-2)*(dim(2)-1)*nsim);
        cov5 = sum(sum(sum(field(1:(dim(1)-1),1:(dim(2)-2),:).*field(2:dim(1),3:dim(2),:))))/...
            ((dim(1)-1)*(dim(2)-2)*nsim);
        cov6 = sum(sum(sum(field(1:(dim(1)-2),1:(dim(2)-2),:).*field(3:dim(1),3:dim(2),:))))/...
            ((dim(1)-2)*(dim(2)-2)*nsim);
        cov7 = sum(sum(sum(field(1:(dim(1)-1),2:dim(2),:).*field(2:dim(1),1:(dim(2)-1),:))))/...
            ((dim(1)-1)*(dim(2)-1)*nsim);
        cov8 = sum(sum(sum(field(1:(dim(1)-1),3:dim(2),:).*field(2:dim(1),1:(dim(2)-2),:))))/...
            ((dim(1)-1)*(dim(2)-2)*nsim);
        cov9 = sum(sum(sum(field(1:(dim(1)-2),2:dim(2),:).*field(3:dim(1),1:(dim(2)-1),:))))/...
            ((dim(1)-2)*(dim(2)-1)*nsim);
        cov10 = sum(sum(sum(field(1:(dim(1)-2),3:dim(2),:).*field(3:dim(1),1:(dim(2)-2),:))))/...
            ((dim(1)-2)*(dim(2)-2)*nsim);
        b1 = [cov0,cov1v,cov2v;cov1v,cov0,cov1v;cov2v,cov1v,cov0]';
        b2 = [cov1h,cov3,cov4;cov7,cov1h,cov3;cov9,cov7,cov1h]';
        b3 = [cov2h,cov5,cov6;cov8,cov2h,cov5;cov10,cov8,cov2h]';
        b4 = b2';
        b5 = b3';
        block = {b1,b2,b3,b4,b5};
        c = 1:3;
        r = [1,4,5];
        covv = cell2mat(block(toeplitz(c,r)));
end

%covv = [cov0 cov1h cov1v cov3 cov7 cov2h cov2v cov4 cov5 cov6 cov8];

end