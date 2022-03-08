function covv = empiricalCov(D, field)
% Compute the covariance function for stationary field empirically. Details
% of the procedure are provided in section 2.2
%--------------------------------------------------------------------------
% ARGUMENTS
% D             the dimension of the field.
% field         the fields that the covariance function calculated for
%--------------------------------------------------------------------------
% OUTPUT
% covv          the 3^D by 3^D covariance matrix
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [50 50]; nsubj = 100; D = length(Dim); FWHM = 5;
% noise = noisegen(Dim, nsubj, FWHM);
% sigma2 = empiricalCov(D, noise);
%
% % performs bad with small sample size, need to ensure covariance spd
% [V,D] = eig(sigma2);
% D = max(10e-10, D);
% sigma2 = V*D*V';
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------

dim = size(field);
nsim = dim(end);

switch D
    case 1
        cov0 = sum(sum(field(1:(dim(1)),:).*field(1:dim(1),:)))/(dim(1)*nsim);
        cov1 = sum(sum(field(1:(dim(1)-1),:).*field(2:dim(1),:)))/((dim(1)-1)*nsim);
        cov2 = sum(sum(field(1:(dim(1)-2),:).*field(3:dim(1),:)))/((dim(1)-2)*nsim);
        c = [cov0,cov1,cov2];
        covv = toeplitz(c); % the covariance function in 1D is a toeplitz matrix
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
    case 3
        sz = [3 3 3];
        [I1, I2, I3] = ind2sub(sz, 1:27);
        f = @(x,y)(x-y);
        diff1 = bsxfun(f,I1,I1');
        diff2 = bsxfun(f,I2,I2');
        diff3 = bsxfun(f,I3,I3');
        covmat = zeros(3^3);
        uniqmat = zeros([5 5 5]); % in each dimension the value varies from -2 to 2
        
        %together 125 possibilities when estimating E(XY) for 3D X,Y
        %first dimension l
        for l = 1:5
            if l == 1
                a1 = 3:dim(1); %index for first dimension of X in E(XY)
                b1 = 1:(dim(1)-2); %index for first dimension of Y in E(XY)
            elseif l == 2
                a1 = 2:dim(1);
                b1 = 1:(dim(1)-1);
            elseif l == 3
                a1 = 1:dim(1);
                b1 = 1:dim(1);
            elseif l == 4
                a1 = 1:(dim(1)-1);
                b1 = 2:dim(1);
            else
                a1 = 1:(dim(1)-2);
                b1 = 3:dim(1);
            end
            %second dimension m
            for m = 1:5
                if m == 1
                    a2 = 3:dim(1); %index for second dimension of X in E(XY)
                    b2 = 1:(dim(1)-2); %index for second dimension of Y in E(XY)
                elseif m == 2
                    a2 = 2:dim(1);
                    b2 = 1:(dim(1)-1);
                elseif m == 3
                    a2 = 1:dim(1);
                    b2 = 1:dim(1);
                elseif m == 4
                    a2 = 1:(dim(1)-1);
                    b2 = 2:dim(1);
                else
                    a2 = 1:(dim(1)-2);
                    b2 = 3:dim(1);
                end
                %third dimension q
                for q = 1:5
                    if q == 1
                        a3 = 3:dim(1); %index for third dimension of X in E(XY)
                        b3 = 1:(dim(1)-2); %index for third dimension of Y in E(XY)
                    elseif q == 2
                        a3 = 2:dim(1);
                        b3 = 1:(dim(1)-1);
                    elseif q == 3
                        a3 = 1:dim(1);
                        b3 = 1:dim(1);
                    elseif q == 4
                        a3 = 1:(dim(1)-1);
                        b3 = 2:dim(1);
                    else
                        a3 = 1:(dim(1)-2);
                        b3 = 3:dim(1);
                    end
                    uniqmat(l,m,q) = sum(sum(sum(sum(field(a1,a2,a3,:).*field(b1,b2,b3,:)))))/...
                    (length(a1)*length(a2)*length(a3)*nsim);
                end
            end
        end

     
        for i = 1:3^3
            for j = 1:3^3
                covmat(i,j) = uniqmat(diff1(i,j)+3, diff2(i,j)+3, diff3(i,j)+3);
            end
        end
        covv = covmat;
end



end