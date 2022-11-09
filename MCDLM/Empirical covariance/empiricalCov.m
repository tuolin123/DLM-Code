function covv = empiricalCov(D, field, mask)
% Compute the covariance function for stationary field empirically. Details
% of the procedure are provided in section 2.2
%--------------------------------------------------------------------------
% ARGUMENTS
% D        The dimension of the data
% field    Data with size Dim * n, where Dim is the image size and n is the
% mask     A mask (size Dim) of the data which is made of 1s and 0s. 1s for where
%          the data is inside the mask and 0s for where the data is ouside
%          the mask. Default is taken to be the mask with 1s everywhere
%--------------------------------------------------------------------------
% OUTPUT
% covv     the 3^D by 3^D covariance matrix
%--------------------------------------------------------------------------
% EXAMPLES
%2D example
% Dim = [91,109];
% nsubj = 100;
% FWHM = 2;
% noise = noisegen(Dim, nsubj, FWHM);
% D = 2;
% mask = zeros(Dim); mask(30:60,40:80) = 1; mask = logical(mask);
% sigma2 = empiricalCov(D, noise, mask);
%
%3D example
% Dim = [91,109,91];
% nsubj = 20;
% FWHM = 2;
% noise = noisegen(Dim, nsubj, FWHM);
% D = 3;
% mask = zeros(Dim); mask(30:60,40:80,30:60) = 1; mask = logical(mask);
% sigma2 = empiricalCov(D, noise, mask);
%
% % performs bad with small sample size, need to ensure covariance spd
% [V,D] = eig(sigma2);
% D = max(10e-10, D);
% sigma2 = V*D*V';
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
dim = size(field);

if ~exist('mask','var' )
    mask = ones(dim(1:D));
end

switch D
    case 1
        field_mask = field.*zero2nan(mask);
        cov0 = mean(field_mask(1:(dim(1)),:).*field_mask(1:dim(1),:), 'all', 'omitnan');
        cov1 = mean(field_mask(1:(dim(1)-1),:).*field_mask(2:dim(1),:), 'all', 'omitnan');
        cov2 = mean(field_mask(1:(dim(1)-2),:).*field_mask(3:dim(1),:), 'all', 'omitnan');
        c = [cov0,cov1,cov2];
        covv = toeplitz(c); % the covariance function in 1D is a toeplitz matrix
    case 2 
        % use the property of block toeplitz to reduce the computation time
        field_mask = field.*zero2nan(mask);
        cov0 = mean(field_mask(1:dim(1),1:dim(2),:).*...
            field_mask(1:dim(1),1:dim(2),:), 'all', 'omitnan');
        cov1v = mean(field_mask(1:(dim(1)-1),1:dim(2),:).*...
            field_mask(2:dim(1),1:dim(2),:), 'all', 'omitnan');
        cov1h = mean(field_mask(1:dim(1),1:(dim(2)-1),:).*...
            field_mask(1:dim(1),2:dim(2),:), 'all', 'omitnan');
        cov2v = mean(field_mask(1:(dim(1)-2),1:dim(2),:).*...
            field_mask(3:dim(1),1:dim(2),:), 'all', 'omitnan');
        cov2h = mean(field_mask(1:dim(1),1:(dim(2)-2),:).*...
            field_mask(1:dim(1),3:dim(2),:), 'all', 'omitnan');
        cov3 = mean(field_mask(1:(dim(1)-1),1:(dim(2)-1),:).*...
            field_mask(2:dim(1),2:dim(2),:), 'all', 'omitnan');
        cov4 = mean(field_mask(1:(dim(1)-2),1:(dim(2)-1),:).*...
            field_mask(3:dim(1),2:dim(2),:), 'all', 'omitnan');
        cov5 = mean(field_mask(1:(dim(1)-1),1:(dim(2)-2),:).*...
            field_mask(2:dim(1),3:dim(2),:), 'all', 'omitnan');
        cov6 = mean(field_mask(1:(dim(1)-2),1:(dim(2)-2),:).*...
            field_mask(3:dim(1),3:dim(2),:), 'all', 'omitnan');
        cov7 = mean(field_mask(1:(dim(1)-1),2:dim(2),:).*...
            field_mask(2:dim(1),1:(dim(2)-1),:), 'all', 'omitnan');
        cov8 = mean(field_mask(1:(dim(1)-1),3:dim(2),:).*...
            field_mask(2:dim(1),1:(dim(2)-2),:), 'all', 'omitnan');
        cov9 = mean(field_mask(1:(dim(1)-2),2:dim(2),:).*...
            field_mask(3:dim(1),1:(dim(2)-1),:), 'all', 'omitnan');
        cov10 = mean(field_mask(1:(dim(1)-2),3:dim(2),:).*...
            field_mask(3:dim(1),1:(dim(2)-2),:), 'all', 'omitnan');
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
        field_mask = field.*zero2nan(mask);
        
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
                    uniqmat(l,m,q) =  mean(field_mask(a1,a2,a3,:).*...
                        field_mask(b1,b2,b3,:), 'all', 'omitnan');
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