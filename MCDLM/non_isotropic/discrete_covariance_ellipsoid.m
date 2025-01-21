
function discov = discrete_covariance_ellipsoid(nu, D)
% Compute the correlation (covariance) function of the stationary field 
% from convolving the white noise with elliptical Gaussian kernel on 
% discrete lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% nu            the standard deviation of the D-dim Gaussian smoothing kernel in
%               each direction, length D vector
% D             the dimension of the stationary field, at least 2
%--------------------------------------------------------------------------
% OUTPUT
% discov        a 3^D*3^D convariance matrix to be used in simulateLocMax_Discrete.m 
%--------------------------------------------------------------------------
% EXAMPLES
% nu = [0.3 5];
% D = 2;
% discrete_covariance_ellipsoid(nu, D)
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------

switch D
    % For 1D case we can use discrete.covariance.m
    case 2
        L = max(max(2,ceil(5*nu))); % how is the support defined
        dL = 1;
        [X, Y] = meshgrid( (-L:dL:L) , ...
                           (-L:dL:L));
        if nu ~= 0
            h = ellipsoid_kernel(D, nu, {X,Y}); 
            h = h / sqrt(sum((h(:).^2)));

            % Convolution shortcut (we only need it near the origin)
            %tic
            rhot2 = convn(h, h, 'same');
            %toc
            % Write out all possible neighboring combinations
            rho0 = rhot2(L+1, L+1);
            rho1v = rhot2(L+2, L+1);
            rho1h = rhot2(L+1, L+2);
            rho2v = rhot2(L+3, L+1);
            rho2h = rhot2(L+1, L+3);
            rho3 = rhot2(L+2, L+2);
            rho4 = rhot2(L+3, L+2);
            rho5 = rhot2(L+2, L+3);
            rho6 = rhot2(L+3, L+3);
            rho7 = rhot2(L, L+2);
            rho8 = rhot2(L, L+3);
            rho9 = rhot2(L-1, L+2);
            rho10 = rhot2(L-1, L+3);
            b1 = [rho0,rho1v,rho2v;rho1v,rho0,rho1v;rho2v,rho1v,rho0]';
            b2 = [rho1h,rho3,rho4;rho7,rho1h,rho3;rho9,rho7,rho1h]';
            b3 = [rho2h,rho5,rho6;rho8,rho2h,rho5;rho10,rho8,rho2h]';
            b4 = b2';
            b5 = b3';
            block = {b1,b2,b3,b4,b5};
            c = 1:3;
            r = [1,4,5];
            discov = cell2mat(block(toeplitz(c,r)));
        else
            discov = eye(3^D); %indepdent case
        end
        
    case 3
        L = max(max(2,ceil(5*nu))); % how is the support defined
        dL = 1;
        [X, Y, Z] = meshgrid( (-L:dL:L) , ...
                              (-L:dL:L ) , ...
                              (-L:dL:L ));  
        if nu ~= 0
            h = ellipsoid_kernel(D, nu, {X,Y,Z}); 
            h = h / sqrt(sum((h(:).^2)));

            % Convolution shortcut (we only need it near the origin)
            %tic
            rhot2 = convn(h, h, 'same');
            
            sz = [3 3 3];
            [I1, I2, I3] = ind2sub(sz, 1:27);
            f = @(x,y)(x-y);
            diff1 = bsxfun(f,I1,I1');
            diff2 = bsxfun(f,I2,I2');
            diff3 = bsxfun(f,I3,I3');
            covmat = zeros(3^3);
            uniqmat = zeros([5 5 5]); % in each dimension the value varies from -2 to 2

            % Together 125 possibilities when estimating E(XY) for 3D X,Y
            % Use for loop to assign all possible neighboring combinations
            %first dimension l
            for l = 1:5
                if l == 1
                    a1 = L-1;
                elseif l == 2
                    a1 = L;
                elseif l == 3
                    a1 = L+1;
                elseif l == 4
                    a1 = L+2;
                else
                    a1 = L+3;
                end
                %second dimension m
                for m = 1:5
                    if m == 1
                        a2 = L-1;
                    elseif m == 2
                        a2 = L;
                    elseif m == 3
                        a2 = L+1;
                    elseif m == 4
                        a2 = L+2;
                    else
                        a2 = L+3;
                    end
                    %third dimension q
                    for q = 1:5
                        if q == 1
                            a3 = L-1;
                        elseif q == 2
                            a3 = L;
                        elseif q == 3
                            a3 = L+1;
                        elseif q == 4
                            a3 = L+2;
                        else
                            a3 = L+3;
                        end
                        uniqmat(l,m,q) = rhot2(a1,a2,a3);
                    end
                end
            end
            for i = 1:3^3
                for j = 1:3^3
                    covmat(i,j) = uniqmat(diff1(i,j)+3, diff2(i,j)+3, diff3(i,j)+3);
                end
            end
            discov = covmat;
            %toc
        else
            discov = eye(3^D); %indepdent case
        end
end


