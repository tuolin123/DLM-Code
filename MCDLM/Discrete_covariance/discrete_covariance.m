function discov = discrete_covariance(nu, D)
% Compute the correlation (covariance) function of the isotropic field from 
% convolving the white noise with isotropic Gaussian kernel on discrete
% lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% nu            the standard deviation of the smoothing kernel
% D             the dimension of the isotropic field
%--------------------------------------------------------------------------
% OUTPUT
% discov        the correlation function to be used in simulateLocMax_Discrete.m 
%--------------------------------------------------------------------------
% EXAMPLES
% nu = 0.3;
% D = 2;
% discrete_covariance(nu, D)
%--------------------------------------------------------------------------
% AUTHOR: Tuo Lin
%--------------------------------------------------------------------------
switch D
    case 1
        L = max(2,ceil(5*nu)); % how is the support defined
        dL = 1;
        [X] = -L:dL:L;
        if nu ~= 0
            h = gaussian_kernel(D, nu, X); % using the gaussian kernel
            h = h / sqrt(sum((h(:).^2)));

            % Convolution shortcut (we only need it near the origin)
            %tic
            rhot2 = convn(h, h, 'same'); % convolution of the gaussian kernel
            %toc
            rho0 = rhot2(L+1); % L+1 is the center of [-L, L]
            rho1 = rhot2(L+2); % get the rho of two adjacent voxels
            rho2 = rhot2(L+3); % get the rho of two voxels with eucledian distance 2
            discov = [rho0 rho1 rho2];
        else
            discov = [1 0 0]; % independent case
        end
        
    case 2
        L = max(2,ceil(5*nu)); % how is the support defined
        dL = 1;
        [X, Y] = meshgrid( (-L:dL:L) , ...
                           (-L:dL:L ));
        if nu ~= 0
            h = gaussian_kernel(D, nu, sqrt((X.^2 + Y.^2))); 
            h = h / sqrt(sum((h(:).^2)));

            % Convolution shortcut (we only need it near the origin)
            %tic
            rhot2 = convn(h, h, 'same');
            %toc
            rho0 = rhot2(L+1, L+1);
            rho1 = rhot2(L+1, L+2);
            rho2 = rhot2(L+2, L+2);
            rho4 = rhot2(L+1, L+3);
            rho5 = rhot2(L+2, L+3);
            rho8 = rhot2(L+3, L+3);
            discov = [rho0 rho1 rho2 rho4 rho5 rho8];
        else
            discov = [1 0 0 0 0 0]; %indepdent case
        end
        
    case 3
        L = max(2,ceil(5*nu)); % how is the support defined
        dL = 1;
        [X, Y, Z] = meshgrid( (-L:dL:L) , ...
                              (-L:dL:L ) , ...
                              (-L:dL:L ));  
        if nu ~= 0
            h = gaussian_kernel(D, nu, sqrt((X.^2 + Y.^2 + Z.^2))); 
            h = h / sqrt(sum((h(:).^2)));

            % Convolution shortcut (we only need it near the origin)
            %tic
            rhot2 = convn(h, h, 'same');
            %toc
            rho0 = rhot2(L+1, L+1, L+1);
            rho1 = rhot2(L+1, L+2, L+1);
            rho2 = rhot2(L+1, L+2, L+2);
            rho3 = rhot2(L+2, L+2, L+2);
            rho4 = rhot2(L+1, L+3, L+1);
            rho5 = rhot2(L+1, L+2, L+3);
            rho6 = rhot2(L+2, L+2, L+3);
            rho8 = rhot2(L+1, L+3, L+3);
            rho9 = rhot2(L+2, L+3, L+3);
            rho12 = rhot2(L+3, L+3, L+3);
            discov = [rho0 rho1 rho2 rho3 rho4 rho5 rho6 rho8 rho9 rho12];
        else
            discov = [1 0 0 0 0 0 0 0 0 0];
        end
end


