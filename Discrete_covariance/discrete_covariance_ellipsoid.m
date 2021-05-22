function discov = discrete_covariance_ellipsoid(nu, D)
%__________________________________________________________________________
% Computes discrete version of the covariance function.
% 
% Input:
%   nu   -    Standard deviation of the D-dim gaussian smoothing kernel in
%             each direction, length D vector
%   D    -    The dimension of the kernel, at least 2.
%
% Output:
%   discov:   -    Array of dimension D convariance matrix
%__________________________________________________________________________

switch D
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
            %discov = [rho0, rho1v, rho2v, rho1h, rho3, rho4, rho2h, rho5,...
            %          rho6, rho7, rho8, rho9, rho10];
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
            %discov = [1 0 0 0 0 0]; %indepdent case
            discov = eye(3^D);
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

            %together 125 possibilities when estimating E(XY) for 3D X,Y
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
%             t = num2cell(rhot2((L+1):(L+3), L+1, L+1));
%             [rho0, rho1, rho2] = deal(t{:});
%             t = num2cell(rhot2(L+1, (L+2):(L+3), L+1));
%             [rho3, rho6] = deal(t{:});
%             t = num2cell(rhot2(L+2, (L+2):(L+3), L+1));
%             [rho4, rho7] = deal(t{:});
%             t = num2cell(rhot2(L+3, (L+2):(L+3), L+1));
%             [rho5, rho8] = deal(t{:});
%             t = num2cell(rhot2(L+1, L+1, (L+2):(L+3)));
%             [rho9, rho18] = deal(t{:});
%             t = num2cell(rhot2(L+2, L+1, (L+2):(L+3)));
%             [rho10, rho19] = deal(t{:});
%             t = num2cell(rhot2(L+3, L+1, (L+2):(L+3)));
%             [rho11, rho20] = deal(t{:});
%             t = num2cell(rhot2(L+1, L+2, (L+2):(L+3)));
%             [rho12, rho21] = deal(t{:});
%             t = num2cell(rhot2(L+2, L+2, (L+2):(L+3)));
%             [rho13, rho22] = deal(t{:});
%             t = num2cell(rhot2(L+3, L+2, (L+2):(L+3)));
%             [rho14, rho23] = deal(t{:});
%             t = num2cell(rhot2(L+1, L+3, (L+2):(L+3)));
%             [rho15, rho24] = deal(t{:});
%             t = num2cell(rhot2(L+2, L+3, (L+2):(L+3)));
%             [rho16, rho25] = deal(t{:});
%             t = num2cell(rhot2(L+3, L+3, (L+2):(L+3)));
%             [rho17, rho26] = deal(t{:});
%             t = num2cell(rhot2(L, (L+2):(L+3), L+1));
%             [rho27, rho28] = deal(t{:});
%             t = num2cell(rhot2(L, L+1, (L+2):(L+3)));
%             [rho29, rho32] = deal(t{:});
%             t = num2cell(rhot2(L, L+2, (L+2):(L+3)));
%             [rho30, rho33] = deal(t{:});
%             t = num2cell(rhot2(L, L+3, (L+2):(L+3)));
%             [rho31, rho34] = deal(t{:});
%             t = num2cell(rhot2(L-1, (L+2):(L+3), L+1));
%             [rho35, rho36] = deal(t{:});
%             t = num2cell(rhot2(L-1, L+1, (L+2):(L+3)));
%             [rho37, rho40] = deal(t{:});
%             t = num2cell(rhot2(L-1, L+2, (L+2):(L+3)));
%             [rho38, rho41] = deal(t{:});
%             t = num2cell(rhot2(L-1, L+3, (L+2):(L+3)));
%             [rho39, rho42] = deal(t{:});
%             t = num2cell(rhot2(L+1, L, (L+2):(L+3)));
%             [rho43, rho46] = deal(t{:});
%             t = num2cell(rhot2(L+2, L, (L+2):(L+3)));
%             [rho44, rho47] = deal(t{:});
%             t = num2cell(rhot2(L+3, L, (L+2):(L+3)));
%             [rho45, rho48] = deal(t{:});
%             t = num2cell(rhot2(L, L, (L+2):(L+3)));
%             [rho49, rho50] = deal(t{:});
%             t = num2cell(rhot2(L-1, L, (L+2):(L+3)));
%             [rho51, rho52] = deal(t{:});
%             t = num2cell(rhot2(L+1, L-1, (L+2):(L+3)));
%             [rho53, rho56] = deal(t{:});
%             t = num2cell(rhot2(L+1, L-1, (L+2):(L+3)));
%             [rho54, rho57] = deal(t{:});
%             t = num2cell(rhot2(L+2, L-1, (L+2):(L+3)));
%             [rho55, rho58] = deal(t{:});
%             t = num2cell(rhot2(L, L-1, (L+2):(L+3)));
%             [rho59, rho60] = deal(t{:});
%             t = num2cell(rhot2(L+3, L-1, (L+2):(L+3)));
%             [rho61, rho62] = deal(t{:});
%             
%             a1 = toeplitz([rho0, rho1, rho2]);
%             a2 = toeplitz([rho3, rho4, rho5], [rho3, rho27, rho35]);
%             a3 = toeplitz([rho6, rho7, rho8], [rho6, rho28, rho36]);
%             a4 = a2';
%             a5 = a3';
%             block_a = {a1,a2,a3,a4,a5};
%             c = 1:3;
%             r = [1,4,5];
%             a = cell2mat(block_a(toeplitz(c,r)));
%             b1 = toeplitz([rho9, rho10, rho11], [rho9, rho29, rho37]);
%             b2 = toeplitz([rho12, rho13, rho14], [rho12, rho30, rho38]);
%             b3 = toeplitz([rho15, rho16, rho17], [rho15, rho31, rho39]);
%             b4 = toeplitz([rho43, rho44, rho45], [rho43, rho49, rho51]);
%             b5 = toeplitz([rho53, rho54, rho55], [rho53, rho59, rho61]);
%             block_b = {b1,b2,b3,b4,b5};
%             b = cell2mat(block_b(toeplitz(c,r)));
%             d1 = toeplitz([rho18, rho19, rho20], [rho18, rho32, rho40]);
%             d2 = toeplitz([rho21, rho22, rho23], [rho21, rho33, rho41]);
%             d3 = toeplitz([rho24, rho25, rho26], [rho24, rho34, rho42]);
%             d4 = toeplitz([rho46, rho47, rho48], [rho46, rho50, rho52]);
%             d5 = toeplitz([rho56, rho57, rho58], [rho56, rho60, rho62]);
%             block_d = {d1,d2,d3,d4,d5};
%             d = cell2mat(block_d(toeplitz(c,r)));
%             block = {a,b,d,b',d'};
%             discov = cell2mat(block(toeplitz(c,r)));
        else
            discov = eye(3^D);
        end
end


