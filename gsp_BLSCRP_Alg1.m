function y=gsp_BLSCRP_Alg1(Un,mu,x0,y0,param)
%   Reconstruction of graph signal for BLSCRP by Algorithm 1
%   Usage:  y=gsp_BLSCRP_Alg1(Un,mu,x0,y0,param)
%            
%   Input parameters :
%         Un         : Fourier basis (u2,u3,...,un), a N*(n-1) matrix
%         mu         : eigenvalues (\lambda_2,...,\lambda_n), a (n-1)*1 vector
%         x0         : Locations of known labeled points
%         y0         : Labels of known labeled points
%         param      : Optional parameters
%   Output parameters:
%         y          : Reconstructed graph signal
%
%   Example:
%         % Give N,n, generate a graph and calculate the eigendecomposition.
%         N = 500;
%         n = 0.1*N;
%         G = gsp_sensor(N);
%         G = gsp_compute_fourier_basis(G);
% 
%         % Construct an original signal with the original bandwidth 0.1*N
%         cutoffcoeff = rand(n,1);
%         cutoffcoeff = sort(cutoffcoeff,'descend');
%         f = G.U(:,1:n)*cutoffcoeff;  
% 
%        % Locations and labels of known labeled points, parameters
%        p = randperm (G.N);
%        labs = fix(0.2*G.N);
%        x0 = p(1:labs);
%        y0 = f(x0);
%        param.gamma = 0.1;
%        Un=G.U(:,2:n);
%        mu=G.e(2:n);
% 
%        % Reconstruction and mean absolute error.
%        y = gsp_BLSCRP_Alg1(G,x0,y0,param);
%        err = mean(abs(f-y));        
%       
%
%   Additional parameters
%   ---------------------
%    param.gamma  : The regularization parameter used to balance the weight 
%                   between the error and the quality index
% 
%   References:
%     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
%     Algorithm for Recovery of Bandlimited Graph Signals Based on the 
%     Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
%     Processing, revised vesion, T-SP-25570-2019, 2020.7.     
%
% Copyright (C) 2020 Qian Zhang and Lihua Yang. All rights reserved.
%
% This function is free software for scientific research, you can 
% redistribute it and/or modify it. But not for commercial use.
%
% This function is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   

% This function is suggested to use in the context of the GSPBOX toolbox,
% when you use this algorithm, please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
%     and 
%     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
%     Algorithm for Recovery of Bandlimited Graph Signals Based on the 
%     Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
%     Processing, 2020. 
%     and if needed
%     J. Paratte and L. Martin. Fast eigenspace approximation using random
%     signals. arXiv preprint arXiv:1611.00938, 2016.

% Author:   Qian Zhang and Lihua Yang
% Date:     1 July 2020
% Version:  2020.7.01
ell=numel(x0);
 
B=Un(x0,:);
 
temp=repmat(1./mu,1,size(B,1));

temp=temp.*B';

Tl=B*temp;

temp2=Tl*ones(ell,1);

GG=param.gamma*Tl+Tl*Tl-1/ell*(temp2*temp2');

d=Tl*y0-1/ell*temp2*(ones(1,ell)*y0);

xi=inv(GG)*d;
 
g=Un*(temp*xi);
 
temp3=- sum(g(x0)-y0)/ell;

y=  temp3+g;
