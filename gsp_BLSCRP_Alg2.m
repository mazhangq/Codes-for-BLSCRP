function [Un,D]=gsp_BLSCRP_Alg2(G,n,param)
%   Approximate first n eigenvectors and eigenvalues for BLSCRP by Algorithm 2
%   Usage:  [U,D]=gsp_BLSCRP_Alg2(G,n,param)
%            
%   Input parameters :
%         G          : Graph structure
%         n          : bandwidth
%         param      : Optional parameters
%   Output parameters:
%         U          : a N*n matrix consists of first n eigenvectors
%         D          : a n*1 vector consists of  first n eigenvalues

%   Additional parameters
%   ---------------------
%    param.order  : The regularization parameter used to balance the weight 
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

Bn = gsp_eigenspace_estimation_simple(G,n,param);
 
An=Bn'*(G.L*Bn);

[S,D]=eig(An);

[D,index] = sort(diag(D),'ascend');

S=S(:,index);

Un=Bn*S;
