function y=gsp_BLSCRP_Analytic(Un,mu,x0,y0,param)
%   Reconstruction of graph signal for BLSCRP by analytic way
%   Usage:  y=gsp_BLSCRP_Analytic(Un,mu,x0,y0,param)
%            
%   Input parameters :
%         Un         : Fourier basis (u1,u2,u3,...,un), a N*n matrix
%         mu         : eigenvalues (\lambda_1,\lambda_2,...,\lambda_n), a n*1 vector
%         x0         : Locations of known labeled points
%         y0         : Labels of known labeled points
%         param      : Optional parameters
%   Output parameters:
%         y          : Reconstructed graph signal
 
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


B=Un(x0,:);
 
C=B'*B+param.gamma*diag(mu);

xi=B'*y0;

x=inv(C)*xi;
 
y=Un*x;
