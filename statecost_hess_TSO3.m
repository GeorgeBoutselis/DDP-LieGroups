function [f_x, f_xx] = statecost_hess_TSO3(g, xi, K_g, K_xi, g_des, xi_des, par_dyn)
% Gradients and hessians of state cost on TSO(3)
% Matrix-vector representation

% Inputs:
%   g:          rotation matrix (3, 3)
%   xi:         body-fixed velocity (3, 1)
%   K_g, K_xi:  p.d. weighting matrices
%   g_des:      desired rotation matrix
%   xi_des:     desired body-fixed velocity
%   par_dyn:    struct with dynamics parameters

% Outputs:
%   f_x:        gradient wrt state (6, 1)
%   f_xx:       hessian wrt state (6, 6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% author: George Boutselis, July 2020 - MATLAB version: R2016a

% If you are using this code, please consider citing our paper:

% George I. Boutselis and Evangelos Theodorou,
% "Discrete-time Differential Dynamic Programming on Lie Groups: 
%  Derivation, Convergence Analysis and Numerical Results", 
% Transactions on Automatic Control (TAC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Copyright (c) 2020 George Boutselis
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_gl = par_dyn.n_gl;
n_x = par_dyn.n_x;
f_xx = zeros(n_x);

f_gskew = (K_g*g_des'*g - g'*g_des*K_g);
f_g = vmap_so3(f_gskew);
f_xi = K_xi*(xi - xi_des);
f_x = [f_g;f_xi];

f_gg = trace(K_g*g_des'*g)*eye(3) - K_g*g_des'*g; f_gg = 0.5*(f_gg + f_gg');
f_xixi = K_xi;
f_xx(1:n_gl,1:n_gl) = f_gg;
f_xx(n_gl + 1:end,n_gl + 1:end) = f_xixi;

end