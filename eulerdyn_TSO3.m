function [g_kplus1, xi_kplus1] = eulerdyn_TSO3(g_k, xi_k, u_k, par_dyn)
% Discrete Euler dynamics on TSO(3)

% Inputs:
%   g_k:        rotation matrix at t_k (3, 3)
%   xi_k:       body-fixed velocity at t_k (3, 1)
%   u_k:        control vector at t_k (3, 1)
%   par_dyn:    struct with dynamics parameters

% Outputs:
%   g_kplus1:   rotation matrix at t_kplus1 (3, 3)
%   xi_kplus1:  body-fixed velocity at t_kplus1 (3, 1)

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

dt = par_dyn.dt;
I = par_dyn.I;

g_kplus1 = g_k * expm(dt * hat_so3(xi_k));
xi_kplus1 = xi_k + dt * (I \ (cross(I * xi_k, xi_k) + u_k));

end