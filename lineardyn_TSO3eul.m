function [Phi, B, fxx, fxu, fuu] = lineardyn_TSO3eul(g, xi, u, order_lin, par_dyn)
% linearized dynamics for system on TSO(3)

% Inputs:
%   g:          pose (n_gl, n_gl)
%   xi:         body-fixed velocity (n_xi, 1)
%   u:          control vector (n_u, 1)
%   order_lin:  order of linearization matrices for outputs
%   par_dyn:    struct with dynamics parameters

% Outputs:
%   Phi:        jacobian wrt state (n_gl+n_xi, n_gl+n_xi)
%   B:          jacobian wrt control (n_gl+n_xi, n_u)
%   fxx:        hessian wrt state (n_gl+n_xi, n_gl+n_xi, n_gl+n_xi)
%   fxu:        hessian wrt state/control (n_gl+n_xi, n_u, n_gl+n_xi)
%   fuu:        hessian wrt control (n_u, n_u, n_gl+n_xi)

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

I = par_dyn.I;
dt = par_dyn.dt;
n_gl = par_dyn.n_gl;
n_xi = par_dyn.n_xi;
n_x = n_gl + n_xi;
n_u = par_dyn.n_u;

I1 = I(1,1); I2 = I(2,2); I3 = I(3,3);
xi1 = xi(1); xi2 = xi(2); xi3 = xi(3);
Phi = zeros(n_x);
B = zeros(n_x, n_u);

Phi(1:n_gl,1:n_gl) = Ad_m(expm(-dt*hat_so3(xi)));
Phi(1:n_gl,n_gl + 1:end) = dt*dexp_m(-dt*xi);
Phi(n_gl + 1:end,n_gl + 1:end) = [1,dt.*I1.^(-1).*(I2 + (-1).*I3).*xi3,dt.*I1.^(-1).*(I2 + (-1).* ...
  I3).*xi2;dt.*I2.^(-1).*((-1).*I1 + I3).*xi3,1,dt.*I2.^(-1).*((-1) ...
  .*I1 + I3).*xi1;dt.*(I1 + (-1).*I2).*I3.^(-1).*xi2,dt.*(I1 + (-1) ...
  .*I2).*I3.^(-1).*xi1,1];

B(n_gl + 1:end,:) = [dt.*I1.^(-1),0,0;0,dt.*I2.^(-1),0;0,0,dt.*I3.^(-1)];

fxx = zeros(n_x,n_x,n_x);
fxu = zeros(n_x,n_u,n_x);
fuu = zeros(n_u,n_u,n_x);

if order_lin==2
    [jac_zz, jac_zxi, jac_xixi] = secBCH_exp(xi, dt);
    jac_xiz = zeros(3,3,3);
    for k=1:3, jac_xiz(:,:,k) = jac_zxi(:,:,k)'; end
    G_sec1 = cat(2,jac_zz,jac_zxi); G_sec2 = cat(2,jac_xiz,jac_xixi); G_sec = cat(1,G_sec1,G_sec2);
    
    fxx(:,:,1:n_gl) = G_sec;
    fxx(n_gl + 1:end,n_gl + 1:end,n_gl + 1) = [0,0,0;0,0,dt.*I1.^(-1).*(I2 + (-1).*I3);0,dt.*I1.^(-1).*(I2 + ( ...
        -1).*I3),0];
    fxx(n_gl + 1:end,n_gl + 1:end,n_gl + 2) = [0,0,dt.*I2.^(-1).*((-1).*I1 + I3);0,0,0;dt.*I2.^(-1).*((-1).*I1 + ...
        I3),0,0];
    fxx(n_gl + 1:end,n_gl + 1:end,n_gl + 3) = [0,dt.*(I1 + (-1).*I2).*I3.^(-1),0;dt.*(I1 + (-1).*I2).*I3.^(-1), ...
        0,0;0,0,0];
end

end

function Z = Ad_m(X)
% "Ad()" matrix representation for SO(3)
% X is rotation matrix

Z = X;

end

function z = ad_m(x)
% "ad()" matrix representation for so(3)
% x is (3, 1)

z = hat_so3(x);

end

function z = dexp_m(x)
% matrix representation for "dexp()" for so(3)
% x is (3, 1)
    
z = eye(3) + 0.5 * ad_m(x) + 0.1667 * ad_m(x) ^ 2;

end