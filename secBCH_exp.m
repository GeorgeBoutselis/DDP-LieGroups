function [jac_zz, jac_zu, jac_uu] = secBCH_exp(u, dt)
% Second order terms of BCH expansion for euler kinematics on TSO(3)

% Inputs:
%   u:          vector (3, 1)
%   dt:         time step of discrete dynamics

% Outputs:
%   jac_zz:     hessian wrt state (3, 3, 3)
%   jac_zu:     hessian wrt state/control (3, 3, 3)
%   jac_uu:     hessian wrt control (3, 3, 3)

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

uhat = hat_so3(u);
e = eye(3);
E = zeros(3,3,3);  % basis for so(3) lie algebra
for i=1:3, E(:,:,i) = hat_so3(e(:,i)); end

uhat_bar = uhat;
L_z = inv_dexp_m(dt*uhat_bar);
L_u = e*dt;

L_zz = zeros(3,3,3);
L_zu = zeros(3,3,3);
for i=1:3
    for j=1:3
        L_zz(i,j,:) = -1/12*dt*ad_m(E(:,:,i))*ad_m(uhat_bar)*e(:,j) + ...
            1/24*dt^2*ad_m(E(:,:,i))*ad_m(uhat_bar)^2*e(:,j);
        L_zu(i,j,:) = 1/2*dt*ad_m(E(:,:,i))*e(:,j) - ...
            1/12*dt^2*ad_m(uhat_bar)*ad_m(E(:,:,i))*e(:,j) + ...
            1/12*dt^2*ad_m(E(:,:,j))*ad_m(uhat_bar)*e(:,i);
    end
end
for k=1:3, L_zz(:,:,k) = L_zz(:,:,k) + L_zz(:,:,k)'; end

G_zz = zeros(3,3,3); G_zu = zeros(3,3,3); G_zu1 = zeros(3,3,3);
G_zz1 = zeros(3,3,3); G_uu1 = zeros(3,3,3);
for i=1:3
    for j=1:3
        G_zz(i,j,:) = dexp_m(-dt*uhat_bar)*reshape(L_zz(i,j,:),3,1,1); 
        G_zu(i,j,:) = dexp_m(-dt*uhat_bar)*reshape(L_zu(i,j,:),3,1,1);
        G_zz1(i,j,:) = (1/12*e - 1/24*ad_m(dt*uhat_bar))*ad_m(hat_so3(L_z(:,i)))*ad_m(dt*uhat_bar)*L_z(:,j);
        G_zu1(i,j,:) = (1/12*e - 1/24*ad_m(dt*uhat_bar))*(ad_m(hat_so3(L_z(:,i)))*ad_m(dt*uhat_bar)*L_u(:,j) +...
                        ad_m(hat_so3(L_u(:,j)))*ad_m(dt*uhat_bar)*L_z(:,i));
        G_uu1(i,j,:) = (1/12*e - 1/24*ad_m(dt*uhat_bar))*ad_m(hat_so3(L_u(:,i)))*ad_m(dt*uhat_bar)*L_u(:,j);
    end
end
for k=1:3
    G_zz1(:,:,k) = G_zz1(:,:,k) + G_zz1(:,:,k)';
    G_uu1(:,:,k) = G_uu1(:,:,k) + G_uu1(:,:,k)';
end

jac_zz = G_zz + G_zz1;
jac_zu = G_zu + G_zu1;
jac_uu = G_uu1;

end


function k = ad_m(X)
% "ad()" matrix representation for so(3)
% x is skew symmetric

k = X;

end

function z = dexp_m(x)
% matrix representation for "dexp()" for so(3)
% x is (3, 1)
    
z = eye(3) + 0.5 * ad_m(x) + 0.1667 * ad_m(x) ^ 2;

end

function z = inv_dexp_m(x)
% matrix representation for inverse of "dexp()" for so(3)
% x is (3, 1)

z = eye(3) - 0.5 * ad_m(x) + 0.0833 * ad_m(x) ^ 2;

end