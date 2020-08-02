% Loads controls and feedback gains from DDP outputs
% Tests on noisy trajectories with/without feedback gains

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

% load results from DDP
load('results.mat')
g0 = results.g0;
xi0 = results.xi0;
u_ddp = results.u_ddp;
S_ddp = results.S_ddp;
g_des = results.g_des;
xi_des = results.xi_des;
dt = results.dt;
horizon = results.horizon;
t = results.t;
g_ddp = results.g_ddp;
xi_ddp = results.xi_ddp;

I = diag([10; 11.1; 13]);  % moment of inertia matrix
n_gl = 3;  % dimension of GL group
n_xi = 3;  % dimension of the Lie algebra
n_x = n_gl + n_xi;
n_u = 3;  % dimension of controls

par_dyn = struct('dt', dt, 'I', I, 'n_gl', n_gl, 'n_xi', n_xi, 'n_u', n_u, ...
                 'n_x', n_gl + n_xi, 'vmap_fun', @(z) vmap_so3(z));

f_dyn = @(g, xi, u) eulerdyn_TSO3(g, xi, u, par_dyn);

%% without gains

t0 = t(1); tf = t(end); quat_des = rotm2quat(g_des);
magn = 0.001; %magnitude for noise variance
N = 90;  %number of samples

ge = zeros(n_gl,n_gl,horizon,N); %open loop
xie = zeros(n_xi,horizon,N);
quate = zeros(4,horizon,N);

for i=1:N
    ge_k = g0; ge(:,:,1,i) = ge_k; quate(:,1,i) = rotm2quat(ge(:,:,1,i))';
    xie_k = xi0; xie(:,1,i) = xie_k;
    for k=1:horizon-1
        u_k = u_ddp(:,k);
        [ge_k,xie_k] = f_dyn(ge_k,xie_k,u_k);
        xie_k = xie_k + magn*randn(n_xi,1);
        ge(:,:,k + 1,i) = ge_k;
        quate(:,k + 1,i) = rotm2quat(ge(:,:,k + 1,i))';
        xie(:,k + 1,i) = xie_k;
    end   
    figure(1); hold on
    plot(t,xie(1,:,i),'g',t,xie(2,:,i),'k',t,xie(3,:,i),'r','linewidth',1.9)
    figure(3); hold on
    plot(t,quate(1,:,i),'g',t,quate(2,:,i),'c',t,quate(3,:,i),'k',t,quate(4,:,i),'r','linewidth',1.9)
end

figure(1); plot(t(end),xi_des(3),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Sampled body-fixed velocities','fontsize',13,'fontweight','bold');
title('Open loop')

figure(3); plot(t(end),quat_des(1),'go',t(end),quat_des(2),'co',t(end),quat_des(3),'ko',t(end),quat_des(4),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Sampled unit quaternions','fontsize',13,'fontweight','bold');
title('Open loop')

%% with gains

gc = zeros(n_gl,n_gl,horizon,N); %closed loop
xic = zeros(n_xi,horizon,N);
quatc = zeros(4,horizon,N);

for i=1:N
    gc_k = g0; gc(:,:,1,i) = gc_k; quatc(:,1,i) = rotm2quat(gc(:,:,1,i))';
    xic_k = xi0; xic(:,1,i) = xic_k;
    deltax_k = zeros(n_x,1);
    for k=1:horizon-1
        u_k = u_ddp(:,k) + S_ddp(:,:,k)*deltax_k;
        [gc_k,xic_k] = f_dyn(gc_k,xic_k,u_k);
        xic_k = xic_k + magn*randn(n_xi,1);
        gc(:,:,k + 1,i) = gc_k;
        quatc(:,k + 1,i) = rotm2quat(gc(:,:,k + 1,i))';
        xic(:,k + 1,i) = xic_k;
        deltax_k(1:n_gl) = vmap_so3(real(logm(g_ddp(:,:,k + 1)'*gc(:,:,k + 1,i))));
        deltax_k(n_gl + 1:end) = xic(:,k + 1,i) - xi_ddp(:,k + 1);
    end   
    figure(2); hold on
    plot(t,xic(1,:,i),'g',t,xic(2,:,i),'k',t,xic(3,:,i),'r','linewidth',1.9)
    figure(4); hold on
    plot(t,quatc(1,:,i),'g',t,quatc(2,:,i),'c',t,quatc(3,:,i),'k',t,quatc(4,:,i),'r','linewidth',1.9)
end

figure(2); plot(t(end),xi_des(3),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Sampled body-fixed velocities','fontsize',13,'fontweight','bold');
title('Closed loop')

figure(4); plot(t(end),quat_des(1),'go',t(end),quat_des(2),'co',t(end),quat_des(3),'ko',t(end),quat_des(4),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Sampled unit quaternions','fontsize',13,'fontweight','bold');
title('Closed loop')