% Discrete DDP on Lie groups
% Application on simple model on TSO(3)
% Runs algorithm for different orders of linearization

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

%% parameters

dt = 0.01;  % time step
t0 = 0;  % initial time
tf = 3;  % end time
t = t0:dt:tf;
horizon = length(t);

I = diag([10; 11.1; 13]);  % moment of inertia matrix
n_gl = 3;  % dimension of GL group
n_xi = 3;  % dimension of the Lie algebra
g0 = eye(n_gl);  % initial rotation matrix
g_des = rot_axis(30, 'x') * rot_axis(70, 'z');  % desired rotation matrix
xi0 = zeros(n_xi, 1);  % initial body-fixed velocity
xi_des = zeros(n_xi, 1);  % desired body-fixed velocity
n_u = 3;  % dimension of controls

par_dyn = struct('dt', dt, 'I', I, 'n_gl', n_gl, 'n_xi', n_xi, 'n_u', n_u, ...
                 'n_x', n_gl + n_xi, 'vmap_fun', @(z) vmap_so3(z));

Lambda_u = 0.1 * eye(n_u) * dt;
K_g = 10000 * eye(n_gl);
K_xi = 10000 * eye(n_xi);

f_dyn = @(g, xi, u) eulerdyn_TSO3(g, xi, u, par_dyn);
hess_dyn = @(g, xi, u, order) lineardyn_TSO3eul(g, xi, u, order, par_dyn);
L_cost = @(u) Lcost_quad_u(u, Lambda_u);
Lhess_cost = @(u) LcostHess_quad_u(u, Lambda_u, par_dyn);
F_final_cost = @(g, xi) statecost_TSO3(g, xi, K_g, K_xi, g_des, xi_des);
Fhess_final_cost = @(g, xi) statecost_hess_TSO3(g, xi, K_g, K_xi, g_des, xi_des, par_dyn);
                 
par_ddp = struct('horizon', horizon, 'g0', g0, 'xi0', xi0, 'g_des', g_des, ...
                 'xi_des', xi_des, 'L_cost', L_cost, 'F_final_cost', F_final_cost, ...
                 'Lhess_cost', Lhess_cost, 'Fhess_final_cost', Fhess_final_cost, ...
                 'f_dyn', f_dyn, 'hess_dyn', hess_dyn);

%% ddp

disp('1st + 2nd order dynamics approximations:')

options_d.iter = 300;  % maximum number of iterations
options_d.gamma0 = 1;  % initial line search parameter
options_d.gammachange = 1/3;  % decrease ratio for line search condition
options_d.convtol = 1e-8;  % convergence criterion for change in cost function
options_d.k = 0.001;  % used for accepting a new control
options_d.rho0 = 0;  % initial value for regularization of Hessian
options_d.rhochange0 = 1.83;  % initial ratio for changing rho
options_d.rho_max = 1e10;
options_d.rho_min = 1e-6;
options_d.order = 1;  % initial order
options_d.orderchange = 1;  % when 1, order will change
options_d.orderchange_tol = 0.1;  % tolerance for switching orders

u_bar = zeros(n_u, horizon - 1);  % nominal control sequence

[g_ddp, xi_ddp, u_ddp, S_ddp, u_bars, J] = ...
                            ddpLie_fun(u_bar,par_dyn,par_ddp,options_d);
quat_ddp =  rotm2quat(g_ddp)';
quat_des = rotm2quat(g_des);
disp('........................................................')

%% only 1st order

disp('1st order dynamics approximations:')

options_d.orderchange = 0;
u_bar = zeros(n_u, horizon - 1);  % nominal control sequence

[g_ddp1, xi_ddp1, u_ddp1, S_ddp1, u_bars1, J1] = ...
                            ddpLie_fun(u_bar,par_dyn,par_ddp,options_d);
quat_ddp1 =  rotm2quat(g_ddp1)';
disp('........................................................')

%% only 2nd order

disp('2nd order dynamics approximations:')

options_d.order = 2;
u_bar = zeros(n_u, horizon - 1);  % nominal control sequence

[g_ddp2, xi_ddp2, u_ddp2, S_ddp2, u_bars2, J2] = ...
                            ddpLie_fun(u_bar,par_dyn,par_ddp,options_d);
quat_ddp2 =  rotm2quat(g_ddp2)';
disp('........................................................')

%% plots

d = zeros(length(J), 1); for i=1:size(u_bars, 3), d(i) = norm(u_bars(:, :, i) - u_ddp); end
d1 = zeros(length(J1), 1); for i=1:size(u_bars1, 3), d1(i) = norm(u_bars1(:, :, i) - u_ddp1); end
d2 = zeros(length(J2), 1); for i=1:size(u_bars2, 3), d2(i) = norm(u_bars2(:, :, i) - u_ddp2); end

figure();
subplot(1,2,1)
semilogy(0:length(d1) - 1,d1,'-go','linewidth',1.9); hold on
semilogy(0:length(d2) - 1,d2,'-ro','linewidth',1.9)
semilogy(0:length(d) - 1,d,'--ko','linewidth',1.9)
set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w'); lg = legend('1^{st} order','2^{nd} order','1^{st} + 2^{nd} order');
lg.FontSize = 11;
xlabel('iteration number','fontsize',13,'fontweight','bold');
ylabel('\hbox{\fontsize{17}{4}\selectfont{$||U-U_*||$}}','interpreter','latex'); grid on

subplot(1,2,2)
loglog(d(1:end-1),d(2:end),'--ko',d1(1:end-1),d1(2:end),'-go',d2(1:end-1),d2(2:end),'-ro','linewidth',1.9)
grid on
xlabel('\hbox{\fontsize{17}{4}\selectfont{$||U_{(i)}-U_*||$}}','interpreter','latex')
ylabel('\hbox{\fontsize{17}{4}\selectfont{$||U_{(i+1)}-U_*||$}}','interpreter','latex')
set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');

figure();
plot(t,quat_ddp(1,:),'g',t,quat_ddp(2,:),'c',t,quat_ddp(3,:),'k',t,quat_ddp(4,:),'r',...
    t(end),quat_des(1),'go',t(end),quat_des(2),'co',t(end),quat_des(3),'ko',t(end),quat_des(4),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Unit quaternions','fontsize',13,'fontweight','bold');

figure();
plot(t,xi_ddp(1,:),'g',t,xi_ddp(2,:),'k',t,xi_ddp(3,:),'r',t(end),xi_des(3),'ro','markersize',8,'linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Body-fixed velocities','fontsize',13,'fontweight','bold');

figure();
plot(t(1:end-1),u_ddp(1,:),'g',t(1:end-1),u_ddp(2,:),'k',t(1:end-1),u_ddp(3,:),'r','linewidth',1.9)
axis([t0 tf -inf inf]); grid on; set(gca,'FontWeight','bold'); set(gcf, 'Color', 'w');
xlabel('time','fontsize',13,'fontweight','bold'); ylabel('Controls','fontsize',13,'fontweight','bold');