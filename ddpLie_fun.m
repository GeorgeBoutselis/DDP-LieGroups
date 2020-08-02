function [g, xi, u, S, u_bars, J] = ddpLie_fun(u_bar, par_dyn, par_ddp, options_d)
% Main function for discrete DDP algorithm on Lie groups

% Inputs:
%   u_bar:      initial control trajectory
%   par_dyn:    struct with dynamics parameters (state dimensionality, etc.)
%   par_ddp:    struct with parameters of the algorithm (cost functions, etc.)
%   options_d:  struct with hyperparameters (initial regularizer, etc.)

% Outputs:
%   g:          obtained pose trajectory (n_gl, n_gl, horizon)
%   xi:         obtained body-fixed velocity trajectory (n_xi, horizon)
%   u:          obtained control sequence (n_u, horizon - 1)
%   S:          obtained feedback gains (n_x, n_x, horizon - 1)
%   u_bars:     control sequences per iteration (n_u, horizon - 1, #iterations)
%   J:          cost per iteration (1, #iterations)

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

g0 = par_ddp.g0;
xi0 = par_ddp.xi0;
horizon = par_ddp.horizon;
f_dyn = par_ddp.f_dyn;
hess_dyn = par_ddp.hess_dyn;
n_gl = par_dyn.n_gl;
n_xi = par_dyn.n_xi;
n_x = par_dyn.n_x;
n_u = par_dyn.n_u;
L_cost = par_ddp.L_cost;
Lhess_cost = par_ddp.Lhess_cost;
F_final_cost = par_ddp.F_final_cost;
Fhess_final_cost = par_ddp.Fhess_final_cost;
vmap_fun = par_dyn.vmap_fun;

J = zeros(1, options_d.iter);  % costs for all iterations
u_bars = zeros(n_u, horizon - 1, options_d.iter);  % control trajectories per iteration

% get nominal trajectory
g = zeros(n_gl, n_gl, horizon); g(:, :, 1) = g0;
xi = zeros(n_xi, horizon); xi(:, 1) = xi0;
c_L = 0;
for k=1:horizon-1
    c_L = c_L + L_cost(u_bar(:, k));
    [g(:, :, k + 1), xi(:, k + 1)] = f_dyn(g(:, :, k), xi(:, k), u_bar(:, k));   
end
c_F = F_final_cost(g(:, :, end), xi(:, end));
cost_traj0 = c_L + c_F;
g_bar = g; xi_bar = xi;  % nominal state trajectories

%%%%%%% Geometric DDP %%%%%%%%
rho = options_d.rho0;  % Hessian regularizer
rho_change = options_d.rhochange0;  % rate of change of regularizer
order_changed = 0;  % if 1, order has changed
small_value = 1e-8;  % to check for nonpositive Quu's

for j=1:options_d.iter
    J(j) = cost_traj0; u_bars(:, :, j) = u_bar;
    fprintf(1, 'Iteration %d out of %d ---- total cost is %f \n', j - 1, options_d.iter - 1, J(j));
    Quu_is_pos = 0;  % if 1, Quu is positive definite for all instances
    
    if options_d.orderchange && ~order_changed && (cost_traj0 <= options_d.orderchange_tol * J(1))
        options_d.order = 2;  % change to second order approximations
        disp('order is changed');
        order_changed = 1;
    end
    
    while ~Quu_is_pos
        [V_x_kplus1, V_xx_kplus1] = Fhess_final_cost(g(:, :, end), xi(:, end));  % trivialized gradient and Hessian of value function at tf
        l = zeros(n_u, horizon - 1);  % feedforward control
        S = zeros(n_u, n_x, horizon - 1);  % feedback control gains
        Q_calc = 1;  % if 1, backpropagation is complete
        
        for k=horizon-1:-1:1         
            [Fi, B, fxx, fxu, fuu] = hess_dyn(g_bar(:, :, k), xi_bar(:, k), u_bar(:, k), options_d.order);  % trivialized linearization matrices
            [L_x, L_u, L_xu, L_uu, L_xx] = Lhess_cost(u_bar(:, k));  % trivialized gradients and Hessians of running cost
            temp0 = zeros(1, 1, n_x);
            temp0(:) = V_x_kplus1;
            
            % backpropagation of value function
            Qx = L_x + Fi' * V_x_kplus1;
            Qu = L_u + B' * V_x_kplus1;
            Qxx = L_xx + Fi' * V_xx_kplus1 * Fi + sum(bsxfun(@times, fxx, temp0), 3);
            Qxu = L_xu + Fi' * V_xx_kplus1 * B + sum(bsxfun(@times, fxu, temp0), 3);
            Qux = Qxu';
            Quu = L_uu + B' * V_xx_kplus1 * B + sum(bsxfun(@times, fuu, temp0), 3);
            Quu_reg = Quu + rho * eye(n_u);
            if min(eig(Quu_reg))<small_value  % Quu not positive definite
                if rho>=options_d.rho_max
                    error('Regularizer has become too large');
                end
                fprintf(1, 'Quu non positive definite at instance %d \n', k);
                rho_change = max(rho_change * options_d.rhochange0, options_d.rhochange0);  % update rate of change
                rho = max(rho * rho_change, options_d.rho_min);  % update regularizer
                Q_calc = 0;
                break            
            end
            l(:, k) = -Quu_reg \ Qu;
            S(:, :, k) = -Quu_reg \ Qux;
            V_x_kplus1 = Qx + S(:, :, k)' * Quu * l(:, k) + S(:, :, k)' *Qu + Qxu * l(:, k);
            V_xx_kplus1 = Qxx + S(:, :, k)' * Quu * S(:, :, k) + S(:, :, k)' * Qux + Qxu * S(:, :, k);
            V_xx_kplus1 = 0.5 * (V_xx_kplus1 + V_xx_kplus1');
        end
        if Q_calc, Quu_is_pos = 1; end
    end
    
    no_cost_reduction = 1;  % if 0, cost has sufficiently reduced for next iteration
    conv_flag = 0;  % if 1, algorithm has converged
    gamma = options_d.gamma0;  % line search parameter
    while no_cost_reduction
        % get new trajectory
        g = zeros(n_gl, n_gl, horizon); g(:, :, 1) = g0;
        xi = zeros(n_xi, horizon); xi(:,1) = xi0;
        u = zeros(n_u, horizon - 1);
        deltax_k = zeros(n_x, 1);  % trivialized state deviations
        c_L = 0;
        for k=1:horizon-1
            deltau_k = gamma * l(:,k) + S(:, :, k) * deltax_k;  % control perturbations
            u(:, k) = u_bar(:,k) + deltau_k;
            c_L = c_L + L_cost(u(:, k));
            [g(:, :, k + 1), xi(:, k + 1)] = f_dyn(g(:, :, k), xi(:, k), u(:, k));
            if (sum(any(isnan(g(:, :, k+1)))) || sum(any(isnan(xi(:, k+1)))))
                c_L = inf;
                break
            end
            deltax_k(1:n_gl) =  vmap_fun(real(logm(g_bar(:, :, k + 1)' * g(:, :, k + 1))));
            deltax_k(n_gl+1:end) = xi(:, k + 1) - xi_bar(:, k + 1);
        end
        c_F = F_final_cost(g(:, :, end), xi(:, end));
        cost_traj1 = c_L + c_F;  % cost of new trajectory
    
        if abs(cost_traj0 - cost_traj1) <= options_d.convtol  % convergence
            conv_flag = 1;
            no_cost_reduction = 0;
        elseif cost_traj1-cost_traj0 <= 0  % enough cost improvement
            no_cost_reduction = 0;
            % relax regularizer
            rho_change = min(rho_change / options_d.rhochange0, 1 / options_d.rhochange0);
            rho = options_d.rho0;
        else  % not enough cost improvement
            gamma = gamma * options_d.gammachange;
        end
    end
     if conv_flag 
        disp('Convergence');
        break
     end
    % new trajectory will be treated as nominal at the next iteration
    cost_traj0 = cost_traj1;
    g_bar = g;
    xi_bar = xi;
    u_bar = u;
end
u_bars = u_bars(:, :, 1:j);
J = J(1:j);

end