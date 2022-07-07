%%
clear; close all; clc;
% Import FD methods
addpath('functions')
addpath('NOT_MY_CODE')
% Initialization mode
init = 'noz';
% Adjustments to the Far Field Conditions
M_inf = 1.5;
if strcmp(init,'flat')
    P0_adjust = 1;
    num = 25;
elseif strcmp(init,'noz')
    load('noz_init.mat')
    % Add 1 more row for ghost nodes
    init_u = [u; u(end,:)];
    init_v = [v; v(end,:)];
    init_P = [P; P(end,:)];
    init_T = [T; T(end,:)];
elseif strcmp(init,'full')
    load('temp.mat')
    init_u = u;
    init_v = v;
    init_P = P;
    init_T = T;
end
% Useful globals for functions
global U_inf P_inf T_inf Ut Pt Tt
% Properties of standard air
T_inf = 288.15; % K
P_inf = 101300; % N/m2
R = 287;        % J/kg/K
Pr = 0.71;      % Peclet number
Cp = 1005;      % J/kg/K
Cv = 718;       % J/kg/K
gamma = 1.4;    % gamma
rho_inf = 1.225;% kg/m3
U_inf = M_inf*sqrt(gamma*R*T_inf); % airspeed of Mach 4
% Domain Properties
BC_cond = 'adiabatic'; % 'adiabatic'
SF = 1.5;           % Safety Factor
dt_min = 2.35e-11;  % Minimum time step
% Pre-allocating space for blocks
num_block_rows = 3; 
num_block_cols = 2;
empty = cell(num_block_rows,num_block_cols);
xx = empty; yy = empty;
dx = empty; dy = empty;
nx = empty; ny = empty;
% Nozzle Properties
P_0 = 1.2e6;    % Combustion Chamber pressure
T_0 = 2000;     % Combustion Chamber temperature
At = 1e-5;      % Throat height
num = num;      % Number of characteristic lines ~ nx1,ny1
% Nozzle Geometry (block 2)
[xx{2},yy{2}] = genmesh(T_0,P_0,At,num,R,gamma,P_inf);
xx{2}(end+1,:) = xx{2}(end,:) + (xx{2}(end,:)-xx{2}(end-1,:));
yy{2}(end+1,:) = yy{2}(end,:);
nx{2} = size(xx{2},1); ny{2} = size(xx{2},2);
dx{2} = xx{2}(end,1) - xx{2}(end-1,1); 
dy{2} = yy{2}(end,2) - yy{2}(end,1);
% Geometry / Conditions
Ae = 2*yy{2}(end,end);
L = xx{2}(end,1);
% Defining grids of each block (Blocks 4 5 6)
for blk = 4:6
    nx{blk} = nx{2}; 
    ny{blk} = ny{2};
    if blk == 4
        lower_y = yy{2}(end,end)-dy{2};
        upper_y = 3*yy{2}(end,end)-dy{2};
    elseif blk == 5
        lower_y = yy{2}(end,1);
        upper_y = yy{2}(end,end);
    elseif blk == 6
        lower_y = 3*yy{2}(end,1)+dy{2};
        upper_y = yy{2}(end,1)+dy{2};
    end
    [xx{blk},yy{blk}] = ndgrid(linspace(L-dx{2},2*L-dx{2},nx{blk}),linspace(lower_y,upper_y,ny{blk}));
    dx{blk} = xx{blk}(2,1)-xx{blk}(1,1); 
    dy{blk} = yy{blk}(1,2)-yy{blk}(1,1);
end
% Parametric geometry
[xi,eta] = ndgrid(linspace(0,L,nx{2}),linspace(0,Ae,ny{2}));
dxi = xi(2,1)-xi(1,1); deta = eta(1,2)-eta(1,1);
% Nozzle Metrics, Jacobian
x_xi = ddx_fwd(xx{2},dxi);
y_xi = ddx_fwd(yy{2},dxi);
x_eta = ddy_fwd(xx{2},deta);
y_eta = ddy_fwd(yy{2},deta);
Jac = x_xi.*y_eta - x_eta.*y_xi;
xi_x =  y_eta./Jac;
xi_y = -x_eta./Jac;
eta_x = -y_xi./Jac;
eta_y =  x_xi./Jac;
% Derivative Transforms
ddx = @(X,ddxi,ddeta) xi_x.*ddxi(X,dxi) + eta_x.*ddeta(X,deta);
ddy = @(X,ddxi,ddeta) xi_y.*ddxi(X,dxi) + eta_y.*ddeta(X,deta);
% Pressure adjustments post-geometry
P_0 = P0_adjust*P_0;
% Nozzle Throat Conditions
Tt = T_0*2/(gamma+1);
Pt = P_0*(2/(gamma+1))^(gamma/(gamma-1));
Ut = sqrt(gamma*R*Tt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Initial Conditions %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-allocating space
u = empty; v = empty;
P = empty; T = empty;
rho = empty; U = empty;
mu = empty; k = empty;
% Nozzle (block 2)
for ii = 1:nx{2}
    u{2}(ii,:) = (Ut*(nx{2}-ii)/nx{2} + U_inf*ii/nx{2})*ones(1,ny{2});
    P{2}(ii,:) = (Pt*(nx{2}-ii)/nx{2} + P_inf*ii/nx{2})*ones(1,ny{2});
    T{2}(ii,:) = (Tt*(nx{2}-ii)/nx{2} + T_inf*ii/nx{2})*ones(1,ny{2});
end
if strcmp(init,'flat')
    u{2} = Ut*ones(nx{2},ny{2});
    P{2} = Pt*ones(nx{2},ny{2});
    T{2} = Tt*ones(nx{2},ny{2});
end
v{2} = zeros(nx{2},ny{2});
if strcmp(init,'noz')
    u{2} = init_u;
    v{2} = init_v;
    P{2} = init_P;
    T{2} = init_T;
end
% Upper Exhaust (block 4)
for j = 1:ny{4}
    u{4}(:,j) = U_inf*ones(nx{blk},1);
    v{4}(:,j) = zeros(nx{blk},1);
    P{4}(:,j) = P_inf*ones(nx{blk},1);
    T{4}(:,j) = T_inf*ones(nx{blk},1);
end
% Center Exhaust (block 5)
for ii = 1:nx{5}
    u{5}(ii,:) = U_inf*ones(1,ny{5});
    v{5}(ii,:) = zeros(1,ny{5});
    P{5}(ii,:) = P_inf*ones(1,ny{5});
    T{5}(ii,:) = T_inf*ones(1,ny{5});
    if strcmp(init,'noz')
        u{5}(ii,:) = init_u(end,:)*(nx{5}-ii)/nx{5} + U_inf*(ii/nx{5})*ones(1,ny{5});
        v{5}(ii,:) = init_v(end,:)*(nx{5}-ii)/nx{5} + 0*(ii/nx{5})*ones(1,ny{5});
        P{5}(ii,:) = init_P(end,:)*(nx{5}-ii)/nx{5} + P_inf*(ii/nx{5})*ones(1,ny{5});
        T{5}(ii,:) = init_T(end,:)*(nx{5}-ii)/nx{5} + T_inf*(ii/nx{5})*ones(1,ny{5});
    end
end
% Lower Exhaust (block 6)
for j = 1:ny{6}
    u{6}(:,j) = U_inf*ones(nx{blk},1);
    v{6}(:,j) = zeros(nx{blk},1);
    P{6}(:,j) = P_inf*ones(nx{blk},1);
    T{6}(:,j) = T_inf*ones(nx{blk},1);
end
if strcmp(init,'full')
    u = init_u;
    v = init_v;
    P = init_P;
    T = init_T;
end
% Enforce BC's
for blk = [2 4 5 6]
    [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
    rho{blk} = P{blk}./(R*T{blk});
    U{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
    mu{blk} = sutherland(T{blk});
    k{blk} = (Cp/Pr)*mu{blk};
end
% Loop Initialization
time = 0;
dE_dx = empty; dF_dy = empty;
E = empty; F = empty;
U_bar = empty; U_next = empty;
Tau_xx = empty; Tau_xy = empty;
Tau_yy = empty; qdot_x = empty;
qdot_y = empty; dU = empty;
Et = empty;
temp = empty; 
Big_E = empty; Big_F = empty;
for blk = [2 4 5 6]
    dE_dx{blk} = zeros(4,nx{blk},ny{blk});
    dF_dy{blk} = zeros(4,nx{blk},ny{blk});
    E{blk} = zeros(4,nx{blk},ny{blk});
    F{blk} = zeros(4,nx{blk},ny{blk});
    U_bar{blk} = zeros(4,nx{blk},ny{blk});
    U_next{blk} = zeros(4,nx{blk},ny{blk});
    Tau_xx{blk} = zeros(nx{blk},ny{blk});
    Tau_xy{blk} = zeros(nx{blk},ny{blk});
    Tau_yy{blk} = zeros(nx{blk},ny{blk});
    qdot_x{blk} = zeros(nx{blk},ny{blk});
    qdot_y{blk} = zeros(nx{blk},ny{blk});
    Et{blk} = zeros(nx{blk},ny{blk});
    dU{blk} = zeros(1,4);
    temp{blk} = zeros(4,nx{blk},ny{blk});
    if blk == 2
        Big_E{blk} = zeros(4,nx{blk},ny{blk});
        Big_F{blk} = zeros(4,nx{blk},ny{blk});
    end
end
myfig = figure('Position', [10 10 1700 900]);
visualize(u,v,P,T,rho,Cv,xx,yy)
frame = 1;
for i = 1:6000
    %% Get Time-Step
    dt = dt_min/SF;
    % Save data from current time step
    U_old = U;
    % Iterate over all blocks
    for jj = 1:num_block_cols
        for ii = 1:num_block_rows
            ind = [ii,jj];
            blk = ii + (num_block_rows)*(jj-1);
            if blk==1 || blk==3
                continue
            end
            % Make sure all data is on the same time-level
            U = U_old;
            for blk = [2 4 5 6]
                [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U{blk},R,Cv);
                rho{blk} = P{blk}./(R*T{blk});
                mu{blk} = sutherland(T{blk});
                k{blk} = (Cp/Pr)*mu{blk};
            end
            ind = [ii,jj];
            blk = ii + (num_block_rows)*(jj-1);
            %% PREDICTOR STEP
            % Update Primitives
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U{blk},R,Cv);
            %%%%%%%% E-term %%%%%%%%
            % Define Derivatives
            if blk == 2
                ddx = @(X,dx) xi_x.*ddx_bwd(X,dx) ...
                                  + eta_x.*ddy_central(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_bwd(X,dx{blk}) ...
                                  + eta_y.*ddy_central(X,dy);
            else
                ddx = @ddx_bwd;
                ddy = @ddy_central;
            end
            % Define values
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_x{blk} = -k{blk}.*ddx(T{blk},dx{blk});
            % Update E-term
            E{blk}(1,:,:) = rho{blk}.*u{blk};
            E{blk}(2,:,:) = rho{blk}.*u{blk}.^2 + P{blk} - Tau_xx{blk};
            E{blk}(3,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            E{blk}(4,:,:) = (Et{blk} + P{blk}).*u{blk} - u{blk}.*Tau_xx{blk} ...
                - v{blk}.*Tau_xy{blk} + qdot_x{blk};
            %%%%%%%% F-term %%%%%%%%
            % Define Derivatives
            if blk == 2
                ddx = @(X,dx) xi_x.*ddx_central(X,dx) ...
                                  + eta_x.*ddy_bwd(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_central(X,dx{blk}) ...
                                  + eta_y.*ddy_bwd(X,dy);
            else
                ddx = @ddx_central;
                ddy = @ddy_bwd;
            end
            % Define values
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_y{blk} = -k{blk}.*ddy(T{blk},dy{blk});
            % Update F-term
            F{blk}(1,:,:) = rho{blk}.*v{blk};
            F{blk}(2,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            F{blk}(3,:,:) = rho{blk}.*v{blk}.^2 + P{blk} - Tau_yy{blk};
            F{blk}(4,:,:) = (Et{blk} + P{blk}).*v{blk} - v{blk}.*Tau_yy{blk} ...
                - u{blk}.*Tau_xy{blk} + qdot_y{blk};
            % Predictor step
            if blk == 2
                for j = 1:4
                    % Grid Transform
                    Big_E{blk}(j,:,:) = y_eta.*squeeze(E{blk}(j,:,:)) - x_eta.*squeeze(F{blk}(j,:,:));
                    Big_F{blk}(j,:,:) = - y_xi.*squeeze(E{blk}(j,:,:)) + x_xi.*squeeze(F{blk}(j,:,:));
                    % Forward Derivatives
                    dE_dx{blk}(j,:,:) = ddx_fwd(squeeze(Big_E{blk}(j,:,:)),dxi);
                    dF_dy{blk}(j,:,:) = ddy_fwd(squeeze(Big_F{blk}(j,:,:)),deta);
                    % Solve Predictor step
                    temp{blk}(j,:,:) = (Jac.*squeeze(U{blk}(j,:,:)) ...
                        - dt*squeeze(dE_dx{blk}(j,:,:)) ...
                        - dt*squeeze(dF_dy{blk}(j,:,:)))./Jac;
                end
            else
                % Forward Derivatives
                for j = 1:4
                    dE_dx{blk}(j,:,:) = ddx_fwd(squeeze(E{blk}(j,:,:)),dx{blk});
                    dF_dy{blk}(j,:,:) = ddy_fwd(squeeze(F{blk}(j,:,:)),dy{blk});
                end
                % Solve Predictor step
                temp{blk} = U{blk} - dt*dE_dx{blk} - dt*dF_dy{blk};
            end
            [~,u{blk},v{blk},T{blk},P{blk},~,~] = cons2prim(temp{blk},R,Cv);
            % Enforce boundary conditions & update variables
            [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output U_bar
            U_bar{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
            
            %% CORRECTOR STEP
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_bar{blk},R,Cv);
            %%%%%%%% E-term %%%%%%%%
            % Define Derivatives
            if blk == 2
                ddx = @(X,dx) xi_x.*ddx_fwd(X,dx) ...
                                  + eta_x.*ddy_central(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_fwd(X,dx{blk}) ...
                                  + eta_y.*ddy_central(X,dy);
            else
                ddx = @ddx_fwd;
                ddy = @ddy_central;
            end
            % Define values
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_x{blk} = -k{blk}.*ddx(T{blk},dx{blk});
            % Update E-term
            E{blk}(1,:,:) = rho{blk}.*u{blk};
            E{blk}(2,:,:) = rho{blk}.*u{blk}.^2 + P{blk} - Tau_xx{blk};
            E{blk}(3,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            E{blk}(4,:,:) = (Et{blk} + P{blk}).*u{blk} - u{blk}.*Tau_xx{blk} ...
                - v{blk}.*Tau_xy{blk} + qdot_x{blk};
            %%%%%%%% F-term %%%%%%%%
            % Define Derivatives
            if blk == 2
                ddx = @(X,dx) xi_x.*ddx_central(X,dx) ...
                                  + eta_x.*ddy_fwd(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_central(X,dx{blk}) ...
                                  + eta_y.*ddy_fwd(X,dy);
            else
                ddx = @ddx_central;
                ddy = @ddy_fwd;
            end
            % Define values
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_y{blk} = -k{blk}.*ddy(T{blk},dy{blk});
            % Update F-term
            F{blk}(1,:,:) = rho{blk}.*v{blk};
            F{blk}(2,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            F{blk}(3,:,:) = rho{blk}.*v{blk}.^2 + P{blk} - Tau_yy{blk};
            F{blk}(4,:,:) = (Et{blk} + P{blk}).*v{blk} - v{blk}.*Tau_yy{blk} ...
                - u{blk}.*Tau_xy{blk} + qdot_y{blk};
            % Corrector Step
            if blk == 2
                for j = 1:4
                    % Grid Transform
                    Big_E{blk}(j,:,:) = y_eta.*squeeze(E{blk}(j,:,:)) - x_eta.*squeeze(F{blk}(j,:,:));
                    Big_F{blk}(j,:,:) = - y_xi.*squeeze(E{blk}(j,:,:)) + x_xi.*squeeze(F{blk}(j,:,:));
                    % Backward Derivatives
                    dE_dx{blk}(j,:,:) = ddx_bwd(squeeze(Big_E{blk}(j,:,:)),dxi);
                    dF_dy{blk}(j,:,:) = ddy_bwd(squeeze(Big_F{blk}(j,:,:)),deta);
                    % Solve Corrector Step
                    temp{blk}(j,:,:) = (0.5*(Jac.*squeeze(U{blk}(j,:,:))+Jac.*squeeze(U_bar{blk}(j,:,:))) ...
                           - 0.5*dt*squeeze(dE_dx{blk}(j,:,:)) ...
                           - 0.5*dt*squeeze(dF_dy{blk}(j,:,:)))./Jac;
                end
            else
                % Backward Derivatives
                for j = 1:4
                    dE_dx{blk}(j,:,:) = ddx_bwd(squeeze(E{blk}(j,:,:)),dx{blk});
                    dF_dy{blk}(j,:,:) = ddy_bwd(squeeze(F{blk}(j,:,:)),dy{blk});
                end
                % Solve Corrector Step
                temp{blk} = 0.5*(U{blk}+U_bar{blk}) - 0.5*dt*dE_dx{blk} - 0.5*dt*dF_dy{blk};
            end
            [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(temp{blk},R,Cv);
            % Enforce boundary conditions & update variables
            [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output next step
            U_next{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);

        end
    end
    % Update variables for plotting
    for blk = [2 4 5 6]
        [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_next{blk},R,Cv);
        [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
        rho{blk} = P{blk}./(R*T{blk});
        mu{blk} = sutherland(T{blk});
        k{blk} = (Cp/Pr)*mu{blk};
    end
    for blk = [2 4 5 6]
        U_next{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
        % Convergence Criteria
        for j = 1:4
            dU{blk}(i,j) = norm(squeeze(U_next{blk}(j,:,:)) - squeeze(U{blk}(j,:,:)));
        end
        if mod(i,200)==0
            fprintf(sprintf('Convergence (blk %1d): [%1.1d %1.1d %1.1d %1.1d] \n',blk,dU{blk}(i,:)./dU{blk}(1,:)))
        end
    end
    
    % Visualize every 100th point
    if mod(i,500)==0 || i==1
        fprintf(sprintf('----- iter %d ----- \n',i))
        visualize(u,v,P,T,rho,Cv,xx,yy)
        i_frame = getframe(myfig);
        im = frame2im(i_frame);
        [imind,cm] = rgb2ind(im,256);
        if frame == 1
            imwrite(imind,cm,'Final.gif','gif', 'Loopcount',inf,...
            'DelayTime',0.2);
        else
            imwrite(imind,cm,'Final.gif','gif','WriteMode','append',...
            'DelayTime',0.2);
        end
        frame = frame + 1;
    end
    
%     check_connections(U)
    % Update loop
    U = U_next;
end
% Visualize last step
visualize(u,v,P,T,rho,Cv,xx,yy)
fprintf('Finished\n')

%% Define useful functions
% Check block connections
function check_connections(U)
    for blk = [2 4 5 6]
        [~,temp_u{blk},temp_v{blk},temp_T{blk},temp_P{blk},~,temp_Et{blk}] = cons2prim(U{blk},287,718);
    end
    err25 = norm(temp_u{2}(end,:) - temp_u{5}(2,:))+norm(temp_v{2}(end,:) - temp_v{5}(2,:))...
        + norm(temp_P{2}(end,:) - temp_P{5}(2,:)) + norm(temp_T{2}(end,:) - temp_T{5}(2,:));
    err45 = norm(temp_u{4}(:,2) - temp_u{5}(:,end))...
     + norm(temp_v{4}(:,2) - temp_v{5}(:,end))...
     + norm(temp_P{4}(:,2) - temp_P{5}(:,end))...
     + norm(temp_T{4}(:,2) - temp_T{5}(:,end));
    err56 = norm(temp_u{5}(:,2) - temp_u{6}(:,end))...
     + norm(temp_v{5}(:,2) - temp_v{6}(:,end))...
     + norm(temp_P{5}(:,2) - temp_P{6}(:,end))...
     + norm(temp_T{5}(:,2) - temp_T{6}(:,end));
    errs = [err25 err45 err56];
    err25 = norm(temp_u{2}(end-1,:) - temp_u{5}(1,:))+norm(temp_v{2}(end-1,:) - temp_v{5}(1,:))...
        + norm(temp_P{2}(end-1,:) - temp_P{5}(1,:)) + norm(temp_T{2}(end-1,:) - temp_T{5}(1,:));
    err45 = norm(temp_u{4}(:,1) - temp_u{5}(:,end-1))...
     + norm(temp_v{4}(:,1) - temp_v{5}(:,end-1))...
     + norm(temp_P{4}(:,1) - temp_P{5}(:,end-1))...
     + norm(temp_T{4}(:,1) - temp_T{5}(:,end-1));
    err56 = norm(temp_u{5}(:,1) - temp_u{6}(:,end-1))...
     + norm(temp_v{5}(:,1) - temp_v{6}(:,end-1))...
     + norm(temp_P{5}(:,1) - temp_P{6}(:,end-1))...
     + norm(temp_T{5}(:,1) - temp_T{6}(:,end-1));
    errs = [errs; err25 err45 err56];
    disp(errs)
end
% Enforce the boundary contions on primitive variables
function [u,v,P,T] = enforce_bcs(u,v,P,T,mode,blk)
    global U_inf P_inf T_inf Ut Pt Tt
    if blk == 2
        % Top /// Nozzle Wall
        u{blk}(:,end) = 0;
        v{blk}(:,end) = 0;
        P{blk}(2:end,end) = 2*P{blk}(2:end,end-1) - P{blk}(2:end,end-2);
        switch mode
            case 'standard'
                T{blk}(1:(end-1),end) = Tt;
            case 'adiabatic'
                T{blk}(1:end,end) = (4*T{blk}(1:end,end-1) - T{blk}(1:end,end-2))/3;
            case 'test'
                T{blk}(1:end,end) = T{blk}(1:end,end-1)/2 + T_inf/2;
        end
        % Bottom /// Nozzle Wall
        u{blk}(:,1) = 0;
        v{blk}(:,1) = 0;
        P{blk}(2:end,1) = 2*P{blk}(2:end,2) - P{blk}(2:end,3);
        switch mode
            case 'standard'
                T{blk}(1:end,1) = Tt;
            case 'adiabatic'
                T{blk}(1:end,1) = (4*T{blk}(1:end,2) - T{blk}(1:end,3))/3;
            case 'test'
                T{blk}(1:end,1) = T{blk}(1:end,2)/2 + T_inf/2;
        end
        % Left /// Choked Bc's
        u{blk}(1,:) = Ut;
        v{blk}(1,:) = 0;
        P{blk}(1,:) = Pt;
        T{blk}(1,:) = Tt;
        % Connectivity (block 5 rightward)
        u{blk}(end,:) = u{5}(2,:);
        v{blk}(end,:) = v{5}(2,:);
        P{blk}(end,:) = P{5}(2,:);
        T{blk}(end,:) = T{5}(2,:);
    elseif blk == 4
        % Nozzle Edge
        u{blk}(1,1) = 0;
        T{blk}(1,1) = T_inf;
        % Top /// Extrapolation
        u{blk}(:,end) = 2*u{blk}(:,end-1) - u{blk}(:,end-2);
        v{blk}(:,end) = 2*v{blk}(:,end-1) - v{blk}(:,end-2);
        P{blk}(:,end) = 2*P{blk}(:,end-1) - P{blk}(:,end-2);
        T{blk}(:,end) = 2*T{blk}(:,end-1) - T{blk}(:,end-2);
        % Left /// Far Field
        u{blk}(1,:) = U_inf;
        v{blk}(1,:) = 0;
        P{blk}(1,2:end) = P_inf;
        T{blk}(1,:) = T_inf;
        % Connectivity
        u{blk}(:,1) = u{5}(:,end-1);
        v{blk}(:,1) = v{5}(:,end-1);
        P{blk}(:,1) = P{5}(:,end-1);
        T{blk}(:,1) = T{5}(:,end-1);
    elseif blk == 5
        % Nozzle edges
        u{blk}(1,1) = 0;
        u{blk}(1,end) = 0;
        v{blk}(1,1) = 0;
        v{blk}(1,end) = 0;
        P{blk}(1,1) = P_inf;
        P{blk}(1,end) = P_inf;
        T{blk}(1,1) = T_inf;
        T{blk}(1,end) = T_inf;
        % Connectivity (block 2 leftward)
        u{blk}(1,:) = u{2}(end-1,:);
        v{blk}(1,:) = v{2}(end-1,:);
        P{blk}(1,:) = P{2}(end-1,:);
        T{blk}(1,:) = T{2}(end-1,:);
        % Connectivity (block 4 above)
        u{blk}(:,end) = u{4}(:,2);
        v{blk}(:,end) = v{4}(:,2);
        P{blk}(:,end) = P{4}(:,2);
        T{blk}(:,end) = T{4}(:,2);
        % Connectivity (block 6 below)
        u{blk}(:,1) = u{6}(:,end-1);
        v{blk}(:,1) = v{6}(:,end-1);
        P{blk}(:,1) = P{6}(:,end-1);
        T{blk}(:,1) = T{6}(:,end-1);
    elseif blk == 6
        % Nozzle Edge
        u{blk}(1,end) = 0;
        T{blk}(1,end) = T_inf;
        % Bottom /// Extrapolation
        u{blk}(:,1) = 2*u{blk}(:,2)-u{blk}(:,3);
        v{blk}(:,1) = 2*v{blk}(:,2)-v{blk}(:,3);
        P{blk}(:,1) = 2*P{blk}(:,2)-P{blk}(:,3);
        T{blk}(:,1) = 2*T{blk}(:,2)-T{blk}(:,3);
        % Left /// Far Field
        u{blk}(1,:) = U_inf;
        v{blk}(1,:) = 0;
        P{blk}(1,:) = P_inf;
        T{blk}(1,:) = T_inf;
        % Connectivity (block 5 above)
        u{blk}(:,end) = u{5}(:,2);
        v{blk}(:,end) = v{5}(:,2);
        P{blk}(:,end) = P{5}(:,2);
        T{blk}(:,end) = T{5}(:,2);
    end
    if blk ~= 2
        % Right /// Extrapolation
        u{blk}(end,:) = 2*u{blk}(end-1,:) - u{blk}(end-2,:);
        v{blk}(end,:) = 2*v{blk}(end-1,:) - v{blk}(end-2,:);
        P{blk}(end,:) = 2*P{blk}(end-1,:) - P{blk}(end-2,:);
        T{blk}(end,:) = 2*T{blk}(end-1,:) - T{blk}(end-2,:);
    end
end
% Visualize the primitive variables
function visualize(u,v,P,T,rho,Cv,xx,yy)
    pltfmt = {'interpreter','latex','FontSize',16};
    for blk = 1:numel(u)
        if blk == 1 || blk == 3
            continue
        end
        % density
        subplot(2,3,1)
        hold on;
        plot = pcolor(xx{blk},yy{blk},rho{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'$\rho$ [$kg/m^3$]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('$\rho$ [$kg/m^3$]',pltfmt{:})
        axis equal tight
        shading interp
        % u-velocity
        subplot(2,3,2)
        hold on;
        plot = pcolor(xx{blk},yy{blk},u{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'u [m/s]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('u [m/s]',pltfmt{:})
        axis equal tight
        shading interp
        % v-velocity
        subplot(2,3,3)
        hold on;
        plot = pcolor(xx{blk},yy{blk},v{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'v [m/s]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('v [m/s]',pltfmt{:})
        axis equal tight
        shading interp
        % Energy
        subplot(2,3,4)
        hold on;
        plot = pcolor(xx{blk},yy{blk},Cv*T{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'e [J/kg]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('e [J/kg]',pltfmt{:})
        axis equal tight
        shading interp
        % Pressure
        subplot(2,3,5)
        hold on;
        plot = pcolor(xx{blk},yy{blk},P{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'P [N/$m^2$]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('P [$N/m^2$]',pltfmt{:})
        axis equal tight
        shading interp
        % Temperature
        subplot(2,3,6)
        hold on;
        plot = pcolor(xx{blk},yy{blk},T{blk});
        set(plot, 'EdgeColor', 'none');
        colormap('turbo')
        cb = colorbar; ylabel(cb,'T [K]',pltfmt{:});
        ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
        title('T [K]',pltfmt{:})
        axis equal tight
        shading interp
        % Colormap colors
        colormap jet
    end
    % Update plot
    drawnow
end