%%
clear all; close all; clc;
% Import FD methods
addpath('functions')
addpath('NOT_MY_CODE')
% Adjust the pressure post-geometry
P0_adjust = 3;
% Useful globals for functions
global U_inf P_inf T_inf Ut Pt Tt
% Properties of standard air
T_inf = 288.15; % K
P_inf = 101300; % N/m2
M_inf = 1.5;    % Mach
R = 287;        % J/kg/K
Pr = 0.71;      % Peclet number
Cp = 1005;      % J/kg/K
Cv = 718;       % J/kg/K
gamma = 1.4;    % gamma
rho_inf = 1.225;% kg/m3
U_inf = M_inf*sqrt(gamma*R*T_inf); % airspeed of Mach 4
% Domain Properties
SF = 1.2;             % Safety Factor
dt_min = 2.35e-11;  % Minimum time step
% Pre-allocating space for blocks
num_block_rows = 1; 
num_block_cols = 2;
empty = cell(num_block_rows,num_block_cols);
xx = empty; yy = empty;
dx = empty; dy = empty;
u = empty; v = empty;
P = empty; T = empty;
rho = empty; U = empty;
mu = empty; k = empty;
nx = empty; ny = empty;
% Nozzle Properties
P_0 = 1.2e6;    % Combustion Chamber pressure
T_0 = 2000;     % Combustion Chamber temperature
At = 1e-5;      % Throat height
num = 25;       % Number of characteristic lines ~ nx1,ny1
% Nozzle Geometry (block 2)
[xx{1},yy{1}] = genmesh(T_0,P_0,At,num,R,gamma,P_inf);
xx{1}(end+1,:) = xx{1}(end,:) + (xx{1}(end,:)-xx{1}(end-1,:));
yy{1}(end+1,:) = yy{1}(end,:);
nx{1} = size(xx{1},1); ny{1} = size(xx{1},2);
dx{1} = xx{1}(end,1) - xx{1}(end-1,1); 
dy{1} = yy{1}(end,2) - yy{1}(end,1);
% Pressure adjustments post-geometry
P_0 = P0_adjust*P_0;
% Geometry / Conditions
Ae = 2*yy{1}(end,end);
L = xx{1}(end,1);
% Parametric geometry
[xi,eta] = ndgrid(linspace(0,L,nx{1}),linspace(0,Ae,ny{1}));
dxi = xi(2,1)-xi(1,1); deta = eta(1,2)-eta(1,1);
% Nozzle Throat Conditions
Tt = T_0*2/(gamma+1);
Pt = P_0*(2/(gamma+1))^(gamma/(gamma-1));
Ut = sqrt(gamma*R*Tt);
% Defining grids of each block (Block 2)
for blk = 2
    nx{blk} = nx{1}; 
    ny{blk} = ny{1};
    if blk == 2
        lower = yy{1}(end,1);
        upper = yy{1}(end,end);
    end
    [xx{blk},yy{blk}] = ndgrid(linspace(L-dx{1},2*L-dx{1},nx{blk}),linspace(lower,upper,ny{blk}));
    dx{blk} = xx{blk}(2,1)-xx{blk}(1,1); 
    dy{blk} = yy{blk}(1,2)-yy{blk}(1,1);
end
% Nozzle Metrics, Jacobian
x_xi = ddx_fwd(xx{1},dxi);
y_xi = ddx_fwd(yy{1},dxi);
x_eta = ddy_fwd(xx{1},deta);
y_eta = ddy_fwd(yy{1},deta);
Jac = x_xi.*y_eta - x_eta.*y_xi;
xi_x =  y_eta./Jac;
xi_y = -x_eta./Jac;
eta_x = -y_xi./Jac;
eta_y =  x_xi./Jac;
% Derivative Transforms
ddx = @(X,ddxi,ddeta) xi_x.*ddxi(X,dxi) + eta_x.*ddeta(X,deta);
ddy = @(X,ddxi,ddeta) xi_y.*ddxi(X,dxi) + eta_y.*ddeta(X,deta);
% Initial Conditions
% Nozzle (block 1)
u{1} = Ut*ones(nx{1},ny{1});
P{1} = Pt*ones(nx{1},ny{1});
T{1} = Tt*ones(nx{1},ny{1});
for ii = 1:nx{1}
    u{1}(ii,:) = (Ut*(1-(ii/nx{1})) + U_inf*(ii/nx{1})^2)*ones(1,ny{1});
    P{1}(ii,:) = (Pt*(1-(ii/nx{1})) + P_inf*(ii/nx{1})^2)*ones(1,ny{1});
    T{1}(ii,:) = (Tt*(1-(ii/nx{1})) + T_inf*(ii/nx{1})^2)*ones(1,ny{1});
end
% u{1} = U_inf*ones(nx{1},ny{1});
% P{1} = P_inf*ones(nx{1},ny{1});
% T{1} = T_inf*ones(nx{1},ny{1});
v{1} = zeros(nx{1},ny{1});
% Exhaust (block 2)
for ii = 1:nx{2}
    u{2}(ii,:) = U_inf*ones(1,ny{2});      % Initial u velocity
    v{2}(ii,:) = zeros(1,ny{2});            % Initial u velocity
    P{2}(ii,:) = P_inf*ones(1,ny{2});      % Initial Pressure
    T{2}(ii,:) = T_inf*ones(1,ny{2});      % Initial Temperature
end
% Enforce BC's
for blk = 1:2
    [u,v,P,T] = enforce_bcs(u,v,P,T,blk);
    rho{blk} = P{blk}./(R*T{blk});
    U{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);% Initial Constituitive
    mu{blk} = sutherland(T{blk});         % Initial viscosity
    k{blk} = (Cp/Pr)*mu{blk};             % Initial conductivity
end
% Loop Initialization
time = 0;
dE_dx = empty; dF_dy = empty;
E = empty; F = empty;
U_bar = empty; U_next = empty;
Tau_xx = empty; Tau_xy = empty;
Tau_yy = empty; qdot_x = empty;
qdot_y = empty; dU = empty;
Et = empty; e = empty;
temp = empty; U_old = empty;
Big_E = empty; Big_F = empty;
for blk = 1:2
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
myfig = figure('Position', [10 10 1800 500]);
visualize(u,v,P,T,rho,Cv,xx,yy)
frame = 1;
for i = 1:10000
    %% Get Time-Step
    dt = dt_min/SF;
    U_old = U;
    % Iterate over all blocks
    for jj = 1:num_block_cols
        for ii = 1:num_block_rows
            % Make sure all data is on the same time-level
            U = U_old;
            for blk = [1 2]
                [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U{blk},R,Cv);
                rho{blk} = P{blk}./(R*T{blk});
            end
            % updata index
            ind = [ii,jj];
            blk = ii + (num_block_rows)*(jj-1);
%             check_connections(U)
            %% PREDICTOR STEP
            % Update Primitives
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U{blk},R,Cv);
            %%%%%%%% E-term %%%%%%%%
            % Define Derivatives
            if blk == 1
                ddx = @(X,dx) xi_x.*ddx_bwd(X,dx) ...
                                  + eta_x.*ddy_central(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_bwd(X,dx{blk}) ...
                                  + eta_y.*ddy_central(X,dy);
            elseif blk == 2
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
            if blk == 1
                ddx = @(X,dx) xi_x.*ddx_central(X,dx) ...
                                  + eta_x.*ddy_bwd(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_central(X,dx{blk}) ...
                                  + eta_y.*ddy_bwd(X,dy);
            elseif blk == 2
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
            if blk == 1
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
            [u,v,P,T] = enforce_bcs(u,v,P,T,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output U_bar
            U_bar{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
            
            %% CORRECTOR STEP
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_bar{blk},R,Cv);
            %%%%%%%% E-term %%%%%%%%
            % Define Derivatives
            if blk == 1
                ddx = @(X,dx) xi_x.*ddx_fwd(X,dx) ...
                                  + eta_x.*ddy_central(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_fwd(X,dx{blk}) ...
                                  + eta_y.*ddy_central(X,dy);
            elseif blk == 2
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
            if blk == 1
                ddx = @(X,dx) xi_x.*ddx_central(X,dx) ...
                                  + eta_x.*ddy_fwd(X,dy{blk});
                ddy = @(X,dy) xi_y.*ddx_central(X,dx{blk}) ...
                                  + eta_y.*ddy_fwd(X,dy);
            elseif blk == 2
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
            if blk == 1
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
            [u,v,P,T] = enforce_bcs(u,v,P,T,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output next step
            U_next{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
            
        end
    end
    % Update variables for plotting
    for blk = 1:2
        [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_next{blk},R,Cv);
        [u,v,P,T] = enforce_bcs(u,v,P,T,blk);
        rho{blk} = P{blk}./(R*T{blk});
        mu{blk} = sutherland(T{blk});
        k{blk} = (Cp/Pr)*mu{blk};
    end
    for blk = 1:2
        U_next{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
        for j = 1:4
            dU{blk}(i,j) = norm(squeeze(U_next{blk}(j,:,:)) - squeeze(U{blk}(j,:,:)));
        end
        if mod(i,150)==0
            fprintf(sprintf('Convergence (blk %1d): [%1.1d %1.1d %1.1d %1.1d] \n',blk,dU{blk}(i,:)))
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
            'DelayTime',0.1);
        else
            imwrite(imind,cm,'Final.gif','gif','WriteMode','append',...
            'DelayTime',0.1);
        end
        frame = frame + 1;
    end
    % Update loop
    U = U_next;
end
% Visualize last step
visualize(u,v,P,T,rho,Cv,xx,yy)
fprintf('Finished\n')

%% Define useful functions
% Check block connections
function check_connections(U)
    for blk = 1:2
        [~,temp_u{blk},temp_v{blk},temp_T{blk},temp_P{blk},~,temp_Et{blk}] = cons2prim(U{blk},287,718);
    end
    err25 = norm(temp_u{1}(end,:) - temp_u{2}(2,:))+norm(temp_v{1}(end,:) - temp_v{2}(2,:))...
        + norm(temp_P{1}(end,:) - temp_P{2}(2,:)) + norm(temp_T{1}(end,:) - temp_T{2}(2,:));
    disp(err25)
    err25 = norm(temp_u{1}(end-1,:) - temp_u{2}(1,:))+norm(temp_v{1}(end-1,:) - temp_v{2}(1,:))...
        + norm(temp_P{1}(end-1,:) - temp_P{2}(1,:)) + norm(temp_T{1}(end-1,:) - temp_T{2}(1,:));
    disp(err25)
end
% Enforce the boundary contions on primitive variables
function [u,v,P,T] = enforce_bcs(u,v,P,T,blk)
    global U_inf P_inf T_inf Ut Pt Tt
    if blk == 1
        % Left /// Choked Bc's
        u{blk}(1,:) = Ut;
        v{blk}(1,:) = 0;
        P{blk}(1,:) = Pt;
        T{blk}(1,:) = Tt;
        % Top /// Nozzle Wall
        u{blk}(:,end) = 0;
        v{blk}(:,end) = 0;
        P{blk}(2:end,end) = 2*P{blk}(2:end,end-1) - P{blk}(2:end,end-2);
        T{blk}(2:end,end) = Tt;
        % Bottom /// Nozzle Wall
        u{blk}(:,1) = 0;
        v{blk}(:,1) = 0;
        P{blk}(2:end,1) = 2*P{blk}(2:end,2) - P{blk}(2:end,3);
        T{blk}(2:end,1) = Tt;
        % Nozzle Edge
        u{blk}(end,end) = 0;
        v{blk}(end,end) = 0;
        T{blk}(end,end) = T_inf;
        u{blk}(end,1) = 0;
        v{blk}(end,1) = 0;
        T{blk}(end,1) = T_inf;
        % Connectivity
        u{blk}(end,:) = u{2}(2,:);
        v{blk}(end,:) = v{2}(2,:);
        P{blk}(end,:) = P{2}(2,:);
        T{blk}(end,:) = T{2}(2,:);
    elseif blk == 2
        % Top /// Far Field
        u{blk}(:,end) = U_inf;
        v{blk}(:,end) = 0;
        P{blk}(:,end) = P_inf;
        T{blk}(:,end) = T_inf;
        % Top /// Extrapolation
%         u{blk}(:,end) = u{blk}(:,end-1);
%         v{blk}(:,end) = v{blk}(:,end-1);
%         P{blk}(:,end) = P{blk}(:,end-1);
%         T{blk}(:,end) = T{blk}(:,end-1);
        % Bottom /// Far Field
        u{blk}(:,1) = U_inf;
        v{blk}(:,1) = 0;
        P{blk}(:,1) = P_inf;
        T{blk}(:,1) = T_inf;
        % Bottom /// Extrapolation
%         u{blk}(:,1) = u{blk}(:,2);
%         v{blk}(:,1) = v{blk}(:,2);
%         P{blk}(:,1) = P{blk}(:,2);
%         T{blk}(:,1) = T{blk}(:,2);

        % Right /// Extrapolation
        u{blk}(end,:) = u{blk}(end-1,:);
        v{blk}(end,:) = v{blk}(end-1,:);
        P{blk}(end,:) = P{blk}(end-1,:);
        T{blk}(end,:) = T{blk}(end-1,:);
        % Nozzle Edge
        u{blk}(1,end) = 0;
        v{blk}(1,end) = 0;
        T{blk}(1,end) = T_inf;
        u{blk}(1,1) = 0;
        v{blk}(1,1) = 0;
        T{blk}(1,1) = T_inf;
        % Connectivity
        u{blk}(1,:) = u{1}(end-1,:);
        v{blk}(1,:) = v{1}(end-1,:);
        P{blk}(1,:) = P{1}(end-1,:);
        T{blk}(1,:) = T{1}(end-1,:);
    end
end
% Visualize the primitive variables
function visualize(u,v,P,T,rho,Cv,xx,yy)
    for blk = 1:numel(u)
        % density
        subplot(2,3,1)
        hold on;
        plot = pcolor(xx{blk},yy{blk},rho{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'\rho [kg/m^3]');
        ylabel('y'); xlabel('x');
        title('\rho [kg/m^3]')
        axis equal tight
        shading interp
        % u-velocity
        subplot(2,3,2)
        hold on;
        plot = pcolor(xx{blk},yy{blk},u{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'u [m/s]');
        ylabel('y'); xlabel('x');
        title('u [m/s]')
        axis equal tight
        shading interp
        % v-velocity
        subplot(2,3,3)
        hold on;
        plot = pcolor(xx{blk},yy{blk},v{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'v [m/s]');
        ylabel('y'); xlabel('x');
        title('v [m/s]')
        axis equal tight
        shading interp
        % Energy
        subplot(2,3,4)
        hold on;
        plot = pcolor(xx{blk},yy{blk},Cv*T{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'e [J/kg]');
        ylabel('y'); xlabel('x');
        title('e [J/kg]')
        axis equal tight
        shading interp
        % Pressure
        subplot(2,3,5)
        hold on;
        plot = pcolor(xx{blk},yy{blk},P{blk});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'P [N/m^2]');
        ylabel('y'); xlabel('x');
        title('P [N/m^2]')
        axis equal tight
        shading interp
        % Temperature
        subplot(2,3,6)
        hold on;
        plot = pcolor(xx{blk},yy{blk},T{blk});
        set(plot, 'EdgeColor', 'none');
        colormap('turbo')
        cb = colorbar; ylabel(cb,'T [K]');
        ylabel('y'); xlabel('x');
        title('T [K]')
        axis equal tight
        shading interp
        % Colormap colors
        colormap jet
    end
    % Update plot
    drawnow
end