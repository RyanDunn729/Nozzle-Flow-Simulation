%%
clear; close all; clc;
% Import FD method & Olson's code
addpath('functions')
% Useful globals for functions
global Ut Pt Tt
% Properties of standard air
T_atm = 288.15; % K
P_atm = 101300; % N/m2
R = 287;        % J/kg/K
Pr = 0.71;      % Peclet number
Cp = 1005;      % J/kg/K
Cv = 718;       % J/kg/K
gamma = 1.4;    % gamma
rho_inf = 1.225;% kg/m3
% Domain Properties
BC_cond = 'standard';
% BC_cond = 'adiabatic';
SF = 1.2;           % Safety Factor
dt_min = 2.35e-11;  % Minimum time step
% Geometry / Conditions
L = 1e-5;
Ae = 8e-6;
nx1 = 75; ny1 = 80;
% Spatial Coordinates
[xx,yy] = ndgrid(linspace(0,L,nx1),linspace(0,Ae,ny1));
% Parametric geometry
[xi,eta] = ndgrid(linspace(0,L/2,nx1),linspace(0,Ae/3,ny1));
dxi = xi(2,1)-xi(1,1); deta = eta(1,2)-eta(1,1);
% Throat Conditions
Tt = 288.15;
Pt = 101300;
Ut = sqrt(gamma*R*Tt);
% Initial Conditions
sos = sqrt(gamma*R*Tt);     % speed of sound
u = 4*Ut*ones(nx1,ny1);     % Initial u velocity
v = zeros(nx1,ny1);         % Initial v velocity
P = Pt*ones(nx1,ny1);       % Initial Pressure
T = Tt*ones(nx1,ny1);       % Initial Temperature
[u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond); % Enforce Boundary Conditions
rho = P./(R*T);             % Initial Density
U = prim2cons(rho,u,v,T,Cv);% Initial Constituitive
mu = sutherland(T);         % Initial viscosity
k = (Cp/Pr)*mu;             % Initial conductivity
% Metrics, Jacobian
x_xi = ddx_fwd(xx,dxi);
y_xi = ddx_fwd(yy,dxi);
x_eta = ddy_fwd(xx,deta);
y_eta = ddy_fwd(yy,deta);
Jac = x_xi.*y_eta - x_eta.*y_xi;
xi_x =  y_eta./Jac;
xi_y = -x_eta./Jac;
eta_x = -y_xi./Jac;
eta_y =  x_xi./Jac;
% Derivative Transforms
ddx = @(X,ddxi,ddeta) xi_x.*ddxi(X,dxi) + eta_x.*ddeta(X,deta);
ddy = @(X,ddxi,ddeta) xi_y.*ddxi(X,dxi) + eta_y.*ddeta(X,deta);
% Loop Initialization
time = 0;
dE_dx = zeros(4,nx1,ny1);
dF_dy = zeros(4,nx1,ny1);
E = zeros(4,nx1,ny1);
F = zeros(4,nx1,ny1);
temp = zeros(4,nx1,ny1);
Big_E = zeros(4,nx1,ny1);
Big_F = zeros(4,nx1,ny1);
U_bar = zeros(4,nx1,ny1);
U_next = zeros(4,nx1,ny1);
Tau_xx = zeros(nx1,ny1);
Tau_xy = zeros(nx1,ny1);
Tau_yy = zeros(nx1,ny1);
qdot_x = zeros(nx1,ny1);
qdot_y = zeros(nx1,ny1);
dU = zeros(1,4);
% Visualize the grids
viz = figure('Position',[10 10 700 500]);
plot(xi(1,1),eta(1,1),'r.','MarkerSize',8)
hold on
plot(xx,yy,'k.','MarkerSize',8)
plot(xi,eta,'r.','MarkerSize',8)
title('Grid Transformation','Interpreter','latex','FontSize',14)
l = legend('Parametric','Spatial','Interpreter','latex','FontSize',14,'Location','NorthEastOutside');
xlabel('x/\xi','FontSize',14)
ylabel('y/\eta','FontSize',14)
axis tight equal
% Initialize Figure
spatial = figure('Position',[10 10 1200 500]);
visualize(u,v,P,T,rho,Cv,xx,yy)
frame = 1;
for i = 1:1500
    %% Time-Step set to minimum
    dt = dt_min/SF;
    %% PREDICTOR STEP
    % Update Primitives
    [rho,u,v,T,P,~,Et] = cons2prim(U,R,Cv);
    %%%%%%%% E-term %%%%%%%%
    % Derivatives
    ddxi  = @ddx_bwd;     % xi-derivatives
    ddeta = @ddy_central; % eta-derivatives
    % values
    Tau_xx = 2*mu.*(ddx(u,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_yy = 2*mu.*(ddy(v,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_xy = mu.*(ddy(u,ddxi,ddeta)+ddx(v,ddxi,ddeta));
    qdot_x = -k.*ddx(T,ddxi,ddeta);
    % update E-term
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + P - Tau_xx;
    E(3,:,:) = rho.*u.*v - Tau_xy;
    E(4,:,:) = (Et + P).*u - u.*Tau_xx - v.*Tau_xy + qdot_x;
    %%%%%%%% F-term %%%%%%%%
    % Derivatives
    ddxi = @ddx_central;  % xi-derivatives
    ddeta = @ddy_bwd;     % eta-derivatives
    % values
    Tau_xx = 2*mu.*(ddx(u,ddxi,ddeta)  - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_yy = 2*mu.*(ddy(v,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_xy = mu.*(ddy(u,ddxi,ddeta)+ddx(v,ddxi,ddeta));
    qdot_y = -k.*ddy(T,ddxi,ddeta);
    % update F-term
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v - Tau_xy;
    F(3,:,:) = rho.*v.^2 + P - Tau_yy;
    F(4,:,:) = (Et + P).*v - v.*Tau_yy - u.*Tau_xy + qdot_y;
    % Adjust Coordinate frame
    for j = 1:4
        % Grid Transform
        Big_E(j,:,:) = y_eta.*squeeze(E(j,:,:)) - x_eta.*squeeze(F(j,:,:));
        Big_F(j,:,:) = -y_xi.*squeeze(E(j,:,:)) + x_xi.*squeeze(F(j,:,:));
        % Predictor derivatives (forward)
        dE_dx(j,:,:) = ddx_fwd(squeeze(Big_E(j,:,:)),dxi );
        dF_dy(j,:,:) = ddy_fwd(squeeze(Big_F(j,:,:)),deta);
        % Solve Predictor step
        temp(j,:,:) = (Jac.*squeeze(U(j,:,:)) - dt*squeeze(dE_dx(j,:,:)) - dt*squeeze(dF_dy(j,:,:)))./Jac;
    end
    [~,u,v,T,P,~,~] = cons2prim(temp,R,Cv);
    % Enforce boundary conditions & update variables
    [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond);
    rho = P./(R*T);
    e = Cv*T;
    mu = sutherland(T);
    k = (Cp/Pr)*mu;
    % Output BigU_bar
    U_bar = prim2cons(rho,u,v,T,Cv);
    
    %% CORRECTOR STEP
    [rho,u,v,T,P,~,Et] = cons2prim(U_bar,R,Cv);
    % E-term
    ddxi = @ddx_fwd;      % xi-derivatives
    ddeta = @ddy_central; % eta-derivatives
    Tau_xx = 2*mu.*(ddx(u,ddxi,ddeta)  - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_yy = 2*mu.*(ddy(v,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_xy = mu.*(ddy(u,ddxi,ddeta)+ddx(v,ddxi,ddeta));
    qdot_x = -k.*ddx(T,ddxi,ddeta);
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + P - Tau_xx;
    E(3,:,:) = rho.*u.*v - Tau_xy;
    E(4,:,:) = (Et + P).*u - u.*Tau_xx - v.*Tau_xy + qdot_x;
    % F-term
    ddxi = @ddx_central;  % xi-derivatives
    ddeta = @ddy_fwd;     % eta-derivatives
    Tau_xx = 2*mu.*(ddx(u,ddxi,ddeta)  - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_yy = 2*mu.*(ddy(v,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_xy = mu.*(ddy(u,ddxi,ddeta)+ddx(v,ddxi,ddeta));
    qdot_y = -k.*ddy(T,ddxi,ddeta);
    F(1,:,:) = rho.*v;
    F(2,:,:) = rho.*u.*v - Tau_xy;
    F(3,:,:) = rho.*v.^2 + P - Tau_yy;
    F(4,:,:) = (Et + P).*v - v.*Tau_yy - u.*Tau_xy + qdot_y;
    % Adjust Coordinate frame
    for j = 1:4
        % Grid Transform
        Big_E(j,:,:) = y_eta.*squeeze(E(j,:,:)) - x_eta.*squeeze(F(j,:,:));
        Big_F(j,:,:) = - y_xi.*squeeze(E(j,:,:)) + x_xi.*squeeze(F(j,:,:));
        % Correct derivatives (backward)
        dE_dx(j,:,:) = ddx_bwd(squeeze(Big_E(j,:,:)),dxi);
        dF_dy(j,:,:) = ddy_bwd(squeeze(Big_F(j,:,:)),deta);
        % Solve Corrector Step
        temp(j,:,:) = (0.5*(Jac.*squeeze(U(j,:,:)) + Jac.*squeeze(U_bar(j,:,:))) ...
            - 0.5*dt*squeeze(dE_dx(j,:,:)) ...
            - 0.5*dt*squeeze(dF_dy(j,:,:)))./Jac;
    end
    % Solve Corrector Step
    [~,u,v,T,P,~,Et] = cons2prim(temp,R,Cv);
    % Enforce boundary conditions & update variables
    [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond);
    rho = P./(R*T);
    e = Cv*T;
    mu = sutherland(T);
    k = (Cp/Pr)*mu;
    % Output next step
    U_next = prim2cons(rho,u,v,T,Cv);

    % My convergence criteria (dU = U_n+1 - U_n)
    for j = 1:4
        dU(i,j) = norm(squeeze(U_next(j,:,:) - U(j,:,:)));
    end
    fprintf(sprintf('Convergence criteria: [%1.1d %1.1d %1.1d %1.1d] \n',dU(i,:)))
    
    % Visualize every 100th point
    if mod(i,50)==0 || i==1
        spatial;
        visualize(u,v,P,T,rho,Cv,xx,yy)
        i_frame = getframe(spatial);
        im = frame2im(i_frame);
        [imind,cm] = rgb2ind(im,256);
        if frame == 1
            imwrite(imind,cm,'Grid_Transform.gif','gif', 'Loopcount',inf,...
            'DelayTime',0.1);
        else
            imwrite(imind,cm,'Grid_Transform.gif','gif','WriteMode','append',...
            'DelayTime',0.1);
        end
        frame = frame + 1;
    end

    % Update loop
    U = U_next;

    % Break if criteria's standard deviation 
    % falls below machine tolerance. This means
    % that the solution is not changing, a.k.a dU_dt = 0
    if i>20 && all(std((dU(i-10:i,:))./dU(1,:)) < eps)
        fprintf('Convergence tolerance reached\n')
        break
    end
end

% Visualized last step
visualize(u,v,P,T,rho,Cv,xx,yy)
fprintf('Finished\n')

%% Useful functions
function R = find_imag(P,xx,yy)
if nnz(imag(P)) ~= 0
    plot(xx(imag(P)~=0),yy(imag(P)~=0),'r.','markersize',20)
    hold on
    plot(xx,yy,'k.')
end
end
% Enforce the boundary contions on primitive variables
function [u,v,P,T] = enforce_bcs(u,v,P,T,mode)
    global Ut Pt Tt
    % Left /// Choked Flow
    u(1,:) = 4*Ut;
    v(1,:) = 0;
    P(1,:) = Pt;
    T(1,:) = Tt;
    % Right /// Extrapolation
    u(end,:) = 2*u(end-1,:) - u(end-2,:);
    v(end,:) = 2*v(end-1,:) - v(end-2,:);
    P(end,:) = 2*P(end-1,:) - P(end-2,:);
    T(end,:) = 2*T(end-1,:) - T(end-2,:);
    % Top /// Nozzle Wall
    u(:,end) = 4*Ut;
    v(:,end) = 0;
    P(:,end) = Pt;
    T(:,end) = Tt;
    % Bottom /// Nozzle Wall
    u(:,1) = 0;
    v(:,1) = 0;
    P(2:end,1) = 2*P(2:end,2) - P(2:end,3); % Extrapolation
    P(1,1) = Pt;                         % Bottom left corner
    switch mode % Different wall boundary conditions
        case 'standard' % Constant temperature
            T(:,1) = Tt;
        case 'adiabatic' % Adiabatic (dT_dy|wall = 0)
            T(:,1) = (4*T(:,2) - T(:,3))/3; % 2nd order accurate forward difference method
    end
end
% Visualize the primitive variables
function visualize(u,v,P,T,rho,Cv,xx,yy)
    % density
    subplot(2,3,1)
    plot = pcolor(xx,yy,rho);
    set(plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'\rho [kg/m^3]');
    ylabel('y'); xlabel('x');
    title('\rho [kg/m^3]')
    axis equal tight
    shading interp
    % u-velocity
    subplot(2,3,2)
    plot = pcolor(xx,yy,u);
    set(plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'u [m/s]');
    ylabel('y'); xlabel('x');
    title('u [m/s]')
    axis equal tight
    shading interp
    % v-velocity
    subplot(2,3,3)
    plot = pcolor(xx,yy,v);
    set(plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'v [m/s]');
    ylabel('y'); xlabel('x');
    title('v [m/s]')
    axis equal tight
    shading interp
    % Energy
    subplot(2,3,4)
    plot = pcolor(xx,yy,Cv*T);
    set(plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'e [J/kg]');
    ylabel('y'); xlabel('x');
    title('e [J/kg]')
    axis equal tight
    shading interp
    % Pressure
    subplot(2,3,5)
    plot = pcolor(xx,yy,P);
    set(plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'P [N/m^2]');
    ylabel('y'); xlabel('x');
    title('P [N/m^2]')
    axis equal tight
    shading interp
    % Temperature
    subplot(2,3,6)
    plot = pcolor(xx,yy,T);
    set(plot, 'EdgeColor', 'none');
    colormap('turbo')
    cb = colorbar; ylabel(cb,'T [K]');
    ylabel('y'); xlabel('x');
    title('T [K]')
    axis equal tight
    shading interp
    colormap jet

    % Update plot
    drawnow
end