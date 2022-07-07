%%
clear; close all; clc;
% Adjustments to the Pressure and Grid Resolution
P0_adjust = 1.7;
num = 25;
% Import FD method & Olson's code
addpath('functions')
addpath('NOT_MY_CODE')
% Useful globals for functions
global Ut Pt Tt P_ex M_ex T_atm
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
% BC_cond = 'standard';
BC_cond = 'standard';
SF = 2;           % Safety Factor
dt_min = 2.35e-11;  % Minimum time step
% Chamber Conditions & Throat size
P_0 = 1.2e6;    % Combustion Chamber pressure
T_0 = 2000;     % Combustion Chamber temperature
At = 1e-5;      % Throat height
num = num;      % Number of characteristic lines ~ nx1,ny1
% Geometry from Olson's Code
[xx,yy] = genmesh(T_0,P_0,At,num,R,gamma,P_atm);
nx1 = size(xx,1); ny1 = size(xx,2);
disp(nx1); disp(ny1);
dx = max(abs(diff(xx)),[],'all'); 
dy = max(abs(diff(yy)),[],'all');
% Pressure adjustments post-geometry
P_0 = P0_adjust*P_0;
% Theoretic solution
A_At = yy(1,end)./yy(:,end);
func = @(M) M.*(2/(gamma+1)*(1+(gamma-1)/2*M.^2))^(-(gamma+1)/(2*gamma-2));
M_ex = zeros(nx1,1);
P_ex = zeros(nx1,1);
for ii = 1:nx1
    M_ex(ii) = fzero(@(x) func(x)-A_At(ii),[1 3]);
    P_ex(ii) = P_0*((1+(gamma/2-1/2)*M_ex(ii)*M_ex(ii)))^(-gamma/(gamma-1));
end
% Geometry / Conditions
Ae = 2*yy(end,end);
L = xx(end,1);
% Parametric geometry
[xi,eta] = ndgrid(linspace(0,L,nx1),linspace(0,Ae,ny1));
dxi = xi(2,1)-xi(1,1); deta = eta(1,2)-eta(1,1);
% Throat Conditions
Tt = T_0/(1+(gamma-1)/2);
Pt = P_0*(2/(gamma+1))^(gamma/(gamma-1));
Ut = sqrt(gamma*R*Tt);
% Initial Conditions
sos = sqrt(gamma*R*Tt);
u = zeros(nx1,ny1);
for ii = 1:nx1
    % Exact Initialization
    u(ii,:) = M_ex(ii)*sqrt(gamma*R*Tt)*ones(1,ny1);
    P(ii,:) = P_ex(ii)*ones(1,ny1);
    % Linear Initialization
    u(ii,:) = sqrt(gamma*R*Tt)*ones(1,ny1);
    P(ii,:) = P_atm*ones(1,ny1);
end
v = zeros(nx1,ny1);
T = Tt*ones(nx1,ny1);
[u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond);
rho = P./(R*T);
U = prim2cons(rho,u,v,T,Cv);
mu = sutherland(T);
k = (Cp/Pr)*mu;
% Initialize figures
spatial = figure('Position',[10 10 1200 725]);
visualize(u,v,P,T,rho,Cv,xx,yy)
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
Big_U = zeros(4,nx1,ny1);
Big_E = zeros(4,nx1,ny1);
Big_F = zeros(4,nx1,ny1);
U_bar = zeros(4,nx1,ny1);
U_next = zeros(4,nx1,ny1);
Big_U_bar = zeros(4,nx1,ny1);
Tau_xx = zeros(nx1,ny1);
Tau_xy = zeros(nx1,ny1);
Tau_yy = zeros(nx1,ny1);
qdot_x = zeros(nx1,ny1);
qdot_y = zeros(nx1,ny1);
dU = zeros(1,4);
frame = 1;
for i = 1:5000
    %% Time-Step set to minimum
    dt = dt_min/SF;
    %% PREDICTOR STEP
    % Update Primitives
    [rho,u,v,T,P,~,Et] = cons2prim(U,R,Cv);
    % E-term
    ddxi  = @ddx_bwd;     % xi-derivatives
    ddeta = @ddy_central; % eta-derivatives
    Tau_xx = 2*mu.*(ddx(u,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_yy = 2*mu.*(ddy(v,ddxi,ddeta) - 1/3*(ddx(u,ddxi,ddeta)+ddy(v,ddxi,ddeta)));
    Tau_xy = mu.*(ddy(u,ddxi,ddeta)+ddx(v,ddxi,ddeta));
    qdot_x = -k.*ddx(T,ddxi,ddeta);
    E(1,:,:) = rho.*u;
    E(2,:,:) = rho.*u.^2 + P - Tau_xx;
    E(3,:,:) = rho.*u.*v - Tau_xy;
    E(4,:,:) = (Et + P).*u - u.*Tau_xx - v.*Tau_xy + qdot_x;
    % F-term
    ddxi = @ddx_central;  % xi-derivatives
    ddeta = @ddy_bwd;     % eta-derivatives
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
        Big_E(j,:,:) = y_eta.*squeeze(E(j,:,:)) - x_eta.*squeeze(F(j,:,:));
        Big_F(j,:,:) = -y_xi.*squeeze(E(j,:,:)) + x_xi.*squeeze(F(j,:,:));
    end
    % Predictor derivatives (forward)
    for j = 1:4
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
        Big_E(j,:,:) = y_eta.*squeeze(E(j,:,:)) - x_eta.*squeeze(F(j,:,:));
        Big_F(j,:,:) = - y_xi.*squeeze(E(j,:,:)) + x_xi.*squeeze(F(j,:,:));
    end
    % Correct derivatives (backward)
    for j = 1:4
        dE_dx(j,:,:) = ddx_bwd(squeeze(Big_E(j,:,:)),dxi);
        dF_dy(j,:,:) = ddy_bwd(squeeze(Big_F(j,:,:)),deta);
        % Solve Corrector Step
        temp(j,:,:) = (0.5*(Jac.*squeeze(U(j,:,:)) +  Jac.*squeeze(U_bar(j,:,:))) ...
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
    if mod(i,50)==0
        fprintf(sprintf('Convergence criteria: [%1.1d %1.1d %1.1d %1.1d] \n',dU(i,:)./dU(1,:)))
    end
    % Visualize data every so often
    if mod(i,500)==0 || i==1
        fprintf(sprintf('----- iter %d ----- \n',i))
        figure(1)
        visualize(u,v,P,T,rho,Cv,xx,yy)
        i_frame = getframe(spatial);
        im = frame2im(i_frame);
        [imind,cm] = rgb2ind(im,256);
        if frame == 1
            imwrite(imind,cm,'nozzle_sim.gif','gif', 'Loopcount',inf,...
            'DelayTime',0.1);
        else
            imwrite(imind,cm,'nozzle_sim.gif','gif','WriteMode','append',...
            'DelayTime',0.1);
        end
        frame = frame + 1;
    end

    % Update loop
    U = U_next;
end
close
save('noz_init.mat')
disp('done')


%% Useful functions
% Enforce the boundary contions on primitive variables
function [u,v,P,T] = enforce_bcs(u,v,P,T,mode)
    global Ut Pt Tt T_atm
    % Left /// Choked Flow
    u(1,:) = Ut;
    v(1,:) = 0;
    P(1,:) = Pt;
    T(1,:) = Tt;
    % Right /// Extrapolation
    u(end,2:end-1) = 2*u(end-1,2:end-1) - u(end-2,2:end-1);
    v(end,2:end-1) = 2*v(end-1,2:end-1) - v(end-2,2:end-1);
    P(end,2:end-1) = 2*P(end-1,2:end-1) - P(end-2,2:end-1);
    T(end,2:end-1) = 2*T(end-1,2:end-1) - T(end-2,2:end-1);
    % Top /// Nozzle Wall
    u(:,end) = 0;
    v(:,end) = 0;
    P(2:end,end) = 2*P(2:end,end-1) - P(2:end,end-2); % Extrapolation
    T(2:end,end) = 2*T(2:end,end-1) - T(2:end,end-2);
    switch mode % Different wall boundary conditions
        case 'standard' % Constant temperature
            T(:,end) = Tt;
        case 'adiabatic' % Adiabatic (dT_dy|wall = 0)
            T(:,end) = (4*T(:,end-1) - T(:,end-2))/3; % 2nd order accurate forward difference method
        case 'test'
            T(:,end) = T(:,end-1)/2 + T_atm/2;
            T(:,end) = P(:,end).*T(:,end-1)./P(:,end-1);
    end
    % Bottom /// Nozzle Wall
    u(:,1) = 0;
    v(:,1) = 0;
    P(2:end,1) = 2*P(2:end,2) - P(2:end,3); % Extrapolation
    T(2:end,1) = 2*T(2:end,2) - T(2:end,3);
    switch mode % Different wall boundary conditions
        case 'standard' % Constant temperature
            T(:,1) = Tt;
        case 'adiabatic' % Adiabatic (dT_dy|wall = 0)
            T(:,1) = (4*T(:,2) - T(:,3))/3; % 2nd order accurate forward difference method
        case 'test'
            T(:,1) = T(:,2)/2 + T_atm/2;
            T(:,1) = P(:,1).*T(:,2)./P(:,2);
    end
end
% Visualize the primitive variables
function visualize(u,v,P,T,rho,Cv,xx,yy)
    global M_ex P_ex
    clf
    pltfmt = {'interpreter','latex',...
        'fontsize',16};
    % density
    subplot(3,3,1)
    temp_plot = pcolor(xx,yy,rho);
    set(temp_plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'\rho [kg/m^3]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('$\rho [kg/m^3]$',pltfmt{:})
    axis equal tight
    shading interp
    % u-velocity
    subplot(3,3,2)
    temp_plot = pcolor(xx,yy,u);
    set(temp_plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'u [m/s]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('u [m/s]',pltfmt{:})
    axis equal tight
    shading interp
    % v-velocity
    subplot(3,3,3)
    temp_plot = pcolor(xx,yy,v);
    set(temp_plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'v [m/s]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('v [m/s]',pltfmt{:})
    axis equal tight
    shading interp
    % Energy
    subplot(3,3,4)
    temp_plot = pcolor(xx,yy,Cv*T);
    set(temp_plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'e [J/kg]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('e [J/kg]',pltfmt{:})
    axis equal tight
    shading interp
    % Pressure
    subplot(3,3,5)
    temp_plot = pcolor(xx,yy,P);
    set(temp_plot, 'EdgeColor', 'none');
    cb = colorbar; ylabel(cb,'P [N/m^2]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('$P [N/m^2]$',pltfmt{:})
    axis equal tight
    shading interp
    % Temperature
    subplot(3,3,6)
    temp_plot = pcolor(xx,yy,T);
    set(temp_plot, 'EdgeColor', 'none');
    colormap('turbo')
    cb = colorbar; ylabel(cb,'T [K]');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('T [K]',pltfmt{:})
    axis equal tight
    shading interp
    colormap jet
    % Mach
    subplot(3,3,7)
    Mach = sqrt(u.^2+v.^2)./sqrt(1.4*287*T);
    temp_plot = pcolor(xx,yy,Mach);
    set(temp_plot, 'EdgeColor', 'none');
    colormap('turbo')
    cb = colorbar; ylabel(cb,'Mach');
    ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
    title('Mach Number',pltfmt{:})
    axis equal tight
    shading interp
    colormap jet
    % Centerline
    subplot(3,3,8:9)
    ind = round(size(u,2)/2);
    plot(xx(:,ind),Mach(:,ind),'-','linewidth',1.5,'color','#0072BD');
    hold on;
    plot(xx(:,1),M_ex,'--','linewidth',1.5,'color','#0072BD')
    plot(xx(:,ind),P(:,ind)/101300,'-','linewidth',1.5,'color','#D95319');
    plot(xx(:,1),P_ex/101300,'--','linewidth',1.5,'color','#D95319')
    xlabel('x',pltfmt{:});
    title('Centerline Flow',pltfmt{:})
    legend('Mach Number','Theory','$P/P_{atm}$','Theory','location','northeastoutside',pltfmt{:})
    axis([0 max(xx(:,ind)) 0 max([max(P_ex/101300),max(M_ex)])])
    yticks(0:round(max(P(:,ind)/101300)))
    grid on
    % Update plot
    drawnow
end