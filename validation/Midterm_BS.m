%%
clear all; close all; clc;
% Import FD methods
addpath('functions')
% Useful globals for functions
global U_inf P_inf T_inf
% Properties of standard air
T_inf = 288.15; % K
P_inf = 101300; % N/m2
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
Mach = 4;           % Inlet airspeed
L = 1e-5; H = 8e-6; % Domain Dimensions
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
% Defining grids of each block
nx = {35,35}; 
ny = {80,80}; % equal heights
[xx{1},yy{1}] = ndgrid(linspace(0,L/2,nx{1}),linspace(0,H,ny{1}));
dx{1} = xx{1}(2,1)-xx{1}(1,1);
[xx{2},yy{2}] = ndgrid(linspace(L/2-dx{1},L-dx{1},nx{1}),linspace(0,H,ny{1}));
for ii = 1:num_block_rows
    for jj = 1:num_block_cols
        ind = ii + (num_block_rows)*(jj-1);
        dx{ind} = xx{ind}(2,1)-xx{ind}(1,1); 
        dy{ind} = yy{ind}(1,2)-yy{ind}(1,1);
    end
end
% Initial Conditions
sos = sqrt(gamma*R*T_inf);  % speed of sound
U_inf = Mach*sos;           % airspeed of Mach 4
for blk = 1:2
    u{blk} = U_inf*ones(nx{blk},ny{blk});      % Initial u velocity
    v{blk} = zeros(nx{blk},ny{blk});           % Initial v velocity
    P{blk} = P_inf*ones(nx{blk},ny{blk});      % Initial Pressure
    T{blk} = T_inf*ones(nx{blk},ny{blk});      % Initial Temperature
end
for blk = 1:2
    [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk); % Enforce Boundary Conditions
    rho{blk} = P{blk}./(R*T{blk});             % Initial Density
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
U_old = empty; Et = empty;
for ii = 1:num_block_rows
    for jj = 1:num_block_cols
        ind = ii-1 + (num_block_rows)*jj;
        dE_dx{ii} = zeros(4,nx{ii},ny{ii});
        dF_dy{ii} = zeros(4,nx{ii},ny{ii});
        E{ii} = zeros(4,nx{ii},ny{ii});
        F{ii} = zeros(4,nx{ii},ny{ii});
        U_bar{ii} = zeros(4,nx{ii},ny{ii});
        U_next{ii} = zeros(4,nx{ii},ny{ii});
        Tau_xx{ii} = zeros(nx{ii},ny{ii});
        Tau_xy{ii} = zeros(nx{ii},ny{ii});
        Tau_yy{ii} = zeros(nx{ii},ny{ii});
        qdot_x{ii} = zeros(nx{ii},ny{ii});
        qdot_y{ii} = zeros(nx{ii},ny{ii});
        dU{ii} = zeros(1,4);
    end
end
% Visualize Blocks
viz1 = figure('Position',[10 10 600 400]);
plot(xx{1}(1,1),yy{1}(1,1),'r.','MarkerSize',8)
hold on
plot(xx{2},yy{2},'k.','MarkerSize',8)
plot(xx{1},yy{1},'r.','MarkerSize',8)
title('Block Structure','Interpreter','latex','FontSize',14)
l = legend('Block 1','Block 2','Interpreter','latex','FontSize',14,'Location','NorthEastOutside');
xlabel('x','Interpreter','latex','FontSize',14)
ylabel('y','Interpreter','latex','FontSize',14)
axis tight equal
viz2 = figure('Position',[10 10 600 400]);
plot(xx{1}(end,1)/10,yy{1}(end,1),'ro','MarkerSize',12)
hold on
plot(xx{2}(1:3,1:12)/10,yy{2}(1:3,1:12),'k.','MarkerSize',16)
plot(xx{1}(end-2:end,1:12)/10,yy{1}(end-2:end,1:12),'ro','MarkerSize',12)
title('Block Structure Overlap','Interpreter','latex','FontSize',14)
l = legend('Block 1','Block 2','Interpreter','latex','FontSize',14,'Location','NorthEastOutside');
xlabel('x','Interpreter','latex','FontSize',14)
ylabel('y','Interpreter','latex','FontSize',14)
axis tight
exportgraphics(viz1,'block_structure.pdf')
exportgraphics(viz2,'block_overlap.pdf')
close all;
myfig = figure('Position', [10 10 1200 500]);
frame = 1;
for i = 1:1500
    %% Get Time-Step
    dt = dt_min/SF;
    U_old = U;
    % Iterate over all blocks
    for jj = 1:num_block_cols
        for ii = 1:num_block_rows
            ind = [ii,jj];
            blk = ii + (num_block_rows)*(jj-1);
            U = U_old;
            %% PREDICTOR STEP
            % Update Primitives
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U{blk},R,Cv);
            % E-term
            if blk == 1
                ddx = @ddx_bwd;
                ddy = @ddy_central;
            elseif blk == 2
                ddx = @ddx_bwd;
                ddy = @ddy_central;
            end
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_x{blk} = -k{blk}.*ddx(T{blk},dx{blk});
            E{blk}(1,:,:) = rho{blk}.*u{blk};
            E{blk}(2,:,:) = rho{blk}.*u{blk}.^2 + P{blk} - Tau_xx{blk};
            E{blk}(3,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            E{blk}(4,:,:) = (Et{blk} + P{blk}).*u{blk} - u{blk}.*Tau_xx{blk} ...
                - v{blk}.*Tau_xy{blk} + qdot_x{blk};
            % F-term
            if blk == 1
                ddx = @ddx_central;
                ddy = @ddy_bwd;
            elseif blk == 2
                ddx = @ddx_central;
                ddy = @ddy_bwd;
            end
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_y{blk} = -k{blk}.*ddy(T{blk},dy{blk});
            F{blk}(1,:,:) = rho{blk}.*v{blk};
            F{blk}(2,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            F{blk}(3,:,:) = rho{blk}.*v{blk}.^2 + P{blk} - Tau_yy{blk};
            F{blk}(4,:,:) = (Et{blk} + P{blk}).*v{blk} - v{blk}.*Tau_yy{blk} ...
                - u{blk}.*Tau_xy{blk} + qdot_y{blk};
            % Predictor derivatives (forward)
            for j = 1:4
                if blk == 1
                    dE_dx{blk}(j,:,:) = ddx_fwd(squeeze(E{blk}(j,:,:)),dx{blk});
                elseif blk == 2
                    dE_dx{blk}(j,:,:) = ddx_fwd(squeeze(E{blk}(j,:,:)),dx{blk});
                end
                dF_dy{blk}(j,:,:) = ddy_fwd(squeeze(F{blk}(j,:,:)),dy{blk});
            end
            % Solve Predictor step
            temp = U{blk} - dt*dE_dx{blk} - dt*dF_dy{blk};
            [~,u{blk},v{blk},T{blk},P{blk},~,~] = cons2prim(temp,R,Cv);
            % Enforce boundary conditions & update variables
            [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output U_bar
            U_bar{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
            
            %% CORRECTOR STEP
            [rho{blk},u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_bar{blk},R,Cv);
            % E-term
            if blk == 1
                ddx = @ddx_fwd;
                ddy = @ddy_central;
            elseif blk == 2
                ddx = @ddx_fwd;
                ddy = @ddy_central;
            end
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_x{blk} = -k{blk}.*ddx(T{blk},dx{blk});
            E{blk}(1,:,:) = rho{blk}.*u{blk};
            E{blk}(2,:,:) = rho{blk}.*u{blk}.^2 + P{blk} - Tau_xx{blk};
            E{blk}(3,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            E{blk}(4,:,:) = (Et{blk} + P{blk}).*u{blk} - u{blk}.*Tau_xx{blk} ...
                - v{blk}.*Tau_xy{blk} + qdot_x{blk};
            % F-term
            if blk == 1
                ddx = @ddx_central;
                ddy = @ddy_fwd;
            elseif blk == 2
                ddx = @ddx_central;
                ddy = @ddy_fwd;
            end
            Tau_xx{blk} = 2*mu{blk}.*(ddx(u{blk},dx{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_yy{blk} = 2*mu{blk}.*(ddy(v{blk},dy{blk}) ...
                - 1/3*(ddx(u{blk},dx{blk})+ddy(v{blk},dy{blk})));
            Tau_xy{blk} = mu{blk}.*(ddy(u{blk},dy{blk})+ddx(v{blk},dx{blk}));
            qdot_y{blk} = -k{blk}.*ddy(T{blk},dy{blk});
            F{blk}(1,:,:) = rho{blk}.*v{blk};
            F{blk}(2,:,:) = rho{blk}.*u{blk}.*v{blk} - Tau_xy{blk};
            F{blk}(3,:,:) = rho{blk}.*v{blk}.^2 + P{blk} - Tau_yy{blk};
            F{blk}(4,:,:) = (Et{blk} + P{blk}).*v{blk} - v{blk}.*Tau_yy{blk} ...
                - u{blk}.*Tau_xy{blk} + qdot_y{blk};
            % Correct derivatives (backward)
            for j = 1:4
                dE_dx{blk}(j,:,:) = ddx_bwd(squeeze(E{blk}(j,:,:)),dx{blk});
                dF_dy{blk}(j,:,:) = ddy_bwd(squeeze(F{blk}(j,:,:)),dy{blk});
            end
            % Solve Corrector Step
            temp = 0.5*(U{blk} + U_bar{blk}) - 0.5*dt*dE_dx{blk} - 0.5*dt*dF_dy{blk};
            [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(temp,R,Cv);
            % Enforce boundary conditions & update variables
            [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
            rho{blk} = P{blk}./(R*T{blk});
            mu{blk} = sutherland(T{blk});
            k{blk} = (Cp/Pr)*mu{blk};
            % Output next step
            U_next{blk} = prim2cons(rho{blk},u{blk},v{blk},T{blk},Cv);
        end
    end
    
    for blk = 1:2
        [~,u{blk},v{blk},T{blk},P{blk},~,Et{blk}] = cons2prim(U_next{blk},R,Cv);
        [u,v,P,T] = enforce_bcs(u,v,P,T,BC_cond,blk);
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
    if mod(i,50)==0 || i==1
        visualize(u,v,P,T,rho,Cv,xx,yy)
        i_frame = getframe(myfig);
        im = frame2im(i_frame);
        [imind,cm] = rgb2ind(im,256);
        if frame == 1
            imwrite(imind,cm,'Block_Structure.gif','gif', 'Loopcount',inf,...
            'DelayTime',0.1);
        else
            imwrite(imind,cm,'Block_Structure.gif','gif','WriteMode','append',...
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
% Enforce the boundary contions on primitive variables
function [u,v,P,T] = enforce_bcs(u,v,P,T,mode,blk)
    global U_inf P_inf T_inf
    if blk == 1
        % Left /// Far Field
        u{blk}(1,:) = U_inf;
        v{blk}(1,:) = 0;
        P{blk}(1,:) = P_inf;
        T{blk}(1,:) = T_inf;
        % Connectivity
        u{blk}(end,:) = u{2}(2,:);
        v{blk}(end,:) = v{2}(2,:);
        P{blk}(end,:) = P{2}(2,:);
        T{blk}(end,:) = T{2}(2,:);
    elseif blk == 2
        % Right /// Extrapolation
        u{blk}(end,:) = 2*u{blk}(end-1,:) - u{blk}(end-2,:);
        v{blk}(end,:) = 2*v{blk}(end-1,:) - v{blk}(end-2,:);
        P{blk}(end,:) = 2*P{blk}(end-1,:) - P{blk}(end-2,:);
        T{blk}(end,:) = 2*T{blk}(end-1,:) - T{blk}(end-2,:);
        % Connectivity
        u{blk}(1,:) = u{1}(end-1,:);
        v{blk}(1,:) = v{1}(end-1,:);
        P{blk}(1,:) = P{1}(end-1,:);
        T{blk}(1,:) = T{1}(end-1,:);
    end
    % Top /// Far Field
    u{blk}(:,end) = U_inf;
    v{blk}(:,end) = 0;
    P{blk}(:,end) = P_inf;
    T{blk}(:,end) = T_inf;
    % Bottom /// Wall
    u{blk}(:,1) = 0;
    v{blk}(:,1) = 0;
    P{blk}(2:end,1) = 2*P{blk}(2:end,2) - P{blk}(2:end,3); % Extrapolation
    P{blk}(1,1) = P_inf;                         % Bottom left corner
    switch mode % Different wall boundary conditions
        case 'standard' % Constant temperature
            T{blk}(:,1) = T_inf;
        case 'adiabatic' % Adiabatic (dT_dy|wall = 0)
            T{blk}(:,1) = (4*T(:,2) - T(:,3))/3; % 2nd order accurate forward difference method
    end
end
% Visualize the primitive variables
function visualize(u,v,P,T,rho,Cv,xx,yy)
    clf
    for ii = 1:numel(u)
        % density
        subplot(2,3,1)
        hold on;
        plot = pcolor(xx{ii},yy{ii},rho{ii});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'\rho [kg/m^3]');
        ylabel('y'); xlabel('x');
        title('\rho [kg/m^3]')
        axis equal tight
        shading interp
        % u-velocity
        subplot(2,3,2)
        hold on;
        plot = pcolor(xx{ii},yy{ii},u{ii});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'u [m/s]');
        ylabel('y'); xlabel('x');
        title('u [m/s]')
        axis equal tight
        shading interp
        % v-velocity
        subplot(2,3,3)
        hold on;
        plot = pcolor(xx{ii},yy{ii},v{ii});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'v [m/s]');
        ylabel('y'); xlabel('x');
        title('v [m/s]')
        axis equal tight
        shading interp
        % Energy
        subplot(2,3,4)
        hold on;
        plot = pcolor(xx{ii},yy{ii},Cv*T{ii});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'e [J/kg]');
        ylabel('y'); xlabel('x');
        title('e [J/kg]')
        axis equal tight
        shading interp
        % Pressure
        subplot(2,3,5)
        hold on;
        plot = pcolor(xx{ii},yy{ii},P{ii});
        set(plot, 'EdgeColor', 'none');
        cb = colorbar; ylabel(cb,'P [N/m^2]');
        ylabel('y'); xlabel('x');
        title('P [N/m^2]')
        axis equal tight
        shading interp
        % Temperature
        subplot(2,3,6)
        hold on;
        plot = pcolor(xx{ii},yy{ii},T{ii});
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