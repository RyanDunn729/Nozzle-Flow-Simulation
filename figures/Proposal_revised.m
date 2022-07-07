%   Britton Jeffrey Olson
%   Ph.D. Candidate
%   Stanford Univeristy
%   Department of Aero/Astro

%%%%%%%%    Introduction and Background     %%%%%%%%
%   This program gives the ideal nozzle geometry using the method of
%   characteristics for a Quasi-2D Diverging Nozzle.  Assume gas is
%   exhausting from a combustion chamber that has no mass flow rate in.  
%   Using 2D nozzle flow relations, an optimal throat area is found that
%   will produce the max amount of thrust for the given ambient pressure
%   and combustion chamber parameters.  This Area is automatically set and
%   fed into the method of charactertistics portion of that code.  The
%   method of characteristics also uses the exit Mach number that
%   corresponds to the ideal exit area.  

%%%%%%%%    Directions for running the program      %%%%%%%%
%   Run default or modify the problem parameters in nozzle.m and noz_cfd.m
%   1)  Run the nozzle.m program-  This will design the nozzle
%   2)  Run noz_mesh.m-  This will set up a grid for the newly designed
%   nozzle
%   3)  Run noz_cfd.m-  This will solve the 2-D Euler equations on the
%   nozzle mesh and plot the result.
%   
%   ***Note*** 
%       You must run them in this order without clearing variables.
%   	Each script depends on the one before it.  Also, solver has the
%   	functions used by noz_cfd.m and needs to be in the same directory.
%
%   Figure(1): Static thrust as function of exit area
%   Figure(2): Nozzle design and plots of Mach number & Pressure vs Length
%   Figure(3): CFD simulation of designed nozzle

clear all;
clc;

%   Problem parameters
T_c = 2000;     % Temperature in the combustion chamber (K)
P_c = 1.2e6;    % Pressure in the combustion chamber (Pa)
P_amb = 101e3;  % Ambient pressure (Pa)
T_amb = 288.15;    % Ambient temperature (K)
gamma = 1.4;   % Ratio of Specific Heats Cp/Cv (Gamma)
R = 287;
width = 2;     % Nozzle width (meters)
h_th = 1e-5;    % Throat  height (meters)

%   Method of Characteristics
num = 25;       % Number of Characteristic lines
theta_i = .03;  % Initial step in theta 
plotter = 1;    % Set to '1' to plot nozzle

dh = h_th/100;
max_iter = 10000;

%   Part A

%find where P becomes u
h(1) = h_th;
A_star = h_th*width;
M =1;
dM1 = .1;
for i=1: max_iter
    h(i) = h(1) + (i-1)*dh;
    Ae(i) = h(i)*width;
    A_Asq = (Ae(i)/A_star)^2;
    A_ratio(i)=sqrt(A_Asq);
    
    %Newton Rhapson on Eq. 5.20 - Anderson text
    res = 1;
    if i > 1
        M = Ma(i-1);
    end
    
     while res > .001
        M2 = M + dM1;
        funa1 = -A_Asq + (1/M^2)*((2/(gamma+1))*(1+(gamma-1)*M^2/2))^((gamma+1)/(gamma-1));
        funa2 = -A_Asq + (1/M2^2)*((2/(gamma+1))*(1+(gamma-1)*M2^2/2))^((gamma+1)/(gamma-1));
        dv_dm = (funa2-funa1)/dM1;
        
        M = M - funa1/dv_dm;
        res = abs(funa1);
        
    end
    Ma(i) = M;
    
    % Find Pressure
    P(i) = P_c*(1+(gamma-1)*Ma(i)^2/2)^(-gamma/(gamma-1));
    
    % Find thrust for each point
    Te(i) = T_c/(1+(gamma-1)*Ma(i)^2/2);
    Tt(i) = T_c/(1+(gamma-1)/2);
    Ve(i) = Ma(i)*sqrt(Te(i)*gamma*R);
    Vt(i) = sqrt(Tt(i)*gamma*R);
    rhot(i) = P(i)/(R*Te(i));
    mdot(i) = rhot(i)*Ve(i)*Ae(i);
    TT(i) = mdot(i)*Ve(i) + (P(i) - P_amb)*Ae(i);
    
    if P(i) < P_amb
        %break
        %Calculate the pressure if shock wave exists at the exit plane
        P_exit = P(i)*(1+(gamma*2/(gamma+1))*(Ma(i)^2-1));
        
         if P_exit <= P_amb
             P(i) = P_exit;
             break
         else
         end
        
    else
    end
    
end

%   Part B  
%   Determine the nominal exit area of the nozzle 
%   to maximize thrust

[a,b]=max(TT);
%   Over or Underexpand the nozzle
b = b;
A_max = Ae(b);
Max_thrust = TT(b);

%   Part C
%   Method of Characteristics

M_e = Ma(b);       %Mach number at ideal exit 


%Find theta_max by using equation 11.33
theta_max = (180/pi)*(sqrt((gamma+1)/(gamma-1))*atan((sqrt((gamma-1)*(M_e^2-1)/(gamma+1))))-atan(sqrt(M_e^2-1)))/2;

%  D_theta for each char line
del_theta = (theta_max - theta_i)/(num-1);

% Find 

for i=1:num
    %   Initialize mach numeber
    
    for j=1:num
        if i==1
            %Theta for each line (first lines)
            theta(i,j) = theta_i + del_theta*(j-1);
            nu(i,j) = theta(i,j);
            K_m(i,j) = theta(i,j) + nu(i,j);
            K_p(i,j) = theta(i,j) - nu(i,j);       
        
        elseif i > 1
        
            K_p(i,j) = -K_m(1,i);
            
            % Find Thetas
            if j >= i
                theta(i,j) = del_theta*(j-i);
            else
                %theta(i,j) = theta(j,i-1);
                theta(i,j) = theta(j,i);
                
            end
            nu(i,j) = theta(i,j) - K_p(i,j);
            K_m(i,j) = theta(i,j) + nu(i,j);
        end
    
    % Prandtl-Meyer function (using Newton Rhapson)
    dM = .1; % Leave at about .1
    if j == 1
        M_ex(i,j) = 1.00;
    else
        M_ex(i,j) = M_ex(i,j-1);
    end
    M = M_ex(i,j);    
    
    res = 1;
    while res > .01
        M2 = M + dM;
        funv1 = (-nu(i,j)*(pi/180)+(sqrt((gamma+1)/(gamma-1))*atan((sqrt((gamma-1)*(M^2-1)/(gamma+1))))-atan(sqrt(M^2-1))));
        funv2 = (-nu(i,j)*(pi/180)+(sqrt((gamma+1)/(gamma-1))*atan((sqrt((gamma-1)*(M2^2-1)/(gamma+1))))-atan(sqrt(M2^2-1))));
        dv_dm = (funv2-funv1)/dM;
        
        M = M - funv1/dv_dm;
        res = abs(funv1);
        
    end
    M_ex(i,j) = M;
    
    % Find the angle mu
    mu(i,j) = (180/pi)*asin(1/M_ex(i,j));
         
    end
    
    % Add last point to char line
    theta(i,num+1) = theta(i,num);
    nu(i,num+1) = nu(i,num);
    K_m(i,num+1) = K_m(i,num);
    K_p(i,num+1) = K_p(i,num);
end

char = zeros(num,num+1,2);
for i=1:num
    
    for j=1:num+1
        
% Draw points of intersection        
        %   Point 1 of all char lines          
        if j == 1 
            char(i,j,1) = 0;
            char(i,j,2) = h_th/2;
        end
        
        %   Where first line hits the symmetry line
        if i == 1 & j==2            
            char(i,j,1) = (-h_th/2)/tan((pi/180)*(theta(1,j-1)-mu(1,j-1)));
            char(i,j,2) = 0;
        end

        %   Where all other lines hit the symmetry line
        if j == i+1 & j>2            
              char(i,j,1) = -char(i-1,j,2)/tan((pi/180)*(.5*theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)))) + char(i-1,j,1);
              char(i,j,2) = 0;
              test(i,j) = (theta(i,j-2)-.5*(mu(i,j-2)+mu(i,j-1)));
              testpty(i,j) = char(i-1,j,2);
              testptx(i,j) = char(i-1,j,1);                
        end
            
        %   All other data points for char 1 calculated
        if i ==1 & j>2 & j ~= i+1
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(j-1,1)+theta(i,j-1))-.5*(mu(j-1,1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(1,1,2) - char(1,1,1)*C_m;
                char(1,j-1,2) - char(1,j-1,1)*C_p];
                iterm(1,:)=inv(A)*B;
                char(i,j,1) = iterm(1,2);
                char(i,j,2) = iterm(1,1);
        end
        
        %   All other points for all char lines calculated
        if i > 1 & j~=i+1 & j>2        
            C_p = tan((pi/180)*(.5*(theta(i,j-2)+theta(i,j-1))+.5*(mu(i,j-2)+mu(i,j-1))));
            C_m = tan((pi/180)*(.5*(theta(i-1,j-1)+theta(i,j-1))-.5*(mu(i-1,j-1)+mu(i,j-1))));
            A = [1,-C_m;1,-C_p];
            B = [char(i-1,j,2) - char(i-1,j,1)*C_m; char(i,j-1,2) - char(i,j-1,1)*C_p];
                
            iterm(1,:) = inv(A)*B;
            char(i,j,1) = iterm(1,2);
            char(i,j,2) = iterm(1,1);  
        end
    end
end

    
%  Fill in similar points (where char lines share points)
for i = 2:num
    for j=2:num
        char(j,i,1) = char(i-1,j+1,1);
        char(j,i,2) = char(i-1,j+1,2);
    end
end
        
% ******Make the nozzle shape and extend the char lines to wall******

%   Initial start point of the nozzle (at throat)
noz(1,1) = 0;
noz(1,2) = h_th/2;

%   Find all the points of the nozzle
for i = 2 : num
    %   Find different slopes and points to intersect
    m1 = tan((pi/180)*(theta(i-1,num)+mu(i-1,num)));    
    if i ==2
        m2 = (pi/180)*theta_max;
    else
        m2 = ((pi/180)*(theta(i-1,num+1)));
    end
    m3 = ((pi/180)*(theta(i-1,num)));
    m4 = tan((m2+m3)/2);
    
    A = [1,-m4; 1,-m1];
    B = [noz(i-1,2) - noz(i-1,1)*m4; char(i-1,num+1,2) - char(i-1,num+1,1)*m1];
                
    iterm(1,:) = inv(A)*B;
    noz(i,1) = iterm(1,2);
    noz(i,2) = iterm(1,1); 
    
    %   Extend char lines to wall
    char(i-1,num+2,1)= noz(i,1);
    char(i-1,num+2,2)= noz(i,2);
end

%Last line
m1 = tan((pi/180)*(theta(num,num)+ mu(num,num)));
m2 = ((pi/180)*(theta(num-1,num)));
m3 = ((pi/180)*(theta(num,num+1)));
m4 = tan((m2+m3)/2);

A = [1,-m4; 1,-m1];
B = [noz(num,2) - noz(num,1)*m4; char(num,num+1,2) - char(num,num+1,1)*m1];
                
iterm(1,:) = inv(A)*B;
noz(num+1,1) = iterm(1,2);
noz(num+1,2) = iterm(1,1); 

    
%   Extend char lines to wall
char(num,num+2,1)= noz(num+1,1);
char(num,num+2,2)= noz(num+1,2);

%   Find  % errors in A/A* and Mexit
error_Area = 100*(width*2*noz(num,2) - A_max)/(A_max)
error_Mach = 100*(M_e - M_ex(num,num))/M_e


%   Plot Mach Number and pressure through nozzle using the quasi-1D
%   area relations.  (Isentropic expansion through nozzle)
Mnoz(1) = 1.0;  %   Choked Flow
M = Mnoz(1);
for i=1: size(noz,1)
    Ae(i) = 2*noz(i,2)*width;
    A_Asq = (Ae(i)/A_star)^2;
    A_ratio(i)=sqrt(A_Asq);
    
    %Newton Rhapson on Eq. 5.20 - Anderson text
    res = 1;
    if i > 1
        M = Mnoz(i-1);
    
    
     while res > .001
        M2 = M + dM1;
        funa1 = -A_Asq + (1/M^2)*((2/(gamma+1))*(1+(gamma-1)*M^2/2))^((gamma+1)/(gamma-1));
        funa2 = -A_Asq + (1/M2^2)*((2/(gamma+1))*(1+(gamma-1)*M2^2/2))^((gamma+1)/(gamma-1));
        dv_dm = (funa2-funa1)/dM1;
        
        M = M - funa1/dv_dm;
        res = abs(funa1);
        
    end
    Mnoz(i) = M;
    end
    % Find Pressure 
    Pnoz(i) = P_c*(1+(gamma-1)*Mnoz(i)^2/2)^(-gamma/(gamma-1));
end



%   ****** See nozzle.m for instructions ******
%  Program to extrapolate the data points from nozzle design and make a
%  uniform grid spacing in the x-direction
%  Change nothing... simply run this script
%  Find the minimum spacing given by the method of characteristics
%  and set as the dx value
dx = 0;
for i=1: size(noz,1)-1
    len = noz(i+1,1) - noz(i,1);
    if (len < dx || i == 1)
        dx = len;
    end
end

%  Explicitly give the dx value here
dx = max(noz(:,1))/ceil(max(noz(:,1))/dx);
n = max(noz(:,1))/dx;
% len = max(noz(:,1));
% n = 50;   %  Note # of points is actually n+1 
% dx = len/n
%  Pick m points in y as some factor of x points
yfactor = .8;
m = ceil(yfactor*n);
%  Make uniform x-distribution of points
xmax = 0;
i = 1;
while xmax < max(noz(:,1))
    xmax = dx*(i-1);
    x(i,1:m) = xmax;
    i = i+1;
end

%  Make the y-points and extrapolate linearly from closest points to fit
%  the nozzle geometry
%  Initialize and assign last value
y(1:size(x,1),1:size(x,2)) = 0;
y(1,size(y,2)) = noz(1,size(noz,2));
y(size(y,1),size(y,2)) = noz(size(noz,1),size(noz,2));
for i = 1 : size(x,1)-1
    j = 1;
    while x(i,1) >= noz(j,1)
        x1 = noz(j,1);
        x2 = noz(j+1,1);
        y1 = noz(j,2);
        y2 = noz(j+1,2);
        j = j + 1;
    end
    slope = (y2 - y1)/(x2 - x1);
    y(i,size(y,2)) = y1 + slope*(x(i,1)-x1);
    %Fill in mesh
    dy = y(i,size(y,2))/(size(y,2)-1);
    for k = 1: size(y,2)
        y(i,k) = dy*(k-1);
    end
end
%Fill in mesh
dy = y(size(y,1),size(y,2))/(size(y,2)-1);
for k = 1: size(y,2)
    y(size(y,1),k) = dy*(k-1);
end

pltfmt = {'interpreter','latex','fontsize',20};
[nx1,ny1] = size(x);
[xi,eta] = ndgrid(linspace(0,max(x,[],'all'),nx1),linspace(-max(y,[],'all'),max(y,[],'all'),2*ny1));

sim1 = figure('Position',[10 10 1400 500]);
subplot(1,2,1)
plot(xi,eta,'k.','MarkerSize',10)
xlabel('$\xi$',pltfmt{:}); ylabel('$\eta$',pltfmt{:});
title('Parametric Coordinates',pltfmt{:});
axis equal tight
grid on

subplot(1,2,2)
plot(x,-y,'r.','MarkerSize',10)
hold on
plot(noz(:,1),-noz(:,2),'k','LineWidth',3)
plot(noz(:,1),noz(:,2),'k','LineWidth',3)
plot(x,y,'r.','MarkerSize',10)
plot(x,-y,'r.','MarkerSize',10)
xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
title('Spatial Coordinates',pltfmt{:});
axis equal tight
grid on

exportgraphics(sim1, 'Nozzle_coords.pdf');



y = [-y(:,end:-1:2),y];
x = repmat(x(:,1),1,size(y,2));
xx = cell(3,2); yy = cell(3,2);
dx = cell(3,2); dy = cell(3,2);
xx{2} = x;
yy{2} = y;
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

full_sim = figure('Position',[10 10 1000 900]);
hold on
colors = {'#0072BD';
          '#D95319';
          '#EDB120';
          '#7E2F8E';
          '#77AC30'};
blks = [2 4 5 6];
for ii = 1:4
    blk = blks(ii);
    plot(xx{blk}(1,1),yy{blk}(1,1),'.','Color',colors{ii},'MarkerSize',14)
end
for ii = 1:4
    blk = blks(ii);
    plot(xx{blk},yy{blk},'.','Color',colors{ii},'MarkerSize',14)
end
plot(xx{5},yy{5},'.','Color',colors{3},'MarkerSize',14)
xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
title('Discretized Model',pltfmt{:});
legend('Block 1','Block 2','Block 3','Block 4','Location','northwest','interpreter','latex','fontsize',16)
axis equal tight
grid on

exportgraphics(full_sim, 'full_sim_coords.pdf');

return


% full_sim = figure('Position',[10 10 1200 500]);
% nx2 = 16; ny2 = 28;
% x_max = max(x,[],'all');
% y_max = max(y,[],'all');
% [xx,yy] = ndgrid(linspace(x_max,2*x_max,nx2),linspace(-3*y_max,3*y_max,ny2));
% subplot(1,2,1)
% plot(xx,yy,'k.','MarkerSize',10)
% xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
% title({'Sim2: Spatial Coordinates',''},pltfmt{:});
% axis equal tight
% grid on
% 
% subplot(1,2,2)
% plot(x(1),-y(1),'r.','MarkerSize',10)
% hold on
% plot(xx,yy,'k.','MarkerSize',10)
% plot(noz(:,1),-noz(:,2),'k','LineWidth',3)
% plot(noz(:,1),noz(:,2),'k','LineWidth',3)
% plot(x,y,'r.','MarkerSize',10)
% plot(x,-y,'r.','MarkerSize',10)
% xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
% title('Full Simulation',pltfmt{:});
% legend('Sim1 nodes','Sim2 nodes','Location','northwest','interpreter','latex','fontsize',14)
% axis equal tight
% grid on
% 
% exportgraphics(full_sim, 'Sim2_coords.pdf');

% figure
% plot(noz(:,1),-noz(:,2),'k','LineWidth',3)
% hold on;
% plot(a(1),-A_max/width/2,'g*')
% title('Max Thrust (minimum length) Nozzle Design')
% xlabel('Nozzle length (m)')
% ylabel('Nozzle height (m)')
% legend('Nozzle shape','Area_e_x_i_t(predicted)','Char. Lines')
% 
% figure
% plot(noz(:,1),Mnoz,'r*')
% hold on;
% plot(noz(:,1),Pnoz/P_amb,'b*')
% hold on;
% plot(noz(size(noz,1),1),M_e,'go')
% hold on;
% plot(noz(size(noz,1),1),1,'go')
% xlabel('Nozzle length (m)')
% ylabel('Mach number and P/P_c')
% legend('Mach Number','P/P_a_m_b','M_e_x_i_t(predicted)','P_a_m_b/P_a_m_b')

% figure(2)
% plot(Ae,TT)
% title('Thrust curve')
% xlabel('Exit Area (m^2)')
% ylabel('Thrust (N)')

%   Plot the nozzle shape
% figure(1);clf;
% subplot(2,1,1);
% plot(noz(:,1),noz(:,2),'k','LineWidth',3)
% hold on;
% [a,b] = max(noz);
% plot(a(1),A_max/width/2,'g*')
% 
% %   Plot for loop for char lines
% for i = 1 : num
%     figure(1)
%     hold on;
%     plot(char(i,:,1),char(i,:,2))
%     hold on;
%     plot(char(i,:,1),-char(i,:,2))
% end
% hold on;
% plot(A_max,Max_thrust,'r*')
% legend('Thrust Curve','Max Thrust')














