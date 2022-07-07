%%
close all; clear all; clc;
pltfmt = {'interpreter','latex','fontsize',20};

nx1 = 16; ny1 = 14;
[xi,eta] = ndgrid(linspace(0,1,nx1),linspace(0,1,ny1));

% L = 1e-5; 
% H = 8e-6;

Ae = 8e-6;
L = 1e-5;
At = Ae/4;

func_x = @(xi,eta) L*xi;
func_y = @(xi,eta) (At/2 + (Ae-At)/4*(1-cospi(xi))).*(1-2*eta);

figure('Position',[10 10 1200 500])
subplot(1,2,1)
plot(xi,eta,'k.','MarkerSize',10)
xlabel('$\xi$',pltfmt{:}); ylabel('$\eta$',pltfmt{:});
title('Sim1: Parametric Coordinates',pltfmt{:});
axis equal tight
grid on
xticks(0:0.2:1)
yticks(0:0.2:1)

subplot(1,2,2)
fplot(@(x) func_x(x,0), @(x) func_y(x,1),[0,1],'k-','LineWidth',2)
hold on
fplot(@(x) func_x(x,0), @(x) func_y(x,0),[0,1],'k-','LineWidth',2)
for i=linspace(0.1,0.9,4)
    fplot(@(x) func_x(x,i), @(x) func_y(x,i),[0,1],'b-','LineWidth',0.5)
    fplot(@(y) func_x(i,y), @(y) func_y(i,y),[0,1],'b-','LineWidth',0.5)
end
plot(func_x(xi,eta),func_y(xi,eta),'r.','MarkerSize',10)
xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
title('Sim1: Spatial Coordinates',pltfmt{:});
axis equal tight
grid on

figure('Position',[10 10 1200 500])
nx2 = 16; ny2 = 28;
[xx,yy] = ndgrid(linspace(L,2*L,nx2),linspace(-2*Ae,2*Ae,ny2));
subplot(1,2,1)
plot(xx,yy,'k.','MarkerSize',10)
xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
title({'Sim2: Spatial Coordinates',''},pltfmt{:});
axis equal tight
grid on

subplot(1,2,2)
plot(func_x(xi(1),eta(1)),func_y(xi(1),eta(1)),'r.','MarkerSize',10)
hold on
plot(xx,yy,'k.','MarkerSize',10)
plot(func_x(xi,eta),func_y(xi,eta),'r.','MarkerSize',10)
fplot(@(x) func_x(x,0), @(x) func_y(x,1),[0,1],'k-','LineWidth',2)
fplot(@(x) func_x(x,0), @(x) func_y(x,0),[0,1],'k-','LineWidth',2)
xlabel('x',pltfmt{:}); ylabel('y',pltfmt{:});
title('Full Simulation',pltfmt{:});
legend('Sim1 nodes','Sim2 nodes','Location','northwest','interpreter','latex','fontsize',14)
axis equal tight
grid on