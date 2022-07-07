%%
clear; close all; clc;

load('FVM_Code/olson_soln.mat')
u = u(1:end-1,2:end);
v = v(1:end-1,2:end);
p = p(1:end-1,2:end);
rho = rho(1:end-1,2:end);

% Assume solution is symmetric
T = p./(R*rho);
ss = sqrt(abs(gamma*R*T));
Mach = sqrt(u.^2 + v.^2)./ss;
u_exact = [fliplr(u) u(:,2:end)];
v_exact = [fliplr(v) v(:,2:end)];
P_exact = [fliplr(p) p(:,2:end)];
T_exact = [fliplr(T) T];
M_exact = [fliplr(Mach) Mach(:,2:end)];

load('my_soln.mat')
M = sqrt(u.^2+v.^2)./sqrt(1.4*287*T);

P_data = {P/101300, P_exact/101300, abs(P_exact-P)/101300};
M_data = {M, M_exact, abs(M_exact-M)};
pltfmt = {'interpreter','latex','FontSize'};
titles = ["My Soln","Olson's Soln","Error"];
colors = {'#0072BD';
          '#D95319';
          '#EDB120'};

initial = figure('Position',[10 10 1600 700]);
pltfmt = {'interpreter','latex',...
    'fontsize',16};
% density
subplot(2,3,1)
temp_plot = pcolor(xx,yy,rho);
set(temp_plot, 'EdgeColor', 'none');
cb = colorbar; ylabel(cb,'\rho [kg/m^3]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('$\rho [kg/m^3]$',pltfmt{:})
axis equal tight
shading interp
% u-velocity
subplot(2,3,2)
temp_plot = pcolor(xx,yy,u);
set(temp_plot, 'EdgeColor', 'none');
cb = colorbar; ylabel(cb,'u [m/s]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('u [m/s]',pltfmt{:})
axis equal tight
shading interp
% v-velocity
subplot(2,3,3)
temp_plot = pcolor(xx,yy,v);
set(temp_plot, 'EdgeColor', 'none');
cb = colorbar; ylabel(cb,'v [m/s]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('v [m/s]',pltfmt{:})
axis equal tight
shading interp
% Energy
subplot(2,3,4)
temp_plot = pcolor(xx,yy,Cv*T);
set(temp_plot, 'EdgeColor', 'none');
cb = colorbar; ylabel(cb,'e [J/kg]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('e [J/kg]',pltfmt{:})
axis equal tight
shading interp
% Pressure
subplot(2,3,5)
temp_plot = pcolor(xx,yy,P);
set(temp_plot, 'EdgeColor', 'none');
cb = colorbar; ylabel(cb,'P [N/m^2]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('$P [N/m^2]$',pltfmt{:})
axis equal tight
shading interp
% Temperature
subplot(2,3,6)
temp_plot = pcolor(xx,yy,T);
set(temp_plot, 'EdgeColor', 'none');
colormap('turbo')
cb = colorbar; ylabel(cb,'T [K]');
ylabel('y',pltfmt{:}); xlabel('x',pltfmt{:});
title('T [K]',pltfmt{:})
axis equal tight
shading interp
colormap jet
exportgraphics(initial,'initial_cond.pdf')

pltfmt = {'interpreter','latex','FontSize'};
centerline = figure('Position',[10 10 1200 700]);
subplot(2,1,1)
ind = round(size(u,2)/2);
plot(xx(:,ind),P(:,ind)/101300,'-','linewidth',2.5,'color',colors{1});
hold on
plot(xx(:,ind),P_exact(:,ind)/101300,'--','linewidth',2.5,'color',colors{2});
plot(xx(:,1),P_ex/101300,':','linewidth',2.5,'color',colors{3})
xlabel('x',pltfmt{:},14);
ylabel('$P/P_{atm}$',pltfmt{:},14);
title('Centerline Pressure',pltfmt{:},14)
legend('My Soln',"Olson's Soln",'Theory','location','northeastoutside',pltfmt{:},14)
% axis([0 max(xx(:,ind)) 0 max([max(P_ex/101300),max(M_ex)])])
% yticks(0:round(max(P_ex/101300)))
grid on
axis tight
subplot(2,1,2)
ind = round(size(u,2)/2);
plot(xx(:,ind),M(:,ind),'-','linewidth',2.5,'color',colors{1});
hold on
plot(xx(:,ind),M_exact(:,ind),'--','linewidth',2.5,'color',colors{2});
plot(xx(:,1),M_ex,':','linewidth',2.5,'color',colors{3})
xlabel('x',pltfmt{:},14);
ylabel('Mach',pltfmt{:},14);
title('Centerline Mach Number',pltfmt{:},14)
legend('My Soln',"Olson's Soln",'Theory','location','northeastoutside',pltfmt{:},14)
% axis([0 max(xx(:,ind)) 0 max([max(P_ex/101300),max(M_ex)])])
% yticks(0:round(max(M_ex)))
grid on
axis tight
exportgraphics(centerline,'Centerline_validation.pdf')

valid = figure('Position',[10 10 2400 800]);
for ii = 1:3
    subplot(2,4,ii)
    pcolor(xx,yy,P_data{ii})
    cb = colorbar;
    ylabel(cb,'$P/P_{atm}$',pltfmt{:},12)
    xlabel('x',pltfmt{:},14); ylabel('y',pltfmt{:},14)
    title(titles{ii},pltfmt{:},14);
    shading interp
end
subplot(2,4,4)
plot(P_data{1}(end,:),yy(end,:),'Linewidth',2)
hold on
plot(P_data{2}(end,:),yy(end,:),'Linewidth',2)
xline(P_ex(end)/101300,'--','Linewidth',2)
ylabel('y',pltfmt{:},12)
xlabel('$P/P_{atm}$',pltfmt{:},14); ylabel('y',pltfmt{:},14)
title("Exit Pressure",pltfmt{:},14);
legend('My soln','Olson soln','expected',pltfmt{:},14,'Location','bestoutside')
grid on

for ii = 1:3
    subplot(2,4,ii+4)
    pcolor(xx,yy,M_data{ii})
    cb = colorbar;
    ylabel(cb,'Mach',pltfmt{:},12)
    xlabel('x',pltfmt{:},14); ylabel('y',pltfmt{:},14)
    title(titles{ii},pltfmt{:},14);
    shading interp
end

subplot(2,4,8)
plot(M_data{1}(end,:),yy(end,:),'Linewidth',2)
hold on
plot(M_data{2}(end,:),yy(end,:),'Linewidth',2)
xline(M_ex(end),'--','Linewidth',2)

ylabel('y',pltfmt{:},12)
xlabel('Mach',pltfmt{:},14); ylabel('y',pltfmt{:},14)
title("Exit Mach",pltfmt{:},14);
legend('My soln','Olson soln','expected',pltfmt{:},14,'Location','bestoutside')
grid on

exportgraphics(valid,'Grid_validation.pdf')