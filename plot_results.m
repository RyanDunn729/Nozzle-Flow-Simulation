%%
clear; close all; clc;

files = {'data_075.mat';
         'data_10.mat';
         'data_17.mat'};

colors = {'#0072BD';
          '#D95319';
          '#EDB120';
          '#7E2F8E';
          '#77AC30'};


load('data_075.mat')

N = 30;
x_start = zeros(N,1);
y_start = linspace(-5e-6,5e-6,N);

streamline_fig = figure('Position', [10 10 1200 800]);
title('Streamlines','interpreter','latex','FontSize',14)
hold on

for blk=[2 4 5 6]
    u_plot = pcolor(xx{blk},yy{blk},sqrt(u{blk}.^2+v{blk}.^2)./sqrt(1.4*287*T{blk}));
    set(u_plot, 'EdgeColor', 'none');
    x_start = xx{blk}(1,1)*ones(N,1);
    y_start = linspace(yy{blk}(1,1),yy{blk}(1,end),N);
    stream = streamline(transpose(xx{blk}),transpose(yy{blk}),transpose(u{blk}),transpose(v{blk}),x_start,y_start);
    set(stream, 'Color', 'red','LineWidth',1.5);
end
xlabel('x','interpreter','latex','FontSize',14)
ylabel('y','interpreter','latex','FontSize',14)
cb = colorbar;
ylabel(cb,'Mach','interpreter','latex','FontSize',14)
axis equal tight
exportgraphics(streamline_fig,'figures/streamline_fig.pdf')
      
% Plot pressure at exit
Pe_plots = figure('Position',[10 10 600 500]);
hold on
for kk = 1:numel(files)
    load(files{kk})
    plot(P{2}(end-1,:)/101300,yy{2}(end-1,:),'LineWidth',2)
end
xlabel('$P_e / P_{atm}$','interpreter','latex','FontSize',14)
ylabel('y','interpreter','latex','FontSize',14)
title('Exit Presure Distribution','interpreter','latex','FontSize',14)
legend('0.75$P_0$','1.0$P_0$','1.7$P_0$','interpreter','latex','location','northeast','fontsize',14)
grid on
exportgraphics(Pe_plots,'figures/Pe_plots.pdf')

% Expansion / Oblique shock angles
exp_lines_fig = figure('Position',[10 10 1500 400]);
for kk = 1:numel(files)
    load(files{kk})
    x = [xx{4} xx{5} xx{6}];
    y = [yy{6} yy{5} yy{4}];
    M = [sqrt(u{6}.^2+v{6}.^2)./sqrt(1.4*287*T{6}), sqrt(u{5}.^2+v{5}.^2)./sqrt(1.4*287*T{5}),sqrt(u{4}.^2+v{4}.^2)./sqrt(1.4*287*T{4})];
    height = size(x,2);
    len = size(x,1);
    top = zeros(len,1);
    top_ind = zeros(len,1);
    bot = zeros(len,1);
    bot_ind = zeros(len,1);
    for ii = 1:len
        for jj = 1:height
            if M(ii,jj)>0.77*max(max(M))
                bot(ii) = y(ii,jj);
                bot_ind(ii) = jj;
                break
            end
        end
        for jj = height:-1:1
            if M(ii,jj)>0.77*max(max(M))
                top(ii) = abs(y(ii,jj));
                top_ind(ii) = jj;
                break
            end
        end
    end
    subplot(1,numel(files),kk)
    hold on
    plot(x(:,1),top,'r','LineWidth',3)
    for blk = [2 4 5 6]
        pcolor(xx{blk},yy{blk},u{blk}./sqrt(1.4*287*T{blk}))
    end
    cb = colorbar;
    caxis([0 2.5])
    ylabel(cb,'Mach','interpreter','latex','FontSize',12)
    shading interp
    plot(x(:,1),top,'r','LineWidth',3)
    plot(x(:,1),bot,'r','LineWidth',3)
    ylabel('y','interpreter','latex','FontSize',12)
    xlabel('x','interpreter','latex','FontSize',12)
    title(sprintf('$P = %1.2f P_{0}$',P0_adjust),'interpreter','latex','fontsize',14)
    axis tight equal

    ang = atan2d(top-min(top),x(:,1)-x(1,1))/2 + atan2d(-bot+max(bot),x(:,1)-x(1,1))/2;    
    fprintf(sprintf('Average Shock Angle: %2.2f\n',mean(ang(ang>0))))
    Me = sqrt(u{2}(end,:).^2+v{2}(end,:).^2)./sqrt(1.4*287*T{2}(end,:));
    M_inf = sqrt(u{5}(end,:).^2+v{5}(end,:).^2)./sqrt(1.4*287*T{5}(end,:));
    
    M1 = Me(round(ny{2}/2));
    M2 = M_inf(round(ny{2}/2));
    nu1 = sqrt(2.4/0.4)*atand(sqrt(0.4*(M1^2-1)/2.4)) - atand(sqrt(M1^2-1));
    nu2 = sqrt(2.4/0.4)*atand(sqrt(0.4*(M2^2-1)/2.4)) - atand(sqrt(M2^2-1));
    fprintf('via Centerline:\n')
    fprintf(sprintf('Centerline Me: %2.2f\n',M1))
    fprintf(sprintf('Centerline M_inf: %2.2f\n',M2))
    fprintf(sprintf('Expansion Angle: %2.2f\n',nu2-nu1))
    
    M1 = mean(Me);
    M2 = mean(M_inf(M_inf>1));
    nu1 = sqrt(2.4/0.4)*atand(sqrt(0.4*(M1^2-1)/2.4)) - atand(sqrt(M1^2-1));
    nu2 = sqrt(2.4/0.4)*atand(sqrt(0.4*(M2^2-1)/2.4)) - atand(sqrt(M2^2-1));
    fprintf('via Averaging:\n')
    fprintf(sprintf('Average Me: %2.2f\n',M1))
    fprintf(sprintf('Average M_inf: %2.2f\n',M2))
    fprintf(sprintf('Expansion Angle: %2.2f\n',nu2-nu1))
    fprintf('--------\n')
end
exportgraphics(exp_lines_fig,'figures/exp_lines.pdf')

% Centerline Flow Properties 
centerline_fig = figure('Position',[10 10 900 700]);
subplot(2,1,1)
hold on
for kk = 1:numel(files)
    load(files{kk})
    PP = [P{2}(:,round(ny{2}/2));P{5}(:,round(ny{5}/2))]/101300;
    x = [xx{2}(:,1);xx{5}(:,1)];
    plot(x,PP,'Color',colors{kk},'LineWidth',2)
end
xline(xx{2}(end-1,1),'--')
ylabel('$P / P_{atm}$','interpreter','latex','FontSize',12)
xlabel('x','interpreter','latex','FontSize',12)
title('Centerline Pressure','interpreter','latex','fontsize',14)
axis([xx{2}(1,1) xx{5}(end,1) 0 5])
grid on
subplot(2,1,2)
hold on
for kk = 1:numel(files)
    load(files{kk})
    U = [u{2}(:,round(ny{2}/2));u{5}(:,round(ny{5}/2))];
    V = [v{2}(:,round(ny{2}/2));v{5}(:,round(ny{5}/2))];
    T = [T{2}(:,round(ny{2}/2));T{5}(:,round(ny{5}/2))];
    mach = sqrt(U.^2+V.^2)./sqrt(1.4*287*T);
    x = [xx{2}(:,1);xx{5}(:,1)];
    plot(x,mach,'Color',colors{kk},'LineWidth',2)
end
xline(xx{2}(end-1,1),'--')
ylabel('M','interpreter','latex','FontSize',12)
xlabel('x','interpreter','latex','FontSize',12)
title('Centerline Mach Number','interpreter','latex','fontsize',14)
axis tight
grid on

for kk = 1:numel(files)
    load(files{kk})
    A_At = yy{2}(1,end)./yy{2}(:,end);
    func = @(M) M.*(2/(gamma+1)*(1+(gamma-1)/2*M.^2))^(-(gamma+1)/(2*gamma-2));
    M_ex = zeros(nx{2},1);
    P_ex = zeros(nx{2},1);
    for ii = 1:nx{2}
        M_ex(ii) = fzero(@(x) func(x)-A_At(ii),[1 4]);
        P_ex(ii) = P_0*((1+(gamma/2-1/2)*M_ex(ii)*M_ex(ii)))^(-gamma/(gamma-1));
    end
    subplot(2,1,2)
    hold on
    plot(xx{2}(1:end-1,1),M_ex(1:end-1),':','Color',colors{kk},'LineWidth',2)
    subplot(2,1,1)
    hold on
    plot(xx{2}(1:end-1,1),P_ex(1:end-1)/101300,':','Color',colors{kk},'LineWidth',2)
end
subplot(2,1,1)
legend('0.75$P_0$','1.0$P_0$','1.7$P_0$','Exit','interpreter','latex','location','northeastoutside','fontsize',12)
subplot(2,1,2)
legend('0.75$P_0$','1.0$P_0$','1.7$P_0$','Exit','interpreter','latex','location','northeastoutside','fontsize',12)
exportgraphics(centerline_fig,'figures/centerline.pdf')
% surf 3d pressure/mach distribution
% load('soln_data/sub_10.mat')
% subplot(2,1,1)
% for blk = [2 4 5 6]
%     surf(xx{blk},yy{blk},P{blk}/101300)
%     hold on
% end
% subplot(2,1,2)
% for blk = [2 4 5 6]
%     mach = u{blk}./sqrt(1.4*287*T{blk});
%     surf(xx{blk},yy{blk},mach)
%     hold on
% end
% shading interp