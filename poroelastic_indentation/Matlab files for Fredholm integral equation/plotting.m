%% plotting

figure 
hold on
set(gcf,'PaperPositionMode','auto');
set(gca,'fontsize',10,'linewidth',0.5);
box on;

plot(x_0_inf,x12_M,'Color',[0 0 0],'LineStyle','-','LineWidth',1); hold on
plot(x_0_inf,x12_a_storage(10,:),'Color',[0 0 0],'LineStyle','--','LineWidth',1); hold on
plot(x_0_inf,x12_theta1,'Color',[0 0 0],'LineStyle','-.','LineWidth',1); hold on

plot([10^1.4,10^1.4],[10^(-4),10^(-4.7)],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5); hold on
plot([10^1.4,10^1.75],[10^(-4.7),10^(-4.7)],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5); hold on
plot([10^1.4,10^1.75],[10^(-4),10^(-4.7)],'Color',[0 0 0],'LineStyle','-','LineWidth',0.5); hold on

annotation(figure(1),'textbox',...
    [0.647071428571428 0.298571428571429 0.0647142857142857 0.0552380952380952],...
    'String',{'$0.7$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');

annotation(figure(1),'textbox',...
    [0.689214285714285 0.234761904761906 0.0697142857142866 0.0552380952380954],...
    'String','$0.35$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FitBoxToText','off');

xlabel('$x_{*}$','interpreter','latex')

legend({'$\sqrt{x_{*}}M\left(s_{*},x_{*}\right)$','$\sqrt{x_{*}}a_{10}\left(s_{*},x_{*}\right)$','$\sqrt{x_{*}}\theta_{1}\left(s_{*},x_{*}\right)$'},'Interpreter','latex','fontsize',10,'location','southwest');

grid on
set(gca,'Xscale','log')
set(gca,'Yscale','log')
set(gca,'YMinorTick','off')
set(gca,'YMinorGrid','off')
set(gca,'XMinorTick','off')
set(gca,'XMinorGrid','off')

set(gcf, 'PaperPosition', [0 0 6 4.5])
f1 = figure(1); 
f1.PaperSize = [6 4.5]; 
ylim([10^(-6) 1])
