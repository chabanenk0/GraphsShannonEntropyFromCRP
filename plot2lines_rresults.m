function plot2lines_rresults(x_init,y_init,x_mera,y_mera,seriesname,meraname)
figure;
[h,h1,h2]=plotyy(x_init,y_init,x_mera,y_mera);
set(h(1),'FontSize',14);
set(h(2),'FontSize',14);
set(h1,'LineWidth',2);
set(h2,'LineWidth',2);
set(h1,'Color','k');
set(h2,'Color','r');
grid on;

legend(seriesname,meraname);
xlabel('time,days','FontSize',14);
ylabel(meraname,'FontSize',14);

