function tvar_plotirf(y,minmax,linecolor,bandcolor)

set(gcf,'Color','white');
set(gca,'Box','on');
set(gca,'Linewidth',0.7);

cl = y(:,2);
cu = y(:,3);
ym = y(:,1);

ymin = minmax(1);
ymax = minmax(2);

t = size(y,1);
h = (1:t)';

hold on;

x_patch = [h; flipud(h)];
y_patch = [cu; flipud(cl)];
patch(x_patch, y_patch, bandcolor, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.08);

plot(h, ym, 'LineWidth',2, 'Color', linecolor, 'Marker','x','MarkerSize',6);
plot(h, cl, '--', 'LineWidth',1.3, 'Color', bandcolor);
plot(h, cu, '--', 'LineWidth',1.3, 'Color', bandcolor);

grid on;
yline(0, 'k', 'LineWidth', 0.8);
yline(100, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);   % remove if not needed

xlim([1 t]);
ylim([ymin ymax]);

set(gca,'FontSize',10);

end