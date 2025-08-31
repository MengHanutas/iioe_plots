%% scrippt for plot figure 2: profiles of two voyages profiles:
%% load data
load DM63deep_station_position.mat
load IN19deep_station_position.mat
load d63.mat
load ctd19.mat
load v19.mat

%%
% Create a new figure
figure('Units', 'normalized', 'OuterPosition', [0.08 0.08 0.4 0.9]);

% Define the position and size of each subplot
subplot_positions = {
    [0.1, 0.89, 0.36, 0.08], [0.53, 0.89, 0.36, 0.08];
    [0.1, 0.81, 0.36, 0.08], [0.53, 0.81, 0.36, 0.08];% Row 1
    [0.1, 0.70, 0.36, 0.08], [0.53, 0.70, 0.36, 0.08];
    [0.1, 0.62, 0.36, 0.08], [0.53, 0.62, 0.36, 0.08];% Row 2
    [0.1, 0.51, 0.36, 0.08], [0.53, 0.51, 0.36, 0.08];
    [0.1, 0.43, 0.36, 0.08], [0.53, 0.43, 0.36, 0.08]% Row 3
    [0.1, 0.32, 0.36, 0.08], [0.53, 0.32, 0.36, 0.08];
    [0.1, 0.24, 0.36, 0.08], [0.53, 0.24, 0.36, 0.08]% Row 4
    [0.1, 0.13, 0.36, 0.08], [0.53, 0.13, 0.36, 0.08];
    [0.1, 0.05, 0.36, 0.08], [0.53, 0.05, 0.36, 0.08]% Row 5
};
%%
subplot( 'Position', subplot_positions{2, 1});
pcolor(d63.lat,d63.upres,d63.SA);
shading flat
cmocean('haline',15)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
  caxis([34.5 36]) 
  set(gca,'xtickLabel',[])

ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');

[c,h]=contour(d63.lat,d63.upres,d63.SA,[35.8,34.7],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w','Labelspacing',300);
for i=1:length(latitudes)
    lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
    hold on
    plot(lat_point,pressure_new{i},'.','Markersize',4,'color',[0.6 0.6 0.6])
end
ylim([1000  4500])
set(gca, 'TickDir','out', 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])

subplot('Position', subplot_positions{1, 1});
pcolor(d63.lat,d63.upres,d63.SA);
shading flat
cmocean('haline',15)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
  caxis([34.5 36]) 
  set(gca,'xtick',[])
title('(a) Absolute Salinity', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.SA,[35.8,34.7],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w','Labelspacing',300);

for i=1:length(latitudes)
    lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
    hold on
    plot(lat_point,pressure_new{i},'.','Markersize',4,'color',[0.6 0.6 0.6])
end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
set(gca, 'XTick', []);
set(gca, 'TickDir','out')
xlim([-39.5 -11.5])
%%

subplot( 'Position', subplot_positions{1, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.SA);
shading flat
cmocean('haline',15)
title('(f) Absolute Salinity', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
caxis([34.5 36]) 
hold on
% [c,h]=contour(d63.lat,d63.upres,ctd19.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.SA,[35.8,34.7],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w','Labelspacing',400);

% [c,h]=contour(d63.lat,d63.upres,d63.SA,[35.8,34.7],'linest', '--', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');

set(gca,'xtick',[])
set(gca, 'TickDir','out')
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
    
    for i=1:length(latitudes_19)
    lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
    hold on
    plot(lat_point,pressure19_new{i},'.','Markersize',4,'color',[0.6 0.6 0.6])
    end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
xlim([-39.5 -11.5])

subplot('Position', subplot_positions{2, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.SA);
shading flat
cmocean('haline',15)
caxis([34.5 36]) 
hold on
% [c,h]=contour(d63.lat,d63.upres,ctd19.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.SA,[35.8,34.7],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(ctd19.lat,ctd19.upres,ctd19.SA,[35.8,34.7],'linest', '-', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w','Labelspacing',400);
% [c,h]=contour(d63.lat,d63.upres,d63.SA,[35.8,34.7],'linest', '--', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');

%set(gca,'xtick',[])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
    c=colorbar;
    c.Position=[0.91 0.81 0.0133 0.16];
    c.FontSize=9;
    c.Label.String='[g/kg]'
    for i=1:length(latitudes_19)
    lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
    hold on
    plot(lat_point,pressure19_new{i},'.','Markersize',4,'color',[0.6 0.6 0.6])
    end
ylim([1000  4500])
set(gca, 'TickDir','out')
set(gca, 'YTick', [1000:500:4500]);
set(gca, 'XTickLabel', []);
xlim([-39.5 -11.5])

%%


subplot(6, 2, 7);
set(gca, 'Position', subplot_positions{4,1});
pcolor(d63.lat,d63.upres,d63.ct);
shading flat
cmocean('thermal',25)
caxis([0 25]) 
hold on
[c,h]=contour(d63.lat,d63.upres,d63.ct,[0:1:5,5:5:25],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
ylabel('Pressure')
set(gca,'xtickLabel',[])
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'TickDir','out')
set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])

subplot( 'Position', subplot_positions{3, 1});
pcolor(d63.lat,d63.upres,d63.ct);
shading flat
cmocean('thermal',25)
title('(b) Conservative Temperature', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
caxis([0 25]) 
hold on
[c,h]=contour(d63.lat,d63.upres,d63.ct,[0:1:5,5:5:25],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% ylabel('Pressure')
set(gca,'xtick',[])
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'TickDir','out')
set(gca, 'YTick', [0:200:800]);
set(gca, 'TickDir','out')
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
%%

subplot( 'Position', subplot_positions{3, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.ct);
shading flat
cmocean('thermal',25)
% colorbar
caxis([0 25])
title('(g) Conservative Temperature', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');

hold on
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.ct,[0:1:5,5:5:25],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.ct,[0:1:5,5:5:25],'linest', ':', 'color', 'w', 'linewi', 1.5);
% clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
    set(gca, 'TickDir','out')
set(gca,'xtick',[])
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])


subplot('Position', subplot_positions{4, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.ct);
shading flat
cmocean('thermal',25)
c=colorbar;
c.Position=[0.91 0.62 0.0133 0.16]
c.Label.String=['[^{\circ}C]']
caxis([0 25])
hold on
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.ct,[0:1:5,5:5:25],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.ct,[0:1:5,5:5:25],'linest', ':', 'color', 'w', 'linewi', 1.5);
% clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
set(gca,'xtickLabel',[])
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
set(gca, 'TickDir','out')
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
%%

%%
subplot('Position', subplot_positions{6, 1});
pcolor(d63.lat,d63.upres,d63.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
% title('(e) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')
set(gca, 'XTickLabel', []);

subplot( 'Position', subplot_positions{5, 1});
pcolor(d63.lat,d63.upres,d63.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
title('(c) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
 set(gca, 'TickDir','out')
xlim([-39.5 -11.5])

%%
subplot('Position', subplot_positions{6, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
% title('(f) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% colorbar
hold on
% [c,h]=contour(d63.lat,d63.upres,ctd19.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.oxy,[120,245],'linest', '--', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
 
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')
set(gca, 'XTickLabel', []);

subplot( 'Position', subplot_positions{5, 2});
pcolor(ctd19.lat,ctd19.upres,ctd19.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
title('(h) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
c=colorbar;
    c.Position=[0.91 0.43 0.0133 0.16];
    c.Label.String='[\mumol/L]' ;
hold on
% [c,h]=contour(d63.lat,d63.upres,ctd19.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(ctd19.lat,ctd19.upres,ctd19.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.oxy,[120,245],'linest', '--', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')

%%
subplot('Position', subplot_positions{6, 1});
pcolor(d63.lat,d63.upres,d63.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
% title('(e) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')
set(gca, 'XTickLabel', []);

subplot( 'Position', subplot_positions{5, 1});
pcolor(d63.lat,d63.upres,d63.oxy);
shading flat
cmocean('oxy',15)
caxis([100 250])
title('(e) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.oxy,[100,245],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
 set(gca, 'TickDir','out')
xlim([-39.5 -11.5])

%%
subplot('Position', subplot_positions{8, 1});
pcolor(d63.lat,d63.upres,d63.no3);
shading flat
cmocean('delta',15)
caxis([0 40])
% title('(e) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
ylabel('Pressure')
hold on
[c,h]=contour(d63.lat,d63.upres,d63.no3,[10,30,35],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')
set(gca, 'XTickLabel', []);

subplot( 'Position', subplot_positions{7, 1});
pcolor(d63.lat,d63.upres,d63.no3);
shading flat
cmocean('delta',15)
caxis([0 40])
title('(d) Nitrate', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% ylabel('Pressure')
hold on
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
[c,h]=contour(d63.lat,d63.upres,d63.no3,[10,30,35],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
 set(gca, 'TickDir','out')
xlim([-39.5 -11.5])

%%
subplot('Position', subplot_positions{8, 2});
pcolor(ctd19.lat,ctd19.upres,v19.no3);
shading flat
cmocean('delta',15)
caxis([0 40])
% title('(f) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% colorbar
hold on
[c,h]=contour(ctd19.lat,ctd19.upres,v19.no3,[10,30,35],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
 set(gca, 'XTickLabel', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')

subplot( 'Position', subplot_positions{7, 2});
pcolor(ctd19.lat,ctd19.upres,v19.no3);
shading flat
cmocean('delta',15)
caxis([0 40])
title('(i) Nitrate', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
c=colorbar;
    c.Position=[0.91 0.24 0.0133 0.16];
    c.Label.String='[\mumol/L]' ;
hold on
[c,h]=contour(ctd19.lat,ctd19.upres,v19.no3,[10,30,35],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');

set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')

%%
subplot('Position', subplot_positions{10, 1});
pcolor(d63.lat,d63.upres,d63.phos);
shading flat
cmocean('delta',15)
caxis([0 3])
% title('(e) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
ylabel('Pressure')
hold on
[c,h]=contour(d63.lat,d63.upres,d63.phos,[1,2,2.2,2.4,2.6],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
 xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')

subplot( 'Position', subplot_positions{9, 1});
pcolor(d63.lat,d63.upres,d63.phos);
shading flat
cmocean('delta',15)
caxis([0 3])
title('(i) Phosphate', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% ylabel('Pressure')
hold on
[c,h]=contour(d63.lat,d63.upres,d63.phos,[1,2,2.2,2.4,2.6],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
% [c,h]=contour(d63.lat,d63.upres,d63.gm_n,[21:1:25,26.3,26.8,27:0.2:28.2],'linest', '-', 'color', 'k', 'linewi', 1);
% clabel(c,h,'color','k');
% [c,h]=contour(d63.lat,d63.upres,d63.oxy,[120,245],'linest', '-', 'color', 'w', 'linewi', 2);
% clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes)
%     lat_point=repmat(latitudes(i),1,length(pressure_new{i}));
%     hold on
%     plot(lat_point,pressure_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
set(gca, 'XTick', []);
 set(gca, 'TickDir','out')
xlim([-39.5 -11.5])

%%
subplot('Position', subplot_positions{10, 2});
pcolor(ctd19.lat,ctd19.upres,v19.phos);
shading flat
cmocean('delta',15)
caxis([0 3])
% title('(f) Oxygen', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
% colorbar
hold on
[c,h]=contour(ctd19.lat,ctd19.upres,v19.phos,[1,2,2.2,2.4,2.6],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
 
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')

subplot( 'Position', subplot_positions{9, 2});
pcolor(ctd19.lat,ctd19.upres,v19.phos);
shading flat
cmocean('delta',15)
caxis([0 3])
title('(j) Phosphate', 'Units', 'normalized', 'Position',  [0, 1.00, 0], 'HorizontalAlignment', 'left');
c=colorbar;
    c.Position=[0.91 0.05 0.0133 0.16];
    c.Label.String='[\mumol/L]' ;
hold on
[c,h]=contour(ctd19.lat,ctd19.upres,v19.phos,[1,2,2.2,2.4,2.6],'linest', '-', 'color', 'w', 'linewi', 1);
clabel(c,h,'color','w');
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
    set(gca, 'FontSize', 8);
    set(gca, 'LineWidth', 1);
% xlabel('Latitude')
% for i=1:length(latitudes_19)
%     lat_point=repmat(latitudes_19(i),1,length(pressure19_new{i}));
%     hold on
%     plot(lat_point,pressure19_new{i},'.','Markersize',10,'color',[0.6 0.6 0.6])
% end
ylim([0  1000])
set(gca, 'YTick', [0:200:800]);
 set(gca, 'XTick', []);
xlim([-39.5 -11.5])
set(gca, 'TickDir','out')
