%% code for plot figur 1 in manuscript using CARS data and WOD23:
%load data:
load('Figure1.mat') %load the data 
% Define region of interest
lon_min = 90; lon_max = 160;
lat_min = -50; lat_max = -8;

%% plot
figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.8, 0.85]);
%subplot(a) DH surface
subplot('Position', [0.06, 0.52, 0.3, 0.42]);
m_proj('Mercator', 'longitudes', [lon_min lon_max], 'latitudes', [lat_min lat_max]);
m_pcolor(lon_subset, lat_subset, squeeze(DH(:,:,1))');
shading flat;
colormap(parula(12))
hold on
[CS, CH]=m_contour(lon_subset, lat_subset, movmean(movmean(squeeze(DH(:,:,1)),20,1),20,2)', [14:0.5:26],'LineColor', 'k');
clabel(CS, CH, [14:1:20,20:0.5:26],'FontSize', 12, 'Color', 'k'); % Label the contour lines
clear CS CH
% hcb = colorbar('Position', [0.37 0.08 0.015 0.86]); % Adjust position to cover both maps
hcb = colorbar('Position', [0.37 0.52 0.015 0.42]);
hcb.Label.String = 'Dynamic Height';
hcb.FontSize=12;
title('(a) Suface Dynamic Height referenced to 2000m','Units', 'normalized', 'Position', [0, 1, 0], 'HorizontalAlignment', 'left');
% xlabel('Longitude');
ylabel('Latitude','FontSize',12);
xticks([])
m_coast('patch', [0.7 0.7 0.7]);
m_grid('linestyle', 'none', 'box', 'fancy', 'tickdir', 'out', 'xticklabels', []);
 for  i=1:length(latitudes_dm)
    
    m_plot(longitudes_dm(i), latitudes_dm(i),'ko','markersize',6,'linewidth',.5,'markerfacecolor','w');
%         m_text(longitudes(i), latitudes(i),[int2str(stations(i))],...
%         'fontweight','bold','hor','cent','fontsize',4,'color','k')
 end
 %% subplot(b) DH 400m
 subplot('Position', [0.06, 0.07, 0.3, 0.42]);
m_proj('Mercator', 'longitudes', [lon_min lon_max], 'latitudes', [lat_min lat_max]);
m_pcolor(lon_subset, lat_subset, squeeze(DH(:,:,2))');
shading flat;
colormap(parula(12))
hold on
[CS, CH]=m_contour(lon_subset, lat_subset, movmean(movmean(squeeze(DH(:,:,2)),20,1),20,2)',[9:0.5:16], 'LineColor', 'k');
clabel(CS, CH,[12:0.5:16], 'FontSize', 12, 'Color', 'k'); % Label the contour lines
% colorbar;
title('(b) 400m Dynamic Height referenced to 2000m','Units', 'normalized', 'Position', [0, 1, 0], 'HorizontalAlignment', 'left');
xlabel('Longitude','Position',[0.02,-1.08],'FontSize',12);
ylabel('Latitude','FontSize',12);
m_coast('patch', [0.7 0.7 0.7]);
m_grid('linestyle', 'none','box', 'fancy', 'tickdir', 'out');
hcb = colorbar('Position', [0.37 0.07 0.015 0.42]);
hcb.Label.String = 'Dynamic Height';
hcb.FontSize=12;
caxis([9 15])
for l = 1:length(latitudes_in)
    m_plot(longitudes_in(l), latitudes_in(l), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
%     m_text(mat_longitudes(l), mat_latitudes(l),[int2str(mat_stations(l))],...
%         'fontweight','bold','hor','cent','fontsize',4,'color','k')
end

%% draw right hand side profile figures
%SA
subplot('Position', [0.48, 0.805, 0.4, 0.125]);
pcolor(repmat(lat_subset,1,79),Pressure(:,lat_indices)',SA_110)
shading flat
hold on
[CS, CH]=contour(repmat(lat_subset,1,79),Pressure(:,lat_indices)',SA_110,[34.7,34.8,35.8] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
% c=colorbar;
% c.Label.String='[g/kg]';
% c.FontSize=12;
cmocean('haline',16)
% ylabel('Pressure','FontSize',12)
xticks([])
caxis([34.5 36.1])
title('(c) CARS Absolute Salinity','Units', 'normalized', 'Position', [0, 1, 0], 'HorizontalAlignment', 'left','FontSize',12);
ylim([0 1000])
set(gca, 'YTick', [0:200:1000],'FontSize',12);
xlim([-40 -8])
%%
subplot('Position', [0.48, 0.68, 0.4, 0.125]);
pcolor(repmat(lat_subset,1,79),Pressure(:,lat_indices)',SA_110)
shading flat
hold on
[CS, CH]=contour(repmat(lat_subset,1,79),Pressure(:,lat_indices)',SA_110,[34.8,35.8] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
c=colorbar;
c.Label.String='[g/kg]';
c.FontSize=12;
c.Position=[0.89 0.68 0.016 0.25]
cmocean('haline',16)
ylabel('Pressure','FontSize',12)
xticks([])
caxis([34.5 36.1])
ylim([1000 4500])
xlim([-40 -8])
set(gca, 'YTick', [2000:1000:4000,4500],'FontSize',12);
%% PV
subplot('Position', [0.48, 0.495, 0.4, 0.125]);
pcolor(repmat(lat_subset,1,79),Pressure(:,lat_indices)',PV_formula)
shading flat
hold on
[CS, CH]=contour(repmat(lat_subset,1,79),Pressure(:,lat_indices)',PV_formula,[0.6e-10 0.6e-10] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
% c=colorbar;
% c.Label.String='[m^{-1}s^{-1}]';
% c.FontSize=12;
% c.Position=[0.89 0.37 0.016 0.25]
caxis([0 1e-10])
% ylabel('Pressure','FontSize',12)
xticks([])
title('(d) CARS PV','Units', 'normalized', 'Position', [0, 1, 0], 'HorizontalAlignment', 'left');
ylim([0 1000])
set(gca, 'YTick', [0:200:1000],'FontSize',12);
xlim([-40 -8])
cmocean('thermal',10)
%%
subplot('Position', [0.48, 0.37, 0.4, 0.125]);
pcolor(repmat(lat_subset,1,79),Pressure(:,lat_indices)',PV_formula)
shading flat
hold on
[CS, CH]=contour(repmat(lat_subset,1,79),Pressure(:,lat_indices)',PV_formula,[0.6e-10 0.6e-10] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
c=colorbar;
c.Label.String='[m^{-1}s^{-1}]';
c.FontSize=12;
c.Position=[0.89 0.37 0.016 0.25]
caxis([0 1e-10])
ylabel('Pressure','FontSize',12)
xticks([])
ylim([1000 4500])
xlim([-40 -8])
set(gca, 'YTick', [2000:1000:4000,4500],'FontSize',12);
cmocean('thermal',10)

%%
load WOA23_oxy.mat
%remove large value
index=find(WOA.an>1000);
WOA.an(index)=nan;
subplot('Position', [0.48, 0.185, 0.4, 0.125]);
pcolor(WOA.lat,WOA.depth,WOA.an')
shading flat
hold on
[CS, CH]=contour(WOA.lat,WOA.depth,WOA.an',[100,120,200,245] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
% c=colorbar;
% c.Label.String='[\mumol/kg]';
% c.FontSize=12;
caxis([100 250])
xlim([-40 -8])
cmocean('oxy',15)
xticks([])
ylim([0 1000])
set(gca, 'YTick', [0:200:1000],'FontSize',12);
title('(e) WOA23 Oxygen','Units', 'normalized', 'Position', [0, 1, 0], 'HorizontalAlignment', 'left');
%%
subplot('Position', [0.48, 0.06, 0.4, 0.125]);
pcolor(WOA.lat,WOA.depth,WOA.an')
shading flat
hold on
[CS, CH]=contour(WOA.lat,WOA.depth,WOA.an',[100,120,200,245] ,'LineColor', 'k');
clabel(CS, CH, 'FontSize', 12, 'Color', 'k'); % Label the contour lines
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
c=colorbar;
c.Label.String='[\mumol/kg]';
c.FontSize=12;
c.Position=[0.89 0.06 0.016 0.25]
caxis([100 250])
xlim([-40 -8])
cmocean('oxy',15)
xlabel('Latitude','FontSize',12)
ylabel('Pressure','FontSize',12)
ylim([1000 4500])
set(gca, 'YTick', [2000:1000:4000,4500],'FontSize',12);