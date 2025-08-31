%% prepare data
% load voyage_20dbar.mat
% v19=ctdHy19v03Slat;
% d63=ctdHy63v02D500lat_new;
% ctd19=ctd19v03Slat;
%% This is script for FIgure8
load ctd19.mat
load d63.mat
load v19.mat
load bio_rms.mat
neutralgrid =[21.2:0.02:28.2];

lats = ctd19.lat;
nut = neutralgrid;
pres = repmat(d63.upres, 1, length(d63.lat)); % since both data sets were interpolated on the same pressure levels at 5 dbar interval
p63iso = ra_projvar(pres, d63.gm_n, neutralgrid);
p19iso = ra_projvar(pres, ctd19.gm_n, neutralgrid);

% Creating Pm at mid point of two voyages.

presm = (p63iso + p19iso)/2; %351*286)
lat = repmat(lats, length(presm), 1);
clear p19iso p63iso 

%
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
sa_mean=(sa19+sa63)./2;

ct19 = ra_projvar(ctd19.ct, ctd19.gm_n, neutralgrid);
ct63 = ra_projvar(d63.ct, d63.gm_n, neutralgrid);

oxy19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
oxy63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
oxy_mean=(oxy19+oxy63)./2;
oxy_change_rho = oxy19 - oxy63;
oxy_rms_gamma=ra_projvar(oxy_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);
%oxy_rms_lat is calculated from test_bio_sign.m

aou19 = ra_projvar(ctd19.aou, ctd19.gm_n, neutralgrid);
aou63 = ra_projvar(d63.aou, d63.gm_n, neutralgrid);
aou_mean=(aou19+aou63)./2;
aou_change_rho = aou19 - aou63;
aou_rms_gamma=ra_projvar(aou_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);

no319 = ra_projvar(v19.no3, ctd19.gm_n, neutralgrid);
no363 = ra_projvar(d63.no3, d63.gm_n, neutralgrid);
no3_rms_gamma=ra_projvar(no3_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);
no3_mean=(no319+no363)./2;
no3_change_rho = no319 - no363;

phos19 = ra_projvar(v19.phos, ctd19.gm_n, neutralgrid);
phos63 = ra_projvar(d63.phos, d63.gm_n, neutralgrid);
phos_rms_gamma=ra_projvar(phos_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);
phos_mean=(phos19+phos63)./2;
phos_change_rho = phos19 - phos63;
%%
% Create a figure
figure('Units', 'inches', 'Position', [1 1 12 8], 'PaperPositionMode', 'auto');

% Define subplot positions to minimize blank areas
subplotPositions = [
    0.08, 0.75, 0.35, 0.20; 
    0.08, 0.55, 0.35, 0.20; % [left, bottom, width, height] for 1st row, 1st column
    0.56, 0.75, 0.35, 0.20; 
     0.56, 0.55, 0.35, 0.20; % 1st row, 2nd column
    0.08, 0.285, 0.35, 0.20; 
    0.08, 0.085, 0.35, 0.20; % 2nd row, 1st column
    0.56, 0.285, 0.35, 0.20;
      0.56, 0.085, 0.35, 0.20;% 2nd row, 2nd column
    
    ];
%% subplot(2) DO change on neutral density
subplot('Position', subplotPositions(1, :));
% pcolor(lats,nut,oxy_change_rho)
pcolor(lats,presm,oxy_change_rho)

shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
mask=oxy_rms_gamma>=abs(oxy_change_rho);
stipple(lat,presm,mask)
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-40 40])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
title('(a) Oxygen','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(2, :));
pcolor(lats,presm,oxy_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h

stipple(lat,presm,mask)
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);xlim([-32 -11.5])
caxis([-40 40])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
%title('(c)','FontSize', 12)
ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
c=colorbar;
c.Position=[0.44 0.55 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
%% subplot(4) AOU change on neutral density
subplot('Position', subplotPositions(5, :));
pcolor(lats,presm,aou_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
mask=aou_rms_gamma>=abs(aou_change_rho);
stipple(lat,presm,mask)
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-40 40])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
title('(c) AOU','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(6, :));
pcolor(lats,presm,aou_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h

stipple(lat,presm,mask)
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);xlim([-32 -11.5])
c=colorbar;
c.Position=[0.44 0.1000 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
% set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-40 40])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
ylabel('Pressure', 'FontSize', 12);
xlabel('Latitude', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);


%%
subplot('Position', subplotPositions(3, :));
pcolor(lats,presm,no3_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
mask=no3_rms_gamma>=abs(no3_change_rho);
stipple(lat,presm,mask)
ylim([0 1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-10 10])
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
title('(b) Nitrate','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(4, :));
pcolor(lats,presm,no3_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
stipple(lat,presm,mask)
ylim([1000 4500])
set(gca, 'YTick', [1000:500:4500]);
%set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-10 10])
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
%ylabel('Neutral Dneisty', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
c=colorbar;
c.Position=[0.92 0.55 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
%% subplot(4) AOU change on neutral density
subplot('Position', subplotPositions(7, :));
pcolor(lats,presm,phos_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
mask=phos_rms_gamma>=abs(phos_change_rho);
stipple(lat,presm,mask)
ylim([0 1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-1 1])
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
title('(d) Phosphate','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(8, :));
pcolor(lats,presm,phos_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(lat(:,1:196),presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linest','--');
clabel(c, h, 'fontsize', 12, 'color', 'k','Labelspacing',400); clear c h
stipple(lat,presm,mask)
ylim([1000 4500])
set(gca, 'YTick', [1000:500:4500]);
c=colorbar;
c.Position=[0.92 0.1000 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
% set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-1 1])
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
%ylabel('Neutral Density', 'FontSize', 12);
xlabel('Latitude', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);


%%
%%
cmap = cmocean('balance', 22);
cmap=cmap([1:10,13:22],:);
colormap(cmap); % Apply the adjusted colormap
