%% figure for po4* and n*

%% load data
% load voyage_20dbar.mat
% v19=ctdHy19v03Slat;
% d63=ctdHy63v02D500lat_new;
% ctd19=ctd19v03Slat;
load ctd19.mat
load d63.mat
load v19.mat
load bio_rms.mat
%% calculate value on isopycnal and mean value on isopycnals 

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

po4star19 = ra_projvar(v19.po4star, ctd19.gm_n, neutralgrid);
po4star63 = ra_projvar(d63.po4star, d63.gm_n, neutralgrid);
po4star_change_rho = po4star19 -po4star63 ;
po4star_rms_gamma=ra_projvar(po4star_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);

nstar19 = ra_projvar(v19.Nstar, ctd19.gm_n, neutralgrid);
nstar63 = ra_projvar(d63.Nstar, d63.gm_n, neutralgrid);
nstar_change_rho = nstar19 -nstar63 ;
nstar_rms_gamma=ra_projvar(nstar_rms_lat,(ctd19.gm_n+d63.gm_n)./2,neutralgrid);

aou19 = ra_projvar(ctd19.aou, ctd19.gm_n, neutralgrid);
aou63 = ra_projvar(d63.aou, d63.gm_n, neutralgrid);
aou_mean=(aou19+aou63)./2;
aou_change_rho = aou19 - aou63;

cmap = cmocean('balance', 22);
cmap=cmap([1:10,13:22],:);
%% figure plot
%% figure for oxygen and AOU
% Create a figure
figure('Units', 'inches', 'Position', [1 1 12 8], 'PaperPositionMode', 'auto');

% Define subplot positions to minimize blank areas
subplotPositions = [
    0.08, 0.75, 0.25, 0.20; 
    0.08, 0.55, 0.25, 0.20; % [left, bottom, width, height]
    0.34, 0.75, 0.25, 0.20; 
     0.34, 0.55, 0.25, 0.20; % 
     0.67, 0.75, 0.25, 0.20; 
     0.67, 0.55, 0.25, 0.20;
    0.08, 0.3, 0.25, 0.20; 
    0.08, 0.1, 0.25, 0.20; % 
    0.34, 0.3, 0.25, 0.20;
      0.34, 0.1, 0.25, 0.20;
       0.67, 0.3, 0.25, 0.20;
      0.67, 0.1, 0.25, 0.20;% 
    
    ];

% subplot('Position', subplotPositions(1, :));
% subplot('Position', subplotPositions(2, :));
% subplot('Position', subplotPositions(3, :));
% subplot('Position', subplotPositions(4, :));
% subplot('Position', subplotPositions(5, :));
% subplot('Position', subplotPositions(6, :));
% subplot('Position', subplotPositions(7, :));
% subplot('Position', subplotPositions(8, :));
% subplot('Position', subplotPositions(9, :));
% subplot('Position', subplotPositions(10, :));
% subplot('Position', subplotPositions(11, :));
% subplot('Position', subplotPositions(12, :));

%% start subplot

%% subplot for po4star
%19
subplot('Position', subplotPositions(3, :));
pcolor(lat,presm,po4star19)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
 set(gca, 'YTick', [0:200:1000]);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-1.6 1.6])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(b) 2019 PO4*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(4, :));
pcolor(lat,presm,po4star19)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
% set(gca, 'YTick', [0:200:1000]);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-1.6 1.6])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')

set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])
 set(gca, 'YTick', [1000:500:4500]);
% set(gca, 'YTick', []);
% set(gca, 'XTick', []);
xlim([-32 -11.5])

%63
subplot('Position', subplotPositions(1, :));
pcolor(lat,presm,po4star63)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
% set(gca, 'YTick', []);

set(gca, 'XTickLabel', []);

set(gca, 'YTick', [0:200:1000]);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-1.6 1.6])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(a) 1963 PO4*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(2, :));
pcolor(lat,presm,po4star63)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h

ylim([0  1000])
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
% set(gca, 'YTick', []);
% set(gca, 'XTick', []);
xlim([-32 -11.5])
caxis([-1.6 1.6])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
 ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])
set(gca, 'YTick', [1500:500:4500]);

% set(gca, 'XTick', []);
xlim([-32 -11.5])
c=colorbar;
c.Position=[0.6 0.55 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'

%19-63
a1=subplot('Position', subplotPositions(5, :));
pcolor(lat,presm,po4star_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
mask=po4star_rms_gamma>=abs(po4star_change_rho);
stipple(lat,presm,mask)
ylim([0  1000])
% set(gca, 'YTick', [0:200:1000]);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
set(gca, 'YTick', [0:200:1000]);
xlim([-32 -11.5])
caxis([-0.8 0.8])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(c)2019-1963 PO4*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
a2=subplot('Position', subplotPositions(6, :));
pcolor(lat,presm,po4star_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% mask=po4star_rms_gamma>=abs(po4star_change_rho);
stipple(lat,presm,mask)
ylim([0  1000])
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-0.8  0.8])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])
% ylabel('Pressure', 'FontSize', 12);
% set(gca, 'YTick', []);
 set(gca, 'YTick', [1500:500:4500]);
xlim([-32 -11.5])
c=colorbar;
c.Position=[0.935 0.55 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
colormap(a1,cmap);
colormap(a2,cmap);

%% subplot for nstar
%19
subplot('Position', subplotPositions(9, :));
pcolor(lat,presm,nstar19)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-8   8])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(e) 2019 N*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(10, :));
pcolor(lat,presm,nstar19)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])

xlim([-32 -11.5])
caxis([-8   8])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')

set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
set(gca, 'YTickLabel', []);
%set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
% set(gca, 'XTick', []);
xlabel('Latitude','FontSize',12)
xlim([-32 -11.5])

%63
subplot('Position', subplotPositions(7, :));
pcolor(lat,presm,nstar63)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
% set(gca, 'YTickLabel', []);
set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-8   8])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(d) 1963 N*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
subplot('Position', subplotPositions(8, :));
pcolor(lat,presm,nstar63)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
ylim([0  1000])
set(gca, 'YTick', [1500:500:4500]);
ylabel('Pressure', 'FontSize', 12);
xlim([-32 -11.5])
caxis([-8   8])
cmocean('delta',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])
% set(gca, 'YTick', []);
xlabel('Latitude','FontSize',12)
% set(gca, 'XTick', []);
xlim([-32 -11.5])
c=colorbar;
c.Position=[0.6 0.1 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'

%19-63
b1=subplot('Position', subplotPositions(11, :));
pcolor(lat,presm,nstar_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
mask=nstar_rms_gamma>=abs(nstar_change_rho);
stipple(lat,presm,mask)
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'YTickLabel', []);
 set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
xlim([-32 -11.5])
caxis([-8   8])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
title('(f) 2019-1963 N*','FontSize', 12)
set(gca, 'FontSize', 12, 'LineWidth', 1);
%
b2=subplot('Position', subplotPositions(12, :));
pcolor(lat,presm,nstar_change_rho)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h

[c, h] = contour(lat,presm,sa_mean ,[34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat,presm,oxy_mean, [100 245], 'color', 'k', 'linestyle','--','linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h

stipple(lat,presm,mask)
set(gca, 'YTickLabel', []);
%  set(gca, 'XTickLabel', []);
set(gca,'TickDir','out')
set(gca, 'YTick', [1500:500:4500]);
xlim([-32 -11.5])
caxis([-8   8])
cmocean('balance',16)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'LineWidth', 1);
ylim([1000  4500])

xlabel('Latitude','FontSize',12)
% set(gca, 'XTick', []);
xlim([-32 -11.5])
c=colorbar;
c.Position=[0.935 0.1 0.0178 0.4000];
c.Label.String='[{\mu}mol L^{-1}]'
colormap(b1,cmap);
colormap(b2,cmap);