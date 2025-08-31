% This document draw the figure 3 which is the CT and SA different with
% significant stiples from ORAS4 and Argo. To claculate the noise , go the
% the Argo_rms.m and oras4_rms.m where the noise are calculated and
% interpolated into the voyage data form
%% load data
% load voyage_20dbar.mat
% load Argo_data.mat
% load ORAS4.mat
load argo_interp_noise.mat
load oras4_interp_noise_04to17.mat
% v19=ctdHy19v03Slat;
% d63=ctdHy63v02D500lat;
% ctd19=ctd19v03Slat;
load d63.mat
load ctd19.mat
load v19.mat
%%
 %SA mean from voyage
 voyageSA_mean=(d63.SA+ctd19.SA)./2;
%oxygen mean from voyage
voyageOxy_mean=(d63.oxy+ctd19.oxy)./2;
%neutral density mean from voyage
voyagegm_n_mean=(d63.gm_n+ctd19.gm_n)./2;
%% calculate the total noise
total_argo_ct_noise=sqrt((argo_ct_eddy_noise).^2 +(argo_ct_seasonal_noise.^2));
total_argo_sa_noise=sqrt((argo_sa_eddy_noise).^2 +(argo_sa_seasonal_noise.^2));
total_oras4_ct_noise=sqrt((oras4_ct_eddy_noise).^2 +(oras4_ct_seasonal_noise.^2));
total_oras4_sa_noise=sqrt((oras4_sa_eddy_noise).^2 +(oras4_sa_seasonal_noise.^2));


%%
% Create a figure
figure('Units', 'inches', 'Position', [1 1 9 8], 'PaperPositionMode', 'auto');

% Define subplot positions to minimize blank areas
subplotPositions = [
    0.08, 0.75, 0.37, 0.20; 
    0.08, 0.55, 0.37, 0.20; % [left, bottom, width, height] for 1st row, 1st column
    0.53, 0.75, 0.37, 0.20; 
     0.53, 0.55, 0.37, 0.20; % 1st row, 2nd column
    0.08, 0.3, 0.37, 0.20; 
    0.08, 0.1, 0.37, 0.20; % 2nd row, 1st column
    0.53, 0.3, 0.37, 0.20;
      0.53, 0.1, 0.37, 0.20;% 2nd row, 2nd column
    
    ];

% Create subplots
%% Left subplot draw CT difference
% clf
subplot('Position', subplotPositions(1, :));
pcolor(ctd19.lat,ctd19.upres,ctd19.ct-d63.ct)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.7,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',220); clear c h
[c, h] = contour(v19.lat(1:196),v19.upres, voyageSA_mean(:,1:196), [34.7  35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(1:21), voyageSA_mean(1:21,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(41:end), voyageSA_mean(41:end,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [100  245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(iioe_ct_delt_rms2 +iioe_ct_may_rms2)>= abs(v19.ct-d63.ct);
mask=(total_argo_ct_noise)>= abs(v19.ct-d63.ct);
stipple(v19.lat,v19.upres,mask)
mask2=(total_oras4_ct_noise(102:end,:))>= abs(v19.ct(102:end,:)-d63.ct(102:end,:));
stipple(v19.lat,v19.upres(102:end),mask2)
caxis([-2 2])
cmocean('balance',20)
xlim([-32, -11.4]);

% c=colorbar;
% c.Label.String = 'Conservative Temperature [^{\circ}C]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);
title('(a)','FontSize', 14)
% xlabel('Latitude', 'FontSize', 12);
% ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

%%
subplot('Position', subplotPositions(2, :));
pcolor(ctd19.lat,ctd19.upres,ctd19.ct-d63.ct)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',200); clear c h
[c, h] = contour(v19.lat(1:196),v19.upres, voyageSA_mean(:,1:196), [34.7  35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(1:21), voyageSA_mean(1:21,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(41:end), voyageSA_mean(41:end,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [ 100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(iioe_ct_delt_rms2 +iioe_ct_may_rms2)>= abs(v19.ct-d63.ct);
mask=(total_argo_ct_noise)>= abs(v19.ct-d63.ct);
stipple(v19.lat,v19.upres,mask)
mask2=(total_oras4_ct_noise(102:end,:))>= abs(v19.ct(102:end,:)-d63.ct(102:end,:));
stipple(v19.lat,v19.upres(102:end),mask2)
caxis([-2 2])
cmocean('balance',20)
xlim([-32, -11.4]);

% c=colorbar;
% c.Label.String = 'Conservative Temperature [^{\circ}C]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);

% xlabel('Latitude', 'FontSize', 12);
ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);


%%
% Right subplot
subplot('Position', subplotPositions(5, :));
pcolor(v19.lat,v19.upres,ctd19.SA-d63.SA)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','labelspacing',300); clear c h
[c, h] = contour(v19.lat(1:196),v19.upres, voyageSA_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(1:21), voyageSA_mean(1:21,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(41:end), voyageSA_mean(41:end,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [ 100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
% mask=(iioe_sa_delt_rms2 +iioe_sa_may_rms2)>= abs(v19.SA-d63.SA);
mask=(total_argo_sa_noise)>= abs(v19.SA-d63.SA);
stipple(v19.lat,v19.upres,mask)
mask2=(total_oras4_sa_noise(102:end,:))>= abs(v19.SA(102:end,:)-d63.SA(102:end,:));
stipple(v19.lat,v19.upres(102:end),mask2)
caxis([-0.4 0.4])
cmocean('balance',16)
xlim([-32, -11.4]);
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);

% c=colorbar;
% c.Label.String = 'Absolute Salinity [g/kg]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(c)','FontSize', 14)
% xlabel('Latitude', 'FontSize', 12);
% ylabel('Y-axis', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);


%%
subplot('Position', subplotPositions(6, :));
pcolor(v19.lat,v19.upres,ctd19.SA-d63.SA)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:25,26.3,26.9,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w','Labelspacing',300); clear c h
[c, h] = contour(v19.lat(1:196),v19.upres, voyageSA_mean(:,1:196), [34.7  35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(1:21), voyageSA_mean(1:21,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat(210:end),v19.upres(41:end), voyageSA_mean(41:end,210:end), [34.8  34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [ 100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
% mask=(iioe_sa_delt_rms2 +iioe_sa_may_rms2)>= abs(v19.SA-d63.SA);
mask=(total_argo_sa_noise)>= abs(v19.SA-d63.SA);
stipple(v19.lat,v19.upres,mask)
mask2=(total_oras4_sa_noise(102:end,:))>= abs(v19.SA(102:end,:)-d63.SA(102:end,:));
stipple(v19.lat,v19.upres(102:end),mask2)
caxis([-0.4 0.4])
cmocean('balance',16)
xlim([-32, -11.4]);
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);


% c=colorbar;
% c.Label.String = 'Absolute Salinity [g/kg]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
% title('(c)','FontSize', 14)
xlabel('Latitude', 'FontSize', 12);
ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
%%
%transfer T and S to the nuetral density level
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
sa_change_rho = sa19 - sa63;
ct19 = ra_projvar(ctd19.ct, ctd19.gm_n, neutralgrid);
ct63 = ra_projvar(d63.ct, d63.gm_n, neutralgrid);
ct_change_rho = ct19 - ct63;
oxy19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
oxy63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
oxy_mean=(oxy19+oxy63)./2;
%
% load Argo_data.mat
Argo.gamma_mean=squeeze(nanmean(Argo.gamma_n(91,:,:,:),4));
%interp Argo mean neutral density into iioe format
argo_gamma_mean=[];
for i=1:65
argo_gamma_mean(:,i)=interp1(Argo.pressure,squeeze(Argo.gamma_mean(i,:)),ctd19.upres);

end
argo_gamma_mean2=[];
for j=1:291
    argo_gamma_mean2(j,:)=interp1(Argo.latitude,squeeze(argo_gamma_mean(j,:)),ctd19.lat);
end

argo_ct_eddy_noise_rho=ra_projvar(argo_ct_eddy_noise, argo_gamma_mean2 ,neutralgrid);
argo_ct_seasonal_noise_rho=ra_projvar(argo_ct_seasonal_noise, argo_gamma_mean2 ,neutralgrid);
argo_sa_eddy_noise_rho=ra_projvar(argo_sa_eddy_noise, argo_gamma_mean2 ,neutralgrid);
argo_sa_seasonal_noise_rho=ra_projvar(argo_sa_seasonal_noise, argo_gamma_mean2 ,neutralgrid);



%%
%oras4 noise
% load ORAS4.mat
% load oras4_interp_noise.mat

ORAS4.gamma_mean=squeeze(nanmean(ORAS4.gamma_n(:,:,:),3));
%interp Argo mean neutral density into iioe format
oras4_gamma_mean=[];
for i=1:78
oras4_gamma_mean(i,:)=interp1(squeeze(ORAS4.pressure(:,i)),squeeze(ORAS4.gamma_mean(i,:)),ctd19.upres);

end
oras4_gamma_mean2=[];
for j=1:291
    oras4_gamma_mean2(:,j)=interp1(ORAS4.lat,squeeze(oras4_gamma_mean(:,j)),ctd19.lat);
end

oras4_ct_eddy_noise_rho=ra_projvar(oras4_ct_eddy_noise, oras4_gamma_mean2' ,neutralgrid);
oras4_ct_seasonal_noise_rho=ra_projvar(oras4_ct_seasonal_noise, oras4_gamma_mean2' ,neutralgrid);
oras4_sa_eddy_noise_rho=ra_projvar(oras4_sa_eddy_noise, oras4_gamma_mean2' ,neutralgrid);
oras4_sa_seasonal_noise_rho=ra_projvar(oras4_sa_seasonal_noise, oras4_gamma_mean2' ,neutralgrid);



%%

%%
% figure('Units', 'inches', 'Position', [1 1 10 5], 'PaperPositionMode', 'auto');

%% Left subplot draw CT difference
% clf
subplot('Position', subplotPositions(3, :));
pcolor(lats, presm, ct_change_rho)
shading flat
hold on
% [c, h] = contour(lat, presm,presm, [22:1:26,26.5,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
% clabel(c, h, 'fontsize', 12, 'color', 'w'); clear c h
[c, h] = contour(lat(:,1:196), presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(argo_ct_seasonal_noise_rho)>= abs(ct_change_rho);
% stipple(lat, presm,mask)
% mask2=(oras4_ct_seasonal_noise_rho(339:end,:) )>= abs(ct_change_rho(339:end,:));
% stipple(lat(339:end,:), presm(339:end,:),mask2)
caxis([-2 2])
cmocean('balance',20)
xlim([-32, -11.4]);
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);

c=colorbar;
c.Label.String = 'Conservative Temperature [^{\circ}C]';
c.Position=[0.91 0.55 0.0278 0.4000];
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(b)','FontSize', 14)
% xlabel('Latitude', 'FontSize', 12);
% ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

%%
subplot('Position', subplotPositions(4, :));
pcolor(lats, presm, ct_change_rho)
shading flat
hold on
% [c, h] = contour(lat, presm,presm, [22:1:26,26.5,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
% clabel(c, h, 'fontsize', 12, 'color', 'w'); clear c h
[c, h] = contour(lat, presm, sa_mean, [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(argo_ct_seasonal_noise_rho)>= abs(ct_change_rho);
% stipple(lat, presm,mask)
% mask2=(oras4_ct_seasonal_noise_rho(339:end,:) )>= abs(ct_change_rho(339:end,:));
% stipple(lat(339:end,:), presm(339:end,:),mask2)
caxis([-2 2])
cmocean('balance',20)
xlim([-32, -11.4]);
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);

sition=[0.91 0.55 0.0278 0.4000];
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')

set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
%%
% Right subplot
subplot('Position', subplotPositions(7, :));
pcolor(lat, presm,sa_change_rho)
shading flat
hold on
[c, h] = contour(lat(:,1:196), presm(:,1:196), sa_mean(:,1:196), [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(argo_sa_seasonal_noise_rho +argo_sa_eddy_noise_rho)>= abs(sa_change_rho);
% stipple(lat, presm,mask)
% mask2=(oras4_sa_seasonal_noise_rho(339:end,:) )>= abs(sa_change_rho(339:end,:));
% stipple(lat(339:end,:), presm(339:end,:),mask2)
caxis([-0.4 0.4])
cmocean('balance',16)
xlim([-32, -11.4]);
ylim([0  1000])
set(gca, 'YTick', [0:200:1000]);
set(gca, 'XTick', []);

c=colorbar;
c.Position=[0.91 0.1000 0.0278 0.4000];
c.Label.String = 'Absolute Salinity [g/kg]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(d)','FontSize', 14)
% xlabel('Latitude', 'FontSize', 12);
% ylabel('Y-axis', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

%%
%%
% Right subplot
subplot('Position', subplotPositions(8, :));
pcolor(lat, presm,sa_change_rho)
shading flat
hold on
[c, h] = contour(lat, presm, sa_mean, [34.7 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(1:284,210:end), presm(1:284,210:end), sa_mean(1:284,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat(310:end,210:end), presm(310:end,210:end), sa_mean(310:end,210:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [100 245], 'color', 'k', 'linewi', 2,'linestyle','--');
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
% mask=(argo_sa_seasonal_noise_rho +argo_sa_eddy_noise_rho)>= abs(sa_change_rho);
% stipple(lat, presm,mask)
% mask2=(oras4_sa_seasonal_noise_rho(339:end,:) )>= abs(sa_change_rho(339:end,:));
% stipple(lat(339:end,:), presm(339:end,:),mask2)
caxis([-0.4 0.4])
cmocean('balance',16)
xlim([-32, -11.4]);
ylim([1000  4500])
set(gca, 'YTick', [1000:500:4500]);
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
xlabel('Latitude', 'FontSize', 12);
% ylabel('Y-axis', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);





%%
figure('Units', 'inches', 'Position', [1 1 8 8], 'PaperPositionMode', 'auto');

%% Left subplot draw CT difference
clf
subplot(2, 2, 1);
pcolor(ctd19.lat,ctd19.upres,ctd19.ct-d63.ct)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:26,26.5,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageSA_mean, [34.8 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [245 245], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
mask=(oras4_ct_seasonal_noise +oras4_ct_eddy_noise)>= abs(v19.ct-d63.ct);
stipple(v19.lat,v19.upres,mask)
caxis([-2 2])
cmocean('balance')
xlim([-32, -11.4]);
ylim([0 5000]);
c=colorbar;
c.Label.String = 'Conservative Temperature [^{\circ}C]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
set(gca,'xtick',[])
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(a)','FontSize', 14,'Position',[0,0,0])
% xlabel('Latitude', 'FontSize', 12);
ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

% Right subplot
subplot(2, 2, 2);
pcolor(v19.lat,v19.upres,ctd19.SA-d63.SA)
shading flat
hold on
[c, h] = contour(v19.lat,v19.upres, voyagegm_n_mean, [22:1:26,26.5,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
clabel(c, h, 'fontsize', 12, 'color', 'w'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageSA_mean, [34.8 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(v19.lat,v19.upres, voyageOxy_mean, [245 245], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
mask=(oras4_sa_seasonal_noise +oras4_sa_eddy_noise)>= abs(v19.SA-d63.SA);
stipple(v19.lat,v19.upres,mask)
caxis([-0.4 0.4])
cmocean('balance')
xlim([-32, -11.4]);
ylim([0 5000]);
c=colorbar;
c.Label.String = 'Absolute Salinity [g/kg]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(b)','FontSize', 14)
set(gca,'xtick',[])
set(gca,'ytick',[])
% xlabel('Latitude', 'FontSize', 12);
% ylabel('Y-axis', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

% Adjust the spacing
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));

subplot(2, 2, 3);
pcolor(lats, presm, ct_change_rho)
shading flat
hold on
% [c, h] = contour(lat, presm,presm, [22:1:26,26.5,27,27.1,27.4:0.2:28.2], 'color', 'w', 'linewi', 1);
% clabel(c, h, 'fontsize', 12, 'color', 'w'); clear c h
[c, h] = contour(lat, presm, sa_mean, [34.8 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [245 245], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
mask=(oras4_ct_seasonal_noise_rho +oras4_ct_eddy_noise_rho)>= abs(ct_change_rho);
stipple(lat, presm,mask)
caxis([-2 2])
cmocean('balance')
xlim([-32, -11.4]);
ylim([0 5000]);
c=colorbar;
c.Label.String = 'Conservative Temperature [^{\circ}C]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(c)','FontSize', 14)
xlabel('Latitude', 'FontSize', 12);
ylabel('Pressure', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);

% Right subplot
subplot(2, 2, 4);
pcolor(lat, presm,sa_change_rho)
shading flat
hold on
[c, h] = contour(lat, presm, sa_mean, [34.8 35.8], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
[c, h] = contour(lat, presm, oxy_mean, [245 245], 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k'); clear c h
% plot(-23,v19.upres,'.','color',[0.4 0.4 0.4])
%add stipples
mask=(oras4_sa_seasonal_noise_rho +oras4_sa_eddy_noise_rho)>= abs(sa_change_rho);
stipple(lat, presm,mask)
caxis([-0.4 0.4])
cmocean('balance')
xlim([-32, -11.4]);
ylim([0 5000]);
c=colorbar;
c.Label.String = 'Absolute Salinity [g/kg]';
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% set(gca, 'ytick', 0:200:2000)
%  ylim([0 2000]);
title('(d)','FontSize', 14)
xlabel('Latitude', 'FontSize', 12);
% ylabel('Y-axis', 'FontSize', 12);
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
set(gca,'ytick',[])
% Adjust the spacing
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));

%%
cmap = cmocean('balance', 22);
cmap=cmap([1:10,13:22],:);
colormap(cmap); % Apply the adjusted colormap
