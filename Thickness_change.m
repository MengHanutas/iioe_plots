%voayge data depth and thickness change 
%%
% load voyage_20dbar.mat
% v19=ctdHy19v03Slat;
% d63=ctdHy63v02D500lat;
% ctd19=ctd19v03Slat;
load v19.mat
load ctd19.mat
load d63.mat
%% spiciness: change on neutral surfaces
neutralgrid = 21.2:0.1:28.2;
%% 
p19=nan*ones(length(neutralgrid),286);
for i=1:286
    ind=find(~isnan(ctd19.gm_n(:,i)));
    if ~isempty(ind)
p19(:,i) = interp1(ctd19.gm_n(ind,i),ctd19.upres(ind),neutralgrid);
    end
end

p63=nan*ones(length(neutralgrid),286);
for i=1:286
    ind=find(~isnan(d63.gm_n(:,i)));
    if ~isempty(ind)
p63(:,i) = interp1(d63.gm_n(ind,i),d63.upres(ind),neutralgrid);
    end
end
% %transfer pressure to height
H19=gsw_z_from_p(p19,ctd19.lat);
H63=gsw_z_from_p(p63,d63.lat);
% height=gsw_z_from_p(repmat(ctd19.upres,1,286),repmat(ctdHy19v03Slat.lat,291,1));
% %height=squeeze(height(:,1));
% %transfer depth from height
D19=gsw_depth_from_z(H19);
D63=gsw_depth_from_z(H63);

Thickness_19=[];
Thickness_63=[];
for i=3:length(neutralgrid)
    Thickness_19(i-2,:)=D19(i,:)-D19(i-2,:);
    Thickness_63(i-2,:)=D63(i,:)-D63(i-2,:); %thickness 349*286
end
Thickness_change=Thickness_19-Thickness_63;
%thickness change fraction
Thickness_change_fraction=Thickness_change./Thickness_63;
%% making averaged of thickness first and then calculate fraction
 i_lat=find(lats>-32 & lats<-23);
 Th19_av=nanmean(Thickness_19(:,i_lat),2);
 Th63_av=nanmean(Thickness_63(:,i_lat),2);
 Th_change_av=nanmean(Thickness_change(:,i_lat),2);
 Th_change_f_av=Th_change_av./Th63_av;
%%
FX = 19; FY = 20;
% FIG1 = figure('units', 'centimeters', 'Position', [0 0 FX FY]);
FIG1 = figure('units', 'centimeters', 'Position', [2 2 FX FY]);
width = 14; height = 17;
left = 2;
bottom = 1.5;

%%
pcolor(ctd19.lat,neutralgrid(2:end-1),Thickness_change)
shading interp
colorbar
cmocean('balance',30)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
caxis([-30 30])
hold on
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [34.8,35.8], 'linest', ':', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [245 245], 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
xlim([-32, -11.4]);
set(gca, 'ytick',21.4:0.2:28.2)
title('Thickness change of isopycnals [m]')

%%
pcolor(ctd19.lat,neutralgrid(2:end-1),Thickness_change)
shading flat
colorbar
cmocean('balance',30)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
caxis([-30 30])
hold on
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [34.8,35.8], 'linest', ':', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [245 245], 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
xlim([-32, -11.4]);
set(gca, 'ytick',21.4:0.2:28.2)
title('Thickness change of isopycnals [m]')

%%
pcolor(ctd19.lat,neutralgrid(2:end-1),movmean(movmean(Thickness_change_fraction,20,2),10,2))
shading interp
colorbar
cmocean('balance',30)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
caxis([-1 1])
hold on
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [34.8,35.8], 'linest', ':', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [245 245], 'linest', '-', 'color', 'w', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold'); clear c h
xlim([-32, -11.4]);
set(gca, 'ytick',21.4:0.2:28.2)
title('Thickness change of isopycnals [m]')
%%
%making average
%STW 25-26.2
i_lat=find(ctd19.lat>=-32 & ctd19.lat <=-23);
i_iso=find(neutralgrid(2:end-1)>=25 & neutralgrid(2:end-1)<=26);
Th_stw=nanmean(Thickness_change(i_iso,i_lat),'all')

%SAMW 26.8-27
i_iso=find(neutralgrid(2:end-1)>=26.8 & neutralgrid(2:end-1)<=27);
Th_samw=nanmean(Thickness_change(i_iso,i_lat),'all')

%AAIW 27.1-27.4
i_iso=find(neutralgrid(2:end-1)>=27.1 & neutralgrid(2:end-1)<=27.4);
Th_aaiw=nanmean(Thickness_change(i_iso,i_lat),'all')

% all layers
i_lat=find(ctd19.lat>=-32 & ctd19.lat <=-23);
Th_layers=nanmean(Thickness_change(:,i_lat),2);
figure
plot(Th_layers,neutralgrid(2:end-1),'-')
hold on
plot(repmat(0,1,length(neutralgrid(2:end-1))),neutralgrid(2:end-1))
set(gca, 'ytick',21.4:0.2:28.2)
ylim([23.8 28])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
xlabel('Thickness change [m]')
ylabel('neutral density')
plot([-40:1:40],repmat(25,1,81),'k-')
plot([-40:1:40],repmat(26.2,1,81),'k-')
plot([-40:1:40],repmat(26.8,1,81),'k-')
plot([-40:1:40],repmat(27,1,81),'k-')
plot([-40:1:40],repmat(27.1,1,81),'k-')
plot([-40:1:40],repmat(27.4,1,81),'k-')

% all layers fraction
i_lat=find(ctd19.lat>=-32 & ctd19.lat <=-23);
Th_layers_f=nanmean(Thickness_change_fraction(:,i_lat)*100,2);
figure
plot(Th_layers_f,neutralgrid(2:end-1),'-')
hold on
plot(repmat(0,1,length(neutralgrid(2:end-1))),neutralgrid(2:end-1))
set(gca, 'ytick',21.4:0.2:28.2)
ylim([23.8 28])
xlim([-50 50])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
xlabel('Thickness change fraction [%]')
ylabel('neutral density')
plot([-50:1:50],repmat(25,1,101),'k-')
plot([-50:1:50],repmat(26.2,1,101),'k-')
plot([-50:1:50],repmat(26.8,1,101),'k-')
plot([-50:1:50],repmat(27,1,101),'k-')
plot([-50:1:50],repmat(27.1,1,101),'k-')
plot([-50:1:50],repmat(27.4,1,101),'k-')


% all layers fraction (averaging thickness first and then calculate
% fraction
i_lat=find(ctd19.lat>=-32 & ctd19.lat <=-23);
Th_19=nanmean(Thickness_19(:,i_lat),2);
Th_63=nanmean(Thickness_63(:,i_lat),2);
Th_change=Th_19-Th_63;
Th_change_f=Th_change./Th_63.*100;
figure
plot(Th_layers_f,neutralgrid(2:end-1),'-')
hold on
plot(repmat(0,1,length(neutralgrid(2:end-1))),neutralgrid(2:end-1))
set(gca, 'ytick',21.4:0.2:28.2)
ylim([23.8 28])
xlim([-50 50])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
xlabel('Thickness change fraction [%]')
ylabel('neutral density')
plot([-50:1:50],repmat(25,1,101),'k-')
plot([-50:1:50],repmat(26.2,1,101),'k-')
plot([-50:1:50],repmat(26.8,1,101),'k-')
plot([-50:1:50],repmat(27,1,101),'k-')
plot([-50:1:50],repmat(27.1,1,101),'k-')
plot([-50:1:50],repmat(27.4,1,101),'k-')
%%
%add depth change
hold on
D_change=D19-D63;
D_layers=nanmean(D_change(:,i_lat),2);
plot(D_layers,neutralgrid,'b.-')
xlim([-60 60])
plot([-60:1:60],repmat(25,1,121),'k-')
plot([-60:1:60],repmat(26.2,1,121),'k-')
plot([-60:1:60],repmat(26.8,1,121),'k-')
plot([-60:1:60],repmat(27,1,121),'k-')
plot([-60:1:60],repmat(27.1,1,121),'k-')
plot([-60:1:60],repmat(27.4,1,121),'k-')


%% put figures together
figure('Units', 'inches', 'Position', [1 1 14 6], 'PaperPositionMode', 'auto');
subplot('Position', [0.07, 0.1, 0.3, 0.85]);
%plot depth change
pcolor(ctd19.lat,neutralgrid,movmean(D_change,10,2))
shading interp
colorbar
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
caxis([-100 100])
hold on
 sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
           
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
sa_mean=(sa19+sa63)./2;
[c, h] = contour(ctd19.lat(1:196), neutralgrid,sa_mean(:,1:196) , [34.7,35.8], 'linest', '-', 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','Labelspacing',300); clear c h
[c, h] = contour(lats(220:end), neutralgrid, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
plot(ctd19.lat(i_lat),repmat(25,length(ctd19.lat(i_lat)),1),'--w','LineWidth',2)
i_n=find(neutralgrid>=25);
plot(repmat(-23,length(neutralgrid(i_n)),1),neutralgrid(i_n),'--w','LineWidth',2)
xlim([-32, -11.4]);
ylim([22 28.2])
set(gca, 'ytick',22:0.2:28.2)
ylabel('Neutral Density [kg.m^{-3}]')
xlabel('Latitude')
title('(a) Depth change [m]')
set(gca,'FontSize',12)
 set(gca,'TickDir','out')

% subplot('Position', [0.41, 0.1, 0.3, 0.85]);
% pcolor(ctd19.lat,neutralgrid(2:end-1),movmean(Thickness_change_fraction.*100,10,2))
% shading interp
% colorbar
% cmocean('balance',20)
% set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
% caxis([-100 100])
% hold on
% sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
% sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
% sa_mean=(sa19+sa63)./2;
% [c, h] = contour(ctd19.lat(1:196), neutralgrid,sa_mean(:,1:196) , [34.7,35.8], 'linest', '-', 'color', 'k', 'linewi', 2);
% clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','Labelspacing',300); clear c h
% [c, h] = contour(lats(220:end), neutralgrid, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
% sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
% sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
% [c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
% clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
% plot(ctd19.lat(i_lat),repmat(25,length(ctd19.lat(i_lat)),1),'--w','LineWidth',2)
% i_n=find(neutralgrid>=25);
% plot(repmat(-23,length(neutralgrid(i_n)),1),neutralgrid(i_n),'--w','LineWidth',2)
% xlim([-32, -11.4]);
% ylim([22 28.2])
% set(gca, 'ytick',22:0.2:28.2)
% set(gca,'TickDir','out')
% % ylabel('Neutral Density [kg.m^{-3}]')
% xlabel('Latitude')
% title('(b) Thickness change fraction[%]')
% set(gca,'FontSize',12)


subplot('Position', [0.41, 0.1, 0.3, 0.85]);
pcolor(ctd19.lat,neutralgrid(2:end-1),movmean(Thickness_change,10,2))
shading interp
colorbar
cmocean('balance',20)
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
caxis([-30 30])
hold on
 sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
           
sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
sa_mean=(sa19+sa63)./2;
[c, h] = contour(ctd19.lat(1:196), neutralgrid,sa_mean(:,1:196) , [34.7,35.8], 'linest', '-', 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','Labelspacing',300); clear c h
[c, h] = contour(lats(220:end), neutralgrid, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
[c, h] = contour(ctd19.lat, neutralgrid, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
plot(ctd19.lat(i_lat),repmat(25,length(ctd19.lat(i_lat)),1),'--w','LineWidth',2)
i_n=find(neutralgrid>=25);
plot(repmat(-23,length(neutralgrid(i_n)),1),neutralgrid(i_n),'--w','LineWidth',2)
xlim([-32, -11.4]);
ylim([22 28.2])
set(gca, 'ytick',22:0.2:28.2)
set(gca,'TickDir','out')
% ylabel('Neutral Density [kg.m^{-3}]')
xlabel('Latitude')
title('(b) Thickness change [m]')
set(gca,'FontSize',12)

subplot('Position', [0.75, 0.1, 0.23, 0.85]);

i_lat=find(ctd19.lat>=-32 & ctd19.lat <=-23);
%add depth change
D_change=D19-D63;
D_layers=nanmean(D_change(:,i_lat),2);
plot(D_layers,neutralgrid,'Color',[0.9294 0.4745 0.1020],'LineWidth',2)
hold on
Th_layers=nanmean(Thickness_change(:,i_lat),2);

% plot(Th_layers,neutralgrid(2:end-1),'Color', [0.0784    0.3333    0.4588],'LineWidth',2,'LineStyle','-')
hold on
Th_f_layers=nanmean(Thickness_change_fraction(:,i_lat),2).*100;

plot(Th_f_layers,neutralgrid(2:end-1),'Color', [0.0784    0.3333    0.4588],'LineWidth',2,'LineStyle','-')
% plot(Th_change_f_av*100,neutralgrid(2:end-1),'k','LineWidth',2)
plot(Th_layers,neutralgrid(2:end-1),'Color', 'k','LineWidth',2,'LineStyle','-')
set(gca, 'ytick',25:0.2:28.2)
ylim([25 28.2])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
xlabel('Depth[m]/ Thickness[%] change')
ylabel('neutral density')
title('(C) Averaged change')

xlim([-150 150])
plot(repmat(0,1,length(neutralgrid(2:end))),neutralgrid(2:end),'k')
plot([-150:1:150],repmat(25,1,301),'k--')
plot([--150:1:150],repmat(26.3,1,301),'k--')
plot([-150:1:150],repmat(26.7,1,301),'k--')
plot([-150:1:150],repmat(27,1,301),'k--')
plot([-150:1:150],repmat(27.1,1,301),'k--')
plot([-150:1:150],repmat(27.6,1,301),'k--')
set(gca,'FontSize',12)
 set(gca,'TickDir','out')

% calculate and add the errorbar
for i=1:length(neutralgrid)
    D_error(i)=1.697.*std(squeeze(D_change(i,i_lat)),'omitnan')./sqrt(30);
    
end
for i=1:length(neutralgrid)-2
    Th_error(i)=1.697.*std(squeeze(Thickness_change(i,i_lat)),'omitnan')./sqrt(30);
end
for i=1:length(neutralgrid)-2
    Th_f_error(i)=1.697.*std(squeeze(Thickness_change_fraction(i,i_lat)).*100,'omitnan')./sqrt(30);
end
ind=~isnan(squeeze(nanmean(D_change(:,i_lat),2))');
 patch([squeeze(nanmean(D_change(ind,i_lat),2))'+D_error(ind) flip(squeeze(nanmean(D_change(ind,i_lat),2))'-D_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9294 0.4745 0.1020],'FaceAlpha',0.25, 'EdgeColor','none')
%  ind=~isnan(squeeze(nanmean(Thickness_change(:,i_lat),2))');
%  Th_neutralgrid=neutralgrid(2:end-1);
%  patch([squeeze(nanmean(Thickness_change(ind,i_lat),2))'+Th_error(ind) flip(squeeze(nanmean(Thickness_change(ind,i_lat),2))'-Th_error(ind))],[Th_neutralgrid(ind) flip(Th_neutralgrid(ind))],[0.0784    0.3333    0.4588],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(Thickness_change_fraction(:,i_lat),2))');
 Th_neutralgrid=neutralgrid(2:end-1);
 patch([squeeze(nanmean(Thickness_change_fraction(ind,i_lat),2))'.*100+Th_f_error(ind) flip(squeeze(nanmean(Thickness_change_fraction(ind,i_lat),2))'.*100-Th_f_error(ind))],[Th_neutralgrid(ind) flip(Th_neutralgrid(ind))],[0.0784    0.3333    0.4588],'FaceAlpha',0.25, 'EdgeColor','none')

 legend('Depth [m]', 'Thickness [%]')