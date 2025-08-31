%% scrippt for calculate heave spiciness and draw figure 4

%% load data:
load d63.mat
load v19.mat
load ctd19.mat
%need to add ra_projvar.m to path for interpolate
%% calculate
neutralgrid = 21.2:0.02:28.2;
%% 
%calculate spiciness 
%transform the tempearture and salinity into the neutral density form
T19 = ra_projvar(ctd19.ct, ctd19.gm_n, neutralgrid);
T63 = ra_projvar(d63.ct, d63.gm_n, neutralgrid);
S19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
S63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
p19=nan*ones(351,286);
for i=1:286
    ind=find(~isnan(ctd19.gm_n(:,i)));
    if ~isempty(ind)
p19(:,i) = interp1(ctd19.gm_n(ind,i),ctd19.upres(ind),neutralgrid);
    end
end

p63=nan*ones(351,286);
for i=1:286
    ind=find(~isnan(d63.gm_n(:,i)));
    if ~isempty(ind)
p63(:,i) = interp1(d63.gm_n(ind,i),d63.upres(ind),neutralgrid);
    end
end

%calculate alpha and betta on the neutral density level
alpha_19_n=gsw_alpha_CT_exact(S19,T19,p19);
alpha_63_n=gsw_alpha_CT_exact(S63,T63,p63);
alpha_mean=(alpha_19_n+alpha_63_n)./2;
betta_19_n=gsw_beta_CT_exact(S19,T19,p19);
betta_63_n=gsw_beta_CT_exact(S63,T63,p63);
betta_mean=(betta_63_n+betta_19_n)./2;

temp_n = T19 - T63; %this is temperature change on each neutral surfaces
salt_n = S19 - S63; %this is salinity change on each neutral surfaces

%temperature and salinity change on neutral surfaces with density unit:
%kg/m3
% alpha_temp_n=alpha_19_n.*T19-alpha_63_n.*T63;
% betta_salt_n=betta_19_n.*S19-betta_63_n.*S63;
%use mean alpha betta 
alpha_temp_n=alpha_mean.*(T19-T63);
betta_salt_n=betta_mean.*(S19-S63);

%% heave

%% N prime : depth change of neutral surfaces
% %transfer pressure to height
H19=gsw_z_from_p(p19,ctd19.lat);
H63=gsw_z_from_p(p63,d63.lat);
% height=gsw_z_from_p(repmat(ctd19.upres,1,286),repmat(ctdHy19v03Slat.lat,291,1));
% %height=squeeze(height(:,1));
% %transfer depth from height
D19=gsw_depth_from_z(H19);
D63=gsw_depth_from_z(H63);
% depth=gsw_depth_from_z(height);
% 
% D19 = ra_projvar(depth, ctd19.gm_n, neutralgrid);
% D63 = ra_projvar(depth, d63.gm_n, neutralgrid);
%depth change:
N_prim=D19-D63;

%% T gradient mean
% T gradient 2019 and 1963:
T_gradient_19=diff(T19)./diff(D19); %size 350*286; this is gradient on the middle of each depth levels
T_gradient_63=diff(T63)./diff(D63); 
%need to interp the gradient from mid of two date point to the point:
% mid_point_depth19=diff(D19)/2+D19(1:end-1,:);
% mid_point_depth63=diff(D63)/2+D63(1:end-1,:);
%mid point neutral density
mid_neutral=neutralgrid(1:end-1)+0.01;
T_gradient_19_interp=nan*ones(351,286);
T_gradient_63_interp=nan*ones(351,286);
for i=1:286 % length of latitude
%     index1=~isnan(mid_point_depth19(:,i));
%     index2=~isnan(mid_point_depth63(:,i));
%     if sum(index1)>=2
    T_gradient_19_interp(:,i)=interp1(mid_neutral,squeeze(T_gradient_19(:,i)),neutralgrid);
    T_gradient_63_interp(:,i)=interp1(mid_neutral,squeeze(T_gradient_63(:,i)),neutralgrid);
%     end
%     if sum(index2)>=2
%     T_gradient_63_interp(index2,i)=interp1(squeeze(mid_point_depth63(index2,i)),squeeze(T_gradient_63(index2,i)),squeeze(D63(index2,i)));
%     end
end

T_gradient_mean=(T_gradient_19_interp+T_gradient_63_interp)./2;

%% S `  
% gradient mean
% S gradient 2019 and 1963:
S_gradient_19=diff(S19)./diff(D19); %size 353*286; this is gradient on the middle of each depth levels
S_gradient_63=diff(S63)./diff(D63); 
%need to interp the gradient from mid of two date point to the point:

S_gradient_19_interp=nan*ones(351,286);
S_gradient_63_interp=nan*ones(351,286);
for i=1:286 % length of latitude
   
    S_gradient_19_interp(:,i)=interp1(mid_neutral,squeeze(S_gradient_19(:,i)),neutralgrid);
    S_gradient_63_interp(:,i)=interp1(mid_neutral,squeeze(S_gradient_63(:,i)),neutralgrid);    
end

S_gradient_mean=(S_gradient_19_interp+S_gradient_63_interp)./2;

%% heave temperatur*-+
% +*/
T_heave=N_prim.*T_gradient_mean;
%% heave salinity
S_heave=N_prim.* S_gradient_mean;

%% heave with al[ha and betta
% T_heave_alpha=N_prim.*(T_gradient_19_interp.*alpha_19_n+T_gradient_63_interp.*alpha_63_n)./2;
% S_heave_alpha=N_prim.*(S_gradient_19_interp.*betta_19_n+S_gradient_63_interp.*betta_63_n)./2;

T_heave_alpha=N_prim.*alpha_mean.*(T_gradient_19_interp+T_gradient_63_interp)./2;
S_heave_alpha=N_prim.*betta_mean.*(S_gradient_19_interp+S_gradient_63_interp)./2;

%% Change on depth level
%% Temperature change on depth level
T_z=(ctd19.ct)-(d63.ct);

%% Salinity change on depth level
S_z=(ctd19.SA)-(d63.SA);

%% mean pressure
lats = ctd19.lat;
pres = repmat(ctd19.upres, 1, length(ctd19.lat)); % since both data sets were interpolated on the same pressure levels at 5 dbar interval
% p63iso = ra_projvar(pres, d63.gm_n, neutralgrid);
% p19iso = ra_projvar(pres, ctd19.gm_n, neutralgrid);

% Creating Pm at mid point of two voyages.
presm = (p63 + p19)/2;
% clear p19iso p63iso 

%% interp
for i=1:286
T_z_n(:,i)=interp1(pres(:,i), T_z(:,i),presm(:,i));
S_z_n(:,i)=interp1(pres(:,i), S_z(:,i),presm(:,i));
end

%multiply alpha and betta
T_z_alpha=alpha_mean.*T_z_n;
S_z_betta=betta_mean.*S_z_n;

%% figures4
%%  Below is the code for the figure 4
%plot heave and spiciness figures prepare for the labels 
neutralgrid =[21.2:0.02:28.2];

lats = ctdHy63v02D500lat.lat;
nut = neutralgrid;
pres = repmat(ctdHy63v02D500lat.upres, 1, length(ctdHy63v02D500lat.lat)); % since both data sets were interpolated on the same pressure levels at 5 dbar interval
p63iso = ra_projvar(pres, d63.gm_n, neutralgrid);
p19iso = ra_projvar(pres, ctd19.gm_n, neutralgrid);

% Creating Pm at mid point of two voyages.
presm = (p63iso + p19iso)/2; %351*286)
% clear p19iso p63iso 
lat = repmat(lats, length(presm), 1);



%%

% Create a figure
 figure('Units', 'inches', 'Position', [1 1 11 8], 'PaperPositionMode', 'auto');

% Define subplot positions to minimize blank areas
subplotPositions = [
    0.07, 0.87, 0.27, 0.08;  % [left, bottom, width, height] for 1st row, 1st column
    0.36, 0.87, 0.27, 0.08;  % 1st row, 2nd column
    0.655, 0.87, 0.27, 0.08;  % 1st row, 3rd column
    0.07, 0.4, 0.27, 0.08;  % 2nd row, 1st column
    0.36, 0.4, 0.27, 0.08;  % 2nd row, 2nd column
    0.655, 0.4, 0.27, 0.08   % 2nd row, 3rd column
    ];

% Create subplots
for i = 1:6
    subplot('Position', subplotPositions(i, :));
    switch i
        case 1
         
            pcolor(lats,nut,T_z_n)
            shading flat
            hold on
            sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
           
            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Temperature [^{\circ}C]';
%             c.Label.FontSize=12;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
            title('(a) \theta''|z','FontSize',12)
%             set(gca,'FontSize',12);

            ylabel('Neutral Density','FontSize',12)
            set(gca,'TickDir','out')
            ylim([21.5, 25])
            set(gca, 'ytick',21.5:1:24.5,'FontSize',12)
            %set(gca, 'xtick',[])
        case 2
            pcolor(lats,nut,temp_n)
            shading flat
            hold on
             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Temperature [^{\circ}C]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
            title('(b) \theta''|n (spiciness)','FontSize',12)
            set(gca,'FontSize',12);
            % xlabel('Latitude','FontSize',16)
            % ylabel('Neutral Density','FontSize',14)
            set(gca,'TickDir','out')
            ylim([21.5, 25])
            
            set(gca, 'ytick',[])
            %set(gca, 'xtick',[])
        case 3
            pcolor(lats,nut,-T_heave)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h

            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
            c=colorbar;
            c.Label.String = 'Temperature [^{\circ}C]';
            c.Label.FontSize=12;
            c.Label.Position=[2.6, 0];
           c.Position=[0.93, 0.55, 0.015, 0.4];
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            title('(c) N'' \theta |z (heave)','FontSize',12,'Interpreter','Tex')
            set(gca,'FontSize',12);
            set(gca,'TickDir','out')
             ylim([21.5, 25])
            
            set(gca, 'ytick',[])
            %set(gca, 'xtick',[])
        case 4
            pcolor(lats,nut,S_z_n)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Salinity [PSU]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')

            title('(d) S''|z','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)
            ylabel('Neutral Density','FontSize',14)
            set(gca,'TickDir','out')
             ylim([21.5, 25])
            set(gca, 'ytick',21.5:1:24.5,'FontSize',12)
            
        case 5
            pcolor(lats,nut,salt_n)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Salinity [PSU]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            title('(e) S''|n (spiciness)','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)
            % ylabel('Neutral Density','FontSize',16)
            set(gca,'TickDir','out')
            ylim([21.5, 25])
           
            set(gca, 'ytick',[])
        case 6
            pcolor(lats,nut,-S_heave)
            shading flat
            hold on
             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
            c=colorbar;
            c.Label.String = 'Salinity [PSU]';
            c.Label.FontSize=12;
             c.Label.Position=[2.6, 0];
           c.Position=[0.93, 0.08, 0.015, 0.4];
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
            title('(f) N'' S|z (heave)','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)

            set(gca,'TickDir','out')
             ylim([21.5, 25])
            
            set(gca, 'ytick',[])
    end
    %     xlabel('X-axis', 'FontSize', 12);
    %     ylabel('Y-axis', 'FontSize', 12);
    %     set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
end

% Save the figure as a high-quality image
% print(gcf, 'published_quality_figure_2x3', '-dpng', '-r300');

cmap = cmocean('balance', 22);
cmap=cmap([1:10,13:22],:);
colormap(cmap); % 

%%

% Create a figure
% figure('Units', 'inches', 'Position', [1 1 11 8], 'PaperPositionMode', 'auto');

% Define subplot positions to minimize blank areas
subplotPositions = [
    0.07, 0.55, 0.27, 0.32;  % [left, bottom, width, height] for 1st row, 1st column
    0.36, 0.55, 0.27, 0.32;  % 1st row, 2nd column
    0.655, 0.55, 0.27, 0.32;  % 1st row, 3rd column
    0.07, 0.08, 0.27, 0.32;  % 2nd row, 1st column
    0.36, 0.08, 0.27, 0.32;  % 2nd row, 2nd column
    0.655, 0.08, 0.27, 0.32   % 2nd row, 3rd column
    ];

% Create subplots
for i = 1:6
    subplot('Position', subplotPositions(i, :));
    switch i
        case 1
         
            pcolor(lats,nut,T_z_n)
            shading flat
            hold on
            sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
           
            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Temperature [^{\circ}C]';
%             c.Label.FontSize=12;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
%             title('(a) \theta''|z','FontSize',12)
%             set(gca,'FontSize',12);

            ylabel('Neutral Density','FontSize',12)
            set(gca,'TickDir','out')
            ylim([25, 28.2])
            set(gca, 'ytick',25:0.2:28.2,'FontSize',12)
            %set(gca, 'xtick',[])
        case 2
            pcolor(lats,nut,temp_n)
            shading flat
            hold on
             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Temperature [^{\circ}C]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
%             title('(b) \theta''|n (spiciness)','FontSize',12)
            set(gca,'FontSize',12);
            % xlabel('Latitude','FontSize',16)
            % ylabel('Neutral Density','FontSize',14)
            set(gca,'TickDir','out')
            ylim([25, 28.2])
           
            set(gca, 'ytick',[])
            %set(gca, 'xtick',[])
        case 3
            pcolor(lats,nut,-T_heave)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h

            caxis([-2 2])
            cmocean('balance',20)
            xlim([-32, -11.4]);
            c=colorbar;
            c.Label.String = 'Temperature [^{\circ}C]';
            c.Label.FontSize=12;
            c.Label.Position=[2.6, 0];
           c.Position=[0.93, 0.55, 0.015, 0.4];
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
%             title('(c) N'' \theta |z (heave)','FontSize',12,'Interpreter','Tex')
            set(gca,'FontSize',12);
            set(gca,'TickDir','out')
            ylim([25, 28.2])
            
            set(gca, 'ytick',[])
            %set(gca, 'xtick',[])
        case 4
            pcolor(lats,nut,S_z_n)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Salinity [PSU]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')

%             title('(d) S''|z','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)
            ylabel('Neutral Density','FontSize',14)
            set(gca,'TickDir','out')
            ylim([25, 28.2])
            set(gca, 'ytick',25:0.2:28.2,'FontSize',12)
        case 5
            pcolor(lats,nut,salt_n)
            shading flat
            hold on
              sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
%             c=colorbar;
%             c.Label.String = 'Salinity [PSU]';
%             c.Label.FontSize=14;
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
%             title('(e) S''|n (spiciness)','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)
            % ylabel('Neutral Density','FontSize',16)
            set(gca,'TickDir','out')
           ylim([25, 28.2])
            
            set(gca, 'ytick',[])
        case 6
            pcolor(lats,nut,-S_heave)
            shading flat
            hold on
             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500,1000:1000:4500], 'linest', ':', 'color', 'w', 'linewi', 2);
            clabel(c, h,[100:100:300,500,1000], 'fontsize', 12, 'color', 'w', 'fontweigh', 'bold','Labelspacing',300); clear c h
            sa19 = ra_projvar(ctd19.SA, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.SA, d63.gm_n, neutralgrid);
            sa_mean=(sa19+sa63)./2;
            [c, h] = contour(lats(1:196), nut, sa_mean(:,1:196), [34.7,35.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            [c, h] = contour(lats(220:end), nut, sa_mean(:,220:end), [34.8 34.8], 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
%             [c, h] = contour(lats(210:end), nut(310:end), sa_mean(310:end,210:end), [34.7,35.8], 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            sa19 = ra_projvar(ctd19.oxy, ctd19.gm_n, neutralgrid);
            sa63 = ra_projvar(d63.oxy, d63.gm_n, neutralgrid);
            [c, h] = contour(lats, nut, (sa19+sa63)./2, [100 245], 'linest', '--', 'color', 'k', 'linewi', 2);
            clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold','labelspacing',200); clear c h
%             sa19 = ra_projvar(meshgrid(ctd19.upres,ctd19.lat)', ctd19.gm_n, neutralgrid);
%             sa63 = ra_projvar(meshgrid(d63.upres,d63.lat)', d63.gm_n, neutralgrid);
%             [c, h] = contour(lats, nut, (sa19+sa63)./2, [100:100:300,500:1000:4500], 'linest', '-', 'color', 'k', 'linewi', 2);
%             clabel(c, h, 'fontsize', 12, 'color', 'k', 'fontweigh', 'bold'); clear c h
            caxis([-0.4 0.4])
            cmocean('balance',16)
            xlim([-32, -11.4]);
            c=colorbar;
            c.Label.String = 'Salinity [PSU]';
            c.Label.FontSize=12;
             c.Label.Position=[2.6, 0];
           c.Position=[0.93, 0.08, 0.015, 0.4];
            set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom')
            % set(gca, 'ytick', 0:200:2000)
            %  ylim([0 2000]);
%             title('(f) N'' S|z (heave)','FontSize',12)
            set(gca,'FontSize',12);
            xlabel('Latitude','FontSize',14)

            set(gca,'TickDir','out')
            ylim([25, 28.2])
           
            set(gca, 'ytick',[])
    end
    %     xlabel('X-axis', 'FontSize', 12);
    %     ylabel('Y-axis', 'FontSize', 12);
    %     set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.5);
end

% Save the figure as a high-quality image
% print(gcf, 'published_quality_figure_2x3', '-dpng', '-r300');

cmap = cmocean('balance', 22);
cmap=cmap([1:10,13:22],:);
colormap(cmap); % 

%% figure 5 ==============================
S_z_mean=((ctd19.SA)+(d63.SA))./2;
for i=1:286
S_z_mean_n(:,i)=interp1(pres(:,i), S_z_mean(:,i),presm(:,i));
end
 i_lat=find(lats>-32 & lats<-23);
S_Z_mean_n_avg=nanmean(S_z_mean_n(:,i_lat),2);
index_s_max=find(S_Z_mean_n_avg==max(S_Z_mean_n_avg));
index_s_min=find(S_Z_mean_n_avg==min(S_Z_mean_n_avg));
%%
figure('Units', 'inches', 'Position', [1 1 10 8], 'PaperPositionMode', 'auto');
%%
%subplot('Position', [0.07, 0.4, 0.43, 0.28]);
subplot('Position', [0.07, 0.83, 0.43, 0.12]);

plot(squeeze(nanmean(T_z_alpha(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(alpha_temp_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)
plot(squeeze(nanmean(-T_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot(repmat(0,351,1),neutralgrid,'k')
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00028,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
text(-0.00028,25.4,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])

xlim([-3e-4 3e-4])
ylim([23.5 25])
% set(gca, 'ytick',[21:1:24,25:0.2:28.2])
set(gca, 'ytick',[23.5:0.5:24.5])
set(gca,'XTicklabel',[])
% ylabel('Nuetral Density [kg.m^{-3}]')
% xlabel('Density Anom. Temperature  [kg.m^{-3}]')
title('(a) CT decomposition')
grid on

set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
set(gca,'FontSize',12)
% calculate and add the errorbar
for i=1:length(neutralgrid)
    T_n_error(i)=1.697.*std(squeeze(alpha_temp_n(i,i_lat)),'omitnan')./sqrt(30);
    T_z_error(i)=1.697.*std(squeeze(T_z_alpha(i,i_lat)),'omitnan')./sqrt(30);
    T_heave_error(i)=1.697.*std(squeeze(-T_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(alpha_temp_n(:,i_lat),2))');
 patch([squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'+T_n_error(ind) flip(squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'-T_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_z_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(T_z_alpha(ind,i_lat),2))'+T_z_error(ind) flip(squeeze(nanmean(T_z_alpha(ind,i_lat),2))'-T_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'+T_heave_error(ind) flip(squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'-T_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')
% legend('\alpha \theta \prime |_{z}',' \alpha \theta \prime | _{n}','-N \prime \alpha \theta _{z}')

%%
subplot('Position', [0.07, 0.33, 0.43, 0.5]);
%subplot('Position', [0.07, 0.67, 0.43, 0.27]);

plot(squeeze(nanmean(T_z_alpha(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(alpha_temp_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)
plot(squeeze(nanmean(-T_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot(repmat(0,351,1),neutralgrid,'k')
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00028,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
text(-0.00028,25.4,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])

xlim([-3e-4 3e-4])
ylim([25 28.2])
% set(gca, 'ytick',[21:1:24,25:0.2:28.2])
set(gca, 'ytick',[25:0.2:28.2])

ylabel('Nuetral Density [kg.m^{-3}]')
% xlabel('Density Anom. Temperature  [kg.m^{-3}]')
% title('(a) CT decomposition')
grid on

set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
set(gca,'FontSize',12)
% calculate and add the errorbar
for i=1:length(neutralgrid)
    T_n_error(i)=1.697.*std(squeeze(alpha_temp_n(i,i_lat)),'omitnan')./sqrt(30);
    T_z_error(i)=1.697.*std(squeeze(T_z_alpha(i,i_lat)),'omitnan')./sqrt(30);
    T_heave_error(i)=1.697.*std(squeeze(-T_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(alpha_temp_n(:,i_lat),2))');
 patch([squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'+T_n_error(ind) flip(squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'-T_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_z_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(T_z_alpha(ind,i_lat),2))'+T_z_error(ind) flip(squeeze(nanmean(T_z_alpha(ind,i_lat),2))'-T_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'+T_heave_error(ind) flip(squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'-T_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')
 legend('\alpha \theta \prime |_{z}',' \alpha \theta \prime | _{n}','-N \prime \alpha \theta _{z}')
%%
subplot('Position', [0.07, 0.08, 0.43, 0.2]);

plot(squeeze(nanmean(T_z_alpha(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(alpha_temp_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)
plot(squeeze(nanmean(-T_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot(repmat(0,351,1),neutralgrid,'k')
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
% plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00008,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
% text(-0.00028,25.2,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00008,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])

xlim([-8e-5 8e-5])
ylim([27.1 nut(end)])
set(gca, 'ytick',[25:0.2:28.2])
ylabel('Nuetral Density [kg.m^{-3}]')
xlabel('Density Anom. Temperature  [kg.m^{-3}]')
% title('(a) CT decomposition')
grid on

set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
set(gca,'FontSize',12)
% calculate and add the errorbar
for i=1:length(neutralgrid)
    T_n_error(i)=1.697.*std(squeeze(alpha_temp_n(i,i_lat)),'omitnan')./sqrt(30);
    T_z_error(i)=1.697.*std(squeeze(T_z_alpha(i,i_lat)),'omitnan')./sqrt(30);
    T_heave_error(i)=1.697.*std(squeeze(-T_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(alpha_temp_n(:,i_lat),2))');
 patch([squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'+T_n_error(ind) flip(squeeze(nanmean(alpha_temp_n(ind,i_lat),2))'-T_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_z_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(T_z_alpha(ind,i_lat),2))'+T_z_error(ind) flip(squeeze(nanmean(T_z_alpha(ind,i_lat),2))'-T_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(T_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'+T_heave_error(ind) flip(squeeze(nanmean(-T_heave_alpha(ind,i_lat),2))'-T_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')

%%
% subplot('Position', [0.55, 0.4, 0.43, 0.55]);
% subplot('Position', [0.55, 0.79, 0.43, 0.16]);
subplot('Position', [0.55, 0.83, 0.43, 0.12]);
plot(squeeze(nanmean(S_z_betta(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(betta_salt_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)

plot(squeeze(nanmean(-S_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00028,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
% text(-0.00028,25.2,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% text(-0.00028,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% plot(neutralgrid,squeeze(nanmean(temp_n(:,i_lat),2))-squeeze(nanmean(T_heave(:,i_lat),2)),'.','LineWidth',1.5)
plot(repmat(0,351,1),neutralgrid,'k')
xlim([-3e-4 3e-4])
ylim([23.5 25])
 set(gca, 'ytick',[23.5:0.5:24.5])
set(gca,'FontSize',12)
% ylabel('Nuetral Density [kg.m^{-3}]')
%set(gca,'ytick',[])
% xlabel('Density Anom. Salinity  [kg.m^{-3}]')
title('(b) SA decomposition')
grid on


% calculate and add the errorbar
for i=1:length(neutralgrid)
    S_n_error(i)=1.697.*std(squeeze(betta_salt_n(i,i_lat)),'omitnan')./sqrt(30);
    S_z_error(i)=1.697.*std(squeeze(S_z_betta(i,i_lat)),'omitnan')./sqrt(30);
    S_heave_error(i)=1.697.*std(squeeze(-S_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(betta_salt_n(:,i_lat),2))');
 patch([squeeze(nanmean(betta_salt_n(ind,i_lat),2))'+S_n_error(ind) flip(squeeze(nanmean(betta_salt_n(ind,i_lat),2))'-S_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_z_betta(:,i_lat),2))');
 patch([squeeze(nanmean(S_z_betta(ind,i_lat),2))'+S_z_error(ind) flip(squeeze(nanmean(S_z_betta(ind,i_lat),2))'-S_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'+S_heave_error(ind) flip(squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'-S_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')

% legend('\beta S \prime |_{z}',' \beta S \prime | _{n}','-N \prime \beta S _{z}')
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')

%%
% subplot('Position', [0.55, 0.4, 0.43, 0.39]);
subplot('Position', [0.55, 0.33, 0.43, 0.5]);
plot(squeeze(nanmean(S_z_betta(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(betta_salt_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)

plot(squeeze(nanmean(-S_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00028,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
text(-0.00028,25.2,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00028,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% plot(neutralgrid,squeeze(nanmean(temp_n(:,i_lat),2))-squeeze(nanmean(T_heave(:,i_lat),2)),'.','LineWidth',1.5)
plot(repmat(0,351,1),neutralgrid,'k')
xlim([-3e-4 3e-4])
ylim([25 28.2])
 set(gca, 'ytick',[25:0.2:28.2])
set(gca,'FontSize',12)
% ylabel('Nuetral Density [kg.m^{-3}]')
%set(gca,'ytick',[])
% xlabel('Density Anom. Salinity  [kg.m^{-3}]')
% title('(b) SA decomposition')
grid on


% calculate and add the errorbar
for i=1:length(neutralgrid)
    S_n_error(i)=1.697.*std(squeeze(betta_salt_n(i,i_lat)),'omitnan')./sqrt(30);
    S_z_error(i)=1.697.*std(squeeze(S_z_betta(i,i_lat)),'omitnan')./sqrt(30);
    S_heave_error(i)=1.697.*std(squeeze(-S_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(betta_salt_n(:,i_lat),2))');
 patch([squeeze(nanmean(betta_salt_n(ind,i_lat),2))'+S_n_error(ind) flip(squeeze(nanmean(betta_salt_n(ind,i_lat),2))'-S_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_z_betta(:,i_lat),2))');
 patch([squeeze(nanmean(S_z_betta(ind,i_lat),2))'+S_z_error(ind) flip(squeeze(nanmean(S_z_betta(ind,i_lat),2))'-S_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'+S_heave_error(ind) flip(squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'-S_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')

legend('\beta S \prime |_{z}',' \beta S \prime | _{n}','-N \prime \beta S _{z}')
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')
%%

subplot('Position', [0.55, 0.08, 0.43, 0.2]);
plot(squeeze(nanmean(S_z_betta(:,i_lat),2)),neutralgrid,'LineStyle','-','Marker','none','MarkerSize',6,'LineWidth',2)
hold on
plot(squeeze(nanmean(betta_salt_n(:,i_lat),2)),neutralgrid,'LineStyle','--','Marker','none','MarkerSize',6,'LineWidth',2)

plot(squeeze(nanmean(-S_heave_alpha(:,i_lat),2)),neutralgrid,'LineStyle',':','Marker','none','MarkerSize',6,'LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_min),7,1),'LineStyle','--','Color','k','LineWidth',2)
plot([-3e4:1e4:3e4],repmat(neutralgrid(index_s_max),7,1),'LineStyle','--','Color','k','LineWidth',2)
% plot([-3e4:1e4:3e4],repmat(25,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(26.3,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(26.7,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
% plot([-3e4:1e4:3e4],repmat(27,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.1,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
plot([-3e4:1e4:3e4],repmat(27.6,7,1),'LineStyle',':','Color','k','LineWidth',1.5)
text(-0.00008,27.25,'S min','FontSize',12)
text(-0.00028,25.7,'S max','FontSize',12)
% text(-0.00028,25.2,'STW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% text(-0.00028,26.9,'SAMW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
text(-0.00008,27.4,'AAIW','FontSize',12,'Rotation',0,'BackgroundColor',[0.6 0.6 0.6])
% plot(neutralgrid,squeeze(nanmean(temp_n(:,i_lat),2))-squeeze(nanmean(T_heave(:,i_lat),2)),'.','LineWidth',1.5)
plot(repmat(0,351,1),neutralgrid,'k')
xlim([-8e-5 8e-5])
ylim([27.1 nut(end)])
set(gca, 'ytick',[25:0.2:28.2])
set(gca,'FontSize',12)
% ylabel('Nuetral Density [kg.m^{-3}]')
%set(gca,'ytick',[])
xlabel('Density Anom. Salinity  [kg.m^{-3}]')
% title('(b) SA decomposition')
grid on


% calculate and add the errorbar
for i=1:length(neutralgrid)
    S_n_error(i)=1.697.*std(squeeze(betta_salt_n(i,i_lat)),'omitnan')./sqrt(30);
    S_z_error(i)=1.697.*std(squeeze(S_z_betta(i,i_lat)),'omitnan')./sqrt(30);
    S_heave_error(i)=1.697.*std(squeeze(-S_heave_alpha(i,i_lat)),'omitnan')./sqrt(30);
end

ind=~isnan(squeeze(nanmean(betta_salt_n(:,i_lat),2))');
 patch([squeeze(nanmean(betta_salt_n(ind,i_lat),2))'+S_n_error(ind) flip(squeeze(nanmean(betta_salt_n(ind,i_lat),2))'-S_n_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.8500 0.3250 0.0980],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_z_betta(:,i_lat),2))');
 patch([squeeze(nanmean(S_z_betta(ind,i_lat),2))'+S_z_error(ind) flip(squeeze(nanmean(S_z_betta(ind,i_lat),2))'-S_z_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0 0.4470 0.7410],'FaceAlpha',0.25, 'EdgeColor','none')
 ind=~isnan(squeeze(nanmean(S_heave_alpha(:,i_lat),2))');
 patch([squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'+S_heave_error(ind) flip(squeeze(nanmean(-S_heave_alpha(ind,i_lat),2))'-S_heave_error(ind))],[neutralgrid(ind) flip(neutralgrid(ind))],[0.9290 0.6940 0.1250],'FaceAlpha',0.25, 'EdgeColor','none')

% legend('\beta S \prime |_{z}',' \beta S \prime | _{n}','-N \prime \beta S _{z}')
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom','TickDir','out')


