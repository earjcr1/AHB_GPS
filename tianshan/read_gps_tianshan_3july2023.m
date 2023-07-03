clear all
close all
% clc
addpath(['~/Google Drive/utils/'])
load('colormaps.mat')

addpath(['~/Google Drive/Anatolia/Material/Strain/'])
addpath(['~/Google Drive/Anatolia/Material/Strain/surfacevel2strain/'])

load(['~/Google Drive/Anatolia/Material/Data/anatolia_geog.mat'])
load(['~/Google Drive/Anatolia/Material/faults/gem_plotready.gmt'])
load(['~/Google Drive/AHB/Material/int_borders.mat'])

addpath(['~/Google Drive/AHB/Material/lab2/lab2_files/mbin/gps_ec'])
addpath(['~/Google Drive/AHB/Material/lab2/lab2_files/mbin/stats'])

inputregion = 'tianshan';
thresholddist = 6;
tolh = 1e-3;
tolv = 1e-3;
minocctime = 2.5;
minocc = 3;

lonbounds = [63 96];
latbounds = [36 52];

gps.names = [];
gps.lons = [];
gps.lats = [];
gps.veast_itrf = [];
gps.vnorth_itrf = [];
gps.vup = [];
gps.seast = [];
gps.snorth = [];
gps.sup = [];
gps.occtime = [];
gps.nocc = [];
gps.coven = [];
%gps.coveu = [];
%gps.covnu = [];

studies = {...
    'dzwang20',...
    'ge15',...
    'hao19',...
    'kreemer14',...
    'kufner21',...
    'li22',...
    'metzger20',...
    'midas',...
    'pan16_vert',...
    'pan18',...
    'pan19',...
    'pan20',...
    'pan21_vert',...
    'pan22_vert',...
    'perry19',...
    'rui19',...
    'stevens21',...
    'wang17',...
    'wangshen20',...
    'wu22_vert',...
    'xiong21',...
    'zhao15',...
    'zhao23_vert',...
    'zheng17',...
    'zhou16',...
    };

for studyindex = 1:length(studies)
    gpslocal = load([studies{studyindex} '/' studies{studyindex} '.mat']);
    gps.names = [gps.names; gpslocal.(studies{studyindex}).names];
    gps.lons = [gps.lons; gpslocal.(studies{studyindex}).lons];
    gps.lats = [gps.lats; gpslocal.(studies{studyindex}).lats];
    if isfield(gpslocal.(studies{studyindex}),'vnorth_itrf')
        gps.veast_itrf = [gps.veast_itrf; gpslocal.(studies{studyindex}).veast_itrf];
        gps.vnorth_itrf = [gps.vnorth_itrf; gpslocal.(studies{studyindex}).vnorth_itrf];
        gps.seast = [gps.seast; gpslocal.(studies{studyindex}).seast];
        gps.snorth = [gps.snorth; gpslocal.(studies{studyindex}).snorth];
    else
        gps.veast_itrf = [gps.veast_itrf; gpslocal.(studies{studyindex}).lons*NaN];
        gps.vnorth_itrf = [gps.vnorth_itrf; gpslocal.(studies{studyindex}).lons*NaN];
        gps.seast = [gps.seast; gpslocal.(studies{studyindex}).lons*NaN];
        gps.snorth = [gps.snorth; gpslocal.(studies{studyindex}).lons*NaN];
    end
    gps.occtime = [gps.occtime; gpslocal.(studies{studyindex}).occtime];
    gps.nocc = [gps.nocc; gpslocal.(studies{studyindex}).nocc];
    
    if isfield(gpslocal.(studies{studyindex}),'vup')
        gps.vup = [gps.vup; gpslocal.(studies{studyindex}).vup];
        gps.sup = [gps.sup; gpslocal.(studies{studyindex}).sup];
    else
        gps.vup = [gps.vup; gpslocal.(studies{studyindex}).veast_itrf*NaN];
        gps.sup = [gps.sup; gpslocal.(studies{studyindex}).snorth*NaN];
    end

    if isfield(gpslocal.(studies{studyindex}),'coven')
        gps.coven = [gps.coven; gpslocal.(studies{studyindex}).coven];
    elseif isfield(gpslocal.(studies{studyindex}),'veast_itrf')
        gps.coven = [gps.coven; gpslocal.(studies{studyindex}).veast_itrf*NaN];
    else
        gps.coven = [gps.coven; gpslocal.(studies{studyindex}).vup*NaN];
   end

end
load(['kreemer14/kreemer14.mat']);

stationstoremove = {'HAKK','VAN1','VANA','INEB','JUGU','STOL','KTVL','ZONG','KNYA','KNY1','KMSH','H310','CLUJ','QTIT','QTAG','CLTK','0075','STEF','DNC4','DHAR','CMGZ','NKU4'};

gps.lons(strcmp(gps.names,'OREN')) = 63;
gps.sup(strcmp(gps.names,'KIZ1')) = tolv;

gps_c = gps;
kreemer14_c = kreemer14;

gps_c.in_cahb = and(and(gps_c.lons>=lonbounds(1),gps_c.lons<=lonbounds(2)),and(gps_c.lats>=latbounds(1),gps_c.lats<=latbounds(2)));
kreemer14_c.in_cahb = and(and(kreemer14_c.lons>=lonbounds(1),kreemer14_c.lons<=lonbounds(2)),and(kreemer14_c.lats>=latbounds(1),kreemer14_c.lats<=latbounds(2)));

gps_c.horizsig = sqrt(gps_c.seast.^2 + gps_c.snorth.^2);
gps_c.toremove_h = or(or(gps_c.horizsig>tolh,or(gps_c.occtime<minocctime,gps_c.nocc<minocc)),gps_c.in_cahb==0);
gps_c.toremove_v = or(or(gps_c.sup>tolv,or(gps_c.occtime<minocctime,gps_c.nocc<minocc)),gps_c.in_cahb==0);
gps_c.toremove_h(gps_c.horizsig==1e3) = 0;
%return
kreemer14_c.horizsig = sqrt(kreemer14_c.seast.^2 + kreemer14_c.snorth.^2);
kreemer14_c.toremove = or(or(kreemer14_c.horizsig>tolh,or(kreemer14_c.occtime<minocctime,kreemer14_c.nocc<minocc)),kreemer14_c.in_cahb==0);

for stationindex = 1:length(stationstoremove)
    gps_c.toremove_h = or(gps_c.toremove_h,strcmp(gps_c.names,stationstoremove{stationindex}));
    gps_c.toremove_v = or(gps_c.toremove_v,strcmp(gps_c.names,stationstoremove{stationindex}));
    kreemer14_c.toremove = or(kreemer14_c.toremove,strcmp(kreemer14_c.names,stationstoremove{stationindex}));
end

gps_c.veast_itrf(gps_c.toremove_h) = NaN;
gps_c.vnorth_itrf(gps_c.toremove_h) = NaN;
gps_c.vup(gps_c.toremove_v) = NaN;
gps_c.seast(gps_c.toremove_h) = NaN;
gps_c.snorth(gps_c.toremove_h) = NaN;
gps_c.sup(gps_c.toremove_v) = NaN;
gps_c.sup(gps_c.toremove_v) = NaN;

gps_c.toremove = and(gps_c.toremove_h,gps_c.toremove_v);
gps_c.lons(gps_c.toremove) = [];
gps_c.lats(gps_c.toremove) = [];
gps_c.veast_itrf(gps_c.toremove) = [];
gps_c.vnorth_itrf(gps_c.toremove) = [];
gps_c.seast(gps_c.toremove) = [];
gps_c.snorth(gps_c.toremove) = [];
gps_c.vup(gps_c.toremove) = [];
gps_c.sup(gps_c.toremove) = [];
gps_c.occtime(gps_c.toremove) = [];
gps_c.coven(gps_c.toremove) = [];
gps_c.nocc(gps_c.toremove) = [];
gps_c.names(gps_c.toremove) = [];
gps_c.horizsig(gps_c.toremove) = [];

kreemer14_c.names(kreemer14_c.toremove) = [];
kreemer14_c.lons(kreemer14_c.toremove) = [];
kreemer14_c.lats(kreemer14_c.toremove) = [];
kreemer14_c.veast_itrf(kreemer14_c.toremove) = [];
kreemer14_c.vnorth_itrf(kreemer14_c.toremove) = [];
kreemer14_c.seast(kreemer14_c.toremove) = [];
kreemer14_c.snorth(kreemer14_c.toremove) = [];

%return
%%
gps_m = gps_c;

for gpsindex = 1:length(gps_m.lons)
    gps_m.names{gpsindex} = gps_m.names{gpsindex}(1:4);
end

for gpsindex = 1:length(gps_m.lons)
    stadists = sqrt(((gps_m.lons - gps_m.lons(gpsindex)).*111.1.*cosd(gps_m.lats)).^2 + ((gps_m.lats - gps_m.lats(gpsindex))*111.1).^2);
    samestations = strcmp(gps_m.names{gpsindex},gps_m.names);
    samestations = find(and(samestations,stadists<=111.1));

    for powerindex = 10:-1:3
    highestprecloc = find((round(gps_m.lons(samestations)*10^powerindex)/10^powerindex)~=gps_m.lons(samestations));
    if ~isempty(highestprecloc)
        gps_m.lons(samestations) = mean(gps_m.lons(samestations(highestprecloc)));
        gps_m.lats(samestations) = mean(gps_m.lats(samestations(highestprecloc)));
        break
    end
    end
end
%return
gps_m.numcloseby = gps_m.lons*NaN;
gps_m.numcloseby_withvert = gps_m.lons*NaN;

for gpsindex = 1:length(gps_m.lons)
    stadists = sqrt(((gps_m.lons - gps_m.lons(gpsindex)).*111.1.*cosd(gps_m.lats)).^2 + ((gps_m.lats - gps_m.lats(gpsindex))*111.1).^2);
    stationscloseby = stadists<thresholddist;
    gps_m.numcloseby(gpsindex) = sum(stationscloseby);
end

[~,sortorder] = sort(gps_m.numcloseby);
sortorder = flip(sortorder);
gps_m.names = gps_m.names(sortorder);
gps_m.lons = gps_m.lons(sortorder);
gps_m.lats = gps_m.lats(sortorder);
gps_m.veast_itrf = gps_m.veast_itrf(sortorder);
gps_m.vnorth_itrf = gps_m.vnorth_itrf(sortorder);
gps_m.vup = gps_m.vup(sortorder);
gps_m.seast = gps_m.seast(sortorder);
gps_m.snorth = gps_m.snorth(sortorder);
gps_m.sup = gps_m.sup(sortorder);
gps_m.nocc = gps_m.nocc(sortorder);
gps_m.occtime = gps_m.occtime(sortorder);
gps_m.coven = gps_m.coven(sortorder);
gps_m.numcloseby = gps_m.numcloseby(sortorder);
gps_m.numcloseby_withvert = gps_m.numcloseby_withvert(sortorder);

%%
gps_o.names = [];
gps_o.lons = [];
gps_o.lats = [];
gps_o.veast_itrf = [];
gps_o.vnorth_itrf = [];
gps_o.vup = [];
gps_o.seast = [];
gps_o.snorth = [];
gps_o.sup = [];
gps_o.coven = [];

touse = ones(length(gps_m.lons),1);

for gpsindex = 1:length(gps_m.lons)
    if 1==touse(gpsindex)
        stadists = sqrt(((gps_m.lons - gps_m.lons(gpsindex)).*111.1.*cosd(gps_m.lats)).^2 + ((gps_m.lats - gps_m.lats(gpsindex))*111.1).^2);
        stationscloseby = (stadists<thresholddist);
        stationscloseby_withhoriz = and(stationscloseby,and(and(~isnan(gps_m.veast_itrf),~isnan(gps_m.vnorth_itrf)),and(~isnan(gps_m.seast),~isnan(gps_m.snorth))));
        stationscloseby_withhoriz = and(stationscloseby_withhoriz,touse==1);
        stationscloseby_withhoriz_withcoven = and(stationscloseby_withhoriz,~isnan(gps_m.coven));
        stationscloseby_withvert = and(and(stationscloseby,and(~isnan(gps_m.sup),~isnan(gps_m.vup))),touse==1);

        numcloseby_withhoriz = sum(stationscloseby_withhoriz);
        numcloseby_withvert = sum(stationscloseby_withvert);

        gps_o.lons(end+1) = mean(gps_m.lons(stationscloseby));
        gps_o.lats(end+1) = mean(gps_m.lats(stationscloseby));

        gps_o.veast_itrf(end+1) = mean(gps_m.veast_itrf(stationscloseby_withhoriz));
        seast_taylor = sqrt(sum(gps_m.seast(stationscloseby_withhoriz).^2))/numcloseby_withhoriz;
        seast_small = (max(gps_m.veast_itrf(stationscloseby_withhoriz)) - min(gps_m.veast_itrf(stationscloseby_withhoriz)))/(2*sqrt(numcloseby_withhoriz));
        seast_large = std(gps_m.veast_itrf(stationscloseby_withhoriz))/sqrt(numcloseby_withhoriz);
        gps_o.seast(end+1) = max([seast_taylor seast_small seast_large]);

        gps_o.vnorth_itrf(end+1) = mean(gps_m.vnorth_itrf(stationscloseby_withhoriz));
        snorth_taylor = sqrt(sum(gps_m.snorth(stationscloseby_withhoriz).^2))/numcloseby_withhoriz;
        snorth_small = (max(gps_m.vnorth_itrf(stationscloseby_withhoriz)) - min(gps_m.vnorth_itrf(stationscloseby_withhoriz)))/(2*sqrt(numcloseby_withhoriz));
        snorth_large = std(gps_m.vnorth_itrf(stationscloseby_withhoriz))/sqrt(numcloseby_withhoriz);
        gps_o.snorth(end+1) = max([snorth_taylor snorth_small snorth_large]);

        gps_o.coven(end+1) = mean(gps_m.coven(stationscloseby_withhoriz_withcoven));

        gps_o.vup(end+1) = mean(gps_m.vup(stationscloseby_withvert));
        sup_taylor = sqrt(sum(gps_m.sup(stationscloseby_withvert).^2))/numcloseby_withvert;
        sup_small = (max(gps_m.vup(stationscloseby_withvert)) - min(gps_m.vup(stationscloseby_withvert)))/(2*sqrt(numcloseby_withvert));
        sup_large = std(gps_m.vup(stationscloseby_withvert))/sqrt(numcloseby_withvert);
        gps_o.sup(end+1) = max([sup_taylor sup_small sup_large]);

        gps_o.names{end+1} = gps_m.names{gpsindex};
        touse(stationscloseby) = 0;
    end
end

gps_o.names = gps_o.names';
gps_o.lons = gps_o.lons';
gps_o.lats = gps_o.lats';
gps_o.veast_itrf = gps_o.veast_itrf';
gps_o.vnorth_itrf = gps_o.vnorth_itrf';
gps_o.vup = gps_o.vup';
gps_o.seast = gps_o.seast';
gps_o.snorth = gps_o.snorth';
gps_o.sup = gps_o.sup';
gps_o.coven = gps_o.coven';
gps_o.coven(isnan(gps_o.coven)) = 0;

%%
[x,y,z] = wgs2xyz(kreemer14.lons,kreemer14.lats,zeros(length(kreemer14.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

kreemer14.veast_eu = kreemer14.veast_itrf + Ve;
kreemer14.vnorth_eu = kreemer14.vnorth_itrf + Vn;

%%
[x,y,z] = wgs2xyz(gps.lons,gps.lats,zeros(length(gps.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

gps.veast_eu = gps.veast_itrf + Ve;
gps.vnorth_eu = gps.vnorth_itrf + Vn;

%%
[x,y,z] = wgs2xyz(kreemer14_c.lons,kreemer14_c.lats,zeros(length(kreemer14_c.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

kreemer14_c.veast_eu = kreemer14_c.veast_itrf + Ve;
kreemer14_c.vnorth_eu = kreemer14_c.vnorth_itrf + Vn;

%%
[x,y,z] = wgs2xyz(gps_c.lons,gps_c.lats,zeros(length(gps_c.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

gps_c.veast_eu = gps_c.veast_itrf + Ve;
gps_c.vnorth_eu = gps_c.vnorth_itrf + Vn;

%%
[x,y,z] = wgs2xyz(gps_m.lons,gps_m.lats,zeros(length(gps_m.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

gps_m.veast_eu = gps_m.veast_itrf + Ve;
gps_m.vnorth_eu = gps_m.vnorth_itrf + Vn;

%%
[x,y,z] = wgs2xyz(gps_o.lons,gps_o.lats,zeros(length(gps_o.lats),1));
%E = [-99.8800	55.4700	-0.2540];
%OM = eul2rot(E)';
OM = -0.27778*[-0.085 -0.531 0.770]';
[Vxyz,Venu] = rotate(OM,[0*OM 0*OM 0*OM],[x y z]);
Ve = Venu(:,1)*1e-3;
Vn = Venu(:,2)*1e-3;

gps_o.veast_eu = gps_o.veast_itrf + Ve;
gps_o.vnorth_eu = gps_o.vnorth_itrf + Vn;
%return
%%
gps_o.output_h = and(isnan(gps_o.sup),~isnan(gps_o.snorth));
gps_o.output_v = and(~isnan(gps_o.sup),~isnan(gps_o.snorth));
%%
ffile = fopen(['~/Google Drive/AHB/Material/GPS/' inputregion '_tol' num2str(tolh*1e3) '_minocc' num2str(minocctime) '_dist' num2str(thresholddist) '_edges_2D.dat'],'wt');
h_export = find(gps_o.output_h);
for gpsindex = 1:length(h_export)
fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
    gps_o.lons(h_export(gpsindex)),...
    gps_o.lats(h_export(gpsindex)),...
    gps_o.veast_eu(h_export(gpsindex))*1e3,...
    gps_o.vnorth_eu(h_export(gpsindex))*1e3,...
    gps_o.seast(h_export(gpsindex))*1e3,...
    gps_o.snorth(h_export(gpsindex))*1e3,...
    gps_o.coven(h_export(gpsindex)),...
    gps_o.names{h_export(gpsindex)});
end
fclose(ffile);

ffile = fopen(['~/Google Drive/AHB/Material/GPS/' inputregion '_tol' num2str(tolh*1e3) '_minocc' num2str(minocctime) '_dist' num2str(thresholddist) '_edges_2D_fake3D.dat'],'wt');
h_export = find(gps_o.output_h);
for gpsindex = 1:length(h_export)
fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
    gps_o.lons(h_export(gpsindex)),...
    gps_o.lats(h_export(gpsindex)),...
    gps_o.veast_eu(h_export(gpsindex))*1e3,...
    gps_o.vnorth_eu(h_export(gpsindex))*1e3,...
    gps_o.vnorth_eu(h_export(gpsindex))*0,...
    gps_o.seast(h_export(gpsindex))*1e3,...
    gps_o.snorth(h_export(gpsindex))*1e3,...
    gps_o.snorth(h_export(gpsindex))*0 + 10,...
    gps_o.coven(h_export(gpsindex)),...
    gps_o.coven(h_export(gpsindex))*0,...
    gps_o.coven(h_export(gpsindex))*0,...
    gps_o.names{h_export(gpsindex)});
end
fclose(ffile);

ffile = fopen(['~/Google Drive/AHB/Material/GPS/' inputregion '_tol' num2str(tolv*1e3) '_minocc' num2str(minocctime) '_dist' num2str(thresholddist) '_edges_3D.dat'],'wt');
v_export = find(gps_o.output_v);
for gpsindex = 1:length(v_export)
fprintf(ffile,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s\n',...
    gps_o.lons(v_export(gpsindex)),...
    gps_o.lats(v_export(gpsindex)),...
    gps_o.veast_eu(v_export(gpsindex))*1e3,...
    gps_o.vnorth_eu(v_export(gpsindex))*1e3,...
    gps_o.vup(v_export(gpsindex))*1e3,...
    gps_o.seast(v_export(gpsindex))*1e3,...
    gps_o.snorth(v_export(gpsindex))*1e3,...
    gps_o.sup(v_export(gpsindex))*1e3,...
    gps_o.coven(v_export(gpsindex)),...
    0*gps_o.snorth(v_export(gpsindex)),...
    0*gps_o.sup(v_export(gpsindex)),...
    gps_o.names{v_export(gpsindex)});
end
fclose(ffile);

save(['gpsoutput_' inputregion '_' num2str(tolh*1e3) '_' num2str(minocctime) '_' num2str(thresholddist) '_2D_3D_itrfeu.mat'],'gps','gps_c','gps_m','gps_o')
%return
%%
figure(1); clf; hold on; box on;
sc = 3e1;
scatter(gps.lons,gps.lats,75,1e3*gps.vup,'filled','markeredgecolor','none')
plot(allcoastlines(:,1),allcoastlines(:,2),'color',[0 0.5 0.75],'linewidth',0.5)
plot(borderlons,borderlats,'color',[0.5 0.5 0.5],'linewidth',0.5)
plot(gem_plotready(:,1),gem_plotready(:,2),'color',[1 0.75 0.75],'linewidth',0.5);

gpsplot = quiver(gps.lons,gps.lats,sc/cosd(37.5)*gps.veast_eu,sc*gps.vnorth_eu,0,'color','k','linewidth',0.75);
kreemer14plot = quiver(kreemer14.lons,kreemer14.lats,sc/cosd(37.5)*kreemer14.veast_eu,sc*kreemer14.vnorth_eu,0,'color',[0.625 0.625 0.625],'linewidth',0.75);
gpsplot.MaxHeadSize = 1e-4;
kreemer14plot.MaxHeadSize = 1e-4;

%text(gps.lons,gps.lats,gps.names,'fontsize',8)

quiver(64.4,15.5,sc/cosd(37.5)*0.04,0,0,'k','linewidth',0.75)
text(65,16.25,'40 mm/yr','fontsize',15,'horizontalalignment','center')

colorbar;
colormap(jet4_1000);
caxis([-5 5])

xlim([lonbounds])
ylim([latbounds])

set(gca,'DataAspectRatio',[1/cosd(37.5) 1 111.1]);
set(gcf,'position',[0 0 1400 800])
%legend('','','','','Data since GSRM2.1','GSRM2.1 [Kreemer et al., 2014]','fontsize',14)
ax = gca;
ax.FontSize = 15;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + 0.005;
bottom = outerpos(2) + 0.04;
ax_width = outerpos(3);
ax_height = outerpos(4) - 0.1;
ax.Position = [left bottom ax_width ax_height];
%legend('','','','','Data since GSRM2.1','GSRM2.1 (Kreemer et al., 2014)','fontsize',15)
title(['Compiled GNSS velocities in study region'],'fontsize',25)
print(['tianshan_all_Jul3.png'],'-dpng','-r300');
%return
%%
figure(2); clf; hold on; box on;
set(gca,'DataAspectRatio',[1/cosd(37.5) 1 111.1]);
sc = 3e1;
scatter(gps_c.lons,gps_c.lats,75,1e3*gps_c.vup,'filled','markeredgecolor','none')
plot(allcoastlines(:,1),allcoastlines(:,2),'color',[0 0.5 0.75],'linewidth',0.5)
plot(borderlons,borderlats,'color',[0.5 0.5 0.5],'linewidth',0.5)
plot(gem_plotready(:,1),gem_plotready(:,2),'color',[1 0.75 0.75],'linewidth',0.5);

gpsplot = quiver(gps_c.lons,gps_c.lats,sc/cosd(37.5)*gps_c.veast_eu,sc*gps_c.vnorth_eu,0,'color','k','linewidth',0.75);
kreemer14plot = quiver(kreemer14_c.lons,kreemer14_c.lats,sc/cosd(37.5)*kreemer14_c.veast_eu,sc*kreemer14_c.vnorth_eu,0,'color',[0.625 0.625 0.625],'linewidth',0.75);
gpsplot.MaxHeadSize = 1e-2;
kreemer14plot.MaxHeadSize = 1e-2;

quiver(64.4,15.5,sc/cosd(37.5)*0.04,0,0,'k','linewidth',0.75)
text(65,16.25,'40 mm/yr','fontsize',15,'horizontalalignment','center')

colorbar;
colormap(jet4_1000);
caxis([-5 5])

xlim([lonbounds])
ylim([latbounds])

set(gcf,'position',[0 0 1400 800])
ax = gca;
ax.FontSize = 15;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + 0.005;
bottom = outerpos(2) + 0.04;
ax_width = outerpos(3);
ax_height = outerpos(4) - 0.1;
ax.Position = [left bottom ax_width ax_height];
%legend('','','','','Data since GSRM2.1','GSRM2.1 (Kreemer et al., 2014)','fontsize',15)
title(['GNSS data with {\leq}' num2str(tolh*1e3) ' mm/yr uncertainty, {\geq}' num2str(minocctime) ' years occupation time, {\geq}3 occupations'],'fontsize',25)
%return
print(['tianshan_tolh' num2str(tolh*1e3) '_minocctime' num2str(minocctime) '_all_Jul3.png'],'-dpng','-r300');

%%
figure(3); clf; hold on; box on;
sc = 3e1;
scatter(gps_m.lons,gps_m.lats,75,1e3*gps_m.vup,'filled','markeredgecolor','none')
plot(allcoastlines(:,1),allcoastlines(:,2),'color',[0 0.5 0.75],'linewidth',0.5)
plot(borderlons,borderlats,'color',[0.5 0.5 0.5],'linewidth',0.5)
plot(gem_plotready(:,1),gem_plotready(:,2),'color',[1 0.75 0.75],'linewidth',0.5);

%output_plot2D = or(gps_m.output_h,gps_m.output_v);
%gpsplot = quiver(gps_m.lons(gps_m.output_h),gps_m.lats(gps_m.output_h),sc/cosd(37.5)*gps_m.veast_eu(gps_m.output_h),sc*gps_m.vnorth_eu(gps_m.output_h),0,'color','k','linewidth',0.75);
gpsplot = quiver(gps_m.lons,gps_m.lats,sc/cosd(37.5)*gps_m.veast_eu,sc*gps_m.vnorth_eu,0,'color','k','linewidth',0.75);
gpsplot.MaxHeadSize = 1e-2;
%kreemer14plot.MaxHeadSize = 1e-2;

%text(gps_m.lons,gps_m.lats,lower(gps_m.names))

%quiver(64.4,15.5,sc/cosd(37.5)*0.04,0,0,'k','linewidth',0.75)
%text(65,16.25,'40 mm/yr','fontsize',15,'horizontalalignment','center')

colorbar;
colormap(jet4_1000);
caxis([-5 5])

xlim([lonbounds])
ylim([latbounds])

set(gca,'DataAspectRatio',[1/cosd(37.5) 1 111.1]);
set(gcf,'position',[0 0 1400 800])
ax = gca;
ax.FontSize = 15;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + 0.005;
bottom = outerpos(2) + 0.04;
ax_width = outerpos(3);
ax_height = outerpos(4) - 0.1;
ax.Position = [left bottom ax_width ax_height];
title(['GNSS data with {\leq}' num2str(tolh*1e3) ' mm/yr uncertainty, {\geq}' num2str(minocctime) ' years occupation time, {\geq}3 occupations, grouped by name'],'fontsize',22.5)
print(['tianshan_tolh' num2str(tolh*1e3) '_minocctime' num2str(minocctime) '_group_Jul3.png'],'-dpng','-r300');

%%
figure(4); clf; hold on; box on;
sc = 3e1;
scatter(gps_o.lons(gps_o.output_v),gps_o.lats(gps_o.output_v),75,1e3*gps_o.vup(gps_o.output_v),'filled','markeredgecolor','none')
plot(allcoastlines(:,1),allcoastlines(:,2),'color',[0 0.5 0.75],'linewidth',0.5)
plot(borderlons,borderlats,'color',[0.5 0.5 0.5],'linewidth',0.5)
plot(gem_plotready(:,1),gem_plotready(:,2),'color',[1 0.75 0.75],'linewidth',0.5);

output_plot2D = or(gps_o.output_h,gps_o.output_v);
%gpsplot = quiver(gps_o.lons(gps_o.output_h),gps_o.lats(gps_o.output_h),sc/cosd(37.5)*gps_o.veast_eu(gps_o.output_h),sc*gps_o.vnorth_eu(gps_o.output_h),0,'color','k','linewidth',0.75);
gpsplot = quiver(gps_o.lons(output_plot2D),gps_o.lats(output_plot2D),sc/cosd(37.5)*gps_o.veast_eu(output_plot2D),sc*gps_o.vnorth_eu(output_plot2D),0,'color','k','linewidth',0.75);
gpsplot.MaxHeadSize = 1e-2;
%kreemer14plot.MaxHeadSize = 1e-2;

%text(gps_o.lons(gps_o.output_h),gps_o.lats(gps_o.output_h),lower(gps_o.names(gps_o.output_h)))

%quiver(64.4,15.5,sc/cosd(37.5)*0.04,0,0,'k','linewidth',0.75)
%text(65,16.25,'40 mm/yr','fontsize',15,'horizontalalignment','center')

colorbar;
colormap(jet4_1000);
caxis([-5 5])

xlim([lonbounds])
ylim([latbounds])

set(gca,'DataAspectRatio',[1/cosd(37.5) 1 111.1]);
set(gcf,'position',[0 0 1400 800])
ax = gca;
ax.FontSize = 15;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + 0.005;
bottom = outerpos(2) + 0.04;
ax_width = outerpos(3);
ax_height = outerpos(4) - 0.1;
ax.Position = [left bottom ax_width ax_height];
title(['GNSS data with {\leq}' num2str(tolh*1e3) ' mm/yr uncertainty, {\geq}' num2str(minocctime) ' years occupation time, {\geq}3 occupations, averaged within ' num2str(thresholddist) ' km'],'fontsize',22.5)
print(['tianshan_tolh' num2str(tolh*1e3) '_minocctime' num2str(minocctime) '_median' num2str(thresholddist) '_Jul3.png'],'-dpng','-r300');