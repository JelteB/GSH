clear;
close all;
clc;

addpath('~/TUDelft/Code/GSH/Tools/')

%%
filename = 'megt90n000cb.img';
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
el4 = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

%%

latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

aa = 18;
figure(1)
imagesc(lonT,latT,el4./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

%% GSHA

cs = GSHA(el4,130);
sc = cs2sc(cs);

n = 1:size(sc,1);

D = 200e9*(100e3)^3/(12*(1-0.5^2));
PHI = (1 + (D)/(500*9.81).*(2.*(n+1)./(2*6378000)).^4).^(-1);

sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*PHI';
end

%% GSHS

mapf = GSHS(sc_flex,lonT,90-latT,179);

%%
figure(2)
imagesc(lonT,latT,mapf./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

figure(3)
imagesc(lonT,latT,(el4-mapf)./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)