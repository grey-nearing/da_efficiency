clear all
close all
clc

latlon = load('../data/scan/SCAN_Site_List.txt');
load('data/all_data.mat');
I = find(outdata.Ndata>1200);
lat = latlon(I,2);
lon = latlon(I,3);

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[680,115,1750,1050]);

ax = usamap('conus');
set(ax, 'Visible', 'off')
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
names = {states.Name};
indexHawaii = strcmp('Hawaii',names);
indexAlaska = strcmp('Alaska',names);
indexConus = 1:numel(states);
indexConus(indexHawaii|indexAlaska) = []; 
stateColor = [0.5 1 0.5];

geoshow(ax(1), states(indexConus),  'FaceColor', 0.5*[1,1,1])
%geoshow(ax(2), states(indexAlaska), 'FaceColor', stateColor)
%geoshow(ax(3), states(indexHawaii), 'FaceColor', stateColor)

geoshow(lat,lon,'DisplayType','Point','Marker','s','markeredgecolor','b','markerfacecolor','c');

title('Location of SCAN Sites','fontsize',22);
 
fname = strcat('figures/Figure3_SCAN_Sites');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

