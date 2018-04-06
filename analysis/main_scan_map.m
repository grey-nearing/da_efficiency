clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd))
%% --- Experiment Parameters ----------------------------------------------

% minimum # of data points to use from each site
Nmin = 1200;

%% --- Load Raw Data ------------------------------------------------------

% load all data
latlon = load('../data/scan/SCAN_Site_List.txt');
load('data/site_data_1.mat');
Ndata = outdata.Ndata;

% check that all 10x10 obs/model cov result in the same number of data points
[Nsites,Nobs,Nmod] = size(Ndata);
for s = 1:Nsites
    if isnan(Ndata(s,1,1)); continue; end
    d = Ndata(s,1,1);
    for o = 1:Nobs
        for m = 1:Nmod
            assert(Ndata(s,o,m) == 0 || Ndata(s,o,m) == Ndata(s,1,1));
        end
    end
end
Ndata = Ndata(:,1,1);

% find sites with enough data
I = find(Ndata>Nmin);
lat = latlon(I,2);
lon = latlon(I,3);

%% --- Make Plot ----------------------------------------------------------

% initialize map
figure(1); close(1); figure(1);
set(gcf,'color','w','position',[680,115,1750,1050]);

ax = usamap('conus');
set(ax,'Visible', 'off','fontsize',36)
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
set(gca,'fontsize',36);
title('Location of SCAN Sites','fontsize',22);

fname = strcat('figures/Figure3_SCAN_Sites');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- END SCRIPT ---------------------------------------------------------