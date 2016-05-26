addpath('../tools');
addpath('..');
more off
% close all
clear all
close all
load('../../data/laser.mat')

for i=1:length(laser)
    figure(1)
    xy = robotlaser_as_cartesian(laser(1, i));
    plot(xy(1,:), xy(2,:),'r.'); 
    xlim([-20 20]); ylim([-20 20]);
    grid on;
    drawnow;
    pause(0.01);
end