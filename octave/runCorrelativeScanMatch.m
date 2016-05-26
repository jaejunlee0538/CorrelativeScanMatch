clear all
% addpath('tools');
%% intel dataset
load('../data/laser.mat')
idx1 = 120;
idx2 = 121;

scan1 = robotlaser_as_cartesian(laser(1, idx1));
scan2 = robotlaser_as_cartesian(laser(1, idx2));

T1 = v2t(laser(1,idx1).pose);
T2 = v2t(laser(1,idx2).pose);
T12 = inv(T1) * T2;
T_true = t2v(T12);
%%
tic
[T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:),'verbose', figure(1))
% [T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:))
toc