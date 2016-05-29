clear all
addpath('tools');
%% intel dataset
load('../data/laser.mat')
%%
idx1 = 82;
idx2 = 83;

scan1 = robotlaser_as_cartesian(laser(1, idx1));
scan2 = robotlaser_as_cartesian(laser(1, idx2));

T1 = v2t(laser(1,idx1).pose);
T2 = v2t(laser(1,idx2).pose);
T12 = inv(T1) * T2;
T_true = t2v(T12);
T_init = T_true + [normrnd(0, 0.5,[1,2]), 0.3*rand(1)-0.15]';
[T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:),'verbose', figure(1),'initGuess',T_init);
% [T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:))
%%
T_err_init = ominus(T_true, T_init)'
T_err= ominus(T_true, T)'


