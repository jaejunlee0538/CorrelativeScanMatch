clear all
addpath('tools');
%% intel dataset
load('../data/laser.mat')
T_err = zeros(length(laser)-1, 3);
for i=1:length(laser)-1
    disp(strcat(num2str(i),' vs ', num2str(i+1)))
    scan1 = robotlaser_as_cartesian(laser(1, i));
    scan2 = robotlaser_as_cartesian(laser(1, i+1));

    T1 = v2t(laser(1,i).pose);
    T2 = v2t(laser(1,i+1).pose);
    T12 = inv(T1) * T2;
    T_true = t2v(T12);
    T_init = T_true + [normrnd(0, 0.5,[1,2]), 0.3*rand(1)-0.15]';
    
    
    subplot(2,2,1);
    title(strcat(num2str(i),' vs ', num2str(i+1)),'fontsize',15);
    
%     [T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:),'verbose', figure(1), 'initGuess', T_init);
    [T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:), 'initGuess', T_init);
    
    T_err(i, :) = ominus(T_true, T)';
    % [T, score] = correlativeScanMatchMR(scan1(1:2,:), scan2(1:2,:))
    pause(0.01);
end

