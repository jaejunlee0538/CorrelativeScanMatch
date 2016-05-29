clear all;
hf = figure(1);
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

figure(2);
%%
plot(ax1, 1:100, rand(1,100),'rx');
plot(ax2, 1:100, rand(1,100), 'b.');

axes(ax1);
imagesc(rand(50,50))
