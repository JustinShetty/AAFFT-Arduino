load('vals.mat');
close all

scrsz = get(groot,'ScreenSize');
figure('OuterPosition',[1 0 scrsz(3).*(3./4) scrsz(4)]);

%sp = subplot(2,1,1);

color = [89, 171, 227]./256;
%histogram(vals10,32,'Normalization','probability','FaceColor',color);

[N1,edges] = histcounts(vals10,32,'Normalization','probability');
hold on
xh = [0 2.^15]; yh = [mean(N1) mean(N1)];
%plot(xh,yh,'--k');

title('Pseudorandom integer generation on [0, 2\^15)');
legend('10 runs (150 points)');
xlabel('Value');
ylabel('Relative Frequency');
ylim([0 0.07]);
xlim([0 32768]);
ax = gca;
ax.XAxis.Exponent = 0;
set(ax,'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subplot(2,1,2);

histogram(vals100,32,'Normalization','probability','FaceColor', 'r');

[N2,edges] = histcounts(vals10,32,'Normalization','probability');
hold on
xh = [0 2.^15]; yh = [mean(N2) mean(N2)];
plot(xh,yh,'--k');

legend('1500 samples');
xlabel('Value');
ylabel('Relative Frequency');
ylim([0 0.07]);
xlim([0 32768]);
ax = gca;
ax.XAxis.Exponent = 0;
set(ax,'FontSize',20);
