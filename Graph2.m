figure
plot(0.01:0.01:0.25,3/4*(0.01:0.01:0.25));
hold on
plot(0.01:0.01:0.25,-3/4*log(1-4*(1-BLASTiden(1,2:26))/3));
plot(0.01:0.01:0.25,-3/4*log(1-4*(1-distance(2:26,1))/3));
xlabel('Time (t)');
ylabel('Evolutionary distance (d)');
legend('Evolutionary distance (Input value)', 'Evolutionary distance (BLAST)', 'Evolutionary distance (REminer)');

figure
plot(0.01:0.01:0.25,3/4*(0.01:0.01:0.25)+3/4*log(1-4*(1-BLASTiden(1,2:26))/3));
hold on
plot(0.01:0.01:0.25,3/4*(0.01:0.01:0.25)+3/4*log(1-4*(1-distance(2:26,1))/3)');
ylim([-0.01 0.01]);
xlabel('Time (t)');
ylabel('Difference (d)');
legend('Difference (BLAST)', 'Difference (REminer)');