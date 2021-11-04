BLASTdist = zeros(10,10);
BLASTdiff = zeros(10,10);
REdist = zeros(10,10);
REdiff = zeros(10,10);

for i = 1:10
    title = sprintf('Simulation_%d.mat',i);
    t = 3/4 * (0.1:0.1:1);
    load(title);
    BLASTdist(i,:) = -3/4*log(1-4*(1-BLASTiden(1,2:11))/3);
    BLASTdiff(i,:) = BLASTdist(i,:) - t;
    REdist(i,:) = -3/4*log(1-4*(1-distance(2:11,1))/3);
    REdiff(i,:) = REdist(i,:) - t;
end

figure
plot(0.1:0.1:1,t);
hold on
plot(0.1:0.1:1,mean(BLASTdist));
plot(0.1:0.1:1,mean(REdist));
xlabel('t');
ylabel('evolutionary distance');
legend('Input value','BLAST','REminer');
figure
plot(0.1:0.1:1,BLASTdiff);
hold on
plot(0.1:0.1:1,REdiff);
ylim([-0.01 0.01]);
xlabel('t');
ylabel('Difference');
legend('BLAST','REminer');
figure
plot(0.1:0.1:1,mean(BLASTdiff));
hold on
plot(0.1:0.1:1,mean(REdiff));
ylim([-0.01 0.01]);
xlabel('t');
ylabel('Mean Difference');
legend('BLAST','REminer');
figure
plot(0.1:0.1:1,std(BLASTdiff));
hold on
plot(0.1:0.1:1,std(REdiff));
ylim([-0.01 0.01]);
xlabel('t');
ylabel('Standard deviation');
legend('BLAST','REminer');