clear;

%% plotting
load('data.mat');
load('data2.mat');
load('data3.mat');
load('data4.mat');
load('Kr.mat');
load('Cophylog.mat');
load('Mash.mat');

% plot(0.01:0.01:0.01*genenum,3/4*(0.01:0.01:0.01*genenum));
plot(0:0.1:0.1*genenum,(0:0.1:0.1*genenum));
hold on
y = zeros(genenum,1);
% plot(0.01:0.01:0.01*genenum,y);
for i = 1:1
        y(:,1) = mean(kEvolDist(i,:,:,1),2)*4/3;
%         plot(0.01:0.01:0.01*genenum,y-3/4*(0.01:0.01:0.01*genenum)');
%         plot(0.01:0.01:0.01*genenum,y);
        plot(0:0.1:0.1*genenum,[0 y'],'o-');
        y(:,1) = mean(kEvolDist(i,:,:,2),2)*4/3;
        plot(0:0.1:0.1*genenum,[0 y'],'s-.');
        y(:,1) = mean(kEvolDist(i,:,:,5),2)*4/3;
        plot(0:0.1:0.1*genenum,[0 y'],'x:');
        y(:,1) = mean(kEvolDist(i,:,:,10),2)*4/3;
        plot(0:0.1:0.1*genenum,[0 y'],'d--');
end
plot(0:0.1:0.1*genenum,meanCophylog*4/3,'^-');
plot(0:0.1:0.1*genenum,mean(Kr(1,:,:),3)*4/3,'v-');
plot(0:0.1:0.1*genenum,[0; meanMash*4/3],'p-');
xlabel('Input evolutionary distance');
ylabel('Estimated evolutionary distance');
% legend('Evolutionary distance (Input value)', 'Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
% legend('Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
legend({'Evolutionary distance (Input value)','Proposed method, $\hat{t}$ (N=1)', 'Proposed method, $\hat{t}$ (N=2)', 'Proposed method, $\hat{t}$ (N=5)', 'Proposed method, $\hat{t}$ (N=10)', 'Cophylog', 'Kr', 'Mash'},'Interpreter','latex');

figure
% plot(0.01:0.01:0.01*genenum,3/4*(0.01:0.01:0.01*genenum));
plot(0.7:0.1:0.1*genenum,(0.7:0.1:0.1*genenum));
hold on
y = zeros(4,1);
% plot(0.01:0.01:0.01*genenum,y);
for i = 1:1
        y(:,1) = mean(kEvolDist(i,:,7:10,1),2)*4/3;
%         plot(0.01:0.01:0.01*genenum,y-3/4*(0.01:0.01:0.01*genenum)');
%         plot(0.01:0.01:0.01*genenum,y);
        plot(0.7:0.1:0.1*genenum,y','o-');
        y(:,1) = mean(kEvolDist(i,:,7:10,2),2)*4/3;
        plot(0.7:0.1:0.1*genenum,y','s-.');
        y(:,1) = mean(kEvolDist(i,:,7:10,5),2)*4/3;
        plot(0.7:0.1:0.1*genenum,y','x:');
        y(:,1) = mean(kEvolDist(i,:,7:10,10),2)*4/3;
        plot(0.7:0.1:0.1*genenum,y','d--');
end
plot(0.7:0.1:0.1*genenum,meanCophylog(8:11)*4/3,'^-');
% plot(0.7:0.1:0.1*genenum,mean(Kr(1,7:10,:),3)*4/3,'v-');
% xlabel('Input evolutionary distance');
% ylabel('Estimated evolutionary distance');
set(gca,'XTick',[0.7:0.1:1]);
xlim([0.7, 1]);
ylim([0.7, 1]);
% legend('Evolutionary distance (Input value)', 'Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
% legend('Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
% legend({'Evolutionary distance (Input value)','Proposed method, $\hat{t}$ (N=1)', 'Proposed method, $\hat{t}$ (N=2)', 'Proposed method, $\hat{t}$ (N=5)', 'Proposed method, $\hat{t}$ (N=10)'},'Interpreter','latex');

figure
semilogx([1 2 5 10 20 50 100],FrobT,'o-');
xlabel('Number of sequence segment, N');
ylabel('Matrix similarity from BLAST data');


% %% plotting
% load('data.mat');
% load('Kr.mat');
% load('Cophylog.mat');
% t = 0:0.1:0.1*genenum;
% % plot(0.01:0.01:0.01*genenum,3/4*(0.01:0.01:0.01*genenum));
% plot(0:0.1:0.1*genenum,zeros(1,genenum+1));
% hold on
% y = zeros(genenum,1);
% % plot(0.01:0.01:0.01*genenum,y);
% for i = 1:1
%         y(:,1) = mean(kEvolDist(i,:,:,1),2)*4/3;
% %         plot(0.01:0.01:0.01*genenum,y-3/4*(0.01:0.01:0.01*genenum)');
% %         plot(0.01:0.01:0.01*genenum,y);
%         plot(0:0.1:0.1*genenum,[0 y']-t,'o-');
%         y(:,1) = mean(kEvolDist(i,:,:,2),2)*4/3;
%         plot(0:0.1:0.1*genenum,[0 y']-t,'s-.');
%         y(:,1) = mean(kEvolDist(i,:,:,5),2)*4/3;
%         plot(0:0.1:0.1*genenum,[0 y']-t,'x:');
%         y(:,1) = mean(kEvolDist(i,:,:,10),2)*4/3;
%         plot(0:0.1:0.1*genenum,[0 y']-t,'d--');
% end
% plot(0:0.1:0.1*genenum,mean(Cophylog(1,:,:),3)*4/3-t,'^-');
% plot(0:0.1:0.1*genenum,mean(Kr(1,:,:),3)*4/3-t,'v-');
% ylim([-0.01,0.01]);
% xlabel('Input evolutionary distance');
% ylabel('Estimated evolutionary distance');
% % legend('Evolutionary distance (Input value)', 'Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
% % legend('Evolutionary distance (FFP, k=1, d=1)', 'Evolutionary distance (FFP, k=1, d=2)', 'Evolutionary distance (FFP, k=1, d=3)', 'Evolutionary distance (FFP, k=1, d=4)', 'Evolutionary distance (FFP, k=1, d=5)', 'Evolutionary distance (FFP, k=1, d=6)', 'Evolutionary distance (FFP, k=1, d=7)', 'Evolutionary distance (FFP, k=1, d=8)', 'Evolutionary distance (FFP, k=1, d=9)', 'Evolutionary distance (FFP, k=1, d=10)');
% legend({'Evolutionary distance (Input value)','Proposed method, $\hat{t}$ (N=1)', 'Proposed method, $\hat{t}$ (N=2)', 'Proposed method, $\hat{t}$ (N=5)', 'Proposed method, $\hat{t}$ (N=10)', 'Cophylog', 'Kr'},'Interpreter','latex');
