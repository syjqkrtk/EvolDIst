clear;

% filename = 'simple_nuc_out.ma';
% filename = 'complex_nuc_out.ma';
filename = 'MT10_out.ma';

file = fopen(filename);
text = fscanf(file,'%s');
fclose(file);
% sim = 10;
% genenum = 11;
sim = 1;
genenum = 11;

taxons = strsplit(text,'>');

result = cell(sim,genenum);
count = 1;

% for i = 1:sim
%     for j = 1:genenum
%         count = count + 1;
%         if j == 1
%             result{i,genenum+1-j} = taxons{count}(8:end);
%         else
%             if j == genenum
%                 result{i,genenum+1-j} = taxons{count}(5:end);                
%             else
%                 result{i,genenum+1-j} = taxons{count}(7:end);
%             end
%         end
%     end
% end
for i = 1:sim
    for j = 1:genenum
        count = count + 1;
        count2 = 0;
        flag = 0;
        while(1)
            count2 = count2+1;
            if taxons{count}(count2) < 97
                if flag == 1
                    result{i,j} = taxons{count}(count2:end);
                    break;
                end
                flag = 1;
            end
        end
    end
end

identity = zeros(genenum,genenum,sim);
idenExptGap = zeros(genenum,genenum,sim);
idenGap = zeros(genenum,genenum,sim);
widentity = zeros(genenum,genenum,sim);
EvolDist = zeros(genenum,genenum,sim);
EvolDistExptGap = zeros(genenum,genenum,sim);
EvolDistGap = zeros(genenum,genenum,sim);
wEvolDist = zeros(genenum,genenum,sim);

for i = 1:sim
    for j = 1:genenum
        for k = 1:genenum
            identity(j,k,i) = sum((result{i,j}==result{i,k})&((result{i,j}~='-')|(result{i,k}~='-')))/sum((result{i,j}~='-')|(result{i,k}~='-'));
            idenExptGap(j,k,i) = sum((result{i,j}==result{i,k})&((result{i,j}~='-')&(result{i,k}~='-')))/sum((result{i,j}~='-')&(result{i,k}~='-'));
            idenGap(j,k,i) = sum(((result{i,j}=='-')|(result{i,k}=='-'))&((result{i,j}~='-')|(result{i,k}~='-')))/sum((result{i,j}~='-')|(result{i,k}~='-'));
            EvolDist(j,k,i) = -3/4*log(1-4*(1-identity(j,k,i))/3);
            EvolDistExptGap(j,k,i) = -3/4*log(1-4*(1-idenExptGap(j,k,i))/3);
            EvolDistGap(j,k,i) = -3/4*log(1-4*(idenGap(j,k,i))/3);
        end
    end
end


t = 0:10;

real(7) = 0.140342   + 0.0594382   + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(8) = 0.125684   + 0.0594382   + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(6) =              0.134914    + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(5) =                            0.100796    + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;

real(9) = 0.0167535  + 0.027723                  + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;
real(10)= 0.0232319  + 0.027723                  + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;
real(11)=              0.0432739                 + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;

real(4) =                                                      0.0546632     + 0.027436  + 0.0894384;
real(3) =                                                                      0.0682325 + 0.0894384;
real(2) =                                                                                  0.0894384;
real(1) =                                                                                  0;

x = -2:0.001:2;
for i = 1:4001
    widentity = identity+x(i)*idenGap;
    MSE(i) = sum((real(2:11)+3/4*log(1-4*(1-widentity(1,2:11))/3)).^2);
end
[MMSE, index] = min(MSE);
widentity = identity+x(index)*idenGap;
wEvolDist = -3/4*log(1-4*(1-widentity)/3);

plot(t,real);
% plot(t,t~=0);
hold on
ylim([0 0.5]);
meanEvol = mean(EvolDist,3);
meanEvolExptGap = mean(EvolDistExptGap,3);
meanEvolGap = mean(EvolDistGap,3);
meanwEvol = mean(wEvolDist,3);
plot(t,meanEvol(1,:));
plot(t,meanEvolExptGap(1,:));
% plot(t,meanEvolGap(1,:));
plot(t,meanwEvol(1,:));
legend('Input Value','Evolutionary distance','Evolutionary distance except gap','Evolutionary distance with weighted sum');
title('Guide tree 1');
% title('Guide tree 2');
xlabel('Taxon');
ylabel('Estimated evolutionary distance');