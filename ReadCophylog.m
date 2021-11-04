trial = 100;
genenum = 11;
Cophylog = zeros(11,trial);
meanCophylog = zeros(11,1);
count = zeros(11,1);

for i = 1:trial
    Str = sprintf('Cophylog\\EvolDist\\Simulation_%d.txt',i);
    file = fopen(Str,'r');
    
    for j = 1:genenum
        Cophylog(j,i) = fscanf(file,'%f',1);
    end
    fclose(file);
end

for i = 1:trial
    for j = 1:genenum
        if ~isnan(Cophylog(j,i))
            meanCophylog(j) = meanCophylog(j) + Cophylog(j,i);
            count(j) = count(j) + 1;
        end
    end
end
meanCophylog = meanCophylog./count;
meanCophylog(1) = 0;
genenum = 10;
save('Cophylog.mat');