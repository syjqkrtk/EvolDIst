trial = 100;
genenum = 10;
Mash = zeros(genenum,trial);
meanMash = zeros(genenum,1);
count = zeros(genenum,1);

for i = 1:trial
    Str = sprintf('Mash\\result_%d.txt',i);
    file = fopen(Str,'r');
    
    for j = 1:genenum
        temp1 = fscanf(file,'%s',1);
        temp2 = fscanf(file,'%s',1);
        Mash(j,i) = fscanf(file,'%f',1);
        temp3 = fscanf(file,'%f',1);
        temp4 = fscanf(file,'%s',1);
    end
    fclose(file);
end

for i = 1:trial
    for j = 1:genenum
        if ~isnan(Mash(j,i))
            meanMash(j) = meanMash(j) + min(Mash(j,i),0.75);
            count(j) = count(j) + 1;
        end
    end
end
meanMash = meanMash./count;
% meanMash(1) = 0;
genenum = 10;
save('Mash.mat');