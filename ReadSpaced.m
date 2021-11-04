trial = 100;
spaced = zeros(11,11,trial);

for i = 1:trial
    Str = sprintf('spaced\\EvolDist\\Simulation_%d.txt',i);
    file = fopen(Str,'r');
    
    genenum = fscanf(file,'%s',1);
    
    for j = 1:str2num(genenum)
        temp = fscanf(file,'%s',1);
        for k = 1:str2num(genenum)
            temp = fscanf(file,'%s',1);
            spaced(j,k,i) = str2num(temp);
            disp(str2num(temp));
        end
    end
    fclose(file);
end

genenum = str2num(genenum)-1;
save('spaced.mat');