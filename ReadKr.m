genenum = 11;
trial = 100;

for i = 1:trial
    Str = sprintf('Simulation_%d.fas',i);
    file = fopen(Str,'r');
    
    while(~feof(file))
        count = count + 1;
        temp = fgetl(file);
        [temp3,temp2] = strtok(temp,',');
        temp4 = strfind(temp,'_GBL');
        if temp4
            prot{i,count} = temp3(1:temp4-1);
        else
            prot{i,count} = temp3;
        end
        data(i,count) = str2num(temp2);
    end
end