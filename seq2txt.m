clear;

file = fopen('MT10_out.seq','r');
allseq = fscanf(file,'%s');
fclose(file);

count = 0;

for i = 1:length(allseq)
    if allseq(i)==62
        count = count + 1;
        count2 = 0;
    end
    if allseq(i)>=65 && allseq(i)<=90
        count2 = count2 + 1;
        seq{count}(count2)=allseq(i);
    end
end

for i = 1:count
    seq{i} = seq{i}(2:end);
end

file = fopen('MT10_out.txt','w');

for i = 1:count
    for j = 1:length(seq{i})
        fprintf(file,seq{i}(j));
    end
    for j = 1:5000
        fprintf(file,'N');
    end
end

fclose(file);