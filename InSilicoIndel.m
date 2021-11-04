%% super parameter

clear;
r = 0.04;
z = 1.6;
maxindel = 10;

%% Reference sequence 1

i = 1;
temp = getgenbank(GeneSerial(i));
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');

PDF = [sum(sequence==1),sum(sequence==2),sum(sequence==3),sum(sequence==4)]/sum(sequence~=0);

ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
% ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];

time = 0.01:0.01:0.25;

TransitionProbability = zeros(4,4,25);
IndelProbability = zeros(1,25);
for i = 1:25
    TransitionProbability(:,:,i) = expm(ratematrix*time(i));
    IndelProbability(i) = 1 - exp(-r*time(i)/2);
end


%% point mutation

NewSequence = cell(1,25);

for i = 1:25
    for j = 1:length(sequence)
        Prob = rand();
        if Prob<TransitionProbability(1,sequence(j),i)
            NewSequence{i}(j) = 'T';
        else if Prob<TransitionProbability(1,sequence(j),i)+TransitionProbability(2,sequence(j),i)
                NewSequence{i}(j) = 'C';
            else if Prob<TransitionProbability(1,sequence(j),i)+TransitionProbability(2,sequence(j),i)+TransitionProbability(3,sequence(j),i)
                    NewSequence{i}(j) = 'A';
                else
                    NewSequence{i}(j) = 'G';
                end
            end
        end
    end
    Str = sprintf('For t = %2.2f, point mutation is completed.', time(i));
    disp(Str);
end

%% insertion

NewSequence2 = cell(1,25);

for i = 1:25
    count = 0;
    for j = 1:length(sequence)
        Prob = rand();
        if Prob<IndelProbability(i)
            n = randraw('zipf',1.6,1);
            Str = sprintf('Point %d, insertion is occured. Insertion length is %d', count, n);
            disp(Str);
            for k = 1:n
                count = count + 1;
                Prob = rand();
                if Prob < PDF(1)
                    NewSequence2{i}(count) = 'T';
                else if Prob < PDF(1) + PDF(2)
                    NewSequence2{i}(count) = 'C';
                    else if Prob < PDF(1) + PDF(2) + PDF(3)
                        NewSequence2{i}(count) = 'A';
                        else
                            NewSequence2{i}(count) = 'G';
                        end
                    end
                end
            end
        end
        count = count + 1;
        NewSequence2{i}(count) = NewSequence{i}(j);
    end
    Prob = rand();
    if Prob<IndelProbability(i)
        n = randraw('zipf',1.6,1);
        for k = 1:n
            count = count + 1;
            Prob = rand();
            if Prob < PDF(1)
                NewSequence2{i}(count) = 'T';
            else if Prob < PDF(1) + PDF(2)
                NewSequence2{i}(count) = 'C';
                else if Prob < PDF(1) + PDF(2) + PDF(3)
                    NewSequence2{i}(count) = 'A';
                    else
                        NewSequence2{i}(count) = 'G';
                    end
                end
            end
        end
    end
    Str = sprintf('For t = %2.2f, insertion is completed. Total length is changed from %d to %d', time(i), length(sequence), length(NewSequence2{i}));
    disp(Str);
end

%% deletion

NewSequence3 = cell(1,25);

for i = 1:25
    count = 0;
    j = 0;
    while j < length(NewSequence2{i})
        Prob = rand();
        if Prob<IndelProbability(i)
            n = randraw('zipf',1.6,1);
            Str = sprintf('Point %d, deletion is occured. Deletion length is %d', j, n);
            j = j+n;
            disp(Str);
        end
        count = count + 1;
        j = j+1;
        if j > length(NewSequence2{i})
            break;
        end
        NewSequence3{i}(count) = NewSequence2{i}(j);
    end
    Str = sprintf('For t = %2.2f, deletion is completed. Total length is changed from %d to %d', time(i), length(NewSequence2{i}), length(NewSequence3{i}));
    disp(Str);
end

%% save file

% Str = sprintf('Simulation_%d.fas',k);
Str = sprintf('Simulation_indel.fas',k);
file = fopen(Str,'w');

fprintf(file,'>Human\n%s\n',upper(rawsequence));
for i = 1:25
    fprintf(file,'>%d\n%s\n',i,NewSequence3{i});
end
fclose(file);

% Str = sprintf('Simulation_%d.txt',k);
Str = sprintf('Simulation_indel.txt',k);
file = fopen(Str,'w');

fprintf(file,'%s',upper(rawsequence));
blank = char(ones(1,5000)*'N');
for i = 1:25
    fprintf(file,'%s',blank);
    fprintf(file,'%s',i,NewSequence3{i});
end
fclose(file);