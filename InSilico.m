%% Reference sequence 1

i = 1;
temp = getgenbank(GeneSerial(i));
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');

PDF = [sum(sequence==1),sum(sequence==2),sum(sequence==3),sum(sequence==4)]/sum(sequence~=0);

% ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];

time = 0.01:0.01:0.25;

TransitionProbability = zeros(4,4,25);
for i = 1:25
    TransitionProbability(:,:,i) = expm(ratematrix*time(i));
end

NewSequence = char(zeros(length(sequence),25)');

for i = 1:25
    for j = 1:length(sequence)
        Prob = rand();
        if Prob<TransitionProbability(1,sequence(j),i)
            NewSequence(i,j) = 'T';
        else if Prob<TransitionProbability(1,sequence(j),i)+TransitionProbability(2,sequence(j),i)
                NewSequence(i,j) = 'C';
            else if Prob<TransitionProbability(1,sequence(j),i)+TransitionProbability(2,sequence(j),i)+TransitionProbability(3,sequence(j),i)
                    NewSequence(i,j) = 'A';
                else
                    NewSequence(i,j) = 'G';
                end
            end
        end
    end
end

file = fopen('Simulation.fas','w');

fprintf(file,'>Human\n%s\n',upper(rawsequence));
for i = 1:25
    fprintf(file,'>%d\n%s\n',i,NewSequence(i,:));
end
fclose(file);

file = fopen('Simulation.txt','w');

fprintf(file,'%s',upper(rawsequence));
blank = char(ones(1,5000)*'N');
for i = 1:25
    fprintf(file,'%s',blank);
    fprintf(file,'%s',i,NewSequence(i,:));
end
fclose(file);

%% logm을 통해 t 계산

TransSim = zeros(4,4,25);
RateSim = zeros(4,4,25);
timeSim = zeros(25,1);
ParSim = 0:0.0001:1;
ParSim2 = ones(1,10001);
OptSim = zeros(4,4,10001);

for i = 1:25
    NewSequenceT = NewSequence(i,find(rawsequence=='t'));
    TransSim(1,1,i) = sum(NewSequenceT=='T')/sum(NewSequenceT~=0);
    TransSim(2,1,i) = sum(NewSequenceT=='C')/sum(NewSequenceT~=0);
    TransSim(3,1,i) = sum(NewSequenceT=='A')/sum(NewSequenceT~=0);
    TransSim(4,1,i) = sum(NewSequenceT=='G')/sum(NewSequenceT~=0);
    NewSequenceC = NewSequence(i,find(rawsequence=='c'));
    TransSim(1,2,i) = sum(NewSequenceC=='T')/sum(NewSequenceC~=0);
    TransSim(2,2,i) = sum(NewSequenceC=='C')/sum(NewSequenceC~=0);
    TransSim(3,2,i) = sum(NewSequenceC=='A')/sum(NewSequenceC~=0);
    TransSim(4,2,i) = sum(NewSequenceC=='G')/sum(NewSequenceC~=0);
    NewSequenceA = NewSequence(i,find(rawsequence=='a'));
    TransSim(1,3,i) = sum(NewSequenceA=='T')/sum(NewSequenceA~=0);
    TransSim(2,3,i) = sum(NewSequenceA=='C')/sum(NewSequenceA~=0);
    TransSim(3,3,i) = sum(NewSequenceA=='A')/sum(NewSequenceA~=0);
    TransSim(4,3,i) = sum(NewSequenceA=='G')/sum(NewSequenceA~=0);
    NewSequenceG = NewSequence(i,find(rawsequence=='g'));
    TransSim(1,4,i) = sum(NewSequenceG=='T')/sum(NewSequenceG~=0);
    TransSim(2,4,i) = sum(NewSequenceG=='C')/sum(NewSequenceG~=0);
    TransSim(3,4,i) = sum(NewSequenceG=='A')/sum(NewSequenceG~=0);
    TransSim(4,4,i) = sum(NewSequenceG=='G')/sum(NewSequenceG~=0);
    
    RateSim(:,:,i) = logm(TransSim(:,:,i));
    OptSim = mean(((reshape(ratematrix,[16,1])) * ParSim - reshape(RateSim(:,:,i) ,[16,1])* ParSim2).^2);
    [temp, index] = min(OptSim);
    timeSim(i) = index/10000;
end

%% Alignment를 읽어 저장하는 프로그램
genenum = 26;
name = cell(genenum,1);
identity = zeros(genenum,genenum);
filename = 'Simulation_BLAST.txt';
file = fopen(filename, 'r');
temp2 = zeros(genenum,genenum);

%% 데이터 읽기
for i = 1:genenum
    while(1)
        Str = fgets(file);
        if(strfind(Str,'Query='))
            disp(i);
            disp(Str);
            temp = strsplit(Str);
            name{i} = temp{2};
            break;
        end
    end
end
fclose(file);

file = fopen(filename, 'r');
for i = 1:genenum
    for j = 1:genenum
        while(1)
            Str = fgets(file);
            if(strfind(Str,'>Query_'))
                temp = strsplit(Str);
                k = find(strcmp(name,temp{2}));
                temp2(i,j)=k;
                break;
            end
        end
        while(1)
            Str = fgets(file);
            if(strfind(Str,'Identities'))
                temp = strsplit(Str);
                identity(i,k) = str2num(temp{4});
                break;
            end
        end
    end
    disp(i);
end
fclose(file);

%% plotting
plot(0.01:0.01:0.25,3/4*(0.01:0.01:0.25));
hold on
plot(0.01:0.01:0.25,3/4*timeSim);
plot(0.01:0.01:0.25,-3/4*log(1-4*(1-identity(1,2:26))/3));
xlabel('Time (t)');
ylabel('Evolutionary distance (d)');
legend('Evolutionary distance (Input value)', 'Evolutionary distance (Calculated)', 'Evolutionary distance (Identity)');