clear

%% Reference sequence 1

i = 1;
temp = getgenbank(GeneSerial(i));
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');

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

%% Parameter 설정
k = 1;
num = 4^k;
genenum = 25;
seqlen = length(rawsequence);
numsequence = (rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');
numSequence = (NewSequence=='T')+2*(NewSequence=='C')+3*(NewSequence=='A')+4*(NewSequence=='G');

kTransProb = ones(num,num,genenum);

for t = 1:genenum
    for i = 1:num
        temp = zeros(1,k);
        word = i-1;
        for l = 1:k
            temp(l) = mod(word,4);
            word = floor(word/4);
        end
        for j = 1:num
            temp2 = zeros(1,k);
            word2 = j-1;
            for l = 1:k
                temp2(l) = mod(word2,4);
                word2 = floor(word2/4);
            end
            for m = 1:k
                kTransProb(i,j,t) = kTransProb(i,j,t)*TransitionProbability(temp(m)+1,temp2(m)+1,t);
            end
        end
    end
end

RawPdf = zeros(1,num);
NewPdf = zeros(genenum,num);
temp = zeros(1,k);

%% word pdf 구하기
for j = 1:seqlen-k+1
    temp = numsequence(1,j:j+k-1);
    if temp(:,:) ~= 0
        index = 0;
        for l = 1:k
            index = index + 4^(l-1)*(temp(1,l)-1);
        end
        index = index+1;
        RawPdf(1,index) = RawPdf(1,index)+1;
    end
end
RawPdf(1,:) = RawPdf(1,:)/sum(RawPdf(1,:));

temp = zeros(1,k);
for i = 1:genenum
    for j = 1:seqlen-k+1
        temp = numSequence(i,j:j+k-1);
        if temp(:,:) ~= 0
            index = 0;
            for l = 1:k
                index = index + 4^(l-1)*(temp(1,l)-1);
            end
            index = index+1;
            NewPdf(i,index) = NewPdf(i,index)+1;
        end
    end
    temp = zeros(1,k);
    Str = sprintf('(%d/%d, %2.2f%%) is completed', i, genenum, 100*i/genenum);
    disp(Str);
end
for i = 1:genenum
    NewPdf(i,:) = NewPdf(i,:)/sum(NewPdf(i,:));
end

%% FFP를 통해 evolutionary distance를 estimation 하기 위한 방정식 전개
temp = zeros(num,k);
temp2 = zeros(num,num,k);
for i = 1:num
    word = i-1;
    for l = 1:k
        temp(i,l) = mod(word,4);
        word = floor(word/4);
    end
end
for i = 1:k
    temp2(:,:,i) = 4*repmat(temp(1:num,i),1,num)+repmat(temp(1:num,i)',num,1);
end
    
tic;
syms x;
T = [0.25+0.75*x, 0.25-0.25*x, 0.25-0.25*x, 0.25-0.25*x;0.25-0.25*x, 0.25+0.75*x, 0.25-0.25*x, 0.25-0.25*x;0.25-0.25*x, 0.25-0.25*x, 0.25+0.75*x, 0.25-0.25*x;0.25-0.25*x, 0.25-0.25*x, 0.25-0.25*x, 0.25+0.75*x];
kT = T(temp2(:,:,1)+1);
for i = 2:k
    kT = kT.*T(temp2(:,:,i)+1);
end
toc;
disp('Equation construction is completed');


%% FFP를 통해 evolutionary distance estimation
EucDist = zeros(1,25);
EvolDist = zeros(1,25);
xSol = zeros(1,25);
p = RawPdf*kT;

for n = 1:genenum
    EucDist(n) = norm(NewPdf(n,:)-RawPdf);

    tic;
    eqn = sum((p-RawPdf).^2)==EucDist(n)^2;
    Sol = double(solve(eqn,x));
    index = find(((0<real(Sol))&(real(Sol)<1)&(imag(Sol)==0)),1);
    xSol(n) = Sol(index);
    EvolDist(n) = -3/4*log(Sol(index));
    toc;
    
    Str = sprintf('(%d/%d, %2.2f%%) is completed', n, genenum, 100*n/genenum);
    disp(Str);
end

%% plotting
plot(0.01:0.01:0.25,3/4*(0.01:0.01:0.25));
hold on
plot(0.01:0.01:0.25,3/4*timeSim);
plot(0.01:0.01:0.25,EvolDist);
xlabel('Time (t)');
ylabel('Evolutionary distance (d)');
legend('Evolutionary distance (Input value)', 'Evolutionary distance (Calculated)', 'Evolutionary distance (FFP)');