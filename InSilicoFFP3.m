%% Transition probability matrix를 역산해 구하는 방법

clear

maxk = 1;
maxd = 10;
trial = 100;
genenum = 10;

kEucDist = zeros(maxk,trial,genenum,maxd);
kEvolDist = zeros(maxk,trial,genenum,maxd);
kRawPdf = cell(maxk,trial,maxd);
kNewPdf = cell(maxk,trial,maxd);

%% Reference sequence

i = 1;
temp = getgenbank('V00662');
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');
seqlen = length(sequence);
numsequence = sequence;

% ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];

time = 0.1:0.1:0.1*genenum;
% time = 0.01:0.01:0.01*genenum;

TransitionProbability = zeros(4,4,genenum);
for i = 1:genenum
    TransitionProbability(:,:,i) = expm(ratematrix*time(i));
end


for tr = 1:trial
    %% Sequence generation
    Str = sprintf('%d trial', tr);
    disp(Str);
    NewSequence = char(zeros(seqlen,genenum)');

    for i = 1:genenum
        for j = 1:seqlen
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
    numSequence = (NewSequence=='T')+2*(NewSequence=='C')+3*(NewSequence=='A')+4*(NewSequence=='G');

    Str = sprintf('Simulation_%d.fas',tr);
    file = fopen(Str,'w');
    fprintf(file,'>%d\n',0);
    fprintf(file,'%s\n',upper(rawsequence));
    for i = 1: genenum
        fprintf(file,'>%d\n',i);
        fprintf(file,'%s\n',NewSequence(i,:));
    end
    fclose(file);
    
    %% Evolutionary distance estimation
    for k = 1:maxk
        tic;
        Str = sprintf('k = %d', k);
        disp(Str);
        num = 4^k;
        for d = [1 2 5 10 20 50 100]
            tic;
            Str = sprintf('d = %d', d);
            disp(Str);
            RawPdf = zeros(1,num,d);
            NewPdf = zeros(genenum,num,d);
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
                    RawPdf(1,index,floor((j*d-1)/(seqlen-k+1))+1) = RawPdf(1,index,floor((j*d-1)/(seqlen-k+1))+1)+1;
                end
            end
            for j = 1:d
                RawPdf(1,:,j) = RawPdf(1,:,j)/sum(RawPdf(1,:,j));
            end

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
                        NewPdf(i,index,floor((j*d-1)/(seqlen-k+1))+1) = NewPdf(i,index,floor((j*d-1)/(seqlen-k+1))+1)+1;
                    end
                end
                temp = zeros(1,k);
    %             Str = sprintf('(%d/%d, %2.2f%%) is completed', i, genenum, 100*i/genenum);
    %             disp(Str);
            end
            for i = 1:genenum
                for j = 1:d
                    NewPdf(i,:,j) = NewPdf(i,:,j)/sum(NewPdf(i,:,j));
                end
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

%             ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
            ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];
            T = expm(ratematrix);

            kT = ones(num,num);
            for i = 1:k
                kT = kT.*T(temp2(:,:,i)+1);
            end
%             disp('Equation construction is completed');


            %% FFP를 통해 evolutionary distance estimation
            EvolDist = zeros(1,genenum);
            xSol = zeros(1,genenum);
            p = zeros(d,num);
            x = 0:0.001:1;
            
            for n = 1:genenum
                %% MMSE를 통해서 가장 잘 맞는 x를 찾고 계산
                total = zeros(1001,d);
                for m = 1:1001
                    for l = 1:d
                        p(l,:) = RawPdf(1,:,l)*kT^x(m);
                        total(m,l) = total(m,l) + sum((NewPdf(n,:,l)-double(p(l,:))).^2);
                    end
                end
                [minval, index] = min(total);
                xSol(1,n) = mean(x(index));
    %             Str = sprintf('(%d/%d, %2.2f%%) is completed', n, genenum, 100*n/genenum);
    %             disp(Str);
            end
            EvolDist = 3/4*xSol;
            kRawPdf{k,tr,d} = RawPdf;
            kNewPdf{k,tr,d} = NewPdf;
            kEvolDist(k,tr,:,d) = EvolDist;
            toc;
        end
        toc;
    end
end

save('data.mat');