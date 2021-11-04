%% Transition probability matrix를 역산해 구하는 방법

clear

maxk = 5;
trial = 50;
genenum = 25;

kEucDist = zeros(maxk,trial,genenum);
kEvolDist = zeros(maxk,trial,genenum);
kRawPdf = cell(maxk,trial,1);
kNewPdf = cell(maxk,trial,1);

%% Reference sequence

i = 1;
temp = getgenbank(GeneSerial(i));
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');
seqlen = length(sequence);
numsequence = sequence;

% ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];

% time = 0.1:0.1:0.1*genenum;
time = 0.01:0.01:0.01*genenum;

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

    
    %% Evolutionary distance estimation
    for k = 1:maxk
        tic;
        Str = sprintf('k = %d', k);
        disp(Str);
        num = 4^k;
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
%             Str = sprintf('(%d/%d, %2.2f%%) is completed', i, genenum, 100*i/genenum);
%             disp(Str);
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

        syms x T kT(x) p(x);
        T = [0.25+0.75*x, 0.25-0.25*x, 0.25-0.25*x, 0.25-0.25*x;0.25-0.25*x, 0.25+0.75*x, 0.25-0.25*x, 0.25-0.25*x;0.25-0.25*x, 0.25-0.25*x, 0.25+0.75*x, 0.25-0.25*x;0.25-0.25*x, 0.25-0.25*x, 0.25-0.25*x, 0.25+0.75*x];
        kT(x) = T(temp2(:,:,1)+1);
        for i = 2:k
            kT(x) = kT(x).*T(temp2(:,:,i)+1);
        end
%         disp('Equation construction is completed');


        %% FFP를 통해 evolutionary distance estimation
        EvolDist = zeros(1,genenum);
        xSol = zeros(1,genenum);
        p(x) = RawPdf*kT(x);
        x = 0:0.01:1;
        
        for n = 1:genenum
            %% MMSE를 통해서 가장 잘 맞는 x를 찾고 계산
            total = zeros(1,101);
            for m = 1:101
                total(m) = total(m) + sum((NewPdf(n,:)-double(p(x(m)))).^2);
            end
            [minval, index] = min(total);
            xSol(1,n) = x(index);
%             Str = sprintf('(%d/%d, %2.2f%%) is completed', n, genenum, 100*n/genenum);
%             disp(Str);
        end
        EvolDist = -3/4*log(xSol);
        kRawPdf{k,tr} = RawPdf;
        kNewPdf{k,tr} = NewPdf;
        kEvolDist(k,tr,:) = EvolDist;
        toc;
    end
end

t = zeros(1,1,genenum);
% t(:,:,:) = 0.01:0.01:0.01*genenum;
t(:,:,:) = 0.1:0.1:0.1*genenum;
Diff = kEvolDist-3/4*repmat(t,maxk,trial);
meanPdf = mean(Diff,2);
stdPdf = std(Diff,0,2);
ratioPdf = kEvolDist./(3/4*repmat(t,maxk,trial));
meanRatio = mean(ratioPdf,2);
stdRatio = std(ratioPdf,0,2);
Str = sprintf('%d_%d_%d.mat',maxk, trial, genenum);
save(Str);