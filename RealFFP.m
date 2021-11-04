%% Gene sequence를 읽어 저장하는 프로그램
clear;

% genenum = 30;
% maxlength = 17245;
% name = {'Human';'Commonchimpanzee';'Pigmychimpanzee';'Gorilla';'Orangutan';'Gibbon';'Baboon';'Horse';'Whiterhinoceros';'Harborseal';'Grayseal';'Cat';'Finwhale';'Bluewhale';'Cow';'Rat';'Mouse';'Opossum';'Wallaroo';'Platypus';'Squirrel';'Fatdormouse';'Guineapig';'Donkey';'Indianrhinoceros';'Dog';'Sheep';'Pig';'Hippopotamus';'Rabbit'};
genenum = 47;
maxlength = 7277;
name = {'B1_I';'B2_I';'I3_I';'NP1_I';'P2_I';'Yam67_I';'C1_I';'C2_I';'C3_I';'C4_I';'Chinahebei_I';'P1_I';'I1_I';'Morocco_I';'T3_I';'M1_II';'HEJA10_III';'JKNSap_III';'JMYHAW_III';'SWUS1_III';'US1_III';'US2_III';'JBOAR1HYO04_III';'JDEERHYO03L_III';'JJTKAN_III';'JMOHYO03L_III';'JRA1_III';'JSOHYO03L_III';'JTHHYO03L_III';'JYOHYO03L_III';'SWJ570_III';'KYRGYZ_III';'ARKELL_III';'HEJA1_IV';'HEJK4_IV';'HEJI4_IV ';'JAKSai_IV';'JKKSAP_IV';'JSMSAP95_IV';'JSNSAPFH_IV';'JSNSAPFH02C_IV';'JTSSAP02_IV';'JYWSAP02_IV';'SWJ131_IV';'SWCH25_IV';'T1_IV';'CCC220_IV';};
locusname = zeros(genenum,15);
length = zeros(genenum);
sequence = zeros(genenum,maxlength);
k = 1;
maxn = 10;
kEvolDist = zeros(genenum,genenum,maxn);

%% 데이터 읽기

for i = 1:genenum
    Data = getgenbank(GeneSerial(i));
    Str = sprintf('>%s\n%s\n',char(name(i,:)),char(Data.Sequence-32));
    length = str2num(Data.LocusSequenceLength);
    sequence(i,1:length) = Data.Sequence;
    Str = sprintf('(%d/%d, %2.2f%%) is loaded', i, genenum, 100*i/genenum);
    disp(Str);
end

%% A=1,T=2,C=3,G=4
rawsequence = char(sequence);
sequenceA = (sequence == 'a');
sequenceT = 2*(sequence == 't');
sequenceC = 3*(sequence == 'c');
sequenceG = 4*(sequence == 'g');

sequence = sequenceA+sequenceT+sequenceC+sequenceG;


for n = [1 2 5 10 20 50 100]
    %% Get WordPdf
    num = 4^k;
    WordPdf = zeros(genenum,num,n);
    temp = zeros(1,k);

    for i = 1:genenum
        for j = 1:maxlength-k+1
            temp = sequence(i,j:j+k-1);
            if temp(:,:) ~= 0
                index = 0;
                for l = 1:k
                    index = index + 4^(l-1)*(temp(1,l)-1);
                end
                index = index+1;
                WordPdf(i,index,floor((j*n-1)/(maxlength-k+1))+1) = WordPdf(i,index,floor((j*n-1)/(maxlength-k+1))+1)+1;
            end
        end
        for j = 1:n
            WordPdf(i,:,j) = WordPdf(i,:,j)/sum(WordPdf(i,:,j));
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

    % ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
    ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];
    T = expm(ratematrix);

    kT = ones(num,num);
    for i = 1:k
        kT = kT.*T(temp2(:,:,i)+1);
    end
    % disp('Equation construction is completed');

    %% Get Evolutionary Distance
    EvolDist = zeros(genenum,genenum);

    tic;
    for i = 1:genenum-1
        tic;
        Str = sprintf('i = %d', i);
        disp(Str);
        for j = i+1:genenum
            tic;
            Str = sprintf('j = %d', j);
            disp(Str);

            %% FFP를 통해 evolutionary distance estimation
            x = 0:0.001:1;

            %% MMSE를 통해서 가장 잘 맞는 x를 찾고 계산
            total = zeros(1001,n);
            total2 = zeros(1001,n);
            p = zeros(n,num);
            q = zeros(n,num);

            for m = 1:1001
                for l = 1:n
                    p(l,:) = WordPdf(i,:,l)*kT^x(m);
                    q(l,:) = WordPdf(j,:,l)*kT^x(m);
                    total(m,l) = total(m,l) + sum((WordPdf(j,:,l)-double(p(l,:))).^2);
                    total2(m,l) = total2(m,l) + sum((WordPdf(i,:,l)-double(q(l,:))).^2);
                end
            end

            [minval,index] = min(total);
            [minval2,index2] = min(total2);
            xSol = max(mean(x(index)), mean(x(index2)));
    %         Str = sprintf('(%d/%d, %2.2f%%) is completed', n, genenum, 100*n/genenum);
    %         disp(Str);

            EvolDist(i,j) = 3/4*xSol;
            toc;
        end
        toc;
    end
    toc;

    %% phylogenetic tree 그리기
    Z = linkage(EvolDist,'average','euclidean');
    Phytree = phytree(Z,name);
    % DistanceMatrix = pdist(result, 'euclidean');
    % Repeattree = seqneighjoin(DistanceMatrix, 'equivar', name);
    Str = sprintf('Simulation_%d.dnd', n);
    phytreewrite(Str,Phytree,'BRANCHNAMES',false);
    tr = phytreeread(Str);
    phytreeviewer(tr);
    
    kEvolDist(:,:,n) = EvolDist;
end

save('data2.mat');