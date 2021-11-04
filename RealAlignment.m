clear;

%% Alignment를 읽어 저장하는 프로그램
genenum = 11;
name = cell(genenum,1);
tempiden = zeros(genenum,genenum,2);
tempidenscore = zeros(genenum,genenum,2);
tempgap = zeros(genenum,genenum,2);
tempgapscore = zeros(genenum,genenum,2);
BLASTiden = zeros(genenum,genenum);
BLASTgap = zeros(genenum,genenum);
BLASTidenscore = zeros(genenum,genenum);
BLASTgapscore = zeros(genenum,genenum);
DistBl = zeros(genenum,genenum);

title = ('MT10_out');
str = sprintf('%s_BLAST.txt',title);
file = fopen(str, 'r');

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

file = fopen(str, 'r');

while Str~=-1
    Str = fgets(file);
    if(strfind(Str,'Query='))
        temp = strsplit(Str);
        i = find(strcmp(name,temp{2}));
        disp(Str);
    end
    if(strfind(Str,'>Query_'))
        temp = strsplit(Str);
        j = find(strcmp(name,temp{2}));
        disp(Str);
    end
    if(strfind(Str,'Score'))
        temp = strsplit(Str);
        tempscore = str2num(temp{4});
    end
    if(strfind(Str,'Identities'))
        temp = strsplit(Str);
        temp2 = strsplit(temp{4},'/');
        tempiden(i,j,1) = tempiden(i,j,1) + str2num(temp2{1});
        tempiden(i,j,2) = tempiden(i,j,2) + str2num(temp2{2});
        tempidenscore(i,j,1) = tempidenscore(i,j,1) + tempscore*str2num(temp2{1});
        tempidenscore(i,j,2) = tempidenscore(i,j,2) + tempscore*str2num(temp2{2});
        temp2 = strsplit(temp{8},'/');
        tempgap(i,j,1) = tempgap(i,j,1) + str2num(temp2{1});
        tempgap(i,j,2) = tempgap(i,j,2) + str2num(temp2{2});
        tempgapscore(i,j,1) = tempgapscore(i,j,1) + tempscore*str2num(temp2{1});
        tempgapscore(i,j,2) = tempgapscore(i,j,2) + tempscore*str2num(temp2{2});
        disp(Str);
    end
end
for i = 1:genenum
    for j = 1:genenum
        if tempiden(i,j,2) ~= 0
            BLASTiden(i,j) = tempiden(i,j,1)./tempiden(i,j,2);
            BLASTidenscore(i,j) = tempidenscore(i,j,1)./tempidenscore(i,j,2);
            BLASTgap(i,j) = tempgap(i,j,1)./tempgap(i,j,2);
            BLASTgapscore(i,j) = tempgapscore(i,j,1)./tempgapscore(i,j,2);
            DistBl(i,j) = -3/4*log(1-4*(1-BLASTiden(i,j))/3);
        end
    end
end
fclose(file);
% save('data3.mat');

% %% REminer 읽기
% Str = sprintf('%s\\Forward\\0.0_5000000.5000000\\0_1000000',title);
% Str2 = sprintf('D:\\Download\\program_reminer_용준\\%s.txt',title);
% file = fopen(Str2, 'r');
% raw = fscanf(file,'%s');
% fclose(file);
% 
% genenum = 11;
% seqdiff = 100;
% count = 1;
% state = 0;
% startend = zeros(genenum,2);
% startend(1,1) = 1;
% startend(genenum,2) = size(raw,2);
% 
% for i = 1:size(raw,2)
%     if raw(i) == 'N' && state == 0
%         state = 1;
%         startend(count,2) = i-1;
%         count = count+1;
%     elseif raw(i) ~= 'N' && state == 1
%         state = 0;
%         startend(count,1) = i;
%     end
% end
% 
% list = ls(Str);
% seq1 = zeros(size(list,1)-2,2);
% seq2 = zeros(size(list,1)-2,2);
% score = zeros(size(list,1)-2,1);
% REiden = zeros(size(list,1)-2,3);
% gap = zeros(size(list,1)-2,3);
% 
% for i = 3:size(list,1)
%     filename = sprintf('%s\\%s',Str,list(i,:));
%     file = fopen(filename,'r');
%     temp = fgetl(file);
%     seq1(i-2,:) = sscanf(temp, 'Seq1 : %d - %d');
%     temp = fgetl(file);
%     seq2(i-2,:) = sscanf(temp, 'Seq2 : %d - %d');
%     temp = fgetl(file);
%     score(i-2,:) = sscanf(temp, 'Score : %d');
%     temp = fgetl(file);
%     REiden(i-2,:) = sscanf(temp, 'Identity : %d / %d (%d %%)');
%     temp = fgetl(file);
%     gap(i-2,:) = sscanf(temp, 'Gap : %d / %d (%d %%)');
%     fclose(file);
% end
% disp('Data Reading is completed');
% Str = sprintf('%d alignment file is detected', size(list,1)-2);
% 
% seq = sortrows([seq1 seq2 REiden(:,1:2) score(:,1)]);
% result = zeros(genenum,genenum,8);
% distance = zeros(genenum,genenum);
% distance2 = zeros(genenum,genenum);
% distance3 = zeros(genenum,genenum);
% distance4 = zeros(genenum,genenum);
% 
% for i = 1:genenum
%     for j = 1:genenum
%         temp = (startend(i,1)-seqdiff<=seq(:,1) & seq(:,2)<=startend(i,2)+seqdiff) & (startend(j,1)-seqdiff<=seq(:,3) & seq(:,4)<=startend(j,2)+seqdiff);
%         result(i,j,1) = startend(i,1);
%         result(i,j,2) = startend(i,2);
%         result(i,j,3) = startend(j,1);
%         result(i,j,4) = startend(j,2);
%         result(i,j,5) = sum(seq(temp,5));
%         result(i,j,6) = sum(seq(temp,6));
%         result(i,j,7) = sum(seq(temp,7));
%         num(i,j) = sum(temp);
%         if result(i,j,6) ~= 0
%             distance(i,j) = result(i,j,5)/result(i,j,6);
%             result(i,j,8) = sum(seq(temp,7).*seq(temp,6))/sum(seq(temp,6));
%         end
%     end
% end
% 
% for i = 2:genenum
%     for j = 1:i-1
%         distance2(i,j) = sum((result(i,:,6) - result(j,:,6)).^2);
%         distance3(i,j) = sum((result(i,:,7) - result(j,:,7)).^2);
%         distance4(i,j) = sum((result(i,:,8) - result(j,:,8)).^2);
%     end
% end
% disp('Data summation is completed');
% 
% Str = sprintf('%s.mat',title);
% save(Str, 'BLASTiden', 'distance', 'distance2', 'distance3', 'distance4', 'result');

%     Str = sprintf('java -jar D:\\Download\\MATLAB\\phylonet_v2_4\\phylonet_v2_4\\phylonet_v2_4.jar rf -m D:\\Download\\MATLAB\\REMinerPattern\\Data\\ReferenceTree\\merge_all_file59.dnd -e D:\\Download\\MATLAB\\REMinerPattern\\Data\\Alignmentbased\\MTgenome59\\%s.dnd -o D:\\Download\\MATLAB\\REMinerPattern\\Data\\Alignmentbased\\MTgenome59\\%s.txt',title,title);
%     [status, result] = dos(Str);
%     disp(Str);
%     disp('Phylogenetic tree is constructed');

real(7) = 0.140342   + 0.0594382   + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(8) = 0.125684   + 0.0594382   + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(6) =              0.134914    + 0.0132982   + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;
real(5) =                            0.100796    + 0.0252885 + 0.00369163    + 0.027436  + 0.0894384;

real(9) = 0.0167535  + 0.027723                  + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;
real(10)= 0.0232319  + 0.027723                  + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;
real(11)=              0.0432739                 + 0.010615  + 0.00369163    + 0.027436  + 0.0894384;

real(4) =                                                      0.0546632     + 0.027436  + 0.0894384;
real(3) =                                                                      0.0682325 + 0.0894384;
real(2) =                                                                                  0.0894384;
real(1) =                                                                                  0;

x = -2:0.001:2;
for i = 1:4001
    wBLASTiden = BLASTiden+x(i)*BLASTgap;
    MSE(i) = sum((real(2:11)+3/4*log(1-4*(1-wBLASTiden(1,2:11))/3)).^2);
    wBLASTidenscore = BLASTidenscore+x(i)*BLASTgapscore;
    MSEscore(i) = sum((real(2:11)+3/4*log(1-4*(1-wBLASTidenscore(1,2:11))/3)).^2);
end
[MMSE, index] = min(MSE);
[MMSEscore, indexscore] = min(MSEscore);
wBLASTiden = BLASTiden+x(index)*BLASTgap;
wBLASTidenscore = BLASTidenscore+x(indexscore)*BLASTgapscore;

% x = -2:0.001:2;
% for i = 1:4001
%     wdistance = distance+x(i)*BLASTgap;
%     MSE(i) = sum((real(2:11)+3/4*log(1-4*(1-wBLASTiden(1,2:11))/3)).^2);
% end
% [MMSE, index] = min(MSE);
% wBLASTiden = BLASTiden+x(index)*BLASTgap;



%% plotting
figure
plot(real(2:11));
hold on
plot(-3/4*log(1-4*(1-BLASTiden(1,2:11))/3));
% plot(0.1:0.1:1,-3/4*log(1-4*(wBLASTgap(1,2:11))/3));
% plot(-3/4*log(1-4*(1-wBLASTiden(1,2:11))/3));
plot(-3/4*log(1-4*(1-BLASTidenscore(1,2:11))/3));
plot(-3/4*log(1-4*(1-wBLASTidenscore(1,2:11))/3));
% plot(-3/4*log(1-4*(1-distance(2:11,1))/3));
ylim([0 0.5]);
xlabel('Taxon');
ylabel('Evolutionary distance (d)');
legend('Evolutionary distance (Input value)', 'Evolutionary distance (BLAST)', 'Evolutionary distance (Weighted by score)', 'Evolutionary distance (Weighted by score + gaps)');
% legend('Evolutionary distance (Input value)', 'Evolutionary distance (BLAST)', 'Evolutionary distance (wBLAST)', 'Evolutionary distance (REminer)');

figure
plot(0.1:0.1:1,real(2:11)+3/4*log(1-4*(1-BLASTiden(1,2:11))/3));
hold on
% plot(0.1:0.1:1,real(2:11)+3/4*log(1-4*(1-wBLASTiden(1,2:11))/3));
plot(0.1:0.1:1,real(2:11)+3/4*log(1-4*(1-BLASTidenscore(1,2:11))/3));
plot(0.1:0.1:1,real(2:11)+3/4*log(1-4*(1-wBLASTidenscore(1,2:11))/3));
plot(0.1:0.1:1,zeros(1,10));
% plot(0.1:0.1:1,real(2:11)+3/4*log(1-4*(1-distance(2:11,1))/3)');
xlabel('Time (t)');
ylabel('Difference (d)');
legend('Difference (BLAST)', 'Difference (Weighted by score)', 'Difference (Weighted by score + gaps)');
% legend('Difference (BLAST)', 'Difference (wBLAST)', 'Difference (REminer)');

BLASTphy = BLASTiden(2:genenum,2:genenum);
BLASTphyscore = BLASTidenscore(2:genenum,2:genenum);

for i = 1:genenum-1
    BLASTphy(i,i) = 1;
    BLASTphyscore(i,i) = 1;
end
REMtree = seqlinkage(squareform(1-BLASTphy,'tovector'), 'average', name(2:end));
Str = sprintf('Data\\PhylogeneticTree\\MT10.dnd');
phytreewrite(Str,REMtree,'BRANCHNAMES',false);
phytreeviewer(REMtree);

REMtree = seqlinkage(squareform(1-BLASTphyscore,'tovector'), 'average', name(2:end));
Str = sprintf('Data\\PhylogeneticTree\\MT10_score.dnd');
phytreewrite(Str,REMtree,'BRANCHNAMES',false);
phytreeviewer(REMtree);