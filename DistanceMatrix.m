load('data2.mat');
load('data3.mat');
kEvolDist = kEvolDist+permute(kEvolDist,[2 1 3]);

DistCo = zeros(genenum,genenum);
file = fopen('Cophylog\Real\Hepatitis_47.txt');
for i = 1:genenum
    for j = 1:genenum
        DistCo(i,j)=fscanf(file,'%f',1);
    end
end
fclose(file);

DistKr = zeros(genenum,genenum);
file = fopen('Kr\Real\Hepatitis_47.txt');
temp=fscanf(file,'%s',1);
for i = 1:genenum
    name=fscanf(file,'%s',1);
    for j = 1:genenum
        DistKr(i,j)=fscanf(file,'%f',1);
    end
end
fclose(file);

% A = DistCo;
% B = DistBl;
% % meanA = mean(reshape(A,[genenum*genenum 1]));
% meanA = sqrt(sum(sum(A.^2)));
% % meanB = mean(reshape(B,[genenum*genenum 1]));
% meanB = sqrt(sum(sum(B.^2)));
% A = A/meanA*meanB;
% FrobCo = sqrt(trace((A-B)*(A-B)'));
% 
% A = DistKr;
% B = DistBl;
% % meanA = mean(reshape(A,[genenum*genenum 1]));
% meanA = sqrt(sum(sum(A.^2)));
% % meanB = mean(reshape(B,[genenum*genenum 1]));
% meanB = sqrt(sum(sum(B.^2)));
% A = A/meanA*meanB;
% FrobKr = sqrt(trace((A-B)*(A-B)'));
n = [1 2 5 10 20 50 100];
for i = 1:7
    A = kEvolDist(:,:,n(i));
    B = DistBl;
    % meanA = mean(reshape(A,[genenum*genenum 1]));
    meanA = sqrt(sum(sum(A.^2)));
    % meanB = mean(reshape(B,[genenum*genenum 1]));
    meanB = sqrt(sum(sum(B.^2)));
    A = A/meanA*meanB;
%     FrobT(i) = sqrt(trace((A-B)*(A-B)'));
    FrobT(i) = trace(A'*B)/sqrt(trace(A'*A)*trace(B'*B));
end
save('data4.mat');