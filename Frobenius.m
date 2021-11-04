A = squareform(DistanceMatrix);
B = refdist2;
meanA = mean(reshape(A,[num*num 1]));
stdA = std(reshape(A,[num*num 1]));
meanB = mean(reshape(B,[num*num 1]));
stdB = std(reshape(B,[num*num 1]));
A = (A-meanA)/stdA*stdB+meanB;
result = sqrt(trace((A-B)*(A-B)'));