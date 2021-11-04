    %% Reference sequence 1

i = 1;
temp = getgenbank(GeneSerial(i));
rawsequence = temp.Sequence;
sequence = 1*(rawsequence=='t')+2*(rawsequence=='c')+3*(rawsequence=='a')+4*(rawsequence=='g');

PDF = [sum(sequence==1),sum(sequence==2),sum(sequence==3),sum(sequence==4)]/sum(sequence~=0);

ratematrix = [-1.037,0.856,0.173,0.009;0.905,-1.300,0.351,0.044;0.165,0.317,-0.776,0.295;0.019,0.093,0.694,-0.806];
% ratematrix = [-0.75,0.25,0.25,0.25;0.25,-0.75,0.25,0.25;0.25,0.25,-0.75,0.25;0.25,0.25,0.25,-0.75];

time = 0.1:0.1:1;

TransitionProbability = zeros(4,4,25);
for i = 1:10
    TransitionProbability(:,:,i) = expm(ratematrix*time(i));
end

for k = 1:10
    NewSequence = char(zeros(length(sequence),25)');

    for i = 1:10
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

    Str = sprintf('Simulation_%d.fas',k);
    file = fopen(Str,'w');

    fprintf(file,'>Human\n%s\n',upper(rawsequence));
    for i = 1:10
        fprintf(file,'>%d\n%s\n',i,NewSequence(i,:));
    end
    fclose(file);

    Str = sprintf('Simulation_%d.txt',k);
    file = fopen(Str,'w');

    fprintf(file,'%s',upper(rawsequence));
    blank = char(ones(1,5000)*'N');
    for i = 1:10
        fprintf(file,'%s',blank);
        fprintf(file,'%s',i,NewSequence(i,:));
    end
    fclose(file);
end
