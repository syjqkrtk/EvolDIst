for t = 1:10

    %% Alignment를 읽어 저장하는 프로그램
    title = sprintf('Simulation_%d',t);
    genenum = 11;
    name = cell(genenum,1);
    tempiden = zeros(genenum,genenum,2);
    BLASTiden = zeros(genenum,genenum);
   
    filename = sprintf('%s_BLAST.txt',title);
    file = fopen(filename, 'r');

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
        if(strfind(Str,'Identities'))
            temp = strsplit(Str);
            temp2 = strsplit(temp{4},'/');
            tempiden(i,j,1) = tempiden(i,j,1) + str2num(temp2{1});
            tempiden(i,j,2) = tempiden(i,j,2) + str2num(temp2{2});
            disp(Str);
        end
    end
    for i = 1:genenum
        for j = 1:genenum
            if tempiden(i,j,2) ~= 0
                BLASTiden(i,j) = tempiden(i,j,1)./tempiden(i,j,2);
            end
        end
    end
    fclose(file);

    %% REminer 읽기
    Str = sprintf('%s\\Forward\\0.0_5000000.5000000\\0_1000000',title);
    Str2 = sprintf('D:\\Download\\program_reminer_용준\\%s.txt',title);
    file = fopen(Str2, 'r');
    raw = fscanf(file,'%s');
    fclose(file);

    genenum = 11;
    seqdiff = 100;
    count = 1;
    state = 0;
    startend = zeros(genenum,2);
    startend(1,1) = 1;
    startend(genenum,2) = size(raw,2);

    for i = 1:size(raw,2)
        if raw(i) == 'N' && state == 0
            state = 1;
            startend(count,2) = i-1;
            count = count+1;
        elseif raw(i) ~= 'N' && state == 1
            state = 0;
            startend(count,1) = i;
        end
    end

    list = ls(Str);
    seq1 = zeros(size(list,1)-2,2);
    seq2 = zeros(size(list,1)-2,2);
    score = zeros(size(list,1)-2,1);
    REiden = zeros(size(list,1)-2,3);
    gap = zeros(size(list,1)-2,3);

    for i = 3:size(list,1)
        filename = sprintf('%s\\%s',Str,list(i,:));
        file = fopen(filename,'r');
        temp = fgetl(file);
        seq1(i-2,:) = sscanf(temp, 'Seq1 : %d - %d');
        temp = fgetl(file);
        seq2(i-2,:) = sscanf(temp, 'Seq2 : %d - %d');
        temp = fgetl(file);
        score(i-2,:) = sscanf(temp, 'Score : %d');
        temp = fgetl(file);
        REiden(i-2,:) = sscanf(temp, 'Identity : %d / %d (%d %%)');
        temp = fgetl(file);
        gap(i-2,:) = sscanf(temp, 'Gap : %d / %d (%d %%)');
        fclose(file);
    end
    disp('Data Reading is completed');
    Str = sprintf('%d alignment file is detected', size(list,1)-2);

    seq = sortrows([seq1 seq2 REiden(:,1:2) score(:,1)]);
    result = zeros(genenum,genenum,8);
    distance = zeros(genenum,genenum);
    distance2 = zeros(genenum,genenum);
    distance3 = zeros(genenum,genenum);
    distance4 = zeros(genenum,genenum);

    for i = 1:genenum
        for j = 1:genenum
            temp = (startend(i,1)-seqdiff<=seq(:,1) & seq(:,2)<=startend(i,2)+seqdiff) & (startend(j,1)-seqdiff<=seq(:,3) & seq(:,4)<=startend(j,2)+seqdiff);
            result(i,j,1) = startend(i,1);
            result(i,j,2) = startend(i,2);
            result(i,j,3) = startend(j,1);
            result(i,j,4) = startend(j,2);
            result(i,j,5) = sum(seq(temp,5));
            result(i,j,6) = sum(seq(temp,6));
            result(i,j,7) = sum(seq(temp,7));
            num(i,j) = sum(temp);
            if result(i,j,6) ~= 0
                distance(i,j) = result(i,j,5)/result(i,j,6);
                result(i,j,8) = sum(seq(temp,7).*seq(temp,6))/sum(seq(temp,6));
            end
        end
    end

    for i = 2:genenum
        for j = 1:i-1
            distance2(i,j) = sum((result(i,:,6) - result(j,:,6)).^2);
            distance3(i,j) = sum((result(i,:,7) - result(j,:,7)).^2);
            distance4(i,j) = sum((result(i,:,8) - result(j,:,8)).^2);
        end
    end
    disp('Data summation is completed');

    Str = sprintf('%s.mat',title);
    save(Str, 'BLASTiden', 'distance', 'distance2', 'distance3', 'distance4', 'result');

%     Str = sprintf('java -jar D:\\Download\\MATLAB\\phylonet_v2_4\\phylonet_v2_4\\phylonet_v2_4.jar rf -m D:\\Download\\MATLAB\\REMinerPattern\\Data\\ReferenceTree\\merge_all_file59.dnd -e D:\\Download\\MATLAB\\REMinerPattern\\Data\\Alignmentbased\\MTgenome59\\%s.dnd -o D:\\Download\\MATLAB\\REMinerPattern\\Data\\Alignmentbased\\MTgenome59\\%s.txt',title,title);
%     [status, result] = dos(Str);
%     disp(Str);
%     disp('Phylogenetic tree is constructed');

    %% plotting
%     figure
%     plot(0.1:0.1:1,3/4*(0.1:0.1:1));
%     hold on
%     plot(0.1:0.1:1,-3/4*log(1-4*(1-BLASTiden(1,2:11))/3));
%     plot(0.1:0.1:1,-3/4*log(1-4*(1-distance(2:11,1))/3));
%     xlabel('Time (t)');
%     ylabel('Evolutionary distance (d)');
%     legend('Evolutionary distance (Input value)', 'Evolutionary distance (BLAST)', 'Evolutionary distance (REminer)');
% 
%     figure
%     plot(0.1:0.1:1,3/4*(0.1:0.1:1)+3/4*log(1-4*(1-BLASTiden(1,2:11))/3));
%     hold on
%     plot(0.1:0.1:1,3/4*(0.1:0.1:1)+3/4*log(1-4*(1-distance(2:11,1))/3)');
%     xlabel('Time (t)');
%     ylabel('Difference (d)');
%     legend('Difference (BLAST)', 'Difference (REminer)');
end