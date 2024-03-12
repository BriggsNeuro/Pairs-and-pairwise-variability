%Calculate spike count correlations for known connected pairs

%load PairsList
numRep = length(PairsList);

for i = 1:13  %normally is numRep
    c=[]; d=[]; r=[]; DD=[]; cellA=[]; cellB=[]; spikesA=[]; spikesB=[]; gratchange=[];
    scA=[]; scB=[]; scA2=[]; scAW=[]; scB2=[]; scBW=[]; pcc2=[]; pccW=[];
    datafile = char(PairsList{i,1});
    load(datafile);
    DD = contains(datafile,'Thomas');
    cellA = cell2mat(PairsList{i,2});
    cellB = cell2mat(PairsList{i,3});
    
    if DD == 1 %analyze Thomas data
        sessions = cell2mat({data.num});
        [a,r] = ismember(cell2mat(PairsList{i,4}),sessions);
        c = find(data(r).RTmat(:,3) == 0);
        d = find(data(r).RTmat(:,3) == 1);
        cellAunit = cellA(1);
        cellBunit = cellB(1);
        switch cellAunit
            case 1
                spikesA = data(r).ajj.spikes.elect1;
            case 2
                spikesA = data(r).ajj.spikes.elect2;
            case 3
                spikesA = data(r).ajj.spikes.elect3;
        end
        switch cellBunit
            case 1
                spikesB = data(r).ajj.spikes.elect1;
            case 2
                spikesB = data(r).ajj.spikes.elect2;
            case 3
                spikesB = data(r).ajj.spikes.elect3;
            case 4
                spikesB = data(r).ajj.spikes.elect4;
        end
        gratchange = data(r).contrastchangets;
        
    else %analyze Ethel data
        [a,cellAunit] = ismember(cellA(1),channels);
        [a,cellBunit] = ismember(cellB(1),channels);
        if cellA(2) == 1
            spikesA = data.unitAts(:,cellAunit);
        elseif cellA(2) == 2
            spikesA = data.unitBts(:,cellAunit);
        end
        if cellB(2) == 1
            spikesB = data.unitAts(:,cellBunit);
        elseif cellB(2) == 2
            spikesB = data.unitBts(:,cellBunit);
        end
        
        spikesA(find(spikesA == 0)) = [];
        spikesB(find(spikesB == 0)) = [];
        
        for s = 1:length(data.answerts)
            g = [];
            g = find(data.contrastchangets > (data.answerts(s) - 7) & data.contrastchangets < data.answerts(s));
            gratchange(s) = data.contrastchangets(g(end));
        end
        
        c = find(data.RTmat(:,3) == 1);
        d = find(data.RTmat(:,3) == 2);
    end
    
    for x = 1:length(gratchange)
        aA = []; bB=[];
        aA = find(spikesA < gratchange(x) & spikesA > gratchange(x) - 1.3);
        if length(aA) > 1
            scA(x) = length(aA);
        else scA(x) = 0;
        end
        bB = find(spikesB < gratchange(x) & spikesB > gratchange(x) - 1.3);
        if length(bB) > 1
            scB(x) = length(bB);
        else scB(x) = 0;
        end
    end
    
    scA2 = scA(c);
    scAW = scA(d);
    scB2 = scB(c);
    scBW = scB(d);
    ProbZ(i,1) = length(find(scA2 == 0)) / length(scA2);
    ProbZ(i,2) = length(find(scAW == 0)) / length(scAW);
    ProbZ(i,3) = length(find(scB2 == 0)) / length(scB2);
    ProbZ(i,4) = length(find(scBW == 0)) / length(scBW);

    figure
    subplot(1,2,1), scatter(scA2,scB2)
    hold on
    subplot(1,2,2), scatter(scAW,scBW)
    
    pcc2 = corr(cat(2,scA2',scB2'));
    pccW = corr(cat(2,scAW',scBW'));
    PairsPCC(i,1) = pcc2(2,1);
    PairsPCC(i,2) = pccW(2,1);
    FR(i,1) = sum(scA2) / length(c) / 1.3;
    FR(i,2) = sum(scAW) / length(d) / 1.3;
    FR(i,3) = sum(scB2) / length(c) / 1.3;
    FR(i,4) = sum(scBW) / length(d) / 1.3;
    
    clear data channels filelist GoodUnitAs GoodUnitBs GoodUnitCs
end