%cross-correlating spike trains for Uprobe recordings using Attn structures
%load EthelAttn92713new
%load EthelAttnCon92713new

%[datafile,path] = uigetfile('C:\Users\Farran Briggs\Documents\MATLAB\matlabfiles\');
%load(datafile);

name2 = 'When was this data file recorded?';
prompt2 = {'1 = pre 3/14 and NOT special; 2 = pre 3/14 Special (8/9/13 and 10/10/13); 3 = post 3/14'};
numlines = 1;
defaultparams2 = {'1'};
params2 = inputdlg(prompt2,name2,numlines,defaultparams2);
newfile = str2double(params2{1});

% duration of the grating
xcorrdur = 1;
%edges of our time bins (1ms bins)
edge = (0:0.001:xcorrdur);
win = 50;
c = []; d = []; uAH = []; uBH = []; uCH = []; lengthAs = []; lengthBs = []; lengthCs = [];
uAH2 = []; uBH2 = []; uCH2 = []; lengthAs2 = []; lengthBs2 = []; lengthCs2 = [];

uAHS = []; uBHS = []; uCHS = []; lengthAsS = []; lengthBsS = []; lengthCsS = [];
uAH2S = []; uBH2S = []; uCH2S = []; lengthAs2S = []; lengthBs2S = []; lengthCs2S = [];

%calculate the shuffle shift in msec (1 grating cycle)
gratTf = str2double(data.gratingparams{8,2});
shufwin = 1000/gratTf / 1000;

%loop over answer timestamps, finding all grating start time-stamps (TSs) that are within
%7 secs of each answer TS, saving the TS closest to the answer TS to
%'gratts'. This calculation varies depending on when the file was recorded.

switch newfile
    case 1
        for s = 1:length(data.answerts)
            f = []; a = [];
            f = find(data.contrastchangets > (data.answerts(s) - 4) & data.contrastchangets < data.answerts(s));
            CGts(s) = data.contrastchangets(f(end));
            a = find(data.gratstartts > (data.answerts(s) - 4) & data.gratstartts < data.answerts(s));
            gratts(s) = CGts(s) - rem((CGts(s) - data.gratstartts(a(end))),shufwin) - 1;
            gratdur(s) = CGts(s) - gratts(s);
            startts(s) = data.gratstartts(a(end));
        end
    case 2
        for s = 1:length(data.answerts)
            f = []; a = [];
            f = find(data.contrastchangets > (data.answerts(s) - 7) & data.contrastchangets < data.answerts(s));
            CGts(s) = data.contrastchangets(f(end));
            a = find(data.gratstartts > (data.answerts(s) - 7) & data.gratstartts < data.answerts(s));
            gratts(s) = CGts(s) - rem((CGts(s) - data.gratstartts(a(end))),shufwin) - 1;
            gratdur(s) = CGts(s) - gratts(s);
            startts(s) = data.gratstartts(a(end));
        end
    case 3
        for s = 1:length(data.answerts)
            a = []; f = [];
            f = find(data.contrastchangets > (data.answerts(s) - 7) & data.contrastchangets < data.answerts(s));
            CGts(s) = data.contrastchangets(f(end));
            a = find(data.gratstartts > (CGts(s) - 3) & data.gratstartts < CGts(s));
            gratts(s) = data.gratstartts(a(end)) - 1;
            startts(s) = data.gratstartts(a(1));
            gratdur(s) = CGts(s) - gratts(s);
        end
end


%separate out attend-toward (RTmat 3rd column = 1) and attend-away (2) trials
%Attend toward= A2 versus Attend away= AW
c = find(data.RTmat(:,3) == 1);
d = find(data.RTmat(:,3) == 2);
A2ts = gratts(c);
AWts = gratts(d);

%loop over our A2 TSs making PSTHs per trial and also generated shuffled
%PSTH, which is shifted by 1 grating cycle

for x = 1:length(A2ts)
    %loop over each unitA
    for y = 1:length(data.unitAts(1,:))
        aA = []; uA = []; aAS = []; uAS = [];
        %find all unitsA TSs that occured within 1 sec analysis window
        aA = find(data.unitAts(:,y) > A2ts(x) & data.unitAts(:,y) < A2ts(x) + 1);
        %convert all TSs to be relative to the current A2TS
        uA = data.unitAts(aA,y) - A2ts(x);
        %keep track of the number of TSs per trial
        lengthAs(y,x) = length(uA);
        %make sure there are spikes on a given trial
        if length(uA) > 1
            %make PSTH
            uAH(:,y,x) = histc(uA,edge);
        else
            %if no spikes, PSTH is zeros
            uAH(:,y,x) = zeros(length(edge),1);
        end
        %do the same steps but for the shuffled window
        aAS = find(data.unitAts(:,y) > A2ts(x) - shufwin & data.unitAts(:,y) < A2ts(x) + 1 - shufwin);
        uAS = data.unitAts(aAS,y) - (A2ts(x) - shufwin);
        lengthAsS(y,x) = length(uAS);
        if length(uAS) > 1
            uAHS(:,y,x) = histc(uAS,edge);
        else
            uAHS(:,y,x) = zeros(length(edge),1);
        end
    end
    
    %do the exact same thing for unitBs and unitCs (if there are any)
    if length(data.unitBts(:,1)) > 1
        for y = 1:length(data.unitBts(1,:))
            aB = []; uB = []; aBS = []; uBS = [];
            aB = find(data.unitBts(:,y) > A2ts(x) & data.unitBts(:,y) < A2ts(x) + 1);
            uB = data.unitBts(aB,y) - A2ts(x);
            lengthBs(y,x) = length(uB);
            if length(uB) > 1
                uBH(:,y,x) = histc(uB,edge);
            else
                uBH(:,y,x) = zeros(length(edge),1);
            end
            
            aBS = find(data.unitBts(:,y) > A2ts(x) - shufwin & data.unitBts(:,y) < A2ts(x) + 1 - shufwin);
            uBS = data.unitBts(aB,y) - (A2ts(x) - shufwin);
            lengthBsS(y,x) = length(uBS);
            if length(uBS) > 1
                uBHS(:,y,x) = histc(uBS,edge);
            else
                uBHS(:,y,x) = zeros(length(edge),1);
            end
        end
    end
    
    if length(data.unitCts(:,1)) > 1
        for y = 1:length(data.unitCts(1,:))
            aC = []; uC = []; aCS = []; uCS = [];
            aC = find(data.unitCts(:,y) > A2ts(x) & data.unitCts(:,y) < A2ts(x) + 1);
            uC = data.unitCts(aC,y) - A2ts(x);
            lengthCs(y,x) = length(uC);
            if length(uC) > 1
                uCH(:,y,x) = histc(uC,edge);
            else
                uCH(:,y,x) = zeros(length(edge),1);
            end
            
            aCS = find(data.unitCts(:,y) > A2ts(x) - shufwin & data.unitCts(:,y) < A2ts(x) + 1 - shufwin);
            uCS = data.unitCts(aCS,y) - (A2ts(x) - shufwin);
            lengthCsS(y,x) = length(uCS);
            if length(uCS) > 1
                uCHS(:,y,x) = histc(uC,edge);
            else
                uCHS(:,y,x) = zeros(length(edge),1);
            end
        end
    end
    
end

%Create one big matrix of PSTHs for all our good units
% uAH = uAH(:,GoodUnitAs,:);
Bs = find(ismember(channels,GoodUnitBs));
Cs = find(ismember(channels,GoodUnitCs));
uBH = uBH(:,Bs,:);
uCH = uCH(:,Cs,:);
units = cat(2,uAH,uBH,uCH);

%make vector of spike counts per trial for all good units
%lengthUnits = cat(1,lengthAs(GoodUnitAs,:),lengthBs(GoodUnitBs,:),lengthCs(GoodUnitCs,:));
lengthUnits = cat(1,lengthAs,lengthBs(Bs,:),lengthCs(Cs,:));
%figure out how many units total and then calculate all possible
%combinations for cross-correlation analysis
numUnits = length(units(1,:,1));
Combos = combnk(1:numUnits,2);
numPairs = length(Combos(:,1));
A2Xcorr = zeros(length(units(1,1,:)),(win*2 + 1),numPairs);
Shufcorr2 = zeros(length(units(1,1,:)),(win*2 + 1),numPairs);
A2XcorrShuf = zeros(length(units(1,1,:)),(win*2 + 1),numPairs);
%create one bit matrix of PSTHs for all your good units' shuffled data
% uAHS = uAHS(:,GoodUnitAs,:);
uBHS = uBHS(:,Bs,:);
uCHS = uCHS(:,Cs,:);
unitsS = cat(2,uAHS,uBHS,uCHS);
%calculate cross-correlations for all pairings (regular and shuffled) then
%do shuffle correction
for f = 1:length(A2ts)
    for g = 1:numPairs
        A2Xcorr(f,:,g) = xcorr(units(:,Combos(g,1),f),units(:,Combos(g,2),f),win); % / lengthUnits(Combos(g,1),f);
        Shufcorr2(f,:,g) = xcorr(units(:,Combos(g,1),f),unitsS(:,Combos(g,2),f),win); % / lengthUnits(Combos(g,1),f);
        A2XcorrShuf(f,:,g) = A2Xcorr(f,:,g)- Shufcorr2(f,:,g);
    end
end
%sum over trials and then reshape
sumA2XcorrShuf = nansum(A2XcorrShuf,1);
sumA2XcorrShuf = reshape(sumA2XcorrShuf,(win*2 +1), numPairs);

%repeat entire process for attend-away trial data
for x = 1:length(AWts)
    for y = 1:length(data.unitAts(1,:))
        aA = []; uA = []; aAS = []; uAS = [];
        aA = find(data.unitAts(:,y) > AWts(x) & data.unitAts(:,y) < AWts(x) + 1);
        uA = data.unitAts(aA,y) - AWts(x);
        lengthAs2(y,x) = length(uA);
        if length(uA) > 1
            uAH2(:,y,x) = histc(uA,edge);
        else
            uAH2(:,y,x) = zeros(length(edge),1);
        end
        aAS = find(data.unitAts(:,y) > AWts(x) - shufwin & data.unitAts(:,y) < AWts(x) + 1 - shufwin);
        uAS = data.unitAts(aAS,y) - (AWts(x) - shufwin);
        lengthAs2S(y,x) = length(uAS);
        if length(uAS) > 1
            uAH2S(:,y,x) = histc(uAS,edge);
        else
            uAH2S(:,y,x) = zeros(length(edge),1);
        end
    end
    
    if length(data.unitBts(:,1)) > 1
        for y = 1:length(data.unitBts(1,:))
            aB = []; uB = []; aBS = []; uBS = [];
            aB = find(data.unitBts(:,y) > AWts(x) & data.unitBts(:,y) < AWts(x) + 1);
            uB = data.unitBts(aB,y) - AWts(x);
            lengthBs2(y,x) = length(uB);
            if length(uB) > 1
                uBH2(:,y,x) = histc(uB,edge);
            else
                uBH2(:,y,x) = zeros(length(edge),1);
            end
            
            aBS = find(data.unitBts(:,y) > AWts(x) - shufwin & data.unitBts(:,y) < AWts(x) + 1 - shufwin);
            uBS = data.unitBts(aB,y) - (AWts(x) - shufwin);
            lengthBs2S(y,x) = length(uBS);
            if length(uBS) > 1
                uBH2S(:,y,x) = histc(uBS,edge);
            else
                uBH2S(:,y,x) = zeros(length(edge),1);
            end
        end
    end
    
    if length(data.unitCts(:,1)) > 1
        for y = 1:length(data.unitCts(1,:))
            aC = []; uC = []; aCS = []; uCS = [];
            aC = find(data.unitCts(:,y) > AWts(x) & data.unitCts(:,y) < AWts(x) + 1);
            uC = data.unitCts(aC,y) - AWts(x);
            lengthCs2(y,x) = length(uC);
            if length(uC) > 1
                uCH2(:,y,x) = histc(uC,edge);
            else
                uCH2(:,y,x) = zeros(length(edge),1);
            end
            
            aCS = find(data.unitCts(:,y) > AWts(x) - shufwin & data.unitCts(:,y) < AWts(x) + 1 - shufwin);
            uCS = data.unitCts(aCS,y) - (AWts(x) - shufwin);
            lengthCs2S(y,x) = length(uCS);
            if length(uCS) > 1
                uCH2S(:,y,x) = histc(uC,edge);
            else
                uCH2S(:,y,x) = zeros(length(edge),1);
            end
        end
    end
end

% uAH2 = uAH2(:,GoodUnitAs,:);
uBH2 = uBH2(:,Bs,:);
uCH2 = uCH2(:,Cs,:);
units2 = cat(2,uAH2,uBH2,uCH2);
%lengthUnits2 = cat(1,lengthAs2(GoodUnitAs,:),lengthBs2(GoodUnitBs,:),lengthCs2(GoodUnitCs,:));
lengthUnits2 = cat(1,lengthAs2,lengthBs2(Bs,:),lengthCs2(Cs,:));
numUnits2 = length(units2(1,:,1));
Combos2 = combnk(1:numUnits2,2);
numPairs2 = length(Combos2(:,1));

AWXcorr = zeros(length(units2(1,1,:)),(win*2 + 1),numPairs2);
ShufcorrW = zeros(length(units2(1,1,:)),(win*2 + 1),numPairs2);
AWXcorrShuf = zeros(length(units2(1,1,:)),(win*2 + 1),numPairs2);

% uAH2S = uAH2S(:,GoodUnitAs,:);
uBH2S = uBH2S(:,Bs,:);
uCH2S = uCH2S(:,Cs,:);
units2S = cat(2,uAH2S,uBH2S,uCH2S);

for f = 1:length(AWts)
    for g = 1:numPairs2
        AWXcorr(f,:,g) = xcorr(units2(:,Combos2(g,1),f),units2(:,Combos2(g,2),f),win); % / lengthUnits2(Combos2(g,1),f);
        ShufcorrW(f,:,g) = xcorr(units2(:,Combos2(g,1),f),units2S(:,Combos2(g,2),f),win); % / lengthUnits2(Combos2(g,1),f);
        AWXcorrShuf(f,:,g) = AWXcorr(f,:,g) - ShufcorrW(f,:,g);
    end
end

sumAWXcorrShuf = nansum(AWXcorrShuf,1);
sumAWXcorrShuf = reshape(sumAWXcorrShuf,(win*2 +1), numPairs);


allunits = cat(2,channels,GoodUnitBs*100,GoodUnitCs*1000);
Output = allunits(Combos);
for r = 1:numPairs
    [max1 ind1] = max(sumA2XcorrShuf(26:75,r));
    [max2 ind2] = max(sumAWXcorrShuf(26:75,r));
    if max1 > max2
        indtemp = ind1 + 25;
        else indtemp = ind2 + 25;
        end
        Output(r,3) = indtemp - 50;
        if Output(r,3) >= 0
            sumA2XcorrShuf(:,r) = sumA2XcorrShuf(:,r) / sum(lengthUnits(Combos(r,1),:)) * 100;
            sumAWXcorrShuf(:,r) = sumAWXcorrShuf(:,r) / sum(lengthUnits2(Combos2(r,1),:)) * 100;
        else
            sumA2XcorrShuf(:,r) = sumA2XcorrShuf(:,r) / sum(lengthUnits(Combos(r,2),:)) * 100;
            sumAWXcorrShuf(:,r) = sumAWXcorrShuf(:,r) / sum(lengthUnits2(Combos2(r,2),:)) * 100;
        end
        sigA2XcorrShuf = nanmean(sumA2XcorrShuf,1) + 2*nanstd(sumA2XcorrShuf,0,1);
        sigAWXcorrShuf = nanmean(sumAWXcorrShuf,1) + 2*nanstd(sumAWXcorrShuf,0,1);
    if (max1 > sigA2XcorrShuf(r) || max2 > sigAWXcorrShuf(r)) && (max1 > 1 || max2 > 1)
        Output(r,4) = sum(sumA2XcorrShuf(indtemp-1:indtemp+1,r));
        Output(r,5) = sum(sumAWXcorrShuf(indtemp-1:indtemp+1,r));
    end
end
badpairs = find(Output(:,4)==0);
Output(badpairs,:) = [];
Xcorr2 = sumA2XcorrShuf;
Xcorr2(:,badpairs) = [];
XcorrW = sumAWXcorrShuf;
XcorrW(:,badpairs) = [];
sigA2XcorrShuf(badpairs) = [];
sigAWXcorrShuf(badpairs) = [];
numGoodPairs = length(Output(:,1));
numFigs = floor(numGoodPairs / 9);
if rem(numGoodPairs,9) > 0
    extra = 1;
else extra = 0;
end
numFigs = numFigs + extra;

for i = 1:numFigs
    figure
    for j = 1:9
        if ((i-1)*9)+j <= numGoodPairs
            subplot(3,3,j); plot((-25:1:24),Xcorr2(26:75,((i-1)*9)+j),'r')
            set(gca,'XLim',[-25 24])
            hold on
            subplot(3,3,j); plot((-25:1:24),XcorrW(26:75,((i-1)*9)+j),'b')
            set(gca,'XLim',[-25 24])
            hold on
            subplot(3,3,j); plot((-25:1:24), sigA2XcorrShuf(((i-1)*9)+j)*ones(1,50), '--r')
            hold on
            subplot(3,3,j); plot((-25:1:24), sigAWXcorrShuf(((i-1)*9)+j)*ones(1,50), '--b')
            title(['unit ',num2str(Output(((i-1)*9)+j,1)),' Xcorr unit ',num2str(Output(((i-1)*9)+j,2))])
        end
    end
end

%for plotting individual xcorrs for example pairs
% figure
% z = 16;
% plot((-5:1:5),Xcorr2(45:55,z),'r')
%             set(gca,'XLim',[-5 5])
%             hold on
%             plot((-5:1:5),XcorrW(45:55,z),'b')
%             set(gca,'XLim',[-5 5])
%             hold on
%             plot((-5:1:5), sigA2XcorrShuf(z)*ones(1,11)*-1, '--r')
%             hold on
%             plot((-5:1:5), sigAWXcorrShuf(z)*ones(1,11)*-1, '--b')
%             title(['unit ',num2str(Output(z,1)),' Xcorr unit ',num2str(Output(z,2))])
%             hold on
%             fit1 = fit((-5:1:5)',Xcorr2(45:55,z),"smoothingspline");
%             plot(fit1,'y')

            

