function [subj2get_11, subj2get_12, subj2get_com,...
    subj2get_21, subj2get_22, subj2get_com2] = getLongSubj(volThresh,gaFilt,kidoScore)
serverDir = '/data/mril/users/all/mrdata/research/processed/CNMC/';
pretermDemo = readtable('/home/sdamera/Documents/DBI/fullPretermDemoSheet.csv');
pretermDemo.Project = regexprep(pretermDemo.Project,'crib','cbl_preterm'); % this is how it's saved on the server

%% For comparing Scans <32 to 32<scans<37
relSubj = pretermDemo.ScanPMA<=32 & pretermDemo.Scanner==450 ... 
    & pretermDemo.Kidokoro <kidoScore & ismember(pretermDemo.Scan,'scan_01');

subj2get_11 = pretermDemo(relSubj,:);
nVols_cens = getNVolsFromCensFile(subj2get_11,serverDir);
subj2get_11 = filterTable(subj2get_11,nVols_cens,volThresh);

if gaFilt
    relSubj = pretermDemo.BirthGA<=32&pretermDemo.ScanPMA>32& ...
        pretermDemo.ScanPMA<37 & pretermDemo.Scanner==450 & ...
        pretermDemo.Kidokoro < kidoScore ;
else
    relSubj = pretermDemo.ScanPMA>32& ...
        pretermDemo.ScanPMA<37 & pretermDemo.Scanner==450 & ...
        pretermDemo.Kidokoro < kidoScore ;
end
subj2get_12 = pretermDemo(relSubj,:);
nVols_cens = getNVolsFromCensFile(subj2get_12,serverDir);
subj2get_12 = filterTable(subj2get_12,nVols_cens,volThresh);

[~,idx] = intersect(subj2get_12.ID,subj2get_11.ID);
[~,idx2] = intersect(subj2get_11.ID,subj2get_12.ID);
tmp = subj2get_11(idx2,:);
subj2get_com = subj2get_12(idx,:);
subj2get_com.Scan1PMA = tmp.ScanPMA;
subj2get_com.meanScan1PMA = tmp.meanScanPMA;
subj2get_com.stdScan1PMA = tmp.stdScanPMA;
subj2get_com.mean_GA_scan1 = tmp.mean_GA_scan;
subj2get_com.std_GA_scan1 = tmp.std_GA_scan;
subj2get_com.Scan2PMA = subj2get_com.ScanPMA;
subj2get_com.meanScan2PMA = subj2get_com.meanScanPMA;
subj2get_com.stdScan2PMA = subj2get_com.stdScanPMA;
subj2get_com.mean_GA_scan2 = subj2get_com.mean_GA_scan;
subj2get_com.std_GA_scan2 = subj2get_com.std_GA_scan;
subj2get_com(:,[7 12:15]) = [];
clear std* mean* idx idx2 ia

%% For comparing 32<scans<37 to TEA Scans

if gaFilt
    relSubj = pretermDemo.BirthGA<=32&pretermDemo.ScanPMA>32 &...
        pretermDemo.ScanPMA<37 & pretermDemo.Scanner==450 & ...
        pretermDemo.Kidokoro < kidoScore;
else
    relSubj = pretermDemo.ScanPMA>=32 &...
        pretermDemo.ScanPMA<37 & pretermDemo.Scanner==450 & ...
        pretermDemo.Kidokoro < kidoScore;
end


subj2get_21 = pretermDemo(relSubj,:);
nVols_cens = getNVolsFromCensFile(subj2get_21,serverDir);
subj2get_21 = filterTable(subj2get_21,nVols_cens,volThresh);
if gaFilt
    relSubj = pretermDemo.BirthGA<=32&pretermDemo.ScanPMA>=37 & ...
        pretermDemo.Scanner==750 & ...
        pretermDemo.Kidokoro < kidoScore;
else
    relSubj =  pretermDemo.ScanPMA>=37 & ...
        pretermDemo.Scanner==750 & ...
        pretermDemo.Kidokoro < kidoScore;
end
subj2get_22 = pretermDemo(relSubj,:);
nVols_cens = getNVolsFromCensFile(subj2get_22,serverDir);
subj2get_22 = filterTable(subj2get_22,nVols_cens,volThresh);

[~,idx] = intersect(subj2get_22.ID,subj2get_21.ID);
[~,idx2] = intersect(subj2get_21.ID,subj2get_22.ID);
tmp = subj2get_21(idx2,:);
subj2get_com2 = subj2get_22(idx,:);
subj2get_com2.Scan1PMA = tmp.ScanPMA;
subj2get_com2.meanScan1PMA = tmp.meanScanPMA;
subj2get_com2.stdScan1PMA = tmp.stdScanPMA;
subj2get_com2.mean_GA_scan1 = tmp.mean_GA_scan;
subj2get_com2.std_GA_scan1 = tmp.std_GA_scan;
subj2get_com2.Scan2PMA = subj2get_com2.ScanPMA;
subj2get_com2.meanScan2PMA = subj2get_com2.meanScanPMA;
subj2get_com2.stdScan2PMA = subj2get_com2.stdScanPMA;
subj2get_com2.mean_GA_scan2 = subj2get_com2.mean_GA_scan;
subj2get_com2.std_GA_scan2 = subj2get_com2.std_GA_scan;
subj2get_com2(:,[7 12:15]) = [];

function nVols_cens = getNVolsFromCensFile(subj2get,serverDir)
for s = 1:size(subj2get,1)
    proj = subj2get.Project{s};
    subj = subj2get.ID{s};
    scan = subj2get.Scan{s};
    censFile = fullfile(serverDir,proj,'/neonatal/wip/rs/preproc/', ...
        subj,scan,'func',...
        [subj '_' scan(end-1:end) '_combined_2.1D']);
    if exist(censFile,'file')
        nVols_cens(s) = sum(textread(censFile));
    else
        nVols_cens(s) = NaN;
    end
end

function subj2get = filterTable(subj2get,nVols_cens,volThresh)

subj2get = subj2get(nVols_cens>=volThresh,:);
goodVols = nVols_cens(nVols_cens>=volThresh);
subj2get.NVol = goodVols';
subj2get = subj2get(:,[1:5 7:9 11]);
[~,ia] = unique(subj2get.ID);
subj2get = subj2get(ia,:);
subj2get.meanGA = repmat(mean(subj2get.BirthGA),size(subj2get,1),1);
subj2get.stdGA = repmat(std(subj2get.BirthGA),size(subj2get,1),1);
subj2get.meanScanPMA = repmat(mean(subj2get.ScanPMA),size(subj2get,1),1);
subj2get.stdScanPMA = repmat(std(subj2get.ScanPMA),size(subj2get,1),1);
subj2get.mean_GA_scan = repmat(mean(subj2get.ExUteroWeeks)...
    ,size(subj2get,1),1);
subj2get.std_GA_scan = repmat(std(subj2get.ExUteroWeeks)...
    ,size(subj2get,1),1);