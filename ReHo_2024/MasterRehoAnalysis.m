function MasterRehoAnalysis

%% TEA Setup
[~, ~, ~,subj2get_21, subj2get_22, subj2get_com2] = getLongSubj(150,0,8);

%read in term scan QA data
termData = readtable('/home/sdamera/Documents/DBI/termData.csv');
termData = termData(termData.Var21>149,:); % keep only subjects that have greater than 150 volumes (> 5min)
termDemo = readtable('/home/sdamera/Documents/DBI/TermDemo.csv');
termDemo(104:end,:) = [];
outDir = '/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/srdPreproc/TEA/nbrhood_27/ReHo_WB';
filt1 = '*_term_reho.nii.gz';
filt2 = '*_preterm_reho.nii.gz';
vox1 = 259200;
vox2 = 259200;
[dat1,dat2, dat1Norm, dat2Norm] = normReHo(outDir,outDir,filt1,filt2,vox1,vox2,[],subj2get_22.ID);
for s = 1:length(dat2)
    tmp = strsplit(dat2(s).name,'_');
    subj2{s} = [tmp{1} '_' tmp{2}];
end
for s = 1:length(dat1)
    tmp = strsplit(dat1(s).name,'_');
    subj{s} = [tmp{1} '_' tmp{2}];
end
[dat1Norm_filt,subj_filt,filtIdx,filtIdx2] = matchSubj(dat1Norm,subj,termData.Var3);
termData = termData(filtIdx2,:);
dat1 = dat1(filtIdx);
[dat2Norm_filt,~,filtIdx,filtIdx2] = matchSubj(dat2Norm,subj2,subj2get_22.ID);
subj2get_22 = subj2get_22(filtIdx2,:);
dat2 = dat2(filtIdx);
%clean up term data
[dat1Norm_filt,~,filtIdx,filtIdx2] = matchSubj(dat1Norm_filt,subj_filt,termDemo.SubjectID);
termData = termData(filtIdx,:);
dat1 = dat1(filtIdx);
termDemo = termDemo(filtIdx2,:);
termBirthGA = termDemo.ga_birthw;
badIdx = termBirthGA<37;
dat1(badIdx) = [];
dat1Norm_filt(badIdx,:) = [];
termBirthGA(badIdx) = [];
termData(badIdx,:) = []; 
termDemo(badIdx,:) = []; 
clear subj subj2 tmp vox* idx badIdx filtIdx*
allData = readtable('/home/sdamera/Downloads/Archive/AllData-hemiAvg.csv');
sex = cat(1,allData.Sex(1:86),allData.Sex(315:423));
mot = cat(1,allData.motion(1:86),allData.motion(315:423));
mot = mot - mean(mot);
%% TEA analyses
% test effect of group while controlling for scan PMA
scanPMA = cat(1,termData.Var14,subj2get_22.ScanPMA);
% scanPMA = scanPMA-mean(scanPMA);
birthGA = cat(1,termBirthGA,subj2get_22.BirthGA);
PNA = scanPMA-birthGA;
scanPMA = scanPMA-mean(scanPMA);
% birthGA = birthGA-mean(birthGA);
termStatus = categorical(cat(1,ones(size(dat1Norm_filt,1),1),zeros(size(dat2Norm_filt,1),1)));

%whole brain vs 0 across groups controlling for scan PMA
gmMsk = cosmo_fmri_dataset('/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/congrad/masks/infant-neo-GM-mask_3mm.nii');
gmMskIdx = find(gmMsk.samples);
pVal = nan(2,size(gmMsk.samples(1,:),2));
beta = nan(2,size(gmMsk.samples(1,:),2));
for i = gmMskIdx
    tmp = dat1Norm_filt(:,i);
    tmp(tmp==0) = NaN;
    tmp2 = dat2Norm_filt(:,i);
    tmp2(tmp2==0) = NaN;
    reho = cat(1,nanmean(tmp,2),nanmean(tmp2,2));
    tbl = table(reho,termStatus,birthGA,scanPMA,sex,mot);
%     mdl = fitlm(tbl,'reho~scanPMA+sex+mot');
    mdl = fitlm(scanPMA,cat(1,dat1Norm_filt(:,i),dat2Norm_filt(:,i)));
    pVal(1,i) = mdl.Coefficients.pValue(1);
    beta(1,i) = mdl.Coefficients.Estimate(1);
    pVal(2,i) = mdl.Coefficients.pValue(2);
    beta(2,i) = mdl.Coefficients.Estimate(2);
end
ds = gmMsk;
ds.samples = beta(1,:);
ds.samples(pVal(1,:)>0.01) = 0;
cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/avgReHo_effect-Intercept_group2-all_01.nii']);
ds = gmMsk;
ds.samples = beta(2,:);
ds.samples(pVal(2,:)>0.01) = 0;
cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/avgReHo_effect-scanPMA_group2-all_01.nii']);

%create ROIs
beta_sig = beta(1,:);
beta_sig(pVal(1,:)>0.001) = 0;
nbrhood = cosmo_cluster_neighborhood(ds,'fmri',1);
dataDir = '/home/sdamera/Documents/DBI';
atlasFn = fullfile(dataDir,'chd_abc/neonatal/wip/srd/derivatives/atlases', ...
'infant-neo-aal_3mm.nii');
atlas = cosmo_fmri_dataset(atlasFn);
roiAtlasVals{1} = [57:58]; % ss
roiAtlasVals{2} = [43:48]; %visual
roiAtlasVals{3} = [79:82]; % auditory
roiFns = {'lSomatCtx','rSomatCtx','lVisCtx','rVisCtx','lAudCtx','rAudCtx'};
clusSize = 20;
track = 1;
for r = 1:3
    for h = 1:2
        roiMask = ismember(atlas.samples,roiAtlasVals{r}(h:2:end));
        threshMap{r,h} = thresh2ClusterSize(beta_sig,nbrhood,clusSize,roiMask);
        ds = atlas;
        ds.samples = threshMap{r,h};
        cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/' roiFns{track} '_50vox.nii']);
        track = track+1;
    end
end

%plot top  and bottom 10% of voxels
meanReho = beta(1,gmMskIdx);
[sortVals,sortIdx] = sort(meanReho,'descend');
sortIdx(isnan(sortVals)) = [];
topVox = sortIdx(1:ceil(length(gmMskIdx)*0.1));
[sortVals,sortIdx] = sort((meanReho),'ascend');
sortIdx(isnan(sortVals)) = [];
bottomVox = sortIdx(1:ceil(length(gmMskIdx)*0.1));
ds = gmMsk;
tmp = zeros(size(ds.samples));
tmp(gmMskIdx(topVox)) = meanReho(topVox);
tmp(gmMskIdx(bottomVox)) = meanReho(bottomVox);
ds.samples = tmp;
clear tmp*
cosmo_map2fmri(ds,'/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/topReho_per-10_raw.nii')

% Whole brain differences due to term status after controlling for scan PMA
gmMsk = cosmo_fmri_dataset('/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/congrad/masks/infant-neo-GM-mask_3mm.nii');
gmMskIdx = find(gmMsk.samples);
pVal = nan(2,size(gmMsk.samples(1,:),2));
beta = nan(2,size(gmMsk.samples(1,:),2));
for i = gmMskIdx
    tmp = dat1Norm_filt(:,i);
    tmp(tmp==0) = NaN;
    tmp2 = dat2Norm_filt(:,i);
    tmp2(tmp2==0) = NaN;
    reho = cat(1,nanmean(tmp,2),nanmean(tmp2,2));
    tbl = table(reho,termStatus,birthGA,scanPMA,sex,mot);
    mdl = fitlm(tbl,'reho~termStatus+scanPMA+sex+mot');
    pVal(1,i) = mdl.Coefficients.pValue(2);
    beta(1,i) = mdl.Coefficients.Estimate(2);
    pVal(2,i) = mdl.Coefficients.pValue(3);
    beta(2,i) = mdl.Coefficients.Estimate(3);
end
ds = gmMsk;
ds.samples = beta(1,:);
ds.samples(pVal(1,:)>0.001) = 0;
cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/ReHo_effect-TermStatus_group-all.nii']);
ds = gmMsk;
ds.samples = beta(2,:);
ds.samples(pVal(2,:)>0.001) = 0;
cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/ReHo_effect-scanPMA_group-all.nii']);

% ROI ReHo analysis
rois = {'lSomatCtx','rSomatCtx','lVisCtx','rVisCtx','lAudCtx','rAudCtx'};
pVal = [];
allTbl = [];
for r = 1:2:length(rois)
    ds = cosmo_stack({cosmo_fmri_dataset(['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/' rois{r} '.nii']),...
        cosmo_fmri_dataset(['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/' rois{r+1} '_50vox.nii'])},1);
    ds.samples = sum(ds.samples);
    msk = cat(2,find(ds.samples>0));
    mskSize(r) = length(msk);
    tmp = dat1Norm_filt(:,msk);
    tmp(tmp==0) = NaN;
    tmp2 = dat2Norm_filt(:,msk);
    tmp2(tmp2==0) = NaN;
    reho = cat(1,nanmean(tmp,2),nanmean(tmp2,2));
    ROI = repmat({rois{r}(2:end)},size(reho,1),1);
    subj = cat(1,termDemo.SubjectID,subj2get_22.ID);
    tbl = table(reho,termStatus,birthGA,scanPMA,PNA,sex,mot,ROI,subj);
    allTbl = cat(1,allTbl,tbl);
    disp(rois{r})
end

for r = 1:6
    tic
    gmMsk = cosmo_fmri_dataset('/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/congrad/masks/infant-neo-GM-mask_3mm.nii');
    gmMskIdx = find(gmMsk.samples);
    pVal = nan(2,size(gmMsk.samples(1,:),2));
    beta = nan(2,size(gmMsk.samples(1,:),2));
    for i = gmMskIdx
        reho = tbl_ptm{r}.ReHo;
        
        tmp = conMap_term2(:,i,r);
        tmp(tmp==0) = NaN;
        con = nanmean(tmp,2);
        subj = cat(1,termDemo.SubjectID,subj2get_22.ID);
        subj = tbl_ptm{r}.subj;
        mot= cat(1,mot_term,mot_ptm);
        Hemi = cat(1,repmat({'lh'},32,1),repmat({'rh'},32,1));
        birthGA2 = tbl_ptm{r}.birthGA;
        scanPMA2 = tbl_ptm{r}.scanPMA;
        tbl = table(reho,con,subj,Hemi,birthGA2,scanPMA2);
        mdl = fitlme(tbl,'con~reho+birthGA2+scanPMA2+(reho|subj)');
        pVal(1,i) = mdl.Coefficients.pValue(1);
        beta(1,i) = mdl.Coefficients.Estimate(1);
        pVal(2,i) = mdl.Coefficients.pValue(2);
        beta(2,i) = mdl.Coefficients.Estimate(2);
    end
    
    ds = gmMsk;
    ds.samples = beta(1,:);
    ds.samples(pVal(1,:)>0.001) = 0;
    cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/Con_effect-ReHo_roi-' rois{r} '_group-all.nii']);
    ds = gmMsk;
    ds.samples = beta(2,:);
    ds.samples(pVal(2,:)>0.001) = 0;
    cosmo_map2fmri(ds,['/home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/TEA/Con_effect-birthGA_roi-' rois{r} '_group-all.nii']);
    toc
end


%% T2 Analyses
[~, ~, ~,subj2get_21, subj2get_22, subj2get_com2] = getLongSubj(150,0,8);
outDir1 = '/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/srdPreproc/PTM/nbrhood_27/ReHo_WB';
outDir2 = '/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/srdPreproc/TEA/nbrhood_27/ReHo_WB';
filt1 = '*_t2_reho.nii.gz';
filt2 = '*_preterm_reho.nii.gz';
vox1 = 56304; % for 1.5T scans
vox2 = 259200; % for TEA scans
[~,dat,~,datNorm,~,~,avgDat1,avgDat2] = normReHo(outDir2,outDir1,filt2,filt1,vox2,vox1);
for s = 1:length(dat)
    tmp = strsplit(dat(s).name,'_');
    subj{s} = [tmp{1} '_' tmp{2}];
end
[datNorm_filt,~,filtIdx,filtIdx2] = matchSubj(datNorm,subj,subj2get_21.ID);
subj2get_21 = subj2get_21(filtIdx2,:);
dat = dat(filtIdx);


% Warp Infant-Neo-WithCerebellum.nii.gz used in TEA analyses to each of the
% GA specific templates
ageTemplates = [24:37 40];
for t = 1:length(ageTemplates)
    disp(['auto_warp.py -base /data/mril/users/all/mrdata/research/'...
        'processed/CNMC/chd_abc/neonatal/wip/srd/derivatives/atlases/'...
        'preterm_templates/template-' num2str(ageTemplates(t)) '.nii.gz'...
        '-input /home/sdamera/Documents/DBI/chd_abc/neonatal/wip/srd/'...
        'derivatives/atlases/infant-neo-withCerebellum.nii -output_dir'...
        '/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc/'...
        'neonatal/wip/srd/derivatives/atlases/preterm_templates/warp' 
        num2str(ageTemplates(t))])
end

% Write shell script to warp ROIs to the appropriate GA templates
% ROIs have been created in the TEA Analyses above or specified elsewhere
roiDir = ['/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc' ...
    '/neonatal/wip/srd/derivatives/Analyses/reho/vol_150/srdPreproc/'...
    'TEA/nbrhood_27/ReHoROIs'];
rois = dir(fullfile(roiDir,'*.nii'));
outDirName = 'normTEA';
suffix = '_ReHo_WB';
outDir = ['/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc'...
        '/neonatal/wip/srd/derivatives/atlases/preterm_templates/'];
outFn = fullfile(outDir,'warpHippoROIs.sh'); % EDIT THIS FILENAME
fid = fopen(outFn,'w');
for t = 1:length(ageTemplates)
    tempDir = ['/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc'...
        '/neonatal/wip/srd/derivatives/atlases/preterm_templates/warp' num2str(ageTemplates(t))];
    fprintf(fid,['cd ' tempDir ' \n']);
    for r = 1:length(rois)
        gaOutDir = fullfile(tempDir,outDirName);
        if ~exist(gaOutDir,'dir')
            mkdir(gaOutDir)
        end
        fprintf(fid,['3dNwarpApply -nwarp "' fullfile(tempDir, 'anat.un.aff.qw_WARP.nii')...
            ' anat.un.aff.Xat.1D" -master ' fullfile(tempDir, 'base.nii') ...
            ' -dxyz 3 -source ' fullfile(rois(r).folder, rois(r).name) ...
            ' -overwrite -prefix ' fullfile(gaOutDir, regexprep(rois(r).name,...
            '.nii',['-warp' num2str(ageTemplates(t)) suffix '.nii \n']))]);
    end
end
fclose(fid);

% Iteratively threshold warped ROIs to achieve desired cluster size
clusSize = 20; % 20 voxel clusters are desired
parfor t = 1:length(ageTemplates)
    tempDir = ['/data/mril/users/all/mrdata/research/processed/CNMC/chd_abc'...
        '/neonatal/wip/srd/derivatives/atlases/preterm_templates/warp' num2str(ageTemplates(t))];
    rois = dir(fullfile(tempDir,outDirName,['*' suffix '.nii']));
    for r = 1:length(rois)
        tic
        ds = cosmo_fmri_dataset(fullfile(rois(r).folder,rois(r).name));
        nbrhood = cosmo_cluster_neighborhood(ds,'fmri',1);
        threshMap = thresh2ClusterSize(ds.samples,nbrhood,clusSize);
        ds.samples = threshMap;
        cosmo_map2fmri(ds,fullfile(rois(r).folder,regexprep(rois(r).name,'.nii','_final.nii')))
        toc
    end 
end



function [data_filt,list1_filt,filtIdx,filtIdx2] = matchSubj(data,list1,list2)

[filt,filtIdx2] = ismember(list1,list2);
if ~isempty(data)
    data_filt = data(filt,:);
else
    data_filt = [];
end
filtIdx  = find(filt);
list1_filt = list1(filt);
filtIdx2(filtIdx2==0) = [];
