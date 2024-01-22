function [dat1,dat2, dat1Norm, dat2Norm,goodDat1Subj,goodDat2Subj,avgDat1,avgDat2] = normReHo(outDir1,outDir2,filt1,filt2,vox1,vox2,varargin)
%varargin{1} and varargin{2} will have subject name filters for the
%reho files
%% Normalize Data
% normalize preterm data
dat1 = dir(fullfile(outDir1,filt1));
dat2 = dir(fullfile(outDir2,filt2));

if ~isempty(varargin)
    if ~isempty(varargin{1})
        for s = 1:length(dat1)
            tmp = strsplit(dat1(s).name,'_');
            subj{s} = [tmp{1} '_' tmp{2}];
        end
        filtIdx = ismember(subj,varargin{1});
        dat1 = dat1(filtIdx);
        clear subj
    end
    if ~isempty(varargin{2})
        for s = 1:length(dat2)
            tmp = strsplit(dat2(s).name,'_');
            subj{s} = [tmp{1} '_' tmp{2}];
        end
        filtIdx = ismember(subj,varargin{2});
        dat2 = dat2(filtIdx);
        clear subj
    end
end

track = 1;
for s = 1:length(dat1)
    ds = cosmo_fmri_dataset(fullfile(dat1(s).folder,dat1(s).name));
    numVox = size(ds.samples,2);
    if numVox == vox1
        goodDat1Subj(track) = s;
        tmp = ds.samples;
        tmp2 = ds.samples(ds.samples~=0);
        avgDat1(s) = nanmean(tmp2);
        tmp(ds.samples~=0) = (tmp2-avgDat1(s));
%                 tmp(ds.samples~=0) = tmp2;
        dat1Norm(track,:) = tmp;
        track = track+1;
    end
end

% normalize term data
track = 1;
for s = 1:length(dat2)
    ds = cosmo_fmri_dataset(fullfile(dat2(s).folder,dat2(s).name));
    numVox = size(ds.samples,2);
    if numVox == vox2
        goodDat2Subj(track) = s;
        tmp = ds.samples;
        tmp2 = ds.samples(ds.samples~=0);
        avgDat2(s) = nanmean(tmp2);
        tmp(ds.samples~=0) = (tmp2-avgDat2(s));
%                 tmp(ds.samples~=0) = tmp2;
        dat2Norm(track,:) = tmp;
        track = track+1;
    end
end
