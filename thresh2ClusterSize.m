function threshMap = thresh2ClusterSize(map,nbrhood,clusSize,varargin)

%Iteratively thresholds map till a particular cluster size is achieved -
%works best (for now) if only 1 cluster is present in map
nbrhood_mat=cosmo_convert_neighborhood(nbrhood,'matrix');

%mask data if it is provided (mask is a logical that is the same length as
%map)

if ~isempty(varargin)
    mask = varargin{1};
    map(~mask) = 0;
end

thresh = max(map):-0.0001:0;
track=1;
clusLen = 0;
while clusLen<clusSize && track<=length(thresh)
    smpl = map;
    smpl(smpl<thresh(track)) = 0;
    smpl((smpl~=0)&~isnan(smpl)) = 1;
    clusters=cosmo_clusterize(smpl,nbrhood_mat);
    if ~isempty(clusters)
        [sortVals,sortIdx] = sort(cellfun(@length,clusters),'descend');
        clusLen = sortVals(1);
        maxClus = clusters(sortIdx(1));
    end
    track = track+1;
    
end
if clusLen==0
    error('cluster thresholding failed - desired size was too large')
end

threshMap = zeros(size(map));
threshMap(maxClus{1}) = 1;