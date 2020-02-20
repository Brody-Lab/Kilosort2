function [rez, X] = splitAllClusters(rez, flag)
% i call this algorithm "bimodal pursuit"
% split clusters if they have bimodal projections
% the strategy is to maximize a bimodality score and find a single vector projection
% that maximizes it. If the distribution along that maximal projection crosses a
% bimodality threshold, then the cluster is split along that direction
% it only uses the PC features for each spike, stored in rez.cProjPC

wPCA = gather(rez.ops.wPCA); % use PCA projections to reconstruct templates when we do splits

ccsplit = rez.ops.AUCsplit; % this is the threshold for splits, and is one of the main parameters users can change

NchanNear   = min(rez.ops.Nchan, 32);
Nnearest    = min(rez.ops.Nchan, 32);
sigmaMask   = rez.ops.sigmaMask;

ik = 0;
Nfilt = size(rez.W,2);
nsplits= 0;

[iC, mask, C2C] = getClosestChannels(rez, sigmaMask, NchanNear); % determine what channels each template lives on

rez.ops.nt0min = getOr(rez.ops, 'nt0min', 20); % the waveforms must be aligned to this sample

[~, iW] = max(abs(rez.dWU(rez.ops.nt0min, :, :)), [], 2); % find the peak abs channel for each template
iW = squeeze(int32(iW));

isplit = 1:Nfilt; % keep track of original cluster for each cluster. starts with all clusters being their own origin.
dt = 1/1000;
nccg = 0;
ncorr=0;
split_crit=0;

while ik<Nfilt
    if rem(ik, 100)==0
      % periodically write updates
      if ik>0
        fprintf('Found %d splits, checked %d/%d clusters, nccg %d \n', nsplits, ik, Nfilt, nccg)
      end
    end
    ik = ik+1;

    %
    isp = find(rez.st3(:,2)==ik); % get all spikes from this cluster
    nSpikes = numel(isp);
    if  nSpikes<300
       continue; % do not split if fewer than 300 spikes (we cannot estimate cross-correlograms accurately)
    end

    ss = rez.st3(isp,1)/rez.ops.fs; % convert to seconds

    clp0 = rez.cProjPC(isp, :, :); % get the PC projections for these spikes
    clp0 = gpuArray(clp0(:,:));
    clp = clp0 - mean(clp0,1); % mean center them

    clp = clp - my_conv2(clp, 250, 1); % subtract a running average, because the projections are NOT drift corrected

    % now use two different ways to initialize the bimodal direction
    % the main script calls this function twice, and does both initializations
    if flag
        [u s v] = svdecon(clp');
        w = u(:,1); % initialize with the top PC
    else
        w = mean(clp0, 1)'; % initialize with the mean of NOT drift-corrected trace
        w = w/sum(w.^2)^.5; % unit-normalize
    end

    % initial projections of waveform PCs onto 1D vector
    x = gather(clp * w);
    s1 = var(x(x>mean(x))); % initialize estimates of variance for the first
    s2 = var(x(x<mean(x))); % and second gaussian in the mixture of 1D gaussians

    mu1 = mean(x(x>mean(x))); % initialize the means as well
    mu2 = mean(x(x<mean(x)));
    p  = mean(x>mean(x)); % and the probability that a spike is assigned to the first Gaussian

    logp = zeros(numel(isp), 2); % initialize matrix of log probabilities that each spike is assigned to the first or second cluster

    % do 50 pursuit iteration
    for k = 1:50
        % for each spike, estimate its probability to come from either Gaussian cluster
        logp(:,1) = -1/2*log(s1) - (x-mu1).^2/(2*s1) + log(p);
        logp(:,2) = -1/2*log(s2) - (x-mu2).^2/(2*s2) + log(1-p);

        lMax = max(logp,[],2);
        logp = logp - lMax; % subtract the max for floating point accuracy
        rs = exp(logp); % exponentiate the probabilities

        pval = log(sum(rs,2)) + lMax; % get the normalizer and add back the max
        logP(k) = mean(pval); % this is the cost function: we can monitor its increase

        rs = rs./sum(rs,2); % normalize so that probabilities sum to 1

        p = mean(rs(:,1)); % mean probability to be assigned to Gaussian 1
        mu1 = (rs(:,1)' * x )/sum(rs(:,1)); % new estimate of mean of cluster 1 (weighted by "responsibilities")
        mu2 = (rs(:,2)' * x )/sum(rs(:,2)); % new estimate of mean of cluster 2 (weighted by "responsibilities")

        s1 = (rs(:,1)' * (x-mu1).^2 )/sum(rs(:,1)); % new estimates of variances
        s2 = (rs(:,2)' * (x-mu2).^2 )/sum(rs(:,2));

        if (k>10 && rem(k,2)==1)
            % starting at iteration 10, we start re-estimating the pursuit direction
            % that is, given the Gaussian cluster assignments, and the mean and variances,
            % we re-estimate w
            StS  = clp' * (clp .* (rs(:,1)/s1 + rs(:,2)/s2))/nSpikes; % these equations follow from the model
            StMu = clp' * (rs(:,1)*mu1/s1 + rs(:,2)*mu2/s2)/nSpikes;

            w = StMu'/StS; % this is the new estimate of the best pursuit direection
            w = normc(w'); % which we unit normalize
            x = gather(clp * w);  % the new projections of the data onto this direction
        end
    end
   

    ilow = rs(:,1)>rs(:,2); % these spikes are assigned to cluster 1
%     ps = mean(rs(:,1));
    plow = mean(rs(ilow,1)); % the mean probability of spikes assigned to cluster 1
    phigh = mean(rs(~ilow,2)); % same for cluster 2
    nremove = min(mean(ilow), mean(~ilow)); % the smallest cluster has this proportion of all spikes

    if sum(ilow)<50 || sum(~ilow)<50
        continue
    end
    % if the cross-correlogram of the split is as or more refractory than the
    % parent cluster, abort
    ccg_data = test_ccgs(ss(ilow),ss(~ilow),500,dt,{'merge','xc'});
    if ccg_data.merge.Q>=ccg_data.xc.Q || ccg_data.merge.Rmin>=ccg_data.xc.Rmin % if both metrics are below threshold.
        nccg = nccg+1; % keep track of how many splits were voided by the CCG criterion
        continue;
    end

    % now decide if the split would result in waveforms that are too similar
    c1  = wPCA * reshape(mean(clp0(ilow,:),1), 3, []); %  the reconstructed mean waveforms for putatiev cluster 1
    c2  = wPCA * reshape(mean(clp0(~ilow,:),1), 3, []); %  the reconstructed mean waveforms for putative cluster 2
    max_lag = round(5e-4 * rez.ops.fs); % up to 0.5ms
    [simscore,amp_score,bestlag] = calculate_simscore(gather(cat(3,c1,c2)),max_lag);
    cc = simscore(2);
    % abort if the mean waveforms are highly correlated.
    if cc>.9 && abs(amp_score(2))<0.2
        ncorr = ncorr+1;
        continue;
    end
    
    if min(plow,phigh)<0.95
        [dip, p_value]=HartigansDipSignifTest(x,500);
    end
    
    %time_separated_stat = abs(CalcCP(rez.st3(isp(ilow),1),rez.st3(isp(~ilow),1))-0.5);

    % finaly criteria to continue with the split: if the split piece is more than 5% of all spikes,
    % if the split piece is more than 300 spikes, and if the confidences for assigning spikes to
    % both clusters exceeds a preset criterion ccsplit
    if  min(plow,phigh)>0.95 || p_value<0.05  % when there are low numbers of spikes or one cluster has many more events than the other, hartigan's dip p_value seems overly conservative and the gaussian mixture model gives a more sensible answe
       % one cluster stays, one goes
       Nfilt = Nfilt + 1;

       % the templates for the splits have been estimated from PC coefficients
       rez.dWU(:,iC(:, iW(ik)),Nfilt) = c2;
       rez.dWU(:,iC(:, iW(ik)),ik)    = c1;

       % the temporal components are therefore just the PC waveforms
       rez.W(:,Nfilt,:) = permute(wPCA, [1 3 2]);
       iW(Nfilt) = iW(ik); % copy the best channel from the original template
       isplit(Nfilt) = isplit(ik); % copy the provenance index to keep track of splits

       rez.st3(isp(ilow), 2)    = Nfilt; % overwrite spike indices with the new index       
       rez.simScore(:, Nfilt)   = rez.simScore(:, ik); % copy similarity scores from the original <--- this is temporary. they get recalculated by recompute_clusters.
       rez.simScore(Nfilt, :)   = rez.simScore(ik, :); % copy similarity scores from the original
       rez.simScore(ik, Nfilt) = 1; % set the similarity with original to 1
       rez.simScore(Nfilt, ik) = 1; % set the similarity with original to 1

       rez.iNeigh(:, Nfilt)     = rez.iNeigh(:, ik); % copy neighbor template list from the original
       rez.iNeighPC(:, Nfilt)     = rez.iNeighPC(:, ik); % copy neighbor channel list from the original

       % try this cluster again
       ik = ik-1; % the cluster piece that stays at this index needs to be tested for splits again before proceeding
       % the piece that became a new cluster will be tested again when we get to the end of the list
       nsplits = nsplits + 1; % keep track of how many splits we did
    else
        split_crit = split_crit+1;
    end
end

fprintf('Finished splitting. Found %d splits, checked %d/%d clusters, rejected because ccg dip %d, because correlated %d, because not well split %d \n', nsplits, ik, Nfilt, nccg, ncorr, split_crit)


isplit = rez.simScore ==1;

rez = recompute_clusters(rez);

rez.isplit = isplit; % keep track of origins for each cluster


% figure(1)
% subplot(1,4,1)
% plot(logP(1:k))
%
% subplot(1,4,2)
% [~, isort] = sort(x);
% epval = exp(pval);
% epval = epval/sum(epval);
% plot(x(isort), epval(isort))
%
% subplot(1,4,3)
% ts = linspace(min(x), max(x), 200);
% xbin = hist(x, ts);
% xbin = xbin/sum(xbin);
%
% plot(ts, xbin)
%
% figure(2)
% plotmatrix(v(:,1:4), '.')
%
% drawnow
%
% % compute scores for splits
% ilow = rs(:,1)>rs(:,2);
% ps = mean(rs(:,1));
% [mean(rs(ilow,1)) mean(rs(~ilow,2)) max(ps, 1-ps) min(mean(ilow), mean(~ilow))]
