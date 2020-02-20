
% REMOVE_KS2_DUPLICATE_SPIKES2 Double-counted spikes are hard to avoid with
% Kilosort's template matching algorithm since the overall fit can be
% improved by having multiple templates jointly account for a single variable waveform.
% 
% Unlike REMOVE_KS2_DUPLICATE_SPIKES, this function takes a more targeted
% approach. First, if identifies nearby clusters whose ccg shows a peak at the
% central sub-millisecond bin ONLY, with the adjacent bin not showing a peak. 
% This is very unlikely to be real, and is much more likely to be
% indicative of double-counted spikes.
%
% From these pairs, it identifies the pair with the larger template as
% being the "main" or "reference" cluster and the duplicate spikes from the
% other cluster are removed.
%
% You can control the time interval over which to consider duplication with
% OVERLAP_S (0.6ms by deafault).
%
% You can control the required spatial proximity using CHANNEL_SEPARATION_UM.
%
%=INPUT
%
%   rez structure
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   overlap_s
%       the time interval, in second, within which a sequence of spikes are
%       vetted for duplicates.
%
%   channel_separation_um
%       When the primay channels of two spikes are within this distance, in
%       microns, then the two spikes are vetted for duplicate.
%
%=EXAMPLE
%
%   >> rez = remove_ks2_duplicate_spikes(rez)
function rez = remove_ks2_duplicate_spikes2(rez, varargin)
    input_parser = inputParser;
    addParameter(input_parser, 'overlap_s', 5e-4, @(x) (isnumeric(x)))
    addParameter(input_parser, 'channel_separation_um', 100, @(x) (ischar(x)))
    addParameter(input_parser, 'verbose', false, @(x) (islogical(x)))    
    parse(input_parser, varargin{:});
    P = input_parser.Results;
    nbins=50;
    spike_times = uint64(rez.st3(:,1));
    spike_templates = deal(double(rez.st3(:,2)));
    rez.U=gather(rez.U);
    rez.W = gather(rez.W);
    templates = zeros(rez.ops.Nchan, size(rez.W,1), size(rez.W,2), 'single');
    for iNN = 1:size(templates,3)
       templates(:,:,iNN) = squeeze(rez.U(:,iNN,:)) * squeeze(rez.W(:,iNN,:))';
    end
    templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
    spike_times = double(spike_times);
    %% Make sure that the spike times are sorted
    if ~issorted(spike_times)
        error('rez.st3 must be time ordered.');
    end
    %% deal with cluster 0
    if any(spike_templates==0)
        error('Currently this function can''t deal with existence of cluster 0. It ought to be run first in the post-processing before cluster 0 is created.');
    end
    %% Determine the channel where each spike had that largest amplitude (i.e., the primary) and determine the template amplitude of each cluster
    whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
    whiteningMatrixInv = whiteningMatrix^-1;
    % here we compute the amplitude of every template...
    % unwhiten all the templates
    tempsUnW = zeros(size(templates));
    for t = 1:size(templates,1)
        tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
    end
    % The amplitude on each channel is the positive peak minus the negative
    tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));
    % The template amplitude is the amplitude of its largest channel
    [tempAmpsUnscaled,template_primary] = max(tempChanAmps,[],2);
    %without undoing the whitening
    %template_amplitude = squeeze(max(templates, [], 2) - min(templates, [], 2));
    %[~, template_primary] = max(template_amplitude, [], 2); 
    template_primary = cast(template_primary, class(spike_templates));

    %% Number of samples in the overlap
    n_samples_overlap = round(P.overlap_s * rez.ops.fs);
    n_samples_overlap = cast(n_samples_overlap, class(spike_times));
    %% Distance between each channel
    chan_dist = ((rez.xcoords - rez.xcoords').^2 + (rez.ycoords - rez.ycoords').^2).^0.5;
    ok_count=0;
    suspicious_count=0;
    spike_idx = [1:length(spike_times)]';
    keep = true(size(spike_idx));
    remove_idx=[];
    reference_idx=[];
    flagged=false(length(templates));
    count=1;
    ccg_checked=false(length(templates));
    for i=length(templates):-1:1
        spikeIdxByTemplate{i} = spike_templates==i;
        nspks(i) = sum(spikeIdxByTemplate{i});
        timesByTemplate{i} = spike_times(spikeIdxByTemplate{i});     
    end
    while 1
        this_iter_removed=false;            
        for i=1:length(templates) 
            if nspks(i)
                i_ccg_times = spike_times(spikeIdxByTemplate{i} &keep);
                if isempty(i_ccg_times)
                    continue
                end
            else
                continue
            end
            for j=1:length(templates)           
                if nspks(j) && i<j && ~flagged(i,j) && chan_dist(template_primary(i),template_primary(j))<P.channel_separation_um
                    if ~any(spikeIdxByTemplate{j} & keep)
                        continue
                    end
                    if ~ccg_checked(i,j)
                        [K, Qi, Q00, Q01, Ri, R_central] = ccg(i_ccg_times,spike_times(spikeIdxByTemplate{j} & keep), nbins, round(P.overlap_s*rez.ops.fs*2));
                        ccg_checked(i,j) = true;
                        Rstat=Ri(1);                        
                        % remove duplicates if all of the following criteria are satisfied
                        % 1. the central bin of the CCG is above chance (+/-0.6ms) by default
                        % 2. the adjacent bins, excluding the central bin, are not significantly above chance
                        % 3. the central bin is higher than the adjacent bins
                        if R_central>0.95 && (isnan(Rstat) || Rstat<0.95) && K(nbins+1)>max(K(nbins+[0 2])) % Ri is NaN if there are no spikes in this range of the ccg
                           flagged(i,j)=true;
                           suspicious_count = suspicious_count+1;
                           Rc(suspicious_count) = R_central;
                           Ro(suspicious_count) = Rstat;
                           reference_cluster_idx= double(tempAmpsUnscaled(j)>tempAmpsUnscaled(i));
                           if reference_cluster_idx==0
                               reference_cluster=i;
                               remove_cluster=j;
                           else
                               reference_cluster=j;
                               remove_cluster=i;
                           end 
                           [remove_ij,reference_ij] = remove_duplicates_from_pair(timesByTemplate{i},timesByTemplate{j},reference_cluster_idx,n_samples_overlap,P.verbose);  
                           if ~isempty(remove_ij)
                               ccg_checked(remove_cluster,:)=false;
                               ccg_checked(:,remove_cluster)=false;                               
                               this_iter_removed=true;
                           end
                           if P.verbose
                                fprintf('%g & %g: Removed %g spikes from cluster %g.\n',i,j,length(remove_ij),remove_cluster);
                           end
                           if reference_cluster_idx==1
                               if any(abs(timesByTemplate{i}(remove_ij)-timesByTemplate{j}(reference_ij))>n_samples_overlap)  
                                 error('Something''s gotten confused.')                    
                               else
                                   curr_spike_idx = spike_idx(spike_templates==i);
                                   remove_idx = [remove_idx; curr_spike_idx(remove_ij)];
                                   curr_spike_idx = spike_idx(spike_templates==j);
                                   reference_idx = [reference_idx; curr_spike_idx(reference_ij)];                           
                               end
                           else
                               if any(abs(timesByTemplate{i}(reference_ij)-timesByTemplate{j}(remove_ij))>n_samples_overlap)
                                 error('Something''s gotten confused.')                    
                               else
                                   curr_spike_idx = spike_idx(spike_templates==j);                           
                                   remove_idx = [remove_idx; curr_spike_idx(remove_ij)];
                                   curr_spike_idx = spike_idx(spike_templates==i);
                                   reference_idx = [reference_idx; curr_spike_idx(reference_ij)];                                  
                               end
                           end
                           keep=~ismember(spike_idx,remove_idx);                                       
                           if length(remove_idx)~=length(reference_idx)
                               error('g');
                           end
                        else
                           ok_count = ok_count+1;
                           Rc_ok(ok_count) = R_central;
                           Ro_ok(ok_count) = Rstat;                    
                        end
                    else
                        p=1;
                    end
                end
            end
        end
        if ~this_iter_removed
            break
        end
        count=count+1;     
        if P.verbose
            fprintf('Iteration %g: Recalculating ccgs for unflagged pairs with spikes removed since last ccg check.\n',count);
        end
    end    
    [remove_idx,uidx] = unique(remove_idx);
    reference_idx = reference_idx(uidx);
    if P.verbose
        figure;h(1)=scatter(1-Rc_ok+normrnd(0,0.01,1,length(Rc_ok)),1-Ro_ok+normrnd(0,0.01,1,length(Ro_ok)));hold on;
        h(1).DisplayName = 'Pairs not flagged for duplicate removal';                    
        if exist('Rc','var')
            h(2)=scatter(1-Rc+normrnd(0,0.01,1,length(Rc)),1-Ro+normrnd(0,0.01,1,length(Ro)),'r');
            h(2).DisplayName = 'Pairs flagged for duplicate removal';
        end
        legend(h);        
        xlabel('pvalue for central bin of ccg');
        ylabel('pvalue for adjacent bin');
    end
    rez = remove_spikes(rez,ismember(spike_idx,remove_idx),'duplicate','reference_time',spike_times(reference_idx),'reference_cluster',spike_templates(reference_idx));
end

function [remove_idx,reference_idx] = remove_duplicates_from_pair(times1,times2,reference_cluster,n_samples_overlap,verbose)  
    % identify those spike times in a paired set of times that are closer
    % together than n_samples_overlap.
    % use flag reference_cluster (0 or 1) to pick which of two neighboring
    % times to keep (the reference) and which to remove. 
    % algorithm only compares those times which are adjacent in the
    % temporal sequence (using diff) and does this repeatedly, removing
    % times from the list until no sequential times are within
    % n_samples_overlap. This is much more efficient than finding the times
    % to remove by doing the matrix operation times1-times2'
    n_duplicates=1;
    n1 = numel(times1);
    n2 = numel(times2);
    times1=times1(:);
    times2=times2(:);
    spike_idx = [1:n1+n2]';
    remove_idx = [];
    diff_order=1;
    reference_idx = [];
    [spike_times,sort_idx] = sort([times1;times2]);
    spike_cluster =[zeros(n1,1);ones(n2,1)];    
    spike_cluster = spike_cluster(sort_idx);
    while 1
        if n_duplicates==0
            diff_order=diff_order+1;
        end
        keep_idx = ~ismember(spike_idx,remove_idx);
        current_spike_idx = spike_idx(keep_idx);
        current_spike_cluster = spike_cluster(keep_idx);
        current_spike_times = spike_times(keep_idx); 
        isis=current_spike_times(1+diff_order:end) - current_spike_times(1:end-diff_order);
        simultaneous = isis<n_samples_overlap;
        if ~any(simultaneous)
            if any(same_cluster) && verbose 
                fprintf('Found %g duplicates in this pair that belong to the same cluster. They won''t be removed.\n',sum(same_cluster));
            end            
            break
        end        
        first_duplicate = find(simultaneous); % indexes the first member of the pair
        first_is_bigger =  current_spike_cluster(first_duplicate)==reference_cluster;
        same_cluster = current_spike_cluster(first_duplicate)==current_spike_cluster(first_duplicate+diff_order);
        first_duplicate = first_duplicate(~same_cluster,:);
        first_is_bigger =  first_is_bigger(~same_cluster,:);    
        n_duplicates = length(first_duplicate);
        if ~isempty(first_duplicate)
            remove_idx = [remove_idx ; current_spike_idx([first_duplicate(~first_is_bigger);(first_duplicate(first_is_bigger)+diff_order)])];
            reference_idx = [reference_idx ; current_spike_idx([(first_duplicate(~first_is_bigger)+diff_order);first_duplicate(first_is_bigger)])];
        end
    end
    [remove_idx,idx] = unique(remove_idx);
    reference_idx = reference_idx(idx);    
    [~,unsort_idx] = sort(sort_idx);
    remove_idx = arrayfun(@(x)find(unsort_idx==x,1),remove_idx);
    reference_idx = arrayfun(@(x)find(unsort_idx==x,1),reference_idx);
    if reference_cluster==0
        remove_idx = remove_idx-n1;
    else
        reference_idx = reference_idx-n1;
    end 
end