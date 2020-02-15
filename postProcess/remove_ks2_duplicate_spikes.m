% REMOVE_KS2_DUPLICATE_SPIKES Remove spikes that were counted multiple
% times due to variability in spike waveform that is not yet addressed by
% Kilosort2's template matching algorithm. 
% 
% Any group of spikes within OVERLAP_S and whose peak deflections are within
% CHANNEL_SEPARATION_UM are removed except for the spike with the largest
% amplitude.
%
% The Kilosort2 outputs are modified and saved, and each original output
% file is copied and renamed "*_original.npy." The original files can be
% restored by calling RESTORE_ORIGINAL_KS2_FILEm
%
% The indices of each removed
% spike and of the corresponding overlapping spike that was preserved
% (i.e., the reference spike) is saved in "removed_spikemat."
%
%=INPUT
%
%   rez structure
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   overlap_s
%       the time interval, in second, within which a sequence of spikes are
%       vetted for duplicates
%
%   channel_separation_um
%       When the primay channels of two spikes are within this distance, in
%       microns, then the two spikes are vetted for duplicate The primary
%       channel of each spike is where its large peak-to-peak amplitude
%       occured. This is not the raw peak-to-peak amplitude of the spike,
%       but tather it is the amplitude of a template fitted to the average
%       whitened and filtered waveform of the cluster to which that spike
%       was assigned.
%
%=EXAMPLE
%
%   >> rez = remove_ks2_duplicate_spikes(rez)
function rez = remove_ks2_duplicate_spikes(rez, varargin)
input_parser = inputParser;
addParameter(input_parser, 'overlap_s', 0.000167, @(x) (isnumeric(x)))
addParameter(input_parser, 'channel_separation_um', 50, @(x) (ischar(x)))
parse(input_parser, varargin{:});
P = input_parser.Results;

spike_times = uint64(rez.st3(:,1));
amplitudes = rez.st3(:,3);
[spike_templates,spike_clusters] = deal(uint32(rez.st3(:,2)));
template_features = rez.cProj;

rez.U=gather(rez.U);
rez.W = gather(rez.W);
templates = zeros(rez.ops.Nchan, size(rez.W,1), size(rez.W,2), 'single');
for iNN = 1:size(templates,3)
   templates(:,:,iNN) = squeeze(rez.U(:,iNN,:)) * squeeze(rez.W(:,iNN,:))';
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels

%% Make sure that the spike times are sorted
if ~issorted(spike_times)
    [spike_times, I] = sort(spike_times);
    amplitudes = amplitudes(I);
    spike_clusters = spike_clusters(I);
    spike_templates = spike_templates(I);
    template_features =  template_features(I,:);
end
%% deal with cluster 0
if any(spike_templates==0)
    error('Currently this function can''t deal with existence of cluster 0. Should be OK since it ought to be run first in the post-processing.');
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
spike_primary = template_primary(spike_templates);

%% Number of samples in the overlap
n_samples_overlap = round(P.overlap_s * rez.ops.fs);
n_samples_overlap = cast(n_samples_overlap, class(spike_times));
%% Distance between each channel
chan_dist = ((rez.xcoords - rez.xcoords').^2 + (rez.ycoords - rez.ycoords').^2).^0.5;
%imagesc(chan_dist)

n_spikes=numel(spike_times);
n_duplicates=1;  % set to 1 to initialize while loop
count=0;
remove_idx = [];
reference_idx = [];
spike_idx = [1:n_spikes]';
current_spike_times = spike_times;
current_spike_idx = spike_idx;
current_primaries = spike_primary;
% only check nearest temporal neighbors in the list of spikes times.
% but go recursively until no nearest neighbors are left that are both within the overlap
% period and sufficiently nearby. 
% this means only ever computing a vector operation (i.e. diff(spike_times))
% rather than a matrix one (i.e. spike_times - spike_times').
while n_duplicates>0
    count=count+1;
    keep_idx = ~ismember(spike_idx,remove_idx);
    current_spike_idx = spike_idx(keep_idx);
    current_spike_times = spike_times(keep_idx); 
    current_primaries = spike_primary(keep_idx);
    isis=diff(current_spike_times);
    simultaneous = isis<n_samples_overlap;
    nearby = chan_dist(sub2ind(size(chan_dist),current_primaries(1:end-1),current_primaries(2:end)))<P.channel_separation_um;
    first_duplicate = find(simultaneous & nearby); % indexes the first member of the pair
    n_duplicates = length(first_duplicate);
    if ~isempty(first_duplicate)
        fprintf('On iteration %g, %g duplicate spike pairs were identified.\n',count,n_duplicates);
        first_is_bigger =  diff(tempAmpsUnscaled(current_primaries([first_duplicate first_duplicate(:)+1])),[],2)<=0;
        remove_idx = [remove_idx ; current_spike_idx([first_duplicate(~first_is_bigger);(first_duplicate(first_is_bigger)+1)])];
        reference_idx = [reference_idx ; current_spike_idx([(first_duplicate(~first_is_bigger)+1);first_duplicate(first_is_bigger)])];
        [remove_idx,idx] = unique(remove_idx);
        reference_idx = reference_idx(idx);
    end
end
rez.reference(remove_idx) = reference_idx;
rez.original_cluster(remove_idx) = rez.st3(remove_idx,2);
rez.st3(remove_idx,2)=0; % set duplicate spikes to belong to cluster 0

rez = recompute_clusters(rez);
