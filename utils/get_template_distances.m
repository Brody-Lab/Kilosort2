function template_dist = get_template_distances(rez)

rez.U=gather(rez.U);
rez.W = gather(rez.W);
templates = zeros(rez.ops.Nchan, size(rez.W,1), size(rez.W,2), 'single');
for iNN = 1:size(templates,3)
   templates(:,:,iNN) = squeeze(rez.U(:,iNN,:)) * squeeze(rez.W(:,iNN,:))';
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
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
[~,template_primary] = max(tempChanAmps,[],2);
%without undoing the whitening
%template_amplitude = squeeze(max(templates, [], 2) - min(templates, [], 2));
%[~, template_primary] = max(template_amplitude, [], 2); 


chan_dist = ((rez.xcoords - rez.xcoords').^2 + (rez.ycoords - rez.ycoords').^2).^0.5;



template_dist = chan_dist(template_primary,template_primary);

end