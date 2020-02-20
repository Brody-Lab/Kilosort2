function rez = find_merges(rez, flag)
% this function merges clusters based on template correlation
% however, a merge is veto-ed if refractory period violations are introduced

if ~isfield(rez,'merges')
    rez.merges=[];
end

ops = rez.ops;
dt = 1/1000;


wPCA = gather(rez.ops.wPCA);

distances_um = get_template_distances(rez);
distances_um = distances_um + eye(length(distances_um)).*1000; % remove
distance_score = 1000-distances_um;
distance_score = distance_score./1000;
n_merged_now=0;
Xsim = rez.simScore; % this is the pairwise similarity score
Nk = size(Xsim,1);
Xsim = Xsim - diag(diag(Xsim)); % remove the diagonal of ones

% sort by firing rate first
nspk = zeros(Nk, 1);
for j = Nk:-1:1
    spks{j} = (rez.st3(:,2)==j); % determine total number of spikes in each neuron
    nspk(j) = sum(spks{j});
end
[~, isort] = sort(nspk); % we traverse the set of neurons in ascending order of firing rates

fprintf('Finding merges...\n')

if ~flag
  % if the flag is off, then no merges are performed
  % this function is then just used to compute cross- and auto- correlograms
   rez.R_CCG = Inf * ones(Nk);
   rez.Q_CCG = Inf * ones(Nk);
   rez.K_CCG = {};
end

for j = 1:Nk
    s1 = rez.st3(spks{isort(j)}, 1)/ops.fs; % find all spikes from this cluster
    if isempty(s1)
        continue
    end
    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts') %this is a check for myself to make sure new cluster are combined correctly into bigger clusters
    end
    

    for k = 1:Nk
        if numel(s1)>nspk(k) || Xsim(isort(j),k)<0.5 || distance_score(isort(j),k)<0.9
            continue
        end

        s2 = rez.st3(spks{k}, 1)/ops.fs; % find the spikes of the pair
        if isempty(s2)
            continue
        end
        
        
        dt = 0.001;
        ccg_data = test_ccgs(s1,s2,500,dt,{'merge','unit1','unit2'});
       
        if flag
            if (ccg_data.merge.Q<=max(0.2,min(ccg_data.unit1.Q,ccg_data.unit2.Q)) && ccg_data.merge.Rmin<=max(0.05,min(ccg_data.unit1.Rmin,ccg_data.unit2.Rmin))) 
                %rez.simScore(isort(j),k)>0.95
            %if (ccg_data.merge.Q<max(ccg_data.unit1.Q,ccg_data.unit2.Q) && ccg_data.merge.Q<max(ccg_data.unit1.Rmin,ccg_data.unit2.Rmin))
            %    test=1;
            %end
                % if you don't make the cluster less refractory than the
                % parent clusters by merging OR the waveform correlation is
                % very high, GO AHEAD.
                % now merge j into i and move on
                rez.st3(spks{isort(j)},2) = k; % simply overwrite all the spikes of neuron j with i (i>j by construction)
                spks{k} = spks{k} | spks{isort(j)}; % update spikes    
                spks{isort(j)} = false(size(spks{isort(j)}));                      
                nspk(k) = nspk(k) + nspk(isort(j)); % update number of spikes for cluster i
                nspk(isort(j)) = 0; % update number of spikes for cluster i                
                n_merges = size(rez.merges,1)+1;
                rez.merges(n_merges,:) = [isort(j) k];
                n_merged_now=n_merged_now+1;
                % YOU REALLY SHOULD MAKE SURE THE PC CHANNELS MATCH HERE
                break; % if a pair is found, we don't need to keep going (we'll revisit this cluster when we get to the merged cluster)
            end
        end
    end
end

fprintf('Found %g merges in %s.\n',n_merged_now,timestr(toc));

if flag
    rez = recompute_clusters(rez);
end