function ccg_data = test_ccgs(times1,times2,nbins,dt,ccg_list)

    if nargin<5
        ccg_list = {'merge','unit1','unit2','xc'};
    end

    if isempty(setxor(times1,times2)) && any(ismember(ccg_list,{'merge','xc'}))
        warning('I''m going to compute xcorrs and merge corr even though both lists of times are identical. Highly wasteful.');
    end
    

    ccg_fun = @(t1,t2)test_ccg_internal(t1,t2,nbins,dt);

    %% first get stats of the individual clusters
    if ismember('unit1',ccg_list)
        ccg_data.unit1 = ccg_fun(times1,times1);
    end
    if ismember('unit2',ccg_list)    
        ccg_data.unit2 = ccg_fun(times2,times2);
    end
    
    %% then the cross-correlogram
    if ismember('xc',ccg_list)    
        ccg_data.xc = ccg_fun(times1,times2);
    end
    
    %% then the merged cluster
    if ismember('merge',ccg_list)
        alltimes = union(times1,times2);
        ccg_data.merge = ccg_fun(alltimes,alltimes);
    end   
end


function ccg_data = test_ccg_internal(times1,times2,nbins,dt)

    [K, Qi, Q00, Q01, rir, Rcentral] = ccg(times1, times2, nbins, dt); % compute the cross-correlogram between spikes in the putative new clusters
    Q = min(Qi/max(Q00, Q01)); % refractoriness metric 1
    R_min = min(rir); % refractoriness metric 2
    R_max = max(rir);
    
    ccg_data = struct('K',K,'Qi',Qi,'Q00',Q00,'Q01',Q01,'rir',rir,'Q',Q,'Rmin',R_min,'Rmax',R_max,'Rcentral',Rcentral);
    
end