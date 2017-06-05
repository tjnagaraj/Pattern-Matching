function [value , pos] = singleton_estimator( y_obs,stage_num,bin_num)
%SINGLETON_ESTIMATOR 
%   @params y_obs - checknode observation
%           thres - Energy threshold
%           stage_num - denotes the stage number of the checknode
%           op - '1' if a singleton  '0' if NOT


global num_shifts M Cconn eta;

value=[];
pos=[];
inds=Cconn{stage_num+1,bin_num+1};
vq = [];
for ind=inds 
      aq=steering_vector(ind); 
      vq= [vq (abs(conj(aq)*y_obs))/M*num_shifts];
end

    value = (1-eta)*M;
    pos = inds(find(vq == max(vq))); 
end

