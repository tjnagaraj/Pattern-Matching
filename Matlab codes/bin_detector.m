function [ op, value, pos , pos_flag] = bin_detector( y_obs,stage_num,bin_num)
%ZEROTON_DETECTOR.M    Funtion to determine a checknode as a zeroton or not
% @params   
%   Input:
%           y_obs   - observation vector  
%   Output:
%           op = indicataes if the observation is 
            %zeroton (op=0)
            %singleton (op=1)
            %doubleton (op=2)
            %multiton (op=3)
            
global num_shifts M Cconn N eta
value=-100;
pos=-100;   
pos_flag=0;
y_obs = y_obs/(2900);
            
if (y_obs(1,1) < (1-2*eta)*0.5)
   op =0 ;   % zeroton
elseif(y_obs(1,1) < (3-4*eta)*0.5)
    op =1;   %singleton
    pos_flag=1;
% position identification    
    inds=Cconn{stage_num+1,bin_num+1};
    vq = [];
%     start1 = tic;
    for ind=inds 
      aq=steering_vector(ind); 
      vq= [vq (abs(conj(aq)*y_obs))];
    end
%     toc(start1)   
    energy_thresh = (0.8*num_shifts);
    thresh_inds = find(vq > energy_thresh);
    n_elem = numel(thresh_inds);
    
    if( n_elem == 0)
         op = 0;
    else
    pos1= inds(thresh_inds);
    vq_sub = vq(thresh_inds);
%     value = (1-eta)*M;
      value = 2879;
        if(n_elem == 1)
        pos=pos1;    
        else
        [~ , ind] = max(vq_sub);
        pos = pos1(ind);    
        end
    end        
elseif y_obs(1,1) < (5-6*eta)/2
    op=2;    % doubleton
else
    op=3;    %multiton
end
end

