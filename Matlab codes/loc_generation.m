function [ loc_res ] = loc_generation( N_db, M , num_loc)
%loc_generation Generates valid query locations in the database. 
%   @Params :  N - length of the database
%              M - length of the query
%              num_loc - number of locations where the query needs to match

loc_res =[];
loc = randi(N_db,[5*num_loc 1]);
loc_ascec=sort(loc);
loc_shift =circshift(loc_ascec,-1);
n_locs = nnz((loc_shift - loc_ascec) > M)
loc_new= loc_ascec(find((loc_shift - loc_ascec) > M));

if(n_locs > num_loc)
    
  loc_res = (loc_new(randperm(n_locs,num_loc)))';
end

end

