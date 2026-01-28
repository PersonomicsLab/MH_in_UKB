function [window_corr,window_n,increments] = ...
    time_corr_filter_v2(data1,data2,time,varargin)

rank_corr = @(x,y) corr(x,y,'type','spearman');

short_cut = 14;
long_cut = 93;
def_cutoffs = [min(time(:)), short_cut, long_cut, max(time(:))];

assert(isequal(size(data1),size(data2)),['Input data1, data2, and time '...
    'should have same dimensions.'])
assert(isequal(size(data1),size(time)),['Input data1, data2, and time '...
    'should have same dimensions.'])

p = inputParser;
addParameter(p,'comparison_test',rank_corr)
addParameter(p,'increments',unique(time))
addParameter(p,'cutoffs',def_cutoffs)
parse(p,varargin{:})

cmp_tst = p.Results.comparison_test;
increments = p.Results.increments;
cutoffs = p.Results.cutoffs;


window_corr = NaN(size(increments));
window_n = zeros(size(increments));

for j=1:(length(cutoffs)-1)
    frame_start = cutoffs(j);
    frame_end = cutoffs(j+1);
    
    time_frame = (time >= frame_start) & (time <= frame_end);
    
    tmp_time = time;
    tmp_time(~time_frame) = NaN;
    
    tmp_data1 = data1;
    tmp_data1(~time_frame) = NaN;
    
    tmp_data2 = data2;
    tmp_data2(~time_frame) = NaN;
    
    for i=1:length(increments)
        time_window = (tmp_time <= increments(i));  % assumes 
        
        if any(time_window)
            data1_recent = tmp_data1(time_window);
            data2_recent = tmp_data2(time_window);
            window_n(i) = length(data1_recent);
            if window_n(i) > 3
                window_corr(i) = cmp_tst(data1_recent,data2_recent);
            else
                window_corr(i) = NaN;
            end
        end
    end
end
nan_idx = isnan(window_corr);

window_corr(nan_idx) = [];
window_n(nan_idx) = [];
increments(nan_idx) = [];

end