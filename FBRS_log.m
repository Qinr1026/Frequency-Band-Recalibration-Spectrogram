
function log_S = FBRS_log(S)

[frame_num, fliter_num] = size(S);

for i=1:frame_num
   for j=1:fliter_num
      log_S(i,j)=log(S(i,j));
   end
end

log_S = log_S';

end