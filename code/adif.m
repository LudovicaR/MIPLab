function c=adif(a,b)
% INPUT: a and b (two phases from the BOLD signal. 
% OUTPUT: c, the difference of the two phases (in interval [-pi, pi]).
% This function compute the absolute differences between two angles. 
% The resulting angle phase is bounded to pi.

 if abs(a-b)>pi
  c=2*pi-abs(a-b);
 else
  c=abs(a-b);
 end
