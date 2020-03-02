function c=adif(a,b)
% INPUT: a and b (two phases from the BOLD signal). 
% OUTPUT: c, angle between the two BOLD signals.
% This function compute the absolute differences between two BOLD phases. 
% The resulting angle is bounded to pi.

 if abs(a-b)>pi
  c=2*pi-abs(a-b);
 else
  c=abs(a-b);
 end
