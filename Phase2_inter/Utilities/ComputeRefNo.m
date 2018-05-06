function [ref1, ref2] = ComputeRefNo(cur)


if cur == 1
    ref1 = 2;
    ref2 = 3;
    
    
elseif cur == 9
    ref1 = 7;
    ref2 = 8;
else
    
    ref1 = cur - 1;
    ref2 = cur + 1;
    
end
