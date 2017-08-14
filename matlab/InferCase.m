% function to infer case numbers
function [casenum] = InferCase(a2,b2)
a = max(a2, b2);
b = min(a2, b2);
if ((a >= 2) & (b == 0))
    casenum = 1;
elseif (a == b)
    casenum = 2;
elseif ((a >= 2) & (b == 1))
    casenum = 3;
end
end