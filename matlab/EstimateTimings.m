function [ TEst, t1Est ] = EstimateTimings( data, diploidLength, diploidRegionMutations )

points = length(data.alpha);

C1.cases = find(data.caseNum == 1);
C2.cases = find(data.caseNum == 2);
C3.cases = find(data.caseNum == 3);
% create symbolic variable for T, the ci, the li, the alphai and the betai
syms T a l alph bet;
tvars = sym(zeros(1, points));
avars = sym(zeros(1, points));
lvars = sym(zeros(1, points));
betavars = sym(zeros(1, points));
for i = 1:points
    tvars(i) = sprintf('t%d', i);
    avars(i) = sprintf('a%d', i);
    lvars(i) = sprintf('l%d', i);
    betavars(i) = sprintf('beta%d', i);
end


% the below code only applies when there are Case III regions (i.e. one
% allele is in two copies and the other is one copy). In this case each case
% III CNV defines a critical value (Tcrit) depending on the parameters of
% this CNV. The equation for T obtained from the likelihood maximisation is
% defined piecewise depending on whether T is less than or greater to each
% of these critical values. In the below the critical values are calculated
% and stored in the vector Tcrit. These critical values are then substituted
% into the likelihood equation to determine which region of the likelihood
% equation contains the zero, and is therefore relevant for the solution.
if (length(C3.cases) > 0)
    clear Tcrit
    breakpoint = (a * bet - alph) / ((a - 1) * l);
    for i = 1:length(C3.cases)
        j = C3.cases(i);
        Tcrit(i) = subs(breakpoint, [a l alph bet], [max(data.a2(j), data.b2(j)) data.length(j) data.alpha(j) data.beta(j)]);
    end

    % test the likelihood at each breakpoint, when a value less than zero is
    % found we conclude solution was in the previous interval
    breaks = [double(sort(Tcrit)), max(double(Tcrit)) + 1];
    i = 0;
    result = 1;
    while ((double(result) > 0) & i <= length(breaks))
        i = i+1;
        Teval = breaks(i);
        if (Teval > 0)
            [ equation, solutions ] = GenerateLikelihoodEquations(data, diploidLength, diploidRegionMutations, Tcrit, Teval);
            result = subs(equation,Teval);
        end
    end
else
    [ equation, solutions ] = GenerateLikelihoodEquations(data, diploidLength, diploidRegionMutations, [], 0);
end

% solve
TEst = vpasolve(equation, [1e-9 1]);

% estimate t values
for i= 1:points
    t1Est(i) = double(subs(solutions(i),TEst));
end

end