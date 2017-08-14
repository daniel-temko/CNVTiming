% evaluates the likelihood equation on data at the point T=Teval
function [ equation, solutions ] = GenerateLikelihoodEquations(data, diploidLength, diploidRegionMutations, Tcrit, Teval)

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

clear maximum solutions denoms derivatives sumterms

% defining maximum(k,l) likelihood-maximising solutions for the ti in 
% case k.l
maximum(1,1) = (((1 - a) * l * T + alph + bet) - sqrt(((1 - a) * l * T + alph + bet)^2 - 4 * (1 - a) * l * T * alph)) / (2 * (1 - a) * l);
maximum(2,1) = ((2 * (1 - a) * l * T + alph + bet) - sqrt((2 * (1 - a) * l * T + alph + bet)^2 - 8 * (1 - a) * l * T *alph)) / (4 * (1 - a) * l);
maximum(3,1) = (((1 - a) * (a + 1) * l * T + a * (alph + bet)) - sqrt(((1 - a) * (a + 1) * l * T + a * (alph + bet))^2 - 4 * (1 - a) * a * l * T * (a + 1) *alph)) / (2 * (1 - a) * a * l);
maximum(3,2) = T;


subcaseNum = [];
% find the applicable subcase for Case 3 CNV's, depending on whether Teval
% is less than or greater than Tcrit. 
for i = 1:length(C3.cases)
    if(double(double(Teval) <= Tcrit(i)))
       subcaseNum(C3.cases(i)) = 1;
    else
        subcaseNum(C3.cases(i)) = 2;
    end
end
C3.subcases1 = find(subcaseNum == 1);
C3.subcases2 = find(subcaseNum == 2);

%find the solutions for each ti
for i = 1:points
    if data.caseNum(i) == 3
        solutions(i) = subs(maximum(data.caseNum(i), subcaseNum(i)), [a l alph bet], [max(data.a2(i), data.b2(i)) data.length(i) data.alpha(i) data.beta(i)]);
    else
        solutions(i) = subs(maximum(data.caseNum(i),1), [a l alph bet], [max(data.a2(i), data.b2(i)) data.length(i) data.alpha(i) data.beta(i)]);
    end
end

%find the components of the likelihood function for each CNV
for i = 1:points
    if data.caseNum(i) == 3
        if subcaseNum(i) == 1
            denoms(i) = (((max(data.a2(i), data.b2(i)) + 1) * T) - (max(data.a2(i), data.b2(i)) * tvars(i)));
        else
            denoms(i) = T;
        end
    else
        denoms(i) = (T - tvars(i));
    end
end

for i = 1:points
    if data.caseNum(i) == 3
        if subcaseNum(i) == 1
            nums(i) = (max(data.a2(i), data.b2(i)) + 1) * data.beta(i);
        else
            nums(i) = data.alpha(i) + data.beta(i);
        end
    else
        nums(i) = data.beta(i);
    end
end

for i = 1:points
    sumterms(i) = nums(i) / denoms(i);
end

likelihood = sum(sumterms) + diploidRegionMutations/T -(sum(avars([C1.cases]) .* lvars([C1.cases])) + sum(2 * avars([C2.cases]) .* lvars([C2.cases])) + sum((avars([C3.subcases1]) + 1) .* lvars([C3.subcases1])) + sum(2 * lvars([C3.subcases2])) + 2 * diploidLength);

% substitute into the likelihood derivative for T
if size(data.length,1) > 1
    equation = subs(likelihood, [tvars, avars, lvars, betavars], [solutions, max(data.a2, data.b2)', data.length', data.beta']);
else
    equation = subs(likelihood, [tvars, avars, lvars, betavars], [solutions, max(data.a2, data.b2), data.length, data.beta]);
end

end

