% generates alphaArray and betaArray, each of the columns of these is a possible
% column of alpha and beta values respectively for the input data set.
% These have a column number equal to iterations
function [alphaArray, betaArray, diploidRegionMutationsArray] = GenerateMutations(data, t1ColName, TColName, diploidLength, iterations)

points = length(data.CNV_ID);
cols = fields(data);
tmp = strfind(cols,t1ColName);
fieldName = cols(find(arrayfun(@(x) length(tmp{x}), 1:length(tmp))));
t1 = getfield(data,fieldName{1});
tmp = strfind(cols,TColName);
fieldName = cols(find(arrayfun(@(x) length(tmp{x}), 1:length(tmp))));
T = getfield(data,fieldName{1});

for i = 1:points
    copiedlength = ((data.a2(i) > 1)+(data.b2(i) > 1)) * data.length(i);
    highvafmean = copiedlength * t1(i);
    alphaArray(i, 1:iterations) = poissrnd(highvafmean, 1, iterations);
    earlylowvaflength = ((data.a2(i) == 1) + (data.b2(i) == 1)) * data.length(i);
    latelowvaflength = (data.a2(i) + data.b2(i)) * data.length(i);
    lowvafmean = (earlylowvaflength * t1(i) + latelowvaflength * (T(i) - t1(i)));
    betaArray(i, 1:iterations) = poissrnd(lowvafmean, 1, iterations);
end
normalmean = 2 * diploidLength * T(1);
diploidRegionMutationsArray(1, 1:iterations) = poissrnd(normalmean, 1, iterations);
end