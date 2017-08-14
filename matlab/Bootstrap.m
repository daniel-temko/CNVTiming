function [ TEstArray, t1EstArray  ] = JointBootstrap( data, t1ColName, TColName, diploidLength, iterations )

[alphaArray, betaArray, normalMutationsArray] = GenerateMutations(data, t1ColName, TColName, diploidLength, iterations);

for i = 1:iterations
    dummyData = data;
    dummyData.alpha = alphaArray(:,i);
    dummyData.beta = betaArray(:,i);
    dummyNormalMutations = normalMutationsArray(1,i);
    fprintf('Bootstrap Iteration %d...\n',i)
    [ TEst, t1Est ] = EstimateTimings( dummyData, diploidLength, dummyNormalMutations );
    TEstArray(1,i) = TEst;
    t1EstArray(:,i) = t1Est';
end

end

