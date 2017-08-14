% analyse each set's CNV data

%%%%%%%%%%%%%%%% UPDATE AS APPLICABLE %%%%%%%%%%%%%%%%%%%
root = '/[PATH]/cnv_model/';
setNames = {'Set.01'}
bootstrapIterations = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(setNames)
    set = setNames{j}
    
    % make the directory for results outputs
    dir = strcat(root, 'results/', set);
    mkdir(dir);
    
    % read in mutation information after
    file = strcat(root, 'input_data/', set, '/', set, '.matlab.input.txt');
    formatSpec = '%s%s%f%f%f%f%f%f';
    data = readtable(file,'Delimiter',',', 'Format', formatSpec);
    points = length(data.alpha);
    % read in total length of the diploid regions and number of mutations
    % in these regions
    file = strcat(root, 'input_data/', set, '/', set, '.normal.stats.txt');
    formatSpec = '%f%f';
    diploidStats = readtable(file,'Delimiter',',', 'Format', formatSpec);
    

    % add case numbers to the data, depending on the number of copies of the
    % A allele and B allele
    data.caseNum = repmat(0,points,1);
    for i = 1:points
            data.caseNum(i) = InferCase(data.a2(i), data.b2(i));
    end
    
    [ TEst, t1Est ] = EstimateTimings(data, diploidStats.length, diploidStats.mutations);
    data.t1_est = t1Est'; 
    data.T_est = repmat(double(TEst),points,1);
    
    % bootstrapping
    [ TEstArray, t1EstArray  ] = Bootstrap( data, 't1_est', 'T_est', diploidStats.length, bootstrapIterations );
    clear t1Diff
    for i = 1:bootstrapIterations
        t1Diff(:,i) = t1EstArray(:,i) - data.t1_est;
    end
    
    % plot bootstrap results
    h = figure;
    boxplot(t1Diff')
    title(set)
    ylabel('t1 Error')
    file = strcat(root, 'results/', set, '/', set, '.t1_errors.pdf');
    print(h,file,'-dpdf')
    
    t1SquareErrors = t1Diff .* t1Diff;
    clear MSE
    for i = 1:points
        MSE(i) = mean(t1SquareErrors(i,:));
    end
    data.t1_MSE = MSE';
    TDiff = TEstArray - data.T_est(1);
    
    % plot bootstrap results
    h = figure;
    boxplot(double(TDiff))
    title(set)
    ylabel('T Error')
    file = strcat(root, 'results/', set, '/', set, '.T_errors.pdf');
    print(h,file,'-dpdf')
    
    TSquareErrors = TDiff .* TDiff;
    TMSE = mean(TSquareErrors);
    data.T_MSE = repmat(double(TMSE),points,1);
      
    % write results for analysis
    output = data;
    file = strcat(root, 'results/', set, '/', set, '.matlab.output.txt');
    writetable(output,file,'Delimiter',',');
end