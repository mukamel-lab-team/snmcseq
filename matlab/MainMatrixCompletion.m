% Demo for matrix completion
% Create a random matrix of low rank, remove values and complete.
%



% A = randn(200,100)*randn(100,200); % 200x200 of rank 100
% B = rand(size(A))<0.9; % remove 10% of the entries

f = '/cndd/fangming/integration/processed/snRNA_logtpm.tsv';
AL = dlmread(f, '\t', 1, 1);

[stds, index] = sort(std(AL'));
top_index = index(end-4999:end);

A = AL(top_index,:);
B = double(A>0.01);
nnz(B(:))/numel(B)

gamma = 1;
fprintf('Completion using nuclear norm minimization... \n');

save('MatrixCOmpletionData.mat','As','B')


k = 50;
profile on
[CompletedMat] = MatrixCompletion(A.*B, B, gamma, k);
profile report 


fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(A.*B)));
fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(CompletedMat)));
fprintf('MSE on known entries: %g \n',NormDiff(CompletedMat,A));

