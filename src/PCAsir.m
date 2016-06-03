function [ Wpca,meanMat,k ] = PCAsir( X,per )
% Variance to maintain
%per = 0.99;

% Get mean, and create matrix to optimize calculations
meanVals    = mean(X,1);
stdVals     = std(X,1)+eps;

meanMat     = ones(size(X))*diag(meanVals);
stdMat      = ones(size(X))*diag(stdVals);

% Zero-mean and unit variance
Xm = (X - meanMat);

% Co-variance
sigma = cov(Xm);
 
% Perform SVD
[U, E, V] = svd(sigma);

% Get the sum of the eigenvalues
tra = trace(E);

% Find the smallest number of eigenvalues to maintain specified variance
sm = 0;
for k=1:length(E)
        sm = sm + E(k,k);
        if( sm/tra > per )
            break
        end        
end

% PCA matrix 
Wpca = V(:,1:k);


end

