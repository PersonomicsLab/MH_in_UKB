clear all; clc

X = rand(1000,10);
X = X - mean(X); % Use X = bsxfun(@minus,X,mean(X)) if you have an older version of MATLAB
[coeff,score,latent,~,explained] = pca(X);


[u,s,v] = nets_svds(cov(X),size(X,2));
explained_svds = diag(s)/sum(diag(s));

C = (X(:,1:5)' * X(:,1:5)) + (X(:,6:10)' * X(:,6:10));
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
explained_eig = diag(D)/sum(diag(D));