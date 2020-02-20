function Y = cmdscale(D,p)

% Same as matlab stats function but without the error checking and enforces
% a proper dissimilarity matrix
D=double(D);
D=D+D';
yl=range(D(:));
D=(D-min(D(:)))./yl;
D = D- diag(diag(D));
n=length(D);

if nargin<2
    p=n;
end

D = D.*D; % square elements of D
B = bsxfun(@plus, bsxfun(@minus, bsxfun(@minus, D, sum(D,1)/n),sum(D,2)/n), sum(D(:))/(n^2))*(-0.5);

if (p == n)
    % compute full eigen-decomposition
    [V,E] = eig((B+B')./2); % guard against spurious complex e-vals from roundoff
else
    % compute only p eigenvectors and eigenvalues.
    [V,E] = eigs((B+B')./2,p,'LA'); % guard against spurious complex e-vals from roundoff
end
[e,i] = sort(diag(E)); e = flipud(e); i = flipud(i); % sort descending

% keep only positive e-vals (beyond roundoff)
keep = find(e > max(abs(e)) * eps(class(e))^(3/4));
if isempty(keep)
    Y = zeros(n,1);
else
    % The following line does the same thing as: Y = V(:,i(keep)) * diag(sqrt(e(keep)));
    Y = bsxfun(@times, V(:,i(keep)), sqrt(e(keep))');
end

% Enforce a sign convention on the solution -- the largest element
% in each coordinate will have a positive sign.
[~,maxind] = max(abs(Y),[],1);
d = size(Y,2);
colsign = sign(Y(maxind + (0:n:(d-1)*n)));
Y = bsxfun(@times,Y,colsign);

end