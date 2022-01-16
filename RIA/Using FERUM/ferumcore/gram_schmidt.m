function r = gram_schmidt(V)

% we must define the size of the matrix R, it depends on the number of random variables
n = length(V);

% Define matrix Ro
r0 = eye(n);
r0(n,:)= V;

r(n,:)=V;

% After it should be interesting that A correspond to the vector alpha resulting of FORM
% For that purpose, one should write:
% alpha = formresults.alpha;
% r = eye(n);
% r(n,:)= alpha;

% computation of the rows of R
for k = n-1: -1: 1
   matrice = zeros(n);
   for j = k+1: n
      matrice(j,:)=((r(j,:)*r0(k,:)')/(r(j,:)*r(j,:)'))*r(j,:);
   end
   sum_j = sum(matrice);
   r(k,:)=r0(k,:)-sum_j;
end

% normalization of each rows separately
for k = 1: n-1
   r(k,:) = r(k,:)/norm(r(k,:));
end

   