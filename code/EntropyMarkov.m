function H = EntropyMarkov(P1,P2)

%%% FOR P1
numberStates=size(P1,1);
[V,D] = eig( P1' );    % eigendecomp of P' (P = transition matrix)
[~,ii] = max(diag(D)) ; 
st = V(:,ii); % eigenvector with highest eigenvalue
p1 = abs(st)/sum(abs(st)); % p(i) (stationary!) - normalised

% Remark: you could also get an estimate of p(i) directly from the observations
% p(i) = sum( state == i )/T; 

% Markov entropy:
Hi = zeros(1,numberStates); 
for row = 1:numberStates  
    Hi(row) = -p1(row)*sum( (P1(row,:)+eps).*log2(P1(row,:)+eps) );
end

H1 = sum(Hi)/log2(numberStates); % why this summation and normalisation

%%% FOR P2
numberStates=size(P2,1);
[V,D] = eig( P2' );    % eigendecomp of P' (P=transition matrix)
[~,ii]= max(diag(D)) ; 
st = V(:,ii); 
p2 = abs(st)/sum(abs(st)); % p(i) (stationary!)

% Remark: you could also get an estimate of p(i) directly from the observations
% p(i) = sum( state == i )/T; 

% Markov entropy:
Hi = zeros(1,numberStates); 
for row=1:numberStates  
    Hi(row) = -p2(row)*sum( (P2(row,:)+eps).*log2(P2(row,:)+eps) );
end

H2 = sum(Hi)/log2(numberStates); 

H=abs(H1-H2);