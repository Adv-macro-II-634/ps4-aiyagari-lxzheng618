
clear,clc
tic

% PARAMETERS

alpha = 1/3; 
beta = 0.99; 
sigma = 2; 

% ASSET
delta = 0.025;
rho = 0.5; 
sigma_e = 0.2;
a_l = 0;

n = 5;
[z_grid, P] = rouwenhorst(rho,sigma_e,n); 

a_l = 0;
a_h = 100;
num_a = 500;
a = linspace(a_l,a_h,num_a);

K_min = 1;
K_max = 100;
K_guess = (K_min + K_max) / 2;

K_tol = 1;
while abs(K_tol)>=0.01;
    
    K_guess = (K_min+K_max)/2;
    r = alpha*K_guess^(alpha-1)*N^(1-alpha)+(1-delta); 
    w = (1-alpha)*K_guess^alpha*N^(-alpha); 
    
   	cons = bsxfun(@minus, r*a, a');
    cons = bsxfun(@plus, cons, permute(z_grid*w, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); 
    ret(cons<0)=-Inf;
    
% INITIAL VALUE FUNCTION GUESS
    
    v_guess = zeros(n, num_a);
    
% VALUE FUNCTION ITERATION
v_tol = 1;
while v_tol > 0.0001;
  % CONSTRUCT TOTAL RETURN FUNCTION
  v_mat = ret + beta * ...
       repmat(permute(z_prob * v_guess, [3 2 1]), [num_a 1 1]);
  
  % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
   [vfn, pol_indx] = max(v_mat, [], 2);
   
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(v_guess(:) - vfn(:)));
 
   v_guess = vfn; 
        
end;

% KEEP DECSISION RULE  
pol_indx = permute(pol_indx, [3 1 2]);
pol_fn = a(pol_indx);
    
% SET UP INITITAL DISTRIBUTION
Mu = zeros(n,num_a);
Mu(1, 4) = 1; 

% ITERATE OVER DISTRIBUTIONS
% loop over all non-zeros states
mu_tol = 1;
while mu_tol > 1e-08
    [emp_ind, a_ind] = find(Mu > 0); 
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (P(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
    end
    
    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew ;
    
end
   
   % CHECK AGGREGATE DEMAND
   aggK = sum( pol_fn(:) .* Mu(:) );
if aggK > 0;
  K_min = K_guess;
end
if aggK < 0;
  K_max=K_guess;
end

end