function f = rat_approx1(varsigma,pfe,n_max)

    L = size(varsigma,2);
    N = size(pfe.res);
    r = pfe.residues;
    pi = pfe.poles;
    
    pi_diff = pi(2:end)-pi(1);    
    f1 = r(1)/(varsigma-pi(1));
    
    ratio1 = (r(2:end)./pi_diff).';
    ratio2 =  varsigma-pi(1);
    
    % Replicate vectors to create tensor of terms
    ratio1_rep = repmat(ratio1,[1,L]);
    ratio2_rep = zeros(1,1,L);
    ratio2_rep(1,1,:) = ratio2;
    ratio2_rep = repmat(ratio2_rep,[N-1,n_max+1,1]);
    pi_diff_rep = repmat(pi_diff,n_max+1,L);
    n_power = repmat([0:n_max],[N-1,1,L]);
    
    % Tensor of terms in Sn (summation over n)
    ToT = (ratio2_rep./pi_diff_rep).^n_power; 
    
    % Summation Sn, i.e. summation over n
    S2 = sum(ToT,2);
    if size(S2,1)==1 
        S2 = squeeze(S2).';
    else
        S2 = squeeze(S2);
    end
    
    % Summation Si, i.e. summation over i
    S1 = sum(ratio1_rep.*S2,1);
    
    % Series expansion approximation of f
    f = f1 - S1;
end

function f2 = rat_approx2(varsigma,pfe,n_max)

    L = size(varsigma,2);
    N = size(pfe.res);
    r = pfe.residues;
    pi = pfe.poles;
    
    % Replicate vectors to create tensor of terms
    n_power = repmat([0:n_max],[N,1]);
    r_rep = repmat(r,[1,n+1]);
    pi_rep = repmat(pi,[1,N+1]);
    ratio1 = r_rep./(pi_rep.^(n_power+1));  %r_{i}/(pi_{i}^{n+1})
    
    varsigma_rep = zeros(1,1,L);
    varsigma_rep(1,1,:) = varsigma;
    varsigma_rep = repmat(varsigma_rep,[1,n_max+1,1]);
    n_power = repmat([0:n_max],[1,1,L]);
    
    % Tensor of monomials
    ToM = varsigma_rep.^n_power;
    
    % Terms Gn (summation over i for every n)
    Fn = sum(ratio1,1);
    Fn_rep = repmat(Gn,[1,1,L]);
    
    % Tensor of terms
    ToT = sum(Fn_rep.*ToM);
    
    % Series expansion approximation of f
    f2 = -squeeze(sum(ToT,2));
    f2 = f2.';    % return f as a row vector
end