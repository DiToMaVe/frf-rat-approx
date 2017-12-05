function [f2, auxilaries] = rat_approx2(varsigma,pfe,n_max)
    
    % Compute series expansion approximation to a complex rational function
    % expressed as a partial fraction expansion, for |z|<|pi_min|
    % (Shurman;2017) with pi_min=min{pi_i, for all i from 1 to N}.
    
    % Reference: J. Shurman, “Laurent series and singularities,” 2017,
    % Handouts for Mathematics 311: Complex  Analysis, Lecture Notes
    % [Online]. Available: http://people.reed.edu/jerry/311/laurent.pdf,
    % [Accessed:  21  April2017].
    
    L = size(varsigma,2);
    N = size(pfe.residues,2);
    r = pfe.residues.';
    ppi = pfe.poles.';
    
    % Replicate vectors 
    n_power = repmat([0:n_max],[N,1]);
    r_rep = repmat(r,[1,n_max+1]);
    pi_rep = repmat(ppi,[1,n_max+1]);
    
    ratio1 = r_rep./(pi_rep.^(n_power+1));  %r_{i}/(pi_{i}^{n+1})
   
    varsigma_rep = zeros(1,1,L);
    varsigma_rep(1,1,:) = varsigma;
    varsigma_rep = repmat(varsigma_rep,[1,n_max+1,1]);
    n_power = repmat([0:n_max],[1,1,L]);
    
    % Tensor of monomials
    ToM = varsigma_rep.^n_power;
    
    % Terms Gn (summation over i for every n)
    Fn = sum(ratio1,1);
    Fn_rep = repmat(Fn,[1,1,L]);
    
    % Tensor of terms
    ToT = sum(Fn_rep.*ToM);
    
    % Series expansion approximation of f
    f2 = -squeeze(sum(ToT,2));
    f2 = f2.';    % return f as a row vector
    
    auxilaries.Fn = Fn;
end