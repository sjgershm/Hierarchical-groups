function mu_rel = HGM_infer(X,Z,opts)
    
    % Compute relative value at each node in one individual's tree path.
    %
    % USAGE: mu_rel = HGM_infer(X,Z,opts)
    
    mu_rel = (Z'*Z + (opts.s2/opts.tau2)*eye(size(Z,2)))\(Z'*X);