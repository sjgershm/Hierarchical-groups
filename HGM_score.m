function score = HGM_score(X,Z,parents,opts)
    
    % Bayesian score for tree structure.
    %
    % USAGE: score = HGM_score(X,Z,opts)
    
    N = size(Z,1);
    K = opts.tau2*(Z*Z') + opts.s2*eye(N);
    lik = -0.5*trace(X'*(K\X)) - sum(log(diag(chol(K))));
    score = nCRP_treeprob(Z,parents,opts.alpha) + lik;
    
end

function logp = nCRP_treeprob(Z,parents,alpha)
    
    % Evaluate probability of a tree under the nCRP.
    
    logp = 0;
    for k = 1:size(Z,2)
        ix = parents==k;
        if any(ix)
            m = sum(Z(:,ix));
            m(m==0) = [];
            logp = logp + gammaln(alpha) + length(m)*log(alpha) - gammaln(sum(m)+alpha) + sum(gammaln(m));
        end
    end
    
end