function vec = evenRDistr(N,M) % evenly distributes N items among M groups, and randomly assigns the rest 
    vec = floor(N/M) * ones(1,M);
    excess = rem(N,M);
    if excess > 0
        extras = randperm(M,excess);
        vec(extras) = vec(extras)+1;
    end
end