# create motif tensor
function createtensor3(A,B,C)
    I = zeros(Int32,1,3);
    V = [];
    dim = [size(A,1) size(B,1) size(B,2)];
    # loop over nodes in C which participate in motif
    K, _ = findnz(diag(C*A*B));
    for ix = K
        Ival, Jval, val = findnz(A.*(B[:,ix]*C[ix,:]')');
        I = [I; [Ival Jval repeat([ix],length(val))]];
        V = [V; val];
    end
    df = DataFrame([I[2:end,:] V],["I","J","K","V"]);
    return df
end