"""
Randomized Singular Value Decomposition (SVD)

code taken from http://www.databookuw.com/databook.pdf, chapter 1.8
"""
function rsvd(X,r,q=1)
    # sample column space of X with P matrix
    ny = size(X,2)
    P  = randn(ny,r)
    Z  = X*P
    # power iterations
    for k = 1:q
        Z = X*(transpose(X)*Z)
    end
    Q, R = qr(Z)
    # compute SVD on projected Y=Q'*X
    Y        = transpose(Matrix(Q))*X   # Matrix() required to make Q thin https://fncbook.github.io/fnc/leastsq/demos/qr-qrfact.html
    UY, S, V = svd(Y)
    U        = Q*UY
    return U, S, V
end
