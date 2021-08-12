using LinearAlgebra

function melhor_coord(A,B) #melhor de acordo com a norma 2 (norma euclidiana)
    C=B\A #mínimos quadrados
    return C 
end

function melhor_base(A,C) #melhor de acordo com a norma 2 (norma euclidiana)
    B=C'\A'
    return B'
end

function erro(X,Y)
    return norm(X-Y)
end

function melhor_base_e_coords(A,k)  #alternating least squares #svd,pca
    m,n=size(A)
    B=randn(m,k)
    C=randn(k,n)
    for i=1:20
        C=melhor_coord(A,B)
        B=melhor_base(A,C)
    end
    return B,C
end

function melhor_base_e_coords_posto1(A)
    b,c=melhor_base_e_coords(A,1)
end

function melhor_base_e_coords_versao2(A, k)
    b1,c1 = melhor_base_e_coords_posto1(A)
    B=[b1;]
    C=[c1;]
    A = A-b1*c1
    for i=1:k-1
        b,c = melhor_base_e_coords_posto1(A)
        A = A-b*c
        B = [b B]
        C = [c; C]
    end
    return B, C
end


function projecao_ortogonal(a,b)
    b=b/norm(b)
    projecao= a'*b
    return projecao
end

function melhores_coord(A,b)
    b=b/norm(b)
    c = A'*b
    return c
end

function my_qr(A)
    l,c=size(A)
    q = A[:, 1]
    normQ = norm(q)
    q = q/normQ
    Q = q
    R = [normQ; zeros(l-1)]
    for i=2:c
        proj = zeros(l,1)
        r_vec = [A[:, i]'*Q[:,1]]
        for j=1:i-1
            proj = proj + (A[:, i]'*Q[:, j]).*Q[:, j]
            if j != 1
                r_vec = [r_vec; A[:, i]'*Q[:, j]]
            end
        end
        result = A[:, i] - proj
        result = result/norm(result)
        r_vec = [r_vec; A[:, i]'*result]
        r_vec = [r_vec; zeros(l - i, 1)]
        Q = [Q result]
        R = [R r_vec]
    end
    return Q, R
end



