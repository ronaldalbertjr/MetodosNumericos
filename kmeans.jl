function clusteriza(A,B)
    mA, nA = size(A)
    mB, nB = size(B)
    C = zeros(nB, nA)
    
    for i=1:nA #Iteramos por todas as colunas de A associando cada uma delas a um cluster    
        vC = zeros(nB)
        cluster = 1
        
        #Iteramos por todos os vetores em B procurando pelo que está mais próximo de A e associamos a coluna corrente
        #de A ao cluster de tal vetor.
        for j=2:nB 
            if(norm(A[:, i] - B[:, j]) < norm(A[:, i] - B[:, cluster]))
                cluster = j
            end
        end
        vC[cluster] = 1
        C[:, i] = vC
    end
    
    return C
end

function kmeans(A, k)
    m,n = size(A)
    B = randn(m, k) #Começamos gerando uma base aleatória
    C = clusteriza(A, B) #Encontramos a melhor clusterização para tal base
    norma_anterior = 99999999999999999999 
    while(norm(A - B*C) < norma_anterior) #Realizamos o processo de clusterização até chegar a um mínimo local
        B = (C'\A')'#Achamos a melhor base para as coordenadas encontradas na clusterização
        C = clusteriza(A, B) #Encontramos a melhor clusterização para a base calculada nessa iteração
        norma_anterior = norm(A - B*C) #Calculamos a norma da iteração corrente
    end
    
    return B, C
end
