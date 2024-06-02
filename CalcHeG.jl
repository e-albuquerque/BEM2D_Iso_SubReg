function CalcHeG(nnos, crv, kmat,npgauss = 12)
    dcrvs = map(x -> nrbderiv(x), crv)
    n = length(crv);	# Number of curves
    ncollocpoints = size(collocCoord, 1)

    E =zeros(ncollocpoints,ncollocpoints)
    H=zeros(2*ncollocpoints,2*ncollocpoints)
    G=zeros(2*ncollocpoints,2*ncollocpoints)
    
    qsi, w = gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1:n

        p = crv[i].order - 1
        shapes  = zeros(npgauss, (p + 1));
        derivs  = zeros(npgauss, (p + 1), 2);
        for gp = 1:size(w, 1)
            shapes[gp,:], derivs[gp,:,:] = bernsteinbasis(p, 0, qsi[gp], 0);
        end
        for j = 1:size(crv[i].conn, 1)

            k = 1
            for k1 = 1:n
                for k2 = 1:crv[k1].number
                    xfonte = crv[k1].fontes[k2].coords[1:2]
                
                    if k1 == i && crv[k1].fontes[k2].pts >= crv[i].range[j,1] && crv[k1].fontes[k2].pts <= crv[i].range[j,2]
                       eet = 2 * (crv[k1].fontes[k2].pts - crv[i].range[j,1]) / (crv[i].range[j,2] - crv[i].range[j,1]) - 1
                        h = integra_hsing(crv[k1].fontes[k2], crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat[1],kmat[2], eet); # Integra��o sobre o
                        g = integra_gsing(xfonte, crv[i], qsi, w, crv[i].C[:,:,j], crv[i].conn[j], kmat, eet); # Integra��o sobre o
                    else
                        g, h = integra_elem(xfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
                    end
                    colunas=[2*(crv[i].conn[j] .+ nnos[i]).-1 2*(crv[i].conn[j] .+ nnos[i])]'[:]
                    H[2*k-1:2*k,colunas] += h;
                    G[2*k-1:2*k,colunas] += g;
                    k += 1   
                 
                end
            end
        end
    end
    k = 1
    for k1 = 1:n
        uu = unique(crv[k1].knots)
        for k2 = 1:crv[k1].number
            ind = sum(crv[k1].fontes[k2].pts .> uu)
            E[k,crv[k1].conn[ind] .+ nnos[k1]] += crv[k1].fontes[k2].basis
            k += 1   
        end
    end
    H , G, E
end
function  integra_elem(xfonte, crv, qsi, w, shapes, derivs, C, conn, k)
# Integra��o sobre os elementos (integral I)

    g = zeros(2,2*crv.order)
    h = zeros(2,2*crv.order)
    N=zeros(2,2*crv.order)

    for i = 1:size(w, 1) # Percorre os pontos de integra��o
    
        R, dRdxi = basisfundecomp(shapes[i,:], derivs[i,:,:], C, crv.coefs[4,conn])
        # Funções de forma
        N[1,1:2:end]=R
        N[2,2:2:end]=R
        cf = crv.coefs[1:2,conn] ./ [crv.coefs[4,conn]';crv.coefs[4,conn]']
        p = cf * R  
        dp = cf * dRdxi;

        dgamadqsi = norm(dp)

        nx = dp[2] / dgamadqsi; # Componente x do vetor normal unit�rio
        ny = -dp[1] / dgamadqsi; # Componente y do vetor normal unit�rio

        uast,tast=calc_solfund(p,xfonte,nx,ny,k[1],k[2])
    
        h = h + tast*N * dgamadqsi * w[i]
        g = g + uast*N * dgamadqsi * w[i]
    end
    return g, h
end

function  integra_gsing(xfonte, crv, qsi, w, C, conn, k, eet)
    # Integra��o sobre os elementos (integral I)
    
    
    eta, Jt = telles(qsi, eet)
    
    g = zeros(2,2*crv.order)
    N=zeros(2,2*crv.order)    
    for i = 1:size(w, 1) # Percorre os pontos de integra��o
        shapes, derivs = bernsteinbasis(crv.order - 1, 0, eta[i], 0);
        R, dRdxi = basisfundecomp(shapes, derivs, C, crv.coefs[4,conn])
           # derivadas das fun��es de forma
           N[1,1:2:end]=R
           N[2,2:2:end]=R
        cf = crv.coefs[1:2,conn] ./ [crv.coefs[4,conn]';crv.coefs[4,conn]']
        p = cf * R  
        dp = cf * dRdxi;
    
        dgamadqsi = norm(dp)  
         
           
        nx = dp[2] / dgamadqsi; # Componente x do vetor normal unit�rio
        ny = -dp[1] / dgamadqsi; # Componente y do vetor normal unit�rio
        uast,tast=calc_solfund(p,xfonte,nx,ny,k[1],k[2])

        
        g = g + uast*N * dgamadqsi * w[i] * Jt[i]
    end
    return g
end

function  integra_hsing(fonte, crv, qsi, w, shapes, derivs, C, conn, E,ν,eet)
    # Integra��o sobre os elementos (integral I)
    xfonte=fonte.coords[1:2]
    
    basisrc=fonte.basis
    matbasisrc=zeros(2,2*crv.order)
    matbasisrc[1,1:2:end]=basisrc
    matbasisrc[2,2:2:end]=basisrc

    hterm=[0	-(1-2*ν)/(4*pi*(1-ν))
    (1-2*ν)/(4*pi*(1-ν)) 	0]
    htermMatrix=hterm*matbasisrc


        h = zeros(2,2*crv.order)
        N=zeros(2,2*crv.order)


        for i = 1:size(w, 1) # Percorre os pontos de integra��o
        
            R, dRdxi = basisfundecomp(shapes[i,:], derivs[i,:,:], C, crv.coefs[4,conn])
            # Funções de forma
            N[1,1:2:end]=R
            N[2,2:2:end]=R
            cf = crv.coefs[1:2,conn] ./ [crv.coefs[4,conn]';crv.coefs[4,conn]']
            p = cf * R  
            dp = cf * dRdxi;
    
             dgamadqsi = norm(dp)
    
            nx = dp[2] / dgamadqsi; # Componente x do vetor normal unit�rio
            ny = -dp[1] / dgamadqsi; # Componente y do vetor normal unit�rio
    
            uast,tast=calc_solfund(p,xfonte,nx,ny,E,ν)
        
            h +=(tast*N*dgamadqsi-htermMatrix/(qsi[i]-eet))* w[i]
        end
        
        if eet==1
            #  h+=-htermMatrix*log(abs((1)/(1+eet)));
            beta_m=1/fonte.dgamadqsi
            h+=-htermMatrix*log(abs(2/beta_m))
        elseif eet==-1
            #  h+=htermMatrix*log(abs((1-eet)/(1)));
            beta_m=1/fonte.dgamadqsi
            h+=htermMatrix*log(abs(2/beta_m))
        else
             h+=htermMatrix*log(abs((1-eet)/(1+eet)));
        end
        return h
    
end

function aplica_CDC(G, H, CDC, E)
# Aplica as condições de contorno trocando as colunas das matrizes H e G

ncdc = 2*size(CDC,1); # número de linhas da matriz CDC
A=1*H
B=1*G
valoresconhecidos=zeros(ncdc)
  #deslocamento em x
  A[:,(1:2:ncdc)[CDC[:,2] .== 0]]=  -G[:,(1:2:ncdc)[CDC[:,2] .== 0]]
  B[:,(1:2:ncdc)[CDC[:,2] .== 0]]=  -H[:,(1:2:ncdc)[CDC[:,2] .== 0]]
  #deslocamento em y
  A[:,(2:2:ncdc)[CDC[:,4] .== 0]]=  -G[:,(2:2:ncdc)[CDC[:,4] .== 0]]
  B[:,(2:2:ncdc)[CDC[:,4] .== 0]]=  -H[:,(2:2:ncdc)[CDC[:,4] .== 0]]

valoresconhecidos[1:2:end]=E\CDC[:,3] # Valores das condições de contorno
valoresconhecidos[2:2:end]=E\CDC[:,5] # Valores das condições de contorno
b=B*valoresconhecidos; # vetor b

return A,b
end

function monta(CDC,x)
    # Separa deslocamento e forças de superfície
    
    ncdc = 2*size(CDC,1);
    u = zeros(size(CDC,1),2)
    t = zeros(size(CDC,1),2)
     #deslocamento em x
     u[CDC[:,2] .== 0,1]=  CDC[CDC[:,2] .== 0,3]
     t[CDC[:,2] .== 0,1]=   x[(1:2:ncdc)[CDC[:,2] .== 0]]
     #deslocamento em y
     u[CDC[:,4] .== 0,2]=  CDC[CDC[:,4] .== 0,5]
     t[CDC[:,4] .== 0,2]=   x[(2:2:ncdc)[CDC[:,4] .== 0]]
      #forças de superfície em x
      u[CDC[:,2] .== 1,1]=  x[(1:2:ncdc)[CDC[:,2] .== 1]]
      t[CDC[:,2] .== 1,1]=  CDC[CDC[:,2] .== 1,3]
      #forças de superfcie em y
      u[CDC[:,4] .== 1,2]=  x[(2:2:ncdc)[CDC[:,4] .== 1]]
      t[CDC[:,4] .== 1,2]=  CDC[CDC[:,4] .== 1,5]
    return u,t
    end


function  calc_solfund(x,x0,nx,ny,E,ν)
    #Calcula as solu��es fundamentais
  
    GE=E/(2*(1+ν));
    # Distance of source and field points
    R = x - x0;
    r = norm(R)
  
    # Components of the unity vector in the radial direction
    rd1 = R[1]/r;
    rd2 = R[2]/r;
  
     # Plane elasticity fundamental solutions
    prod1 = 4*pi*(1-ν);
    prod2 = (3-4*ν)*log(1/r);
    prod3 = rd1*nx + rd2*ny;
    
     u11 = (prod2 + rd1^2)/(2*prod1*GE);
    u22 = (prod2 + rd2^2)/(2*prod1*GE);
    u12 = (rd1*rd2)/(2*prod1*GE);
    u21=u12;
  
    t11 = -(prod3*((1-2*ν)+2*rd1^2))/(prod1*r);
    t22 = -(prod3*((1-2*ν)+2*rd2^2))/(prod1*r);
    t12 = -((prod3*2*rd1*rd2)-(1-2*ν)*(rd1*ny-rd2*nx))/(prod1*r); # Checar o
        # sinal de + pois no programa do Brebbia est� diferente.
    t21 = -((prod3*2*rd1*rd2)-(1-2*ν)*(rd2*nx-rd1*ny))/(prod1*r); # Checar o
        # sinal de + pois no programa do Brebbia est� diferente.
        uast = [u11 u12
             u21 u22];
  
         tast = [t11   t12
             t21   t22];
  
  uast,tast
  end