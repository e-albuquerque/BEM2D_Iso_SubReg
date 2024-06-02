# include("Elastico.jl")
include("dad.jl")
include("decomp.jl")
include("nurbs.jl")
include("CalcHeG.jl")
include("telles.jl")
include("formatiso.jl")
using Plots
using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
gr()



# O número de nós (knots), m, o número de pontos de controle, k, e a ordem da
# curva, n , estão relacionados por:
#                           m = k + n + 1

PONTOS,SEGMENTOS,CCSeg,kmat=dad_2() #Arquivo de entrada de dados
# PONTOS,SEGMENTOS,MALHA,CCSeg,kmat=dad_2() #Arquivo de entrada de dados

# NOS,ELEM=format_dad(PONTOS,SEGMENTOS,MALHA)# formata os dados (cria as
crv=format_dad_iso(PONTOS,SEGMENTOS)# formata os dados

# display(mostra_geo(crv))

dcrv=map(x->nrbderiv(x),crv)
n = length(crv);	# N�mero total de elementos

p=0;#refinamento p

for i=1:n
    degree=crv[i].order-1
    coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
	crv[i] = nrbmak(coefs,knots)
end



h=0;#refinamento h
if h>0
for i=1:n
    novosnos=range(0,stop=1,length=h+2)
	degree=crv[i].order-1
	coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
	crv[i] = nrbmak(coefs,knots)
end
end
 

z=0
for k=1:n
    for i=1:crv[k].number
	    global z=z+1
    end
end
numcurva=zeros(Integer,z)
collocPts=zeros(z)
CDC=zeros(z,5)
collocCoord=zeros(z,2)
z=0;
nnos=zeros(Integer,n)
for k=1:n
    p=crv[k].order-1;
    nnos[k]=crv[k].number;

    for i=1:crv[k].number
        global z=z+1;
        # numcurva[z]=k;
        # collocPts[z]=sum(crv[k].knots[(i+1):(i+p)])/p;
        # if(i==2)
        #     collocPts[z-1]=(collocPts[z]+collocPts[z-1])/2;
        # end
        # if(i==nnos[k])
        #     collocPts[z]=(collocPts[z]+collocPts[z-1])/2;
        # end

       CDC[z,:] = [z; CCSeg[k,2:5]];
    end
end
nnos2=cumsum([0 nnos'],dims=2);

# E=zeros(length(collocPts),length(collocPts));
# for i=1:length(collocPts)
#     collocCoord[i,:]=nrbeval(crv[numcurva[i]], collocPts[i]);
#     B, id = nrbbasisfun(crv[numcurva[i]],collocPts[i])
#     E[i,id+nnos2[numcurva[i]]]=B
# end
# plot(collocCoord[:,1],collocCoord[:,2])
# legend('Curva resultante','Polígono de controle','Pontos de controle','Pontos fonte')
H,G,E=CalcHeG(nnos2,crv,kmat,10)
H[1:2:end,1:2:end]+=E/2
H[2:2:end,2:2:end]+=E/2

A,b= aplica_CDC(G,H,CDC,E)
x=A\b # Calcula o vetor x

uc,tc=monta(CDC,x) # Separa deslocamento e força

u=E*uc
t=E*tc

cont=1
for c in crv
for f in c.fontes
    global cont
collocCoord[cont,:]=f.coords[1:2]
collocPts[cont]=f.pts
cont+=1
end
end

plot(collocCoord[:,2],u)

