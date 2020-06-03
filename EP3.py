import numpy as np


n=int(input('Entre com a dimensão do sistema: '))


#pontos da placa
'''
x=list(range(n+1))
y=list(range(n+1))

for i in range (1,n+1):
    x[i]=i/n
    y[i]=i/n
'''

#Encontra o sistema Ax=b
'''
#encontrar valores da temp no bordo da placa
bordo=[]
for i in range(n+1):
    linha=[]
    for j in range(n+1):
        if i==0 or j==0:
            linha.append(-3)
        else:
            if i==n or j==n:
                if i==n:
                    linha.append(6*y[j]-3)
                else:
                    linha.append(6*x[i]-3)
            else:
                linha.append(0) 
    bordo.append(linha)

print(np.matrix(bordo))

#bloco 1 de "B" do sis Ax=B
B1=list(range(n-1))
for j in range (n-1):
    if j==0 :
        B1[j]=bordo[0][1]+bordo[1][0]
    else :
        if j+2==n :
            B1[j]=bordo[0][n-1]+bordo[1][n]
        else :
            B1[j]=bordo[0][j]
#print(B1)

#bloco 2 a n-2 de "B" do sis Ax=B
def meio_b (x):

    B2=list(range(n-1))
    for j in range (n-1):
        if j==0 :
            B2[j]=bordo[x][0]
        else :
            if j+2==n :
                B2[j]=bordo[x][n]
            else :
                B2[j]=0
    return(B2)

#B2=meio_b(2)
#print(K2)

#bloco n-1 de "B" do sis Ax=B
B3=list(range(n-1))
for j in range (n-1):
    if j==0 :
        B3[j]=bordo[n][1]+bordo[n-1][0]
    else :
        if j+2==n :
            B3[j]=bordo[n][n-1]+bordo[n-1][n]
        else :
            B3[j]=bordo[n][j]
#print(B3)

#matriz tridiagonal com diag princial 4 e secundarias -1
D1=[]
for i in range(n-1):
    linha = [ ]
    for j in range(n-1):
        if i==j :
            linha.append(4)
        else :
            if (i+1==j) or (i==j+1) :
                linha.append(-1)
            else :
                linha.append(0)
    D1.append(linha)

D=np.matrix(D1)


N=(-1)*np.eye(n-1) #matriz id mult. por -1
nula= np.zeros((n-1,n-1))  #matriz nula

#matriz A do sistema Ax=b
A = []
for i in range((n-1)):
    linha = [ ]
    for j in range((n-1)):
        if i==j :
            linha.append(D)
        else :
            if (i+1==j) or (i==j+1) :
                linha.append(N)
            else :
                linha.append(nula)
    A.append(linha)
A1=np.block(A)

#print(A1)

#vetor B
B=list(range(n-1))  
for j in range(n-1):
    if j==0 :
        B[j]=B1
    else :
        if j==n-2 :
            B[j]=B3
        else :
            B[j]=meio_b(j)

#print(B)
'''
   

def metodoSOR(x,omega,erro,mat):
    '''calcula o método de SOR para a matriz T'''
    criterioParada = 100   #criterio de parada inicial 
    QtIterações = 0        #contador das iterações 
    T = mat.copy()         #matriz que será alterada
    Tanterior = mat.copy()      #matriz que armazena a iteração x**(n-1)

    while criterioParada > erro and QtIterações < 500:
        QtIterações = QtIterações + 1
        for i in range(1,n):
            for j in range(1,n):
                T[i][j]=((1-omega)*Tanterior[i][j])+((omega/4)*(T[i-1][j]+T[i+1][j]+T[i][j-1]+T[i][j+1]))
        S = T - Tanterior
        criterioParada = abs(max(S.min(), S.max(), key=abs)) #norma infinita de x-x0
        Tanterior = T.copy()
    
    return(QtIterações,T)


def calculoSOR(VetOmega,Matriz):
    '''Faz o calculo do método de SOR para os omegas de 0 a 100'''
    VetorIteracao=np.zeros(101)
    for i in range (101):
        Iteracoes,mat = metodoSOR(n,Omega[i],Erro,Matriz)
        VetorIteracao[i]=Iteracoes

    return(VetorIteracao)


def QuantidadeIteracoes(Iteracoes,W):
    '''Descobre qual omega que gera a menor quantidade de iterações'''
    QuantidadeIteracao=Iteracoes[0]
    PosicaoOmega=0
    OmegaQual=W[0]
    for i in range (100):
        if Iteracoes[i] > Iteracoes[i+1] :
            if QuantidadeIteracao!=Iteracoes[i+1]:
                OmegaQual=W[i+1]
                PosicaoOmega=i+1
                QuantidadeIteracao=Iteracoes[i+1]

    return(QuantidadeIteracao,PosicaoOmega,OmegaQual)


Erro=float(input('Entre com o valor do Erro ')) #pede para o usuario entrar com um erro


Omega = np.zeros(101) #vetor com a variação do parâmetro omega

#definimos os valores da variação do parâmetro omega com k=100
for i in range (101):
    Omega[i]= 1+(i/100)     #sempre varia de 0 a 100

#print(Omega)


'''faz os calculos para o chute inicial 0s'''
#definimos T com os valores do Bordo e chute inicial 0s
T=np.zeros((n+1,n+1))
for i in range(n+1):
    for j in range(n+1):
        if i==0 or j==0:
            T[i][j]=-3
        else:
            if i==n or j==n:
                if i==n:
                    T[i][j]=6*(j/n)-3
                else:
                    T[i][j]=6*(i/n)-3
            else:
                T[i][j]=0


IterZeros = np.zeros(101)        # vetor com a quantidade de iterações com chute inicial 0s

#calculo do SOR com o chute inicial 0s
IterZeros=calculoSOR(Omega,T) 

#print(IterZeros)

#ve qual o menor numero de iterações para o chute inicial com 0s
QTzeros,PosicaoOmegaZeros,OmegaZeros= QuantidadeIteracoes(IterZeros,Omega) 

print('O Omega ',PosicaoOmegaZeros,'resulta no menor número de iterações para o chute inicial com pontos zeros, o número de iterações é:', QTzeros ,'e este Omega é:',OmegaZeros)



'''faz os calculos para o chute inicial aleatorio'''
#redefinimos T com os valores do Bordo e chute inicial aleatórios entre 0 e 1
T1=T.copy()
for i in range(1,n):
    for j in range(1,n):
        T1[i][j]=np.random.random_sample(1)   #altera os valores do chute inicial para números aleatórios entre 0 e 1
#print(T1)

IterAleatorio = np.zeros(101)   # vetor com a quantidade de iterações com chute inicial de valores aleatorios entre 0 e 1

#calculo do SOR com o chute inicial com elementos aleatórios 
IterAleatorio=calculoSOR(Omega,T1)

#print(IterAleatorio)


#ve qual o menor numero de iterações para o chute inicial com numeros aleatorios entre 0 e 1
QTaleatorios,PosicaoOmegaAleatorios,OmegaAleatorios = QuantidadeIteracoes(IterAleatorio,Omega) 

print('O Omega ',PosicaoOmegaAleatorios,'resulta no menor número de iterações para o chute inicial com pontos entre 0 e 1, e o número de iterações é:', QTaleatorios ,'e este Omega é:',OmegaAleatorios)