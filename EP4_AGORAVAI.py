import numpy as np  #biblioteca para trabalhar com vetores e matrizes 
from matplotlib import pyplot as plt    #biblioteca para plotar os graficos 
import tarefa2completa as EP2  #importo o ep2 para resolver a decomposição e resolver o sistema

#rotina que calcula o vetor da diagonal da matriz A
def CalcDiag(h,lamb,n): 
    A1=[]
    for i in range (1,n+1):
        A1.append((2/h)+(lamb*2*h)/3)
    return(A1)

#rotina que calcula o vetor da subdiagonal da matriz A
def CalcSubdiag (h,lamb,n):  
    A2=[]
    for i in range (1,n):
        A2.append((-1/h)+(lamb*h)/6)
    return (A2)

#rotina que faz o calculo da integral xf(x) 
def CalcIntegral1 (x,ex): 
    if ex==1:
        int1=(x**2)/2
    elif ex==2:
        int1= 4*(x**3)-3*(x**4)-(x**2)
    elif ex==3:
        int1= -2*(np.pi)*x*(np.cos(np.pi*x))+2*(np.sin(np.pi*x))
    elif ex==4:
        int1=(x**2)+2*(x**3)-3*(x**4)
    return(int1)

#rotina que faz o calculo da integral de f(x)
def CalcIntegral2 (x,ex): 
    if ex==1:
        int2=x
    elif ex==2:
        int2= 6*(x**2)-4*(x**3)-2*x
    elif ex==3:
        int2= -2*(np.pi)*(np.cos(np.pi*x))
    elif ex==4:
        int2=2*x+3*x**2-4*x**3
    return(int2)

#rotina que calcula a matriz d
def Calc_D(intg1,intg2,n,x,ex,h): 
    D=[]
    for i in range (1,n+1):
        D.append((intg1(x[i],ex)-intg1(x[i-1],ex)-intg1(x[i+1],ex)+intg1(x[i],ex)-x[i-1]*(intg2(x[i],ex)-intg2(x[i-1],ex))+x[i+1]*(intg2(x[i+1],ex)-intg2(x[i],ex)))/h)
    return (D)

#rotina para o calculo da solucao exata
def CalcSolExata(x,ex):
    if ex==1:
        SOL=0.5*x*(1-x)
    elif ex==2:
        SOL=(x**2)*(x-1)**2
    elif ex==3:
        SOL=np.sin(np.pi*x)
    elif ex==4:
        SOL=x-x**2-x**3+x**4
    return(SOL)

#rotina que faz o calculo das phi's para encontrar a solucao aproximada para cada xi
def CalcPhi (k,i,x,h):
    if k >= x[i-1] and k <= x[i] :
        phi=(k-x[i-1])/h
    else:
        if k >= x[i] and k <= x[i+1]:
            phi=(x[i+1]-k)/h
        else: 
            phi=0
    return(phi)

#rotina que faz o calculo do somatorio da solucao aproximada para cada xi/yi
def CalcSolAprox (k,n,c,x,h):
    u2=[]
    for i in range (1,n+1):
        u2.append(c[i-1]*(CalcPhi(k,i,x,h)))
    soma=np.sum(u2)
    return(soma)


# Rotina que faz o calculo do EP para os n
def CalculoEP4(A,B,n,ex):

    h=(B-A)/(n+1)

    x=list(range(n+2)) #calculo do vetor com dos xi 
    for i in range(n+2):
        x[i]=A+i*h

    #print(x)

    yi=[]       #vetor que usarei para o calculo do erro
    sol=[]      #vetor da solucao exata

    #   Calculo dos vetores (matriz) A e d do sistema Ac=d para o exercicio escolhido
    if ex==1:

        Q=n+2    
        yi=x     #altero x para yi para fazer o calculo da solucao aproximada la embaixo

        lambd=0
        a = CalcDiag(h,lambd,n)
        b = CalcSubdiag(h,lambd,n)
        v = Calc_D(CalcIntegral1,CalcIntegral2,n,x,ex,h)

        for i in range (Q):                 #calculo da solucao exata
            sol.append(CalcSolExata(x[i],ex))

    else:    
        Q=(10*n)+1

        for i in range (Q):            #calculo do vetor yi
                yi.append(i/(10*n))

        #print(yi)

        if ex==2:
            lambd=0
            a=CalcDiag(h,lambd,n)
            b=CalcSubdiag(h,lambd,n)
            v=Calc_D(CalcIntegral1,CalcIntegral2,n,x,ex,h)

            for i in range (Q):      #faz o calculo da solucao exata nos nós xi
                sol.append(CalcSolExata(yi[i],ex))

        elif ex==3:
            lambd= np.pi**2
            a=CalcDiag(h,lambd,n)
            b=CalcSubdiag(h,lambd,n)
            v=Calc_D(CalcIntegral1,CalcIntegral2,n,x,ex,h)

            for i in range (Q):      #faz o calculo da solucao exata nos nós yi
                sol.append(CalcSolExata(yi[i],ex))
            #print(yi)

        elif ex==4:
            lambd= 0
            a=CalcDiag(h,lambd,n)
            b=CalcSubdiag(h,lambd,n)
            v=Calc_D(CalcIntegral1,CalcIntegral2,n,x,ex,h)

            for i in range (Q):      #faz o calculo da solucao exata nos nós yi
                sol.append(CalcSolExata(yi[i],ex))

    #print(a)
    #print(b)
    #print(v)

    #   calculo da decomposicao e resolução do sistema Ac=d  (Usando o EP 2)
    L,D = EP2.decomp(a,b,n)
    c = EP2.resolve_sistemas(L,D,v,a,n)
    
    #print("o vetor solução do sistema é: " , c)

    #calculo do vetor com as solucoes aproximada para cada yi
    SolAprox=[]
    for i in range (Q):      
        SolAprox.append(CalcSolAprox(yi[i],n,c,x,h))

    #print('O vetor com as soluções aproximada para cada xi/yi é:')
    #print(SolAprox)
    #print('O vetor com as soluções exata para cada xi/yi é:')
    #print(sol)


    #plotagem do grafico da solucao numerica do ex 2, 3 e 4
        
    if ex==2 or ex==3 or ex==4:
        plt.plot(yi,SolAprox)
        #plt.plot(yi,sol)
        plt.xlabel("Valores de yi")
        plt.ylabel("Valor do erro associado a yi")
        plt.title("Comportamento da Solução Numérica")
        plt.show()

    erro = np.max(abs(np.subtract(sol,SolAprox)))     #calculo do erro

    print('O erro para n =',n,'é: ',erro)

    return(erro)

#defino o intervalo
a=0
b=1 

exercicio=int(input('Digite o numero do exercicio que quer aplicar o EP (Ex. 1, 2, 3 ou 4(Segunda parte do 3)):'))

N=[15,31,63,127,255]   #vetor com os n's que preciso testar o ep

#N=[15]

ERRO=[]      

for i in range (len(N)):      #calculo do erro para cada n
    ERRO.append(CalculoEP4(a,b,N[i],exercicio))

#print(ERRO)


#não consegui usar isso
'''
N=[254,255,256]
er=[]
for i in range (len(N)):
    er.append(CalculoEP4(a,b,N[i],2))  
ordemconvergencia = np.log(er[2]/er[1])/np.log(er[1]/er[0])
print(ordemconvergencia)
'''


if exercicio==2 or exercicio==3 or exercicio==4:

    #verificação da ordem que o erro converge para 0
    for i in range (1,len(ERRO)):
        print('O erro entre n=',N[i-1],'e n=',N[i],'está caindo aproximadamente', ERRO[i-1]/ERRO[i],'vezes')   #não vai ser 4 por conta dos erros de arredondamento
    
    # Plotagem do Grafico do erro para o ex 2 e 3
    plt.plot(N,np.log(ERRO))
    plt.xlabel("Valores do n")
    plt.ylabel("Valor do log do erro associado a n")
    plt.title("Comportamento do Erro")
    plt.show()
