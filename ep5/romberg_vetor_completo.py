import numpy as np

#Parâmetros que serão usados para calcular o método de Romberg

tol=10**(-9)    #tolerância

#print(tol)

NMAX=20        #número maximo de linhas 

#funções auxiliares 

def funcao (x,ex):     #rotina que avalia f(x)
    if ex==1:
        f=x**2
    elif ex==2:
        f=1/(1-x)
    elif ex==3:
        f= x**6
    elif ex==4:
        f=1/(1+25*(x**2))        
    return(f)


def somatorio (I,H,a,ex):   #rotina que faz o calculo do somatorio para usar nos trapezios

    vetsoma=[]
    for j in range (1,2**(I-1)+1):
        #print(j)
        x=(a+((2*j)-1)*H)
        #print(x)
        vetsoma.append(funcao(x,ex))
    s=np.sum(vetsoma)
    #print(s)
    return(s)


def trapz(a,b,T,i,ex):     #rotina que faz o cálculo dos trápezios 
    h=(b-a)/(2**i)
    if i==0:
        t = 0.5 * h * (funcao(a,ex)+funcao(b,ex))   #calculo do 1-trapezio
    else:
        t = 0.5 * T + h * somatorio(i,h,a,ex)      #calculo dos n-trapezios
    return(t)


def romb(a,b,tol,ITMAX,ex):      #rotina que faz o calculo do método de Romberg
    T=[]        #cria a matriz de Romberg
    T.append ([trapz(a,b,T,0,ex)])     #calculo do 1º elemento da 1ª coluna da tabela de Romberg
    
    #print(T[0])    #print da 1ª linha da tabela

    i=1  #auxilia na construção das linhas 

    while (i < ITMAX):      #construção da matriz de Romberg linha por linha, apartir da 2ª linha
        
        if i>0: 
            T.append([trapz(a,b,T[i-1][0],i,ex)])  #calculo dos elementos da 1ª coluna da tabela de Romberg
            #print(T)

            for k in range(1,i+1):
                T[i].append (T[i][k-1]+ (T[i][k-1]-T[i-1][k-1])/((4**k)-1))     #calculo dos demais elementos da tabela de Romberg

            #print(T[i])    #print da i-linha da tabela


        #print(T[i][i] - T[i][i-1])
        #print(abs(T[i][i]-T[i][i-1]) , (tol*abs(T[i][i])))   

        if (abs(T[i][i]-T[i][i-1]) <= (tol*abs(T[i][i]))):       #Verificação do critério de parada
            break

        i += 1      #acrescenta +1 no i para construir uma nova linha caso o critério de parada não seja satisfeito

    if i>=ITMAX :
        print('Número máximo de iteração atingido.')
        exit()
           
    #print(T[i][i])
    #print(i)

    return(T[i][i],i)
    
# Escolha dos exemplos

print('\n Método de Romberg para cálculo de integrais definidas \n')
print('Exemplos de integrais para calcular usando o Método de Romberg')
print('\n 1. ∫ x**2 dx no intervalo [0,1] \n')
print('\n 2. ∫ 1/(1-x) dx no intervalo [0,0.995] \n')
print('\n 3. ∫ x**6 dx no intervalo [0,1] \n')
print('\n 4. ∫ 1/(1+25*(x**2)) dx no intervalo [-1,1] \n')
ex=int(input('Qual integral deseja calcular usando o Método de Romberg?: '))

if ex==1:
    integral,n=romb(0,1,tol,NMAX,ex)      #calculo da integral usando Romberg
    print('∫ x**2 dx ≅',integral,' no intervalo [0,1] com n=',n)
elif ex==2:
    integral,n=romb(0,0.995,tol,NMAX,ex)       #calculo da integral usando Romberg
    print('∫ 1/(1-x) dx ≅',integral,' no intervalo [0,0.995] com n=',n)
elif ex==3:
    integral,n=romb(0,1,tol,NMAX,ex)       #calculo da integral usando Romberg
    print('∫ x**6 dx ≅',integral,' no intervalo [0,1] com n=',n)
elif ex==4:
    integral,n=romb(-1,1,tol,NMAX,ex)
    print('∫ 1/(1+25*(x**2)) dx ≅',integral,' no intervalo [-1,1] com n=',n)
else:
    print('Teste invalido')