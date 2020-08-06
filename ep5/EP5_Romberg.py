import numpy as np
from scipy import integrate 

#Parâmetros que serão usados para calcular o método de Romberg

tol=10**(-9)    #tolerância

#print(tol)

NMAX=20        #número maximo de linhas 

#funções auxiliares 

funcao1=lambda x: x**2
funcao2= lambda x: 1/(1-x)
funcao3=lambda x: 1/(1+(x**2))  
funcao4=lambda x: x**6
funcao5=lambda x: np.exp(x)

def somatorio (I,H,a,f):   #rotina que faz o calculo do somatorio para usar nos trapezios

    vetsoma=[]
    for j in range (1,2**(I-1)+1):
        #print(j)
        x=(a+((2*j)-1)*H)
        #print(x)
        vetsoma.append(f(x))
    s=np.sum(vetsoma)
    #print(s)
    return(s)


def trapz(a,b,T,i,f):     #rotina que faz o cálculo dos trápezios 
    h=(b-a)/(2**i)
    if i==0:
        t = 0.5 * h * (f(a)+f(b))   #calculo do 1-trapezio
    else:
        t = 0.5 * T + h * somatorio(i,h,a,f)      #calculo dos n-trapezios
    return(t)


def romb(f,a,b,tol,ITMAX):      #rotina que faz o calculo do método de Romberg
    T=[]        #cria a matriz de Romberg
    T.append ([trapz(a,b,T,0,f)])     #calculo do 1º elemento da 1ª coluna da tabela de Romberg
    
    #print(T[0])    #print da 1ª linha da tabela

    i=1  #auxilia na construção das linhas 

    while (i < ITMAX):      #construção da matriz de Romberg linha por linha, apartir da 2ª linha
        
        if i>0: 
            T.append([trapz(a,b,T[i-1][0],i,f)])  #calculo dos elementos da 1ª coluna da tabela de Romberg
            #print(T)

            for k in range(1,i+1):
                T[i].append (T[i][k-1]+ (T[i][k-1]-T[i-1][k-1])/((4**k)-1))     #calculo dos demais elementos da tabela de Romberg

            #print(T[i])    #print da i-linha da tabela


        #print(T[i][i] - T[i][i-1])
        #print(abs(T[i][i]-T[i][i-1]) , (tol*abs(T[i][i])))   

        if (abs(T[i][i]-T[i][i-1]) <= (tol*abs(T[i][i]))):       #Verificação do critério de parada
            break

        i += 1      #acrescenta +1 no i para construir uma nova linha caso o critério de parada não seja satisfeito

    if i>=ITMAX : #condição que encerra o programa caso o número máximo de iterações estipulado seja atingidos
        print(' Número máximo de iteração atingido.')
        exit()

    if i<=8:  #impressão da tabela de Romberg caso o n seja menor que 8 e maior que 1
        print(' Tabela de Romberg para o exemplo escolhido: \n')
        for j in range (i+1):
            print(T[j])
           
    #print(T[i][i]) 
    #print(i)

    return(T[i][i],i)
    
# Escolha dos exemplos

print('\n Método de Romberg para integração de funções suaves \n')
print(' Exemplos de integrais para calcular usando o Método de Romberg')
print('\n 1. ∫ x**2 dx no intervalo [0,1] \n')
print('\n 2. ∫ 1/(1-x) dx no intervalo [0,0.995] \n')
print('\n 3. ∫ 1/(1+(x**2)) dx no intervalo [-5,5] \n')
print('\n 4. ∫ x**6 dx no intervalo [0,1] \n')
print('\n 5. ∫ e**x dx no intervalo [0,1] \n')
#print('\n 6. ∫ 1/(1-x) dx no intervalo [-5,5] (exemplo em que a f(x) não é suave) \n')

ex=int(input(' Qual integral deseja calcular usando o Método de Romberg?: '))

if ex==1:
    integral,n=romb(funcao1,0,1,tol,NMAX)      #calculo da integral usando Romberg
    print('\n A aproximação pelo Método de Romberg para ∫ x**2 dx no intervalo [0,1] com n=',n,'é:', integral)
elif ex==2:
    integral,n=romb(funcao2,0,0.995,tol,NMAX)       #calculo da integral usando Romberg
    print('\n A aproximação pelo Método de Romberg para ∫ 1/(1-x) dx no intervalo [0,0.995] com n=',n,'é:', integral)
elif ex==3:
    integral,n=romb(funcao3,-5,5,tol,NMAX)      #calculo da integral usando Romberg  #FUNÇÃO DE RUNGE 
    print('\n A aproximação pelo Método de Romberg para ∫ 1/(1+(x**2)) dx no intervalo [-5,5] com n=',n,'é:', integral)  
elif ex==4:
    integral,n=romb(funcao4,0,1,tol,NMAX)       #calculo da integral usando Romberg
    print('\n A aproximação pelo Método de Romberg para ∫ x**6 dx no intervalo [0,1] com n=',n,'é:', integral)
elif ex==5:
    integral,n=romb(funcao5,0,1,tol,NMAX)
    print('\n A aproximação pelo Método de Romberg para ∫ e**x dx no intervalo [0,1] com n=',n,'é:', integral) #calculo da integral usando Romberg
#elif ex==6:  #Retire os # caso queira fazer o teste para uma função que não seja suave. Obs: não esqueça de retirar tbm do print da escolha dos exemplos.
 #   integral,n=romb(funcao2,-5,5,tol,NMAX)       #calculo da integral usando Romberg
  #  print('\n A aproximação pelo Método de Romberg para ∫ 1/(1-x) dx no intervalo [0,0.995] com n=',n,'é:', integral)
else:
    print('\n Teste invalido')

#inte=integrate.romberg(funcao1,0,1,show=True)   #altere a função e o intervalo de integração para verificar os valores usando o método de Romberg do scipy