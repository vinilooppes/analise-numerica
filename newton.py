# Método de Newton para aproximações de raízes de funções univariadas

import math

def newton(f, derivada, x0, epsilon, max_iter=50):
    if abs(f(x0))<=epsilon:
        return x0
    print ("k\t x0\t\t f(x0)") #\t é o espaçamento
    k=1
    while k<=max_iter:

        x1=x0-(f(x0)/derivada(x0))  #cálculo do método de Newton
        print ('%d\t%e\t%e'%(k,x1,f(x1)))

        if abs(f(x1))<=epsilon:
            return x1

        #recebe novas variáveis se a condição não for satisfeita
        x0=x1 
        k=k+1

    print("ERRO: Número máximo de iterações atingido")
    return x1

if __name__ =="__main__":

    #defina uma função
    def f(x): 
        return 3*x**2+3*x-3

    #defina a derivada da função
    def derivada(x):
        return 6*x+3

#defina um chute inicial e a precisão
raiz = newton(f,derivada,1.0,0.0001) 

print(raiz)

