import math
# definir um intervalo [a,b] e um erro e

print("Defina um intervalo [a,b]")
a=float(input("entre com o valor de a: "))
b=float(input("entre com o valor de b: "))

print("Defina o valor do erro desejado")
e=float(input("entre com o valor do erro: "))

r=float(input("digite o número que você quer calcular a raíz quadrada: "))

#definir uma função
def f(x):
    return x**2-r

#Teorema de Bolzano

if f(a)*f(b)<0:
    while (math.fabs(b-a)/2>e):
        xi=(a+b)/2
        if f(xi)==0:
            print("A raíz é: ",xi)
            break
        else:
            if f(a)*f(xi)<0:
                b=xi
            else:
                a=xi
    print ("O valor da raíz é: ",xi)
        
else:
    print("Não existem raízes neste intervalo")