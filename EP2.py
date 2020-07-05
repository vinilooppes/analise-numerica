n=int(input("Entre com a dimensão da matriz A "))



#Entra com a diagonal matriz A com dimensão N

a=list(range(n))

print("Entre com os elementos da diagonal de A")

for i in a:
    a[i]=float(input())



#Entra com a subdiagonal matriz A com dimensão N-1

b=list(range(n-1))

print("Entre com os elementos da subdiagonal de A")

for i in b:
    b[i]=float(input())


d=list(range(n))         #vetor D com dimensão N
l=list(range(n-1))       #vetor L com dimensão N-1


def decomp (a,b):                     #Rotina que faz o calculo da decomposição 
    
    d[0]=a[0]

    for i in range(1,n):
        l[i-1]=b[i-1]/d[i-1]               #calculo dos elementos do vetor l, da matriz L
        d[i]=a[i]-(l[i-1]**2)*d[i-1]       #calculo dos elementos do vetor d, da matriz D

    return l,d  

    

L,D =decomp (a,b)  

print("O vetor representando L é: ", L)
print("O vetor representando D é: ", D)




#entra com o vetor V do sistema Ax=v (dimensão N)

v=list(range(n))

print("Entre com o valores do vetor b do sistema linear Ax=b")
for i in v:
    v[i]=float(input())                    


#Vetores "incognitas" para resolver o sistema (dimensão N)
y=list(range(n))
z=list(range(n))
x=list(range(n))


#Rotina que resolve sistemas lineares Ax=v com A tridiagonal simétrica
def resolve_sistemas(v):               

    n=len(a)

    y[0]=v[0]

    for i in range(1,n):                #começa no 1 e vai até o anterior ao N
        y[i]=v[i]-l[i-1]*y[i-1]         
        


    for i in range(n):
        z[i]=y[i]/d[i]
    
    

    x[n-1]=z[n-1]

    for i in range(n-2,-1,-1):        #começa no n-1 e vai até 0 
        x[i]=z[i]-l[i]*x[i+1]
     

    return (x)

sol= resolve_sistemas(v)
print("o vetor solução do sistema é: " , sol)

