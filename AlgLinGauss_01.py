#comentarios y sugerencias gmunoz@udistrital.edu.co

#funciona en Python o en Jython
#Los programas requieren python o jython
#
#Descargar los dos archivos AlgLinMat.py y AlgLinGauss.py

#En Python o Jython cambiar el directorio a la carpeta de los archivos
#import os
#os.chdir('/dir1/dir2/')

#Ejecutar
#import AlgLinGauss

#reducciOn por Gauss-Jordan
#\item Vaya a la columna no cero extrema izquierda.
#\item Si el primer renglOn tiene un cero en la columna del paso (1), intercAmbielo con uno que tenga un elemento no cero en la misma columna.
#\item Obtenga ceros abajo del elemento delantero, sumando múltiplos adecuados del renglón superior a los renglones debajo de él.
#\item Cubra el renglón superior y repita el mismo proceso comenzando por el paso (1) aplicado a la sub-matriz restante. Repita este proceso con el resto de los renglones.
#\item Comenzando con el último renglón no cero, avance hacia arriba: para cada renglón obtenga un 1 delantero e introduzca ceros arriba de él, sumando múltiplos adecuados a los renglones correspondientes.
#import sys
#print(sys.path)  
#from sympy import *
import numpy as np
import sympy

def matriz(A):
    return sympy.Matrix(A)
    #return sympy.Matrix(np.matrix(A))


def row_op(A,i,j,k=0):
    #if k==0: R_i <-> Rj
    #elif i=j: kRi -> Ri
    #else: kRj + Ri -> Ri 
    r,c=A.shape
    E=sympy.eye(r)
    if k==0:
        E[i,i]=0
        E[j,j]=0
        E[i,j]=1
        E[j,i]=1        
    else:
        E[i,j]=k
    D=E*A
    sympy.pprint(D)
    return D
    
def str2num(s):
    if s.find('.')!=-1:
        return float(s)
    elif s.find('/')!=-1:
        return sympy.Rational(s)
    elif s.strip()=="-":
        return -1
    else:
        return int(s)
    
    
def op(Amat,s):
    #if s=='nada':
    #    return Amat.ultop
    s=s.lower()
    #print(s)
    def interpretaRenglon(s1):
        s=s1.replace(' ','')
        #print("interpretando "+s)
        d=s.find("r")
        if d==-1:
            raise ValueError("fala la 'r'. Error en "+s)
        elif d==0:
            k=1
        else:
            try:
                k=str2num(s[:d])
            except ValueError as err:
                raise ValueError("antes de 'r' sólo debe haber un entero. Error en "+s+'. '+str(err))
        #a1=s[d+1:]
        #if not a1.isdigit():
        #print ("coef ="+str(k))
        try:
            a=int(s[d+1:])#-1
        except ValueError as err:
            raise ValueError("después de 'r' sélo debe haber un dígito. Error en "+s+'. '+str(err))
        #a=int(a1)
        return (k,a)
            
    #print (str(k)+", "+str(a))
    (size,z)=Amat.shape#size()
    m=sympy.eye(size)#matriz.identidad(size,size,"no_mostrar")
    d=s.find("<->")
    if d!=-1:
        (k1,a1)=interpretaRenglon(s[:d])
        (k2,a2)=interpretaRenglon(s[d+3:])
        if (k1 != 1) or (k2 != 1):
            raise ValueError("al intercambiar renglones no se puede multiplicar. Error en "+s)
        r="R"+str(a1)+" <-> R"+str(a2)
        m[a1,a1]=0    
        m[a2,a2]=0    
        m[a1,a2]=1    
        m[a2,a1]=1    
        #print(str(m))
    else:
        d=s.find("->")
        if d==-1:
            raise ValueError("no se encontrO ni '->' ni '<->'. Error en "+s)
        else:
            (k2,a2)=interpretaRenglon(s[d+2:])
            e=s[:d].find("+")
            f=s[:d].rfind("-")
            if e!=-1:
                (k11,a11)=interpretaRenglon(s[:e])
                (k12,a12)=interpretaRenglon(s[e+1:d])
                if (a11==a2) and (a12==a2):
                    raise ValueError("en este caso es mejor multiplicar el renglOn por "+str(k11+k12)+". Error en "+s)
                elif (a11==a2):
                    #if (k11!=1):
                        #raise ValueError("el renglOn destino no se multiplica por "+str(k11)+". Error en "+s)
                    r=str(k12)+"R"+str(a12)+"+"+str(k11)+"R"+str(a2)+" -> R"+str(a2)
                    m[a2,a12]=k12
                    m[a2,a11]=k11
                    #print(str(m))
                elif (a12==a2):
                    #if (k12!=1):
                        #raise ValueError("el renglOn destino no se multiplica por "+str(k12)+". Error en "+s)
                    r=str(k11)+"R"+str(a11)+"+"+str(k12)+"R"+str(a2)+" -> R"+str(a2)
                    m[a2,a11]=k11
                    m[a2,a12]=k12
                    #print(str(m))
                else:#if (a11!=a2) and (a12!=a2):
                    raise ValueError("el renglón destino y un origen deben coincidir. Error en "+s)
            elif f!=-1:
                c=s[:d].count('r')
                if c==0:
                    raise ValueError("faltan 'R' en "+s[:d])
                elif c==2:
                    (k11,a11)=interpretaRenglon(s[:f])
                    (k12,a12)=interpretaRenglon(s[f+1:d])
                    if (a11==a2) and (a12==a2):
                        raise ValueError("en este caso es mejor multiplicar el renglón. Error en "+s)
                    elif (a11==a2):
                        #if (k11!=1):
                            #raise ValueError("el renglOn destino no se multiplica por "+str(k11)+". Error en "+s)
                            
                        r=str(-k12)+"R"+str(a12)+"+"+str(k11)+"R"+str(a2)+" -> R"+str(a2)
                        m[a2,a12]=-k12
                        m[a2,a11]=k11
                        #print(str(m))
                    elif (a12==a2):
                        #if (k12!=1):
                            #raise ValueError("el renglOn destino no se multiplica por "+str(k11)+". Error en "+s)
                        r=str(-k11)+"R"+str(a11)+"+"+str(k12)+"R"+str(a2)+" -> R"+str(a2)
                        m[a2,a11]=-k11
                        m[a2,a12]=k12
                        #print(str(m))
#                        raise ValueError("el renglOn destino no se multiplica por "+str(-k12)+". Error en "+s)
                    else:#if (a11!=a2) and (a12!=a2):
                        raise ValueError("el renglón destino y un origen deben coincidir. Error en "+s)
                elif c==1:
                    (k1,a1)=interpretaRenglon(s[:d])
                    if a1 != a2 :
                        raise ValueError("el renglón destino y el origen deben coincidir. Error en "+s)
                    r=str(k1)+"R"+str(a1)+" -> R"+str(a2)
                    m[a2,a2]=k1
                    #print(str(m))
                
                else:
                    raise ValueError("se operan máximo dos renglones. Error en "+s)
            else:
                (k1,a1)=interpretaRenglon(s[:d])
                if a1 != a2 :
                    raise ValueError("el renglón destino y el origen deben coincidir. Error en "+s)
                r=str(k1)+"R"+str(a1)+" -> R"+str(a2)
                m[a2,a2]=k1
                #print(str(m))
    #print("La matriz elemental asociada es \n"+str(m))
    print(r)
    #print(m)
    return m*Amat
    #sympy.pprint(m*Amat)
    
def califica_matriz(mat):
    #print(str(mat))
    i=0
    j=0
    (m,n)=mat.shape
    repetir=True
    while repetir:
        cont_ceros=0
        cont_div=0
        #print('i=',i,'j=',j,'mn',m,n)
        delant=mat[i,j]
        ceros_abajo=True
        for i1 in range(i+1,m):
            a=mat[i1,j]
            if bool(a):
                ceros_abajo=False
                #print('i1=',i1,' j=',j,' a=',a)
                try:
                    amod=a%delant
                    if amod==0:
                        cont_div=cont_div+1
                except:
                    pass
                    #print()
            else:
                cont_ceros=cont_ceros+1
        #print('i=',i,'j=',j,'mn',m,n)
        if j==n-1 or i==m-1 or not ceros_abajo:# or not bool(delant):
            repetir=False
        elif ceros_abajo:
            j=j+1
            if bool(delant):
                i=i+1
    #print(mat,delant,cont_ceros,cont_div,ceros_abajo,i,j)    
    return (delant,cont_ceros,cont_div,ceros_abajo,i,j)
        

    print()
    #cuentas cuantas columnas a la izq cumplen e1
    #la primera col tiene elem delantero
    #la primncol que no cumple e1 cuantos ceros tiene abajo
    #la primera col que no cumple cuatos no ceros divide abajo
    
#def ceros_bajo():


print("""#Reducción por Gauss-Jordan
#1 Vaya a la columna no cero extrema izquierda.
#2 Si el primer renglón tiene un cero en la columna del paso (1), intercámbielo con uno que tenga un elemento no cero en la misma columna.
#3 Obtenga ceros abajo del elemento delantero, sumando múltiplos adecuados del renglón superior a los renglones debajo de él.
#4 Cubra el renglón superior y repita el mismo proceso comenzando por el paso (1) aplicado a la sub-matriz restante. Repita este proceso con el resto de los renglones.
#5 Comenzando con el último renglón no cero, avance hacia arriba: para cada renglón obtenga un 1 delantero e introduzca ceros arriba de él, sumando múltiplos adecuados a los renglones correspondientes.""")
print()
print("Ejemplo tomado de la página 19 del libro 'Álgebra Lineal con  aplicaciones' de Nakos-Joyner de la editorial Thomson impreso en 1999.")
A=matriz([[0, 3, -6, -4, -3, -5],[-1, 3, -10, -4, -4, -2],[ 4, -9, 34, 0, 1, -21],[2, -6, 20, 2, 8, -8]])
print('La matriz ingresada fue')
print('A=matriz([[0, 3, -6, -4, -3, -5],[-1, 3, -10, -4, -4, -2],[ 4, -9, 34, 0, 1, -21],[2, -6, 20, 2, 8, -8]])')
#print('A=matriz("0 3 -6 -4 -3 -5;-1 3 -10 -4 -4 -2 ; 4 -9 34 0 1 -21;2 -6 20 2 8 -8")')
try:
    A=eval(input('Puede ingresar otra matriz, si no, solamente oprima "Enter"\n A='))
except:
    print(str(A))

repetir=True
sugerencia=0
B=A
(Am,An)=A.shape
while repetir:
    print('------------------------------------------')
    operacion=input('Escriba la operación entre renglones\n operación = ')
    try:
        C=B
        A=op(A,operacion)
        B=A
        sympy.pprint(B,use_unicode=True)
        sugerencia=0
        D=C-B
        if bool(D):
            msg="la matriz cambió"
            (Bij,Bceros,Bdiv,Babajo,Bi,Bj)=califica_matriz(B)
            (Cij,Cceros,Cdiv,Cabajo,Ci,Cj)=califica_matriz(C)
            if Bi==Am-1 or Bj==An-1:
                repetir=False
                print('Felicitaciones!! la matriz ya está en forma escalón')
            elif Bj>Cj:
                print('Muy bien, avanzó otra columna')
            elif Bj<Cj:
                print('Muy mal, perdió alguna columna')
            else:
                if (not bool(Cij)) and bool(Bij) :
                    print('Super bien, ya no es cero el elemento delantero ' + str(Bi)+','+str(Bj))
                elif (not bool(Bij)) and bool(Bij) :
                    print('Super mal, ahora es cero el elemento delantero ' + str(Bi)+','+str(Bj))
                else:
                    if Bceros>Cceros:
                        print('Bien, obtuvo un cero más')
                    elif Bceros<Cceros:
                        print('Mal, perdió algún cero')
                    else:
                        if Cdiv>0:
                            print('Hubiera podido obtener un cero')
                        elif Bdiv>0:
                            print('Ya obtuvo un divisor en el elemento delantero') 
                        else:
                            print('No entiendo para que realizó esa operación')
            #print(califica_matriz(B)
            #print(califica_matriz(C))
            #calificar las matrices
            #decir si mejora o empeora
            #explicar el porque
        else:
            msg='La matriz no cambió' 
            print(msg)
    except ValueError as e:
        operacion==operacion.replace(' ','').lower()
        if operacion=='salir':
            repetir=False
            print("hasta luego")
        elif operacion=='sugerir' or operacion=='ayuda':
            C=A.ultop
            print(str(C))
            (Cm,Cn)=C.shape
            (Cij,Cceros,Cdiv,Cabajo,Ci,Cj)=califica_matriz(C)
            CiNoCero=Ci+1
            while not bool(C[CiNoCero,Cj]):
                CiNoCero=CiNoCero+1

            if Ci==Cm or Cj==Cn:
                sug1='La matriz ya está en forma escalón'
                sug2=sug1
                sug3=sug2
                repetir=False
            elif not bool(Cij):
                sug1='Nuestro objetivo ahora es obtener un pivote en la columna '+str(Cj)+', pero observe que en el renglón '+str(Ci)+' hay un cero'
                sug2='Debe intercambiar el renglón '+str(Ci)+' con el renglón '+ str(CiNoCero)
                sug3='R'+str(Ci)+' <-> R'+str(CiNoCero)
            else:
                sug1='Ahora hay que obtener ceros bajo el pivite'
                sug2='Hay que sumar un multiplo del renglón '+str(Ci)+' al renglón '+str(CiNoCero)
                sug3=str(frac(-C[CiNoCero,Cj],C[Ci,Cj]))+'R'+str(Ci)+' + R'+str(CiNoCero)+' -> R'+str(CiNoCero)  
            sugerencia=sugerencia+1
            if sugerencia==1:
                print(sug1)
            elif sugerencia==2:
                print(sug2)
            else:
                print(sug3)
                sugerencia=0
                
        else:
            print('No entiendo: '+str(e))
            print('Algunos ejemplos de operaciones son: \n \t r1<->r2 \n \t -2r1+r2->r2 \n \t 3r2->r2')
#print("Según el paso 1, observamos que la primera columna de la matriz no es cero, por lo tanto comenzamos con esta columna, j=1")
#print("Según el paso 2, hay que mirar el primer renglón, i=1, de la columna j.  elemento de esa columna A[] intercambiar el renglón 1 por algún renglón j, que no tenga cero, a[j,1]!=0 , preferiblement si el primer e")
