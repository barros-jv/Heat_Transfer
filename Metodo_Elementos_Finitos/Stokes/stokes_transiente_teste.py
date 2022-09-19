"""
Projeto CFD Elementos Finitos
Problema de Stokes transiente.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.cm as cm
import matplotlib.tri as mtri
import time
import timeit

ini = time.time()
#----------------------------------------------------------
# Malha
# Elemento mini pede malha triangular

Lx = 2.0                        # comprimento da malha em x
Ly = 1.0                        # comprimento da malha em y
nx = 30                          # número de nós em x
ny = 30                         # número de nós em y
npoints = nx*ny                  # número de vértices
nelem = 2*(nx-1)*(ny-1)         # número de espaços triangulares e centróides
nn = npoints + nelem            # número de nós (vértices + centróides)
# print('nn = ', nn)

#----------------------------------------------------------
# Dados do problema
rho = 1.0
ni = 1.0
dt = 0.01
vx_ini = 1.0
vy_ini = 0.0


#----------------------------------------------------------
# criação dos pontos em x e y

xv = np.linspace(0, Lx, nx)
yv = np.linspace(0, Ly, ny)

[X,Y] = np.meshgrid(xv,yv)

X = np.reshape(X,npoints)
Y = np.reshape(Y,npoints)


#----------------------------------------------------------
# Criação da IEN triângular

Tri = matplotlib.tri.Triangulation(X,Y)
IEN = Tri.triangles
ne = IEN.shape[0]
# print(IEN)
# print(ne)

# #----------------------------------------------------------
# # plot da malha de triangulos

# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_title("Malha")

# ax = plt.triplot(X,Y,IEN)
# ax = plt.plot(X,Y,'bo', markersize=2.8, linewidth=1.1)
# for i in range(0,npoints):
#     plt.text(X[i]+0.02,Y[i]+0.03,str(i),color='b')
# plt.show()

#----------------------------------------------------------
# Automatização da definição dos pontos da ICC
# 

# pontos inferiores
iccinf = []
for i in range (0,nx):
    iccinf.append(i)
ICCINF = np.array(iccinf)

# pontos superiores
iccsup = []
for i in range((nx*(ny-1)), npoints):
    iccsup.append(i)
ICCSUP = np.array(iccsup)

# pontos a esquerda
iccesq = []
for i in range(nx, npoints-nx, nx):
    iccesq.append(i)
ICCESQ = np.array(iccesq)

# pontos a direita
iccdir = []
for i in range(2*nx-1, npoints-1, nx):
    iccdir.append(i)
ICCDIR = np.array(iccdir)

# alterar saída a direita
# sempre tem que criar a lista iccout
# metade inferior a direita fica sem fluxo 
iccout = []
wall = npoints//2
# # cria parede na parte de baixo
# for i in range(2*nx-1, wall, nx):
#     iccout.append(i)

# # cria parede na parte de cima
# for i in range(((npoints-nx)-1),wall, -nx):
#     iccout.append(i)

# # coloca orifícios na saída
# for i in range(2*nx-1, npoints-1, 2*nx):
#     iccout.append(i)

# # divide os espaços igualmente em 3 e o meio é a parede
# lim = (ny-2)//3
# iccout = iccdir[lim:-lim]

# corrige os elementos da iccdir e iccout
for i in iccout:
    val_remove = i
    iccdir.remove(val_remove)

icc = iccinf + iccsup + iccesq + iccdir + iccout
ICC = np.array(icc)

#----------------------------------------------------------
# # plot com números das CC's
# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_title("Malha")

# ax = plt.triplot(X,Y,IEN)
# ax = plt.plot(X,Y,'bo', markersize=2.8, linewidth=1.1)
# for i in icc:
#     plt.text(X[i]+0.02,Y[i]+0.03,str(i),color='b')
# plt.show()
# # print(icc)


#----------------------------------------------------------
# Criação do centróide
# Primeiro passo, aumentar as matrizes X e Y para receber os valores dos centróides
# Segundo passo, ajustar a IEN para receber os pontos dos centróides
# Terceiro, criar os pontos e inserir na posição correta
X = X.copy()
Y = Y.copy()
# print('Y = ', len(Y))
X.resize(npoints+nelem,refcheck=False)
Y.resize(npoints+nelem,refcheck=False)
# print(X)
# print(len(X))
# print(Y)
# IEN original com 3 colunas. Adiciona uma coluna com valores 0. IEN fica com 4 colunas
# adiciona   matriz 1, matriz 2 (dimensão da matriz), tipo do argumento.     
# np.hstack((  IEN   , np.zeros((ne      , 1       ), dtype=IEN.dtype)))
IEN = np.hstack((IEN, np.zeros((ne, 1), dtype=IEN.dtype)))
# print(IEN)
malha = time.time()
t_malha = malha - ini

#----------------------------------------------------------
# Definir matrizes K, M e G, velocidade e pressão.
vx = np.zeros( (nn), dtype='float')
vy = np.zeros( (nn), dtype='float')
p = np.zeros( (npoints), dtype='float')

K = np.zeros( (nn,nn), dtype='float')
M = np.zeros( (nn,nn), dtype='float')
# print(K.shape)
# print('K = ', len(K))
# print(M.shape)
# print(len(M))

Gx = np.zeros( (nn, npoints), dtype='float')
Gy = np.zeros( (nn, npoints), dtype='float')
# print('Gx = ', len(Gx))
# print('Gy = ', len(Gy))

Dx = np.zeros( (npoints, nn), dtype='float')
Dy = np.zeros( (npoints, nn), dtype='float')
D = np.hstack((Dx, Dy))
# print('D = ', len(D))
# print('Dshape = ', D.shape)

#----------------------------------------------------------
# Fazer loop para calcular matrizes K, M e G.
for i in range(0, nelem):
    # pega os pontos conforme a posição que aparecem na IEN
    [v1,v2,v3,v4] = IEN[i]
    # O valor de X e Y está na posição do ponto retirado da IEN
    Xv1 = X[v1]
    Xv2 = X[v2]
    Xv3 = X[v3]
    Yv1 = Y[v1]
    Yv2 = Y[v2]
    Yv3 = Y[v3]
    # Calcula o valor em x e y do centróide
    Xc = (Xv1 + Xv2 + Xv3)/3.0
    Yc = (Yv1 + Yv2 + Yv3)/3.0
    X[npoints + i] = Xc
    Y[npoints + i] = Yc
    # Adiciona o ponto centróide na posição correta da IEN
    IEN[i,3] = npoints + i    

    # calcular área do triângulo (a, b, c) e peso z para usar nas matrizes K, M, e G
    a1 = Xv2*Yv3 - Xv3*Yv2
    a2 = Xv3*Yv1 - Xv1*Yv3
    a3 = Xv1*Yv2 - Xv2*Yv1

    b1 = Yv2 - Yv3
    b2 = Yv3 - Yv1
    b3 = Yv1 - Yv2
    
    c1 = Xv3 - Xv2
    c2 = Xv1 - Xv3
    c3 = Xv2 - Xv1

    area = (1.0/2.0)*(Xv2*Yv3 + Xv3*Yv1 + Xv1*Yv2 - Xv2*Yv1 - Xv3*Yv2 - Xv1*Yv3)
    # print('area = ', area)

    # matriz de massa
    melem = (area/840)*np.array([[83, 13, 13, 45],
                                [13, 83, 13, 45],
                                [13, 13, 83, 45],
                                [45, 45, 45, 243]])

    # matriz de rigidez (laplaciano) e matriz dos coeficientes do laplaciano
    B = (1.0/(2.0*area))*np.array([[b1, b2, b3],
                                    [c1, c2, c3]])
    
    BT = B.transpose()

    # peso para ajustar as matrizes mini
    z = (1/(4*area))*(b2**2 + b2*b3 + b3**2 + c2**2 + c2*c3 + c3**2)
    # print('z = ', z)
    k = area*(BT@B)
    # print(k)

    # definição das matrizes mini em função das matrizes lineares
    # matriz k
    kelem_lin = area*(BT@B) + (9.0/10.0)*z
    k_mini = kelem_lin
    k_mini = np.hstack( (k_mini, (-27.0*z/10.0)*np.ones((k_mini.shape[0], 1), dtype=k_mini.dtype)))
    k_mini = np.vstack( (k_mini, (-27.0*z/10.0)*np.ones((k_mini.shape[1]), dtype=k_mini.dtype)))
    k_mini[-1,-1] = 81.0*z/10.0

    # matrizes gx e gy
    gx_elem_lin = (1.0/6.0)* np.array([[b1,b2,b3],
                                    [b1,b2,b3],
                                    [b1,b2,b3]])
    gx_elem_linT = gx_elem_lin.transpose()
    gx_mini = (9.0/20.0)*gx_elem_lin + gx_elem_linT     
    b = np.array([b1, b2, b3])   
    gx_mini = np.vstack( (gx_mini, (-9.0/40.0)*b*np.ones((gx_mini.shape[1]), dtype=gx_mini.dtype)))
    # gx_miniT = gx_mini.transpose()

    gy_elem_lin = (1.0/6.0)* np.array([[c1,c2,c3],
                                    [c1,c2,c3],
                                    [c1,c2,c3]])
    gy_elem_linT = gy_elem_lin.transpose()
    gy_mini = (9.0/20.0)*gy_elem_lin + gy_elem_linT
    c = np.array([c1, c2, c3])
    gy_mini = np.vstack( (gy_mini, (-9.0/40.0)*c*np.ones((gy_mini.shape[1]), dtype=gy_mini.dtype)))
    # gy_miniT = gy_mini.transpose()

    for ilocal in range(0,4):
        iglobal = IEN[i, ilocal]
        for jlocal in range(0,4):
            jglobal = IEN[i, jlocal]

            K[iglobal, jglobal] += (k_mini[ilocal, jlocal])
            M[iglobal, jglobal] += melem[ilocal, jlocal]

        for jlocal in range(0,3):
            jglobal = IEN[i, jlocal]

            Gx[iglobal, jglobal] += gx_mini[ilocal, jlocal]
            Gy[iglobal, jglobal] += gy_mini[ilocal, jlocal]


# Gerar as matrizes Dx e Dy (transpor matrizes Gx e Gy)
Dx = Gx.transpose()
Dy = Gy.transpose()


#----------------------------------------------------------
# Criar a matriz A para solução do problema linear
# A = junção das matrizes M, K, G, D (transposta de G) e nulas para completar
# Método de unir matrizes é mais rápido que fazer o loop para criá-la. 90% mais rápido.
# Sol = np.linalg.solve(A, b)
# Multiplicar os dados físicos (rho, ni, etc na montagem)
Mdt = M/dt
# print('M = ', M)
# print('Mdt = ', Mdt)

# # ini_stack = time.time()
A1 = np.hstack(((Mdt+ni*K),np.zeros((K.shape), dtype=K.dtype), (-1/rho)*Gx))
A2 = np.hstack((np.zeros((K.shape),dtype=K.dtype),(Mdt+ni*K), (-1/rho)*Gy))
A3 = np.hstack((Dx, Dy, np.zeros((npoints,npoints), dtype=K.dtype)))
A = np.vstack((A1, A2, A3))
# print('A = ', A)
# print(A.shape)
V_total = vx + vy



matriz = time.time()
t_matriz = matriz - malha

#----------------------------------------------------------
# Stokes transiente
# Resolver sistema linear A.x = b
# Vetor b e x recebem os valores de vx, vy e p. A posição no vetor é atribuída à variável.
# Condições de Contorno para A, vx, vy e pressão
# npoints = número de vértices
# nelem = número de espaços triangulares e centróides
# nn = npoints + nelem - número de nós (vértices + centróides)
# vx = tamanho nn
# vy = tamanho nn
# p = tamanho npoints


# Condições de contorno na matriz A
# Condições de contorno nas paredes (iccinf e iccsup)

for i in iccsup:
    A[i,:] = 0.0                    # zera as linhas referentes a vx
    A[i,i] = 1.0                    # coloca 1 na diagonal referente a vx
    A[i+nn, :] = 0.0                # zera as linhas referentes a vy
    A[i+nn, i+nn] = 1.0             # coloca 1 na diagonal referente a vy


for i in iccinf:
    A[i,:] = 0.0                    # zera as linhas referentes a vx
    A[i,i] = 1.0                    # coloca 1 na diagonal referente a vx
    A[i+nn, :] = 0.0                # zera as linhas referentes a vy
    A[i+nn, i+nn] = 1.0             # coloca 1 na diagonal referente a vy

# Condições de contorno entrada (iccesq)
for i in iccesq:
    A[i,:] = 0.0                    # zera as linhas referentes a vx
    A[i,i] = 1.0                    # coloca 1 na diagonal referente a vx
    A[i+nn,:] = 0.0                 # zera as linhas referentes a vy
    A[i+nn, i+nn] = 1.0             # coloca 1 na diagonal referente a vy

# Condições de contorno na saída (iccdir)
for i in iccdir:
    A[i+2*nn,:] = 0.0               # zera as linhas referentes a p
    A[i+2*nn, i+2*nn] = 1.0         # coloca 1 na diagonal referente a p

# Colocar barreiras na saída
for i in iccout:
    A[i,:] = 0.0                    # zera as linhas referentes a vx
    A[i,i] = 1.0                    # coloca 1 na diagonal referente a vx
    A[i+nn, :] = 0.0                # zera as linhas referentes a vy
    A[i+nn, i+nn] = 1.0             # coloca 1 na diagonal referente a vy
   
    


# print('A c/ CC = ', A[0:nn])
# print('A pós solução = ', A[0:nn])

sol = np.zeros((2*nn+npoints))
b = np.hstack((vx, vy, p))
# print('b  e sol = ', len(b), len(sol))
# vx = sol[0:nn]
# vy = sol[nn: 2*nn]
# p = sol[2*nn:]

# Condições de contorno para b e resolve o sistema linear
t_total_iter = 0.0
list = []
for t in range(0,31):
    ini_iter = time.time()
    Vx = sol[0:nn]
    Vy = sol[nn: 2*nn]
    b[0::] = 0.0
    b[0:nn] = Mdt@Vx
    b[nn:2*nn] = Mdt@Vy
    
    # print('Vx = ', Vx)
    # print('Vy = ', Vy)
    # print('p = ', p)
    # print('b = ', b)

    # definir condições de contorno para b
    for i in icc:
        # Condições de contorno nas paredes (iccinf e iccsup)
        if i in iccsup:
            b[i] = 0.0                      # velocidade em vx na CC
            b[i+nn] = 0.0                   # velocidade em vy na CC
        if i in iccinf:
            b[i] = 0.0                      # velocidade em vx na CC
            b[i+nn] = 0.0                   # velocidade em vy na CC
        if i in iccesq:
            b[i] = 1.0                      # velocidade em vx na CC
            b[i+nn] = 0.0                   # velocidade em vy na CC
        if i in iccdir:
            b[i+2*nn] = 0.0                 # pressão na CC
        if i in iccout:
            b[i] = 0.0                      # velocidade em vx na CC
            b[i+nn] = 0.0                   # velocidade em vy na CC

    sol = np.linalg.solve(A,b)
    
    fim_iter = time.time()
    t_iter = fim_iter - ini_iter
    t_total_iter += t_iter

    # # print('sol = ', sol)
    # plot 2D cor (triangulo)
    vx = sol[0:npoints]
    S = vx
    Z = S.reshape(nx, ny)
    ax = plt.axes()
    ax.set_aspect('equal')
    ax.triplot(Tri,'ko-', markersize=0.1, linewidth=0.1)
    ax.tricontourf(Tri,vx.reshape(npoints),cmap='jet')
    # ax.tricontourf(Tri,vy.reshape(npoints),cmap='jet')
    # ax.tricontourf(Tri,p.reshape(npoints),cmap='jet')

    surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                      cmap=matplotlib.cm.jet, extent=(X.min(),
                      X.max(), Y.min(), Y.max()))
    plt.colorbar(surf,shrink=1.0, aspect=20)
    # plt.pause(dt)
    # plt.clf()
    plt.savefig('t_30' + str(t))
    plt.close()

    # list.append(max(vx))

    # plt.show()

Vx = sol[0:nn]
vx = sol[0:npoints]

Vy = sol[nn:2*nn]
vy = sol[nn:nn+npoints]

p = sol[2*nn:]

fim = time.time()
t_loop= fim - matriz

T_total = fim - ini
# print(list)
print('Nº de iterações = ', t)
print('dt = ', dt)
print('Vx máxima = ', max(vx))
print('Lx = ', Lx)
print(nx, 'X', ny)
print('Tempo malha = ', t_malha)
print('Montagem matrizes = ', t_matriz)
print('Tempo total das iterações = ', t_total_iter)
print('Solução sistema linear = ', t_loop)
print('Tempo total = ', T_total)

