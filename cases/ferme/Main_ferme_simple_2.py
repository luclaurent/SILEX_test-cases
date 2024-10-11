#############################################################################
#      Import des librairies
#############################################################################
import scipy # librairie de calcul scientifique
import scipy.linalg # librairie de calcul scientifique
import pylab # librairie de trace de courbes et de graphiques

#############################################################################
print("Calcul d'un treillis de barres")
#############################################################################

#############################################################################
#      Maillage
#############################################################################

# tableau de coordonnees des noeuds [m]
nodes=scipy.array([[0.0  ,	   0.0	], # noeud 1
                   [7.0/4.0  , 0.0	], # noeud 2
                   [7.0/2.0   , 0.0	], # ...
                   [3.0*7.0/4  , 0.0],
                   [7.0   , 0.0	],
                   [7.0/4.0   , 0.75],
                   [7.0/2.0    , 1.5],
                   [3.0*7.0/4   , 0.75]])

# table de connectivite des elements (les numeros des noeuds commencent a 1)
elements=scipy.array([[     1 ,	 2], # element 1
                      [     2 ,	 3], # element 2
                      [     3 ,	 4], # ...
                      [     4 ,	 5],
                      [     1 ,	 6],
                      [     6 ,	 7],
                      [     7 ,	 8],
                      [     8 ,	 5],
                      [     6 ,	 2],
                      [     8 ,	 4],
                      [     6 ,	 3],
                      [     7 ,	 3],
                      [     3 ,	 8]])

#############################################################################
#      Materiau et sections des barres
#############################################################################

# Module d'Young
Young  = 2.1e11 # [Pa]

# Aire de la section des barres
Section = 0.000608 # [m^2]

#############################################################################
#      Initialisation de quelques variables
#############################################################################

# nombre de noeuds
#nnodes = 8 # si on le donne directement
nnodes = nodes.shape[0] # pour un calcul automatique a partir du nombre de lignes du tableau nodes

# nombre de degres de liberte
ndof   = nnodes*2 # 2 ddl par noeud

# nombre d'elements
#nelem  = 13 # si on le donne directement
nelem = elements.shape[0] # pour un calcul automatique a partir du nombre de lignes du tableau nodes

# degres de libertes fixes (attention, python commence a 0)
Fixed_Dofs = scipy.array([0, 8, 1, 9])

# degres de libertes libres 
#SolvedDofs = scipy.array([ 2,  3,  4,  5,  6,  7, 10, 11, 12, 13, 14, 15]) # si on le donne directement
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# Initialisation des deplacements aux noeuds
U=scipy.zeros(ndof)

# Initialisation des forces aux noeuds
F=scipy.array([    0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
           0.   ,     0.   ,     0.   ,     0.   ,     0.   , -3403.755,
           0.   , -3403.755,     0.   , -3403.755])

#############################################################################
#      Calcul de la matrice de rigidite K
#############################################################################

# Initialisation d'une matrice nulle de taille ndof X ndof
K=scipy.zeros((ndof,ndof))

# Boucle sur les elements
for e in range(nelem):
    idnode_a = elements[e,0] # numero du noeud "a" de l'element "e"
    idnode_b = elements[e,1] # numero du noeud "b" de l'element "e"
    dofx_a   = (idnode_a-1)*2 # ddl ux du noeud "a" 
    dofy_a   = (idnode_a-1)*2+1 # ddl uy du noeud "a" 
    dofx_b   = (idnode_b-1)*2 # ddl ux du noeud "b" 
    dofy_b   = (idnode_b-1)*2+1 # ddl uy du noeud "b" 
    dofelem = [dofx_a,dofy_a,dofx_b,dofy_b] # liste des 4 ddls de l'element "e"

    xa=nodes[idnode_a-1,0] # coordonnee "x" du noeud "a"
    ya=nodes[idnode_a-1,1] # coordonnee "y" du noeud "a"
    xb=nodes[idnode_b-1,0] # coordonnee "x" du noeud "b"
    yb=nodes[idnode_b-1,1] # coordonnee "y" du noeud "b"

    lx        = xb-xa # difference le long de "x" des coordonnees des noeuds
    ly        = yb-ya # difference le long de "y" des coordonnees des noeuds

    lelem     = scipy.sqrt(lx**2+ly**2) # pythagore pour calculer la longueur de l'element
    cos_theta = lx/lelem # calcul du cosinus de l'angle de l'element avec l'axe "x"
    sin_theta = ly/lelem # calcul du sinus de l'angle de l'element avec l'axe "x"

    cc=(cos_theta*cos_theta)
    cs=(cos_theta*sin_theta)
    ss=(sin_theta*sin_theta)

    # remplissage de la matrice de rigidite ke de l'element "e"
    ke = scipy.mat([[cc,cs,-cc,-cs],
                    [cs,ss,-cs,-ss],
                    [-cc,-cs,cc,cs],
                    [-cs,-ss,cs,ss]])*Young*Section/lelem

    # 2 boucles imbriquees pour assembler la matrice de rigidite ke de l'element "e" dans la matrice de rigidite K du treilli
    for i in range(4): # boucle sur les 4 lignes de ke
        for j in range(4): # boucle sur les 4 colonnes de ke
            K[dofelem[i],dofelem[j]]=K[dofelem[i],dofelem[j]]+ke[i,j]


#############################################################################
#       Resolution du probleme
#############################################################################

# Pour avoir la documentation sur la librairie de resolution taper : print(scipy.linalg.solve.__doc__)
U[SolvedDofs] = scipy.linalg.solve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

#############################################################################
# Solution attendue pour le calcul de la ferme:
#############################################################################
#>>> U
#array([  0.00000000e+00,   0.00000000e+00,   1.96473284e-19,
#        -7.55959310e-04,   3.92946568e-19,  -8.57727111e-04,
#         1.96473284e-19,  -7.55959310e-04,   0.00000000e+00,
#         0.00000000e+00,   1.13706719e-04,  -7.55959310e-04,
#        -2.12884826e-19,  -8.17739387e-04,  -1.13706719e-04,
#        -7.55959310e-04])
