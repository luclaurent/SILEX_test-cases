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
nodes=scipy.array([[0.0  ,	   0.0	], # x,y du noeud 1
                   [  ???? , ????	], # x,y du noeud 2
                   [   ????,???? ], # ...
                   [  ???? , ????],
                   [  ???? ,???? ]])

# table de connectivite des elements (les numeros des noeuds commencent a 1)
elements=scipy.array([[     1 ,	 2],# element 1
                      [   ????   ,????	 ], # element 2
                      [   ????   ,????	 ], # ...
                      [     ???? ,????	 ]])

#############################################################################
#      Materiau et sections des barres
#############################################################################

# Module d'Young
Young  = ???? # [Pa]

# Aire de la section des barres
Section = ???? # [m^2]

#############################################################################
#      Initialisation de quelques variables
#############################################################################

# nombre de noeuds
nnodes = ???? # si on le donne directement
#nnodes = nodes.shape[0] # pour un calcul automatique a partir du nombre de lignes du tableau nodes

# nombre de degres de liberte
ndof   = nnodes*???? # 2 ddl par noeud

# nombre d'elements
nelem  = ???? # si on le donne directement
#nelem = elements.shape[0] # pour un calcul automatique a partir du nombre de lignes du tableau nodes

# degres de libertes fixes (attention, python commence a 0)
Fixed_Dofs = scipy.array([0, ????, ????, ....])

# degres de libertes libres 
SolvedDofs = scipy.array([ ????,  ????, ....]) # si on le donne directement
#SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# Initialisation des deplacements aux noeuds
U=scipy.zeros(ndof)

# Initialisation des forces aux noeuds
F=scipy.array([    ????  ,    ????   ,    ????   , .....])

#############################################################################
#      Calcul de la matrice de rigidite K
#############################################################################

# Initialisation d'une matrice nulle de taille ndof X ndof
K=scipy.zeros((ndof,ndof))

# Boucle sur les elements
for e in range(nelem):
    idnode_a = elements[e,0] # numero du noeud "a" de l'element "e"
    idnode_b = elements[????,????] # numero du noeud "b" de l'element "e"
    dofx_a   = (idnode_a-1)*2 # ddl ux du noeud "a" 
    dofy_a   = (idnode_a-1)*2+1 # ddl uy du noeud "a" 
    dofx_b   = ???? # ddl ux du noeud "b" 
    dofy_b   = ???? # ddl uy du noeud "b" 
    dofelem = [dofx_a,dofy_a,????,????] # liste des 4 ddls de l'element "e"

    xa=nodes[idnode_a-1,0] # coordonnee "x" du noeud "a"
    ya=nodes[????,????] # coordonnee "y" du noeud "a"
    xb=nodes[????,????] # coordonnee "x" du noeud "b"
    yb=nodes[????,????] # coordonnee "y" du noeud "b"

    lx        = xb-xa # difference le long de "x" des coordonnees des noeuds
    ly        = ????-???? # difference le long de "y" des coordonnees des noeuds

    lelem     = scipy.sqrt(????) # pythagore pour calculer la longueur de l'element
    cos_theta = ???? # calcul du cosinus de l'angle de l'element avec l'axe "x"
    sin_theta = ???? # calcul du sinus de l'angle de l'element avec l'axe "x"

    # remplissage de la matrice de rigidite ke de l'element "e"
    ke = scipy.mat([[????,????,????,????],
                    [????,????,????,????],
                    [????,????,????,????],
                    [????,????,????,????]])

    # 2 boucles imbriquees pour assembler la matrice de rigidite ke de l'element "e" dans la matrice de rigidite K du treilli
    for i in range(4): # boucle sur les 4 lignes de ke
        for j in range(4): # boucle sur les 4 colonnes de ke
            K[dofelem[i],dofelem[j]]=K[dofelem[i],dofelem[j]]+ke[????,????]


#############################################################################
#       Resolution du probleme
#############################################################################

# Pour avoir la documentation sur la librairie de resolution taper : print(scipy.linalg.solve.__doc__)
U[SolvedDofs] = scipy.linalg.solve(K[SolvedDofs,:][:,SolvedDofs],????)

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
