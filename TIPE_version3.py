import numpy as np
import xlrd
import matplotlib.pyplot as plt

## Dimensions et discrétisation de la grille
x_max = 60000 #dimension de l'espace sur lequel se propage les rayons
z_max = -1500

dx = 5 # dimensions des cases discrétisées
dz = 5

x_profil1 = 10000 # Pour le tracé du profil de célérité à l'abscisse donnée
x_profil2 = 50000 # Pour le tracé du profil de célérité à l'abscisse donnée

## Réglage du sonar
ouverture_sonar = np.pi/40 # angle d'ouverture du sonar

x_init = 0 # position du sonar
z_init = -40

theta_init = np.pi/30 # Pour le tracé d'un seul rayon
nbre_rayons = 20 # Nombre de rayons tirés par le sonar


## Choix de l'environnement
env1 = "1" # choix du profil de température (récupéré sur les flotteurs ARGO)
env2 = "2"
changement = x_max // 2 # abscisse à laquelle se fait la rupture de profil


# Création de listes contentant les données de température et profondeur sur les fichiers Excel
def creation_profil(env): 
    t="Database/model" + env + ".xls"
    doc=xlrd.open_workbook(t)
    Feuille=doc.sheet_by_index(0)
    rows=Feuille.nrows 
    cols=Feuille.ncols
    Abscisse1 = []
    Ordonnee1 = []
    
    for i in range(1,rows):
        Abscisse1.append(float(Feuille.cell_value(rowx=i, colx=1)))
        Ordonnee1.append(-float(Feuille.cell_value(rowx=i, colx=0)))
        
    # print(Abscisse1[a]) # renvoie tjr la meme valeur
    return Abscisse1, Ordonnee1
    

Abscisse_mois1, Ordonnee_mois1 = creation_profil(env2)
Abscisse_mois2, Ordonnee_mois2 = creation_profil(env1)



## Insertion de la fonction température avec récupération de données
def recherche(z,A): #prend en entrée un z négatif
    a=0
    b=len(A)-1
    m=(a+b)//2
    # print(a,b)
    # if z == 0:
        # return 0
    while a < b :
        # print(a,b,m,A[m],z)
        if A[m] == z:
            return m
        elif A[m]>z :
            a=m+1
        else:
            b=m-1
        m=(a+b)//2
    return m


def temperature2(x,z): # z = profondeur
    if x <= changement:        
        a= recherche(z,Ordonnee_mois1)
        return Abscisse_mois1[a]
    else:
        a = recherche(z,Ordonnee_mois2)
        return Abscisse_mois2[a]
    

## Salinité
def salinite(z,x):
    return 38

## Pression
def pression(z,x):
    return 1+0.0981 * (-z)

## Formule simplifiée calculant la célérité d'une onde dans l'eau de mer
def celerite2(z,x):
    celerite = 1449.2 + 4.6*temperature2(x,z) - 0.055*(temperature2(x,z)**2) + 0.00029 * (temperature2(x,z)**3) + (1.34 - 0.010*temperature2(x,z))*(salinite(z,x) - 35) + 1.58 * 10**(-6)*(pression(z,x)*10**5) #la formule simplifiée choisise prend la pression en Pascals.
    return celerite


## Représentation de la célérité d'une onde en fonction de la profondeur et du profil
def representation_celerite(x):
    Abscisse1 = []
    Ordonnee1 = []
    for i in range(-z_max):
        Abscisse1.append(celerite2(-i,x))
        Ordonnee1.append(-i)
    return Abscisse1, Ordonnee1

Abscisse1, Ordonnee1 = representation_celerite(x_profil1)
Abscisse2, Ordonnee2 = representation_celerite(x_profil2)

plt.plot(Abscisse1,Ordonnee1)
plt.xlabel("célérité de l'onde (m/s)")
plt.ylabel("Profondeur en m")
plt.axis([1500,1550,z_max + 1,0])
plt.show()

plt.plot(Abscisse2,Ordonnee2)
plt.xlabel("célérité de l'onde (m/s)")
plt.ylabel("Profondeur en m")
plt.axis([1500,1550,z_max + 1,0])
plt.show()



## Loi de Snell-Descartes
def theta2(theta1,c1,c2):
    # print((c2/c1)*sin(theta1),theta1,c1,c2)
    if (c2/c1)*np.sin(theta1) >= 1 or (c2/c1)*np.sin(theta1) <= -1:
        return - theta1
    return np.arcsin((c2/c1)*np.sin(theta1))

## Reflexion
def reflexion(theta):
    return - theta

## Tangente a la courbe
def image_tangente(theta,xn,zn,x):
    a = np.tan(theta)
    return a*x + (zn - a*xn)

def antecedent_tangente(theta,xn,zn,z):
    a = np.tan(theta)
    return (z - zn + a*xn)*(1/a)


## Corps du programme : Propagation d'un rayon de case en case
def propagation(theta_init,x_init,z_init):

    abscisse = []
    ordonnee = []

    z1, z2 = z_init, z_init
    x1, x2 = x_init, x_init
    theta = theta_init

    dzn = z_init + (-z_init % dz)
    dxn = x_init - (x_init % dx)


    while x1 < x_max and z1 >= z_max:
        if z1 == 0 and theta > 0:
           theta = reflexion(theta)
        # elif z1 == z_max and theta < 0: # ajoute la reflexion au fond de l'eau
            # theta = reflexion(theta)
        elif x1 % dx == 0 and z1 % dz != 0:
            if image_tangente(theta,x1,z1,dxn + dx) >= dzn:
                z2 = dzn
                x2 = antecedent_tangente(theta,x1,z1,z2)
                c1 = celerite2(dzn,dxn)
                c2 = celerite2(dzn + dz,dxn)
                thetabis = theta2(np.pi/2 - theta,c1,c2)
                if thetabis > 0:
                    theta = np.pi/2 - thetabis
                else:
                    theta = - theta
            elif image_tangente(theta,x1,z1,dxn + dx) <= (dzn - dz):
                z2 = dzn - dz
                x2 = antecedent_tangente(theta,x1,z1,z2)
                c1 = celerite2(dzn,dxn)
                c2 = celerite2(dzn - dz,dxn)
                thetabis = theta2(np.pi/2 - abs(theta),c1,c2)
                if thetabis > 0:
                    theta = -(np.pi/2 - thetabis)
                else:
                    theta = - theta
            else:
                x2 = dxn + dx
                z2 = image_tangente(theta,x1,z1,x2)
                c1 = celerite2(dzn,dxn)
                c2 = celerite2(dzn,dxn + dx)
                theta = theta2(theta,c1,c2)
        elif z1 % dz == 0:
            if theta >= 0:
                if image_tangente(theta,x1,z1,dxn + dx) <= (dzn + dz):
                    x2 = dxn + dx
                    z2 = image_tangente(theta,x1,z1,x2)
                    c1 = celerite2(dzn,dxn)
                    c2 = celerite2(dzn,dxn + dx)
                    theta = theta2(theta,c1,c2)
                else:
                    z2 = dzn + dz
                    x2 = antecedent_tangente(theta,x1,z1,z2)
                    c1 = celerite2(dzn,dxn)
                    c2 = celerite2(dzn + dz,dxn)
                    thetabis = theta2(np.pi/2 - theta,c1,c2)
                    if thetabis > 0:
                        theta = np.pi/2 - thetabis
                    else:
                        theta = - theta
            elif theta < 0:
                if image_tangente(theta,x1,z1,dxn + dx) >= (dzn - dz):
                    x2 = dxn + dx
                    z2 = image_tangente(theta,x1,z1,x2)
                    c1 = celerite2(dzn,dxn)
                    c2 = celerite2(dzn,dxn + dx)
                    theta = theta2(theta,c1,c2)
                else:
                    z2 = dzn - dz
                    x2 = antecedent_tangente(theta,x1,z1,z2)
                    c1 = celerite2(dzn,dxn)
                    c2 = celerite2(dzn - dz,dxn)
                    thetabis = theta2(np.pi/2 - abs(theta),c1,c2)
                    if thetabis > 0:
                        theta = -(np.pi/2 - thetabis)
                    else:
                        theta = - theta
        abscisse.append(x2)
        ordonnee.append(z2)
        x1, z1 = x2, z2
        dzn = z1 + (-z1 % dz)
        dxn = x1 - (x1 % dx)
        # print(x1,z1,theta)

    return abscisse, ordonnee


## Representation du trajet d'un seul rayon
# abscisse, ordonnee = propagation(theta_init,x_init,z_init)

# plt.plot(abscisse,ordonnee)
# plt.axis([0,x_max + 1,z_max + 1,0])
# plt.show()
# print("This is the end")


## Propagation de plusieurs rayons
def propagation2(nbre_rayons,x_init,z_init,ouverture_sonar):
    total_abscisse = []
    total_ordonnee = []
    discretisation_angle = (1/(nbre_rayons+1))*(ouverture_sonar)
    # print("aa",((pi/30)/(nbre_rayons+1)))
    for i in range(1,nbre_rayons + 1):
        # print(pi/60 - i*discretisation_angle)
        abscisse2, ordonnee2 = propagation((ouverture_sonar)/2 - i*discretisation_angle,x_init,z_init)
        total_abscisse.append(abscisse2)
        total_ordonnee.append(ordonnee2)
    return total_abscisse, total_ordonnee

## Tracé de plusieurs rayons
total_abscisse, total_ordonnee = propagation2(nbre_rayons,x_init,z_init,ouverture_sonar)

for j in range(nbre_rayons):
    plt.plot(total_abscisse[j],total_ordonnee[j])
    plt.axis([0,x_max + 1,z_max + 1,0])
plt.xlabel("Distance (m)")
plt.ylabel("Profondeur (m)")
plt.show()


