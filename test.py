import math

## Dumbass CIV102 Project Code ##

import matplotlib.pyplot as plt
import numpy as np

n = 1000
L = 1280
x = 1
E = 4000
MU = 0.2
supportAlocation = 30/2

def buildSFD(xP, P):
    FsupportA = P - P*((xP-supportAlocation)/1060)
    FsupportB = P*(xP-(supportAlocation))/1060
    plotpoints = {}
    for i in np.arange(0, L+x):
        if i < xP:
            if i < supportAlocation:
                SFDvalue = 0
            elif i < supportAlocation + 1060:
                SFDvalue = FsupportA
            else:
                SFDvalue = FsupportA + FsupportB
        else:
            if i < supportAlocation:
                SFDvalue = -1 * P
            elif i < supportAlocation + 1060:
                SFDvalue = FsupportA - P
            else:
                SFDvalue = 0
        plotpoints[i] = SFDvalue
    return plotpoints

def printSFD(xP, P):
    plt.title("Bridge SFD") 
    plt.xlabel("x (mm)") 
    plt.ylabel("shear force (kN)") 
    plt.plot(*zip(*sorted(buildSFD(xP, P).items())))
    plt.show()

def buildBMD(xP, P):
    BMDvalue = 0
    plotpoints = []
    for i in np.arange(0, L):
        BMDvalue += buildSFD(xP, P)[i]
        plotpoints.append(BMDvalue)
    return plotpoints

def printBMD(xP, P):
    plt.title("Bridge BMD") 
    plt.xlabel("x (mm)") 
    plt.ylabel("bending moment (kNm)") 
    # plt.plot(*zip(*sorted(buildBMD(xP, P))))
    plt.plot(buildBMD(xP, P))
    plt.show()    

def ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    if(shape == 0):#shape = 0 during test bridge    
            area1 = widthtop * topthickness #top part of bridge
            area2 = 10 * topthickness #tabs 10 is tab width
            area3 = leftthickness * (height - 2*topthickness)
            area4 = (widthbottom - 2*leftthickness)
            y_bar = (area1*(height-topthickness/2)) + 2*(area2-(topthickness-topthickness/2)) + 2*(area3*((height-2*topthickness)/2)) + area4*(bottomthickness/2)
            y_bar = y_bar / (area1+area2+area3+area4)


def I0(b, h):
    I0 = (b * h**3)/12
    return I0

def second_moment_of_area(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    centroid = ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape)
    area1 = widthtop * topthickness #top part of bridge
    area2 = 10 * topthickness #tabs 10 is tab width
    area3 = leftthickness * (height - 2*topthickness)
    area4 = (widthbottom - 2*leftthickness)


    term1 = I0(widthtop, topthickness) + (area1 * ((height - (topthickness/2) - centroid)**2)) #top
    term2 = I0(10, topthickness) + area2 * ((height - topthickness - topthickness/2 - centroid)**2) #tabs
    term3 = I0(leftthickness, height - 2*topthickness) + area3 * (centroid- ((height-2*topthickness)/2)**2)
    term4 = I0(widthbottom - 2*leftthickness, bottomthickness) + (area4 * (((bottomthickness/2) - centroid)**2))
    second_moment_of_area = term1 + 2* term2 + 2*term3 + term4
    return second_moment_of_area

def first_moment_of_area(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness):
    
    centroid = ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness)
    
    if 0 <= centroid < bottomthickness:
        centroidlocation = "bottom"
    elif bottomthickness <= centroid < height - topthickness:
        centroidlocation = "middle"
    elif height - topthickness <= centroid < height:
        centroidlocation = "top"

    if centroidlocation == "bottom":
        Q = widthbottom * centroid * (centroid/2)
    elif centroidlocation == "top":
        Q = widthtop * (height - centroid) * ((height - centroid)/2)
    elif centroidlocation == "middle":
        term1 = widthbottom * bottomthickness * (centroid - (bottomthickness/2))
        term2 = rightthickness * (centroid - bottomthickness) * ((centroid - bottomthickness)/2)
        term3 = leftthickness * (centroid - bottomthickness) * ((centroid - bottomthickness)/2)
        Q = term1 + term2 + term3
    
    return Q

def get_properties():
    """
    user input to say what these properties are at the locations
    """
    height = [75]*L
    widthTop = [100]*L
    widthBottom = [80]*L
    topThickness = [1.27]*L
    bottomThickness = [1.27]*L
    rightThickness = [1.27]*L
    leftThickness = [1.27]*L
    shape = [0] * L

    return  height, widthTop, widthBottom, topThickness, bottomThickness, rightThickness, leftThickness, shape
   
def y_top(y_bar, height):
    return y_bar - height
def geometric_properties():
    """
    returns I , y_bar, Q, yTop
    """
    height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape = get_properties()
    y_bar = []
    I = []
    Q = []
    yTop = []

    for i in range(L):
        y_bar.append(ybar(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i], shape[i]))
        I.append(second_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i], shape[i]))
        Q.append(first_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i]))
        yTop.append(y_bar[i] - height[i])
    

    return I, y_bar, Q, yTop

        
def V_fail(Tau):
    
    height = get_properties()[0]
    widthtop = get_properties()[1]
    widthbottom = get_properties()[2]
    topthickness = get_properties()[3]
    bottomthickness = get_properties()[4]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]

    V_fail_values = []

    for i in np.arange(0, L, x):

        I = second_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        Q = first_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        centroid = ybar(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        
        if 0 <= centroid < bottomthickness:
            centroidlocation = "bottom"
        elif bottomthickness <= centroid < height - topthickness:
            centroidlocation = "middle"
        elif height - topthickness <= centroid < height:
            centroidlocation = "top"    

        if centroidlocation == "bottom":
            b = widthbottom
        elif centroidlocation == "middle":
            b = rightthickness + leftthickness
        elif centroidlocation == "top":
            b = widthtop
        
        V_fail = I*Tau*b/Q

        V_fail_values.append(V_fail)

    return V_fail_values

def V_buck(E, mu, max_diaphragm_distance):

    height = get_properties()[0]
    widthtop = get_properties()[1]
    widthbottom = get_properties()[2]
    topthickness = get_properties()[3]
    bottomthickness = get_properties()[4]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]

    V_buck_values = []

    for i in np.arange(0, L, x):
        
        I = second_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        Q = first_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        centroid = ybar(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])

        T_crit = ((5*(3.14**2)*E)/(12(1-(mu**2))))*(((rightthickness/max_diaphragm_distance)**2) + ((leftthickness/(height-topthickness))**2))

        if 0 <= centroid < bottomthickness:
            centroidlocation = "bottom"
        elif bottomthickness <= centroid < height - topthickness:
            centroidlocation = "middle"
        elif height - topthickness <= centroid < height:
            centroidlocation = "top"    

        if centroidlocation == "bottom":
            b = widthbottom
        elif centroidlocation == "middle":
            b = rightthickness + leftthickness
        elif centroidlocation == "top":
            b = widthtop

        V_buck = T_crit*I*b/Q

        V_buck_values.append(V_buck)

    return V_buck_values



def M_failMatT(sigT, BMD):
    global L
    """
    calculate moment the moment the causes tension failure at each x
    returns array of the moment required to break bridge at each location
    
    sectional_properties array with I, Q, yBar, Ybot, Ytop
        each indivisual properties should be an array with the value at each location along the bridge
    sigT = sigma tentsion
    BMD 1-d array of the bending moments at each mm of the bridge

    """
    # change with right indices
    height = get_properties()[0]
    widthtop = get_properties()[1]
    widthbottom = get_properties()[2]
    topthickness = get_properties()[3]
    bottomthickness = get_properties()[4]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]

    I = geometric_properties()[0]
    y_bar = geometric_properties()[1]
    Q = geometric_properties()[2]
    y_top = geometric_properties()[3]
    
    # ##########################
    moment_fail = []
    for i in range(L):
        if BMD[i] >0: #positive moment the max tension will be at the bottom
            moment_fail.append(sigT * I[i]/y_bar[i]) #M = sigma*I / y
        elif BMD[i] <0:
            moment_fail.append(-1 * sigT * I[i] /y_bar[i]) 
        else:
            moment_fail.append(0)
        
    return moment_fail


    
def M_failMatC(sigC, BMD):
    global L
    """
    calculate moment the moment the causes tension failure at each x
     returns array of the moment required to break bridge at each location
    sectional_properties array with I, Q, yBar, Ybot, Ytop

    sigC = sigma compression
    BMD 1-d array of the bending moments at each mm of the bridge

    """
    # change with right indices
    I = geometric_properties()[0]
    y_bar = geometric_properties()[1]
    Q = geometric_properties()[2]
    y_top = geometric_properties()[3]
    # ##########################
    
    moment_fail = []
    for i in range(L):
        if BMD[i] >0: #positive moment the max tension will be at the bottom
            moment_fail.append(sigC * I[i]/y_top[i]) #M = sigma*I / y
        elif BMD[i] <0:
            moment_fail.append(-1 * sigC * I[i] /y_bar[i]) 
        else:
            moment_fail.append(0)
    return moment_fail

def MFailBuck(BMD, t):
    global E, MU

    """
    calculated moment needed to break bridge due to buckling at every location along the bridge

     sectional_properties array with I, Q, yBar, Ybot, Ytop
        each indivisual properties should be an array with the value at each location along the bridge
    
        in sectional properties have an array with a value that describes the shape of the bridge so the corresponding correct cases can be used
    global E and mu

    3 cases
        case 1 when two sides restained constant force k = 4
        case 2 when one side restained constant froce k =0.425
        case 3 carries compressive stresses (like the webs) k=6

        sigma = (4*pi**2*E)/(12(1-MU**2))*(t/b)**2
        M = sigma * I /Y
        b= width
    """
    #needed to change depending on sectional properties
    I = geometric_properties()[0]
    y_bar = geometric_properties()[1]
    Q = geometric_properties()[2]
    y_top = geometric_properties()[3]

    case = 1
    sigma = 0
    moments = []
    shape = get_properties()[7]
    for i in range(0, L):
        if shape[i] == 0: #the demo bridge
            sigma1 = case1(i, t, 80)
            sigma2 = case2(i,t,10)
            sigma3 = case3(i, t, 32.02)
            # print(i)
            sigma = min(sigma1, sigma2, sigma3)
            # take min of potential sigmas calculated like if a memeber uses both case1 and 2
                   
        moments.append(sigma*I[i] / y_bar[i]) #idk what Y to reference   
        # moments.append(1)
    return(moments)    


def case1(i,t, b):
    
    n1 = (4*math.pi**2 * E)
    n2 = (12*(1-MU**2))
    n3 = (t/b)**2
    return n1/n2 * n3

def case2(i, t, b):
    
    sigma = (0.425*math.pi**2 * E)/(12*(1-MU**2))*(t/b)**2
    return sigma

def case3(i, t, b):
   
    sigma = (6*math.pi**2 * E)/(12*(1-MU**2))*(t/b)**2
    return sigma

def testFail():
    """
    increment p until the bridge breaks
    if you want to change position through in another loop for position of P

    for every iteration of p
        compare the moment at each time and see if its less then the max
    """
    
    broken = False
    p = 100
    while not broken:
        printSFD(550 + 30/2, p)
        printSFD(1250  + 30/2, p)
        BMD = buildBMD(550, p) #function to get the BMD using p in the 
        BMD = buildBMD(1250, p)
        print("die")
        printBMD(550, p)
        printBMD(1250, p)
        m_buckle = MFailBuck(BMD, 1.27)
        m_tension = M_failMatT(30, BMD)
        m_compression = M_failMatC(6, BMD)

        for i in range(L):
            if m_buckle[i] < BMD[i]:
               broken = True
               print("BUCKLE FAIL at " + str(i) + " at load" + str(p))
            elif m_tension[i] < BMD[i]:
               broken = True
               print("Tension FAIL at " + str(i)+ " at load" + str(p))
            elif m_compression[i] < BMD[i]:
               broken = True
               print("COMPRESSION FAIL at " + str(i)+ " at load" + str(p))
        print("p is " + str(p))
        p+=1
        
               
def deflection(BMD):
    
    """
    assuming a max positive moment and the bridge is deflecting down

    for midspsan delfection:
        use similar triangles with the tangent drawn at the support on the left
        delta_C_a = tangential deviation on the right support to tangent at a

    curvature max becomes the height of the similar trianggles
     delta_c_a / L = (delta_mid_a + deflection)/L/2
    """            
    I = geometric_properties()[0]
    y_bar = geometric_properties()[1]
    Q = geometric_properties()[2]
    y_top = geometric_properties()[3]


    curvatures = []
    for i in range(BMD):
        curvatures.append(BMD[i]/(E*I[i]))
    max_curv = max(curvatures)
    curvature_midpoint_scale = 1 - (curvatures.index(max)/L-1/2) #not sure how this will work for if it is curved before midpoint
    delta_c_a = (1/2*(L/2)*max_curv * (L/2 + (1/3)*(L/2)))+ (1/2*L/2 * max_curv * (L/2*(2/3))) #sum(Ai*di) where di is distance from c to centroid of area
    delta_m_a = (1/2 *(L/2) *max_curv*curvature_midpoint_scale) * (1/3*(L/2))
    deflection = delta_c_a/L * (L/2) -delta_m_a
    return deflection

if __name__ == "__main__":
    testFail()