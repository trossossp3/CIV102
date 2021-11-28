import math

## Dumbass CIV102 Project Code ##
import scipy.integrate as it
import matplotlib.pyplot as plt
import numpy as np

n = 1000
L = 1290+1 #add one so u can do arr[1250]
x = 1
E = 4000
MU = 0.2
supportAlocation = 30/2

sfd = [[0.0, 0.0]]*L
sfd = np.array(sfd)
support_x = 15
support_x_right = 1075
point_loads = []
bmd = np.array([0.0]*L)

def buildSFD(xP, P):
    global sfd
    global point_loads
    temp = [P, xP]
    added = False
    for i in range(len(point_loads)):
        if point_loads[i][1] == temp[1]:
            point_loads[i][0] = temp[0]
            added = True
    if not added:
        point_loads.append(temp)
    sum1  = 0
    sum2 = 0
    # points loads will be a 2d array with [P, xp] as values

    for i in range(len(point_loads)):
        sfd[point_loads[i][1]] = point_loads[i][0]
        sum1 += point_loads[i][0] * (point_loads[i][1] - support_x)
        sum2 += point_loads[i][0]
    # by = (xP-support_x) * P / (support_x_right - support_x)
    by = sum1/(support_x_right-support_x)
    ay = sum2 - by

    sfd[support_x] = [ay, 1]
    sfd[support_x_right] = [by, 1]
    
    sum = 0.0
    for i in range(0, L):
        # print(sfd[i])
        # print(sfd[i][0])
        if (sfd[i][1] == 0):
            pass
        elif (sfd[i][1] ==1):
            sum += sfd[i][0]
        else:
            sum-= sfd[i][0]
        sfd[i][0] = sum

    arr = []
    for i in range(len(sfd)):
        arr.append(sfd[i][0])
    return arr
 


def printSFD():
    plt.title("Bridge SFD") 
    plt.xlabel("x (mm)") 
    plt.ylabel("shear force (kN)") 
    arr = []
    for i in range(len(sfd)):
        arr.append(sfd[i][0])  
    zeroliney = [0, 0]
    zerolinex = [0, L]
    plt.plot(zerolinex, zeroliney)
    plt.plot(arr)
    plt.show()


def buildBMD():
    global bmd
    arr = []
    for i in range(len(sfd)):
        arr.append(sfd[i][0])   

    bmd=it.cumtrapz(arr)


def printBMD():
    plt.title("Bridge BMD") 
    plt.xlabel("x (mm)") 
    plt.ylabel("bending moment (kNm)") 
    zeroliney = [0, 0]
    zerolinex = [0, L]
    plt.plot(zerolinex, zeroliney)
    # plt.plot(*zip(*sorted(buildBMD(xP, P))))
    plt.gca().invert_yaxis()
    plt.plot(bmd)
    plt.show()    



def ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    if(shape == 0):#shape = 0 during test bridge    
            area1 = widthtop * topthickness #top part of bridge
            area2 = 10 * bottomthickness #tabs 10 is tab width
            area3 = leftthickness * (height - topthickness)
            area4 = (widthbottom - 2*leftthickness) * bottomthickness
            d1 = (height-topthickness/2)
            d2 = (height - topthickness - topthickness/2)
            d3 = (height-2*topthickness)/2
            d4 = bottomthickness/2
            y_bar = area1*d1 + 2*(area2*d2) + 2*(area3*d3) + area4*(d4)     
            y_bar = y_bar / (area1+2*area2+2*area3+area4)
            return y_bar

    elif(shape ==1):
        area1 = widthtop * topthickness #top part of bridge
        area2 = 10 * bottomthickness #tabs 10 is tab width
        area3 = leftthickness * (height - topthickness - bottomthickness)
        area4 = (widthbottom - 2*leftthickness) * bottomthickness
        d1 = (height-topthickness/2)
        d2 = (height - topthickness - topthickness/2)
        d3 = (height-2*topthickness)/2
        d4 = bottomthickness/2
        y_bar = area1*d1 + 2*(area2*d2) + 2*(area3*d3) + area4*(d4)     
        y_bar = y_bar / (area1+2*area2+2*area3+area4)
        return y_bar

def I0(b, h):
    I0 = (b * h * h * h)/12
    return I0

def second_moment_of_area(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    centroid = ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape)
    if shape == 0:
        area1 = widthtop * topthickness #top part of bridge
        area2 = 10 * topthickness #tabs 10 is tab width
        area3 = leftthickness * (height - 2*topthickness)
        area4 = (widthbottom - 2*leftthickness) * bottomthickness
        d1L = (height-topthickness/2)
        d2L = (height - topthickness - topthickness/2)
        d3L = (height-2*topthickness)/2
        d4L = bottomthickness/2
        d1C = d1L - centroid
        d2C = d2L - centroid
        d3C = d3L - centroid
        d4C = d4L - centroid

        term1 = I0(widthtop, topthickness) + (area1 * (d1C)**2) #top
        term2 = I0(10, topthickness) + area2 * (d2C)**2 #tabs
        term3 = I0(leftthickness, height - 2*topthickness) + (area3 * (d3C)**2)
        term4 = I0(widthbottom - 2*leftthickness, bottomthickness) + (area4 * (d4C)**2)
        second_moment_of_area = term1 + 2* term2 + 2*term3 + term4
    return second_moment_of_area

""""
def first_moment_of_area(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    
    centroid = ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape)
    
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
"""

def first_moment_of_area(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape):
    
    centroid = ybar(height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape)
    
    if shape == 0:
        area1 = widthbottom*bottomthickness
        area2 = (centroid - bottomthickness) * rightthickness
        area3 = (centroid - bottomthickness) * leftthickness
        d1 = centroid - (bottomthickness/2)
        d2 = (centroid - bottomthickness)/2
        d3 = (centroid - bottomthickness)/2

        Q = area1*d1 + area2*d2 + area3*d3

    return Q

def get_properties():
    """
    user input to say what these properties are at the locations
    """
    height = [75+35]*L
    widthTop = [100]*L
    widthBottom = [80]*L
    topThickness = [1.27*2]*L
    #increases thickness
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
    height, widthtop, widthbottom, topthickness, bottomthickness, rightthickness, leftthickness, shape= get_properties()
    y_bar = []
    I = []
    Q = []
    yTop = []

    for i in range(L):
        y_bar.append(ybar(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i], shape[i]))
        I.append(second_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i], shape[i]))
        Q.append(first_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i], shape[i]))
        yTop.append(height[i]-y_bar[i])
    

    return I, y_bar, Q, yTop

def V_fail(Tau):
  
    I = geometric_properties()[0]
    Q = geometric_properties()[2]
    centroid = geometric_properties()[1]
    bottomthickness = get_properties()[4]
    height = get_properties()[0]
    topthickness = get_properties()[3]
    widthbottom = get_properties()[2]
    widthtop = get_properties()[1]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]
    V_fail_values = []

    for i in range(L):
        if 0 <= centroid[i] < bottomthickness[i]:
            centroidlocation = "bottom"
        elif bottomthickness[i] <= centroid[i] < height[i] - topthickness[i]:
            centroidlocation = "middle"
        elif height[i] - topthickness[i] <= centroid[i] < height[i]:
            centroidlocation = "top"    

        if centroidlocation == "bottom":
            b = widthbottom[i]
        elif centroidlocation == "middle":
            b = rightthickness[i] + leftthickness[i]
        elif centroidlocation == "top":
            b = widthtop[i]

        V_fail = (I[i]*Tau*b)/Q[i]

        V_fail_values.append(V_fail)

    return V_fail_values

def V_buck(max_diaphragm_distance):

    height = get_properties()[0]
    widthtop = get_properties()[1]
    widthbottom = get_properties()[2]
    topthickness = get_properties()[3]
    bottomthickness = get_properties()[4]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]
    I = geometric_properties()[0]
    Q = geometric_properties()[2]
    centroid = geometric_properties()[1]
    V_buck_values = []

    for i in range(L):

        Tcrit = ((5*((math.pi)**2)*E)/(12*(1-(MU**2)))) * (((rightthickness[i]/max_diaphragm_distance)**2) + ((leftthickness[i]/(height[i]- topthickness[i]))**2))

        if 0 <= centroid[i] < bottomthickness[i]:
            centroidlocation = "bottom"
        elif bottomthickness[i] <= centroid[i] < height[i] - topthickness[i]:
            centroidlocation = "middle"
        elif height[i] - topthickness[i] <= centroid[i] < height[i]:
            centroidlocation = "top"    

        if centroidlocation == "bottom":
            b = widthbottom[i]
        elif centroidlocation == "middle":
            b = rightthickness[i] + leftthickness[i]
        elif centroidlocation == "top":
            b = widthtop[i]

        V_buck = Tcrit*I[i]*b/Q[i]

        V_buck_values.append(V_buck)

    return V_buck_values        

def M_failMatT(sigT):
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
    for i in range(L-1):
        if bmd[i] > 0: #positive moment the max tension will be at the bottom
            moment_fail.append((sigT * I[i])/y_bar[i]) #M = sigma*I / y
        elif bmd[i] <0:
            moment_fail.append(-1 * sigT * I[i] /y_top[i]) 
        else:
            moment_fail.append(0)
        
    return moment_fail


    
def M_failMatC(sigC):
    global L
    """
    calculate moment the moment the causes compression failure at each x
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
    for i in range(L-1):
        if bmd[i] >0: #positive moment the max tension will be at the bottom
            moment_fail.append(sigC * I[i]/y_bar[i]) #M = sigma*I / y
        elif bmd[i] <0:
            moment_fail.append(-1 * sigC * I[i] /y_top[i]) 
        else:
            moment_fail.append(0)
    return moment_fail

def MFailBuck(t):
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
    height = get_properties()[0]
    widthtop = get_properties()[1]
    widthbottom = get_properties()[2]
    topthickness = get_properties()[3]
    bottomthickness = get_properties()[4]
    rightthickness = get_properties()[5]
    leftthickness = get_properties()[6]
    case = 1
    sigma = 0
    moments = []
    shape = get_properties()[7]
    for i in range(0, L):
        if shape[i] == 0: #the demo bridge
            sigma1 = case1(i, topthickness[i], 80)
            
            # sigma1 = case1(i, t, 80)
            sigma2 = case2(i, t,10)
            sigma3 = case3(i, t, 32.02)
            # print(i)
            sigma = min(sigma1, sigma2, sigma3)
            # take min of potential sigmas calculated like if a memeber uses both case1 and 2
            n1 =  sigma*I[i] / y_top[i]
            n2 = sigma*I[i] / y_bar[i] 
        moments.append(min(n1, n2)) #idk what Y to reference   
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
    # global sfd
    broken1 = False
    broken2 = False
    broken3 = False
    broken4 = False
    broken5 = False
    p = 550
    # while not broken1 or not broken2 or not broken3:
    while not broken1 and not broken2 and not broken3 and not broken4 and not broken5:
        buildSFD(550, p)
        # printSFD()

        arr = buildSFD(1250,p)
        # print(arr[676])
        with open('your_file.txt', 'w') as f:
            for item in arr:
                f.write("%s\n" % item)
        printSFD()
        buildBMD()
        printBMD()
        m_buckle = MFailBuck(1.27)
        m_tension = M_failMatT(30)  
        m_compression = M_failMatC(6)
        v_fails = V_fail(4)
        v_bucks = V_buck(520)
        for i in range(L-1):
            if abs(m_buckle[i]) < abs(bmd[i]):
               broken1 = True
               print("BUCKLE FAIL at " + str(i) + " at load" + str(p))
               print(bmd[i])
            if abs(m_tension[i]) < abs(bmd[i]):
               broken2 = True
               print("Tension FAIL at " + str(i)+ " at load" + str(p))
               print(bmd[i])
            if abs(m_compression[i]) < abs(bmd[i]):
               broken3 = True
               print("COMPRESSION FAIL at " + str(i)+ " at load" + str(p))
               print(bmd[i])
            if abs(v_fails[i]) < abs(arr[i]):
                # print(arr[i])
                broken4 = True
                print("sheer failure at " + str(i) + " at load " + str(p) + " sheer fail value = " + str(arr[i]))
            if abs(v_bucks[i]) < abs(arr[i]):
                print(arr[i])
                broken5 = True
                print("sheer buckle failure at " + str(i) + " at load " + str(p) + " sheer fail value = " + str(arr[i]))
                pass
        p+=1
def testFailTrain():
     # global sfd
    broken1 = False
    broken2 = False
    broken3 = False
    broken4 = False
    broken5 = False
    p = 42
    # while not broken1 or not broken2 or not broken3:
    while not broken1 and not broken2 and not broken3 and not broken4 and not broken5:
        #for values at the right
        # buildSFD(372, p)
        # buildSFD(548,p)
        # buildSFD(712, p)
        # buildSFD(888,p)
        # buildSFD(1052, p)
        # arr =buildSFD(1228,p)


        # for values at center
        buildSFD(117, p)
        buildSFD(293,p)
        buildSFD(457, p)
        buildSFD(633,p)
        buildSFD(797, p)
        arr =buildSFD(973,p)
        print(arr[676])
        with open('your_file.txt', 'w') as f:
            for item in arr:
                f.write("%s\n" % item)
        # printSFD()
        buildBMD()
        # printBMD()
        m_buckle = MFailBuck(1.27)
        m_tension = M_failMatT(30)  
        m_compression = M_failMatC(6)
        v_fails = V_fail(4)
        v_bucks = V_buck(520)
        for i in range(L-1):
            if abs(m_buckle[i]) < abs(bmd[i]):
               broken1 = True
               print("BUCKLE FAIL at " + str(i) + " at load" + str(p))
               print(bmd[i])
            if abs(m_tension[i]) < abs(bmd[i]):
               broken2 = True
               print("Tension FAIL at " + str(i)+ " at load" + str(p))
               print(bmd[i])
            if abs(m_compression[i]) < abs(bmd[i]):
               broken3 = True
               print("COMPRESSION FAIL at " + str(i)+ " at load" + str(p))
               print(bmd[i])
            if abs(v_fails[i]) < abs(arr[i]):
                print(arr[i])
                broken4 = True
                print("sheer failure at " + str(i) + " at load " + str(p) + " sheer fail value = " + str(arr[i]))
            if abs(v_bucks[i]) < abs(arr[i]):
                print(arr[i])
                broken5 = True
                print("sheer buckle failure at " + str(i) + " at load " + str(p) + " sheer fail value = " + str(arr[i]))
                
        p+=1

def deflection(bmd):
    
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
    for i in range(bmd):
        curvatures.append(bmd[i]/(E*I[i]))
    max_curv = max(curvatures)
    curvature_midpoint_scale = 1 - (curvatures.index(max)/L-1/2) #not sure how this will work for if it is curved before midpoint
    delta_c_a = (1/2*(L/2)*max_curv * (L/2 + (1/3)*(L/2)))+ (1/2*L/2 * max_curv * (L/2*(2/3))) #sum(Ai*di) where di is distance from c to centroid of area
    delta_m_a = (1/2 *(L/2) *max_curv*curvature_midpoint_scale) * (1/3*(L/2))
    deflection = delta_c_a/L * (L/2) -delta_m_a
    return deflection





if __name__ == "__main__":
    # buildSFD(100,50)
    # printSFD()
    # buildBMD()
    # printBMD()
    # init()
    testFail()
    # testFailTrain()

