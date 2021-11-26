import math
L = 1000
E = 4000
MU = 0.2

def M_failMatT(sectional_properties, sigT, BMD):
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
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]

    # ##########################
    moment_fail = []
    for i in range(BMD):
        if i >0: #positive moment the max tension will be at the bottom
            moment_fail.append(sigT * I[i]/yBot[i]) #M = sigma*I / y
        elif i <0:
            moment_fail.append(-1 * sigT * I[i] /yTop[i]) 
        else:
            moment_fail.append(0)


    
def M_failMatC(sectional_properties, sigC, BMD):
    global L
    """
    calculate moment the moment the causes tension failure at each x
     returns array of the moment required to break bridge at each location
    sectional_properties array with I, Q, yBar, Ybot, Ytop

    sigC = sigma compression
    BMD 1-d array of the bending moments at each mm of the bridge

    """
    # change with right indices
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    # ##########################
    moment_fail = []
    for i in range(BMD):
        if i >0: #positive moment the max tension will be at the bottom
            moment_fail.append(sigC * I[i]/yBot[i]) #M = sigma*I / y
        elif i <0:
            moment_fail.append(-1 * sigC * I[i] /yTop[i]) 
        else:
            moment_fail.append(0)

def MFailBuck(sectional_properties,  BMD):
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

        sigma = (4*pi^2*E)/(12(1-MU^2))*(t/b)^2
        M = sigma * I /Y
        b= width
    """
    #needed to change depending on sectional properties
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    t = sectional_properties[3]
    b= sectional_properties[4]
    shape = sectional_properties[5]

    case = 1
    sigma = 0
    moments = []
    for i in range(L):
        if shape(i) == 0:
            sigma = case1(sectional_properties, i)
            # take min of potential sigmas calculated like if a memeber uses both case1 and 2
            pass
        elif shape(i) == 1:
            #do
            pass
        
        moments.append(sigma*I(i) / yTop(i)) #idk what Y to reference        


def case1(sectional_properties, i):
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    t = sectional_properties[3]
    b= sectional_properties[4]
    sigma = (4*math.pi^2 * E)/(12-(1-MU^2))*(t(i)/b(i))^2
    return sigma
def case2(sectional_properties, i):
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    t = sectional_properties[3]
    b= sectional_properties[4]
    sigma = (0.425*math.pi^2 * E)/(12-(1-MU^2))*(t(i)/b(i))^2
    return sigma
def case3(sectional_properties, i):
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    t = sectional_properties[3]
    b= sectional_properties[4]
    sigma = (6*math.pi^2 * E)/(12-(1-MU^2))*(t(i)/b(i))^2
    return sigma

def testFail():
    """
    increment p until the bridge breaks
    if you want to change position through in another loop for position of P

    for every iteration of p
        compare the moment at each time and see if its less then the max
    """
    sectional_properties = sectional_properties()
    broken = False
    p = 1
    while not broken:
        BMD = calcBMD(p) #function to get the BMD using p in the BMD
        m_buckle = MFailBuck(sectional_properties, BMD)
        m_tension = M_failMatT(sectional_properties, BMD)
        m_compression = M_failMatC(sectional_properties, BMD)
        for i in range(L):
            if m_buckle(i) < BMD(i):
               broken = True
               print("BUCKLE FAIL at " + i + " at load" + p)
            elif m_tension(i) < BMD(i):
               broken = True
               print("Tension FAIL at " + i+ " at load" + p)
            elif m_compression(i) < BMD(i):
               broken = True
               print("COMPRESSION FAIL at " + i+ " at load" + p)
               
def deflection(BMD, sectional_properties):
    
    """
    assuming a max positive moment and the bridge is deflecting down

    for midspsan delfection:
        use similar triangles with the tangent drawn at the support on the left
        delta_C_a = tangential deviation on the right support to tangent at a

    curvature max becomes the height of the similar trianggles
     delta_c_a / L = (delta_mid_a + deflection)/L/2
    """            
    I = sectional_properties[0]
    yBot = sectional_properties[1]
    yTop = sectional_properties[2]
    t = sectional_properties[3]
    b= sectional_properties[4]

    curvatures = []
    for i in range(BMD):
        curvatures.append(BMD(i)/(E*I(i)))
    max_curv = max(curvatures)
    curvature_midpoint_scale = 1 - (curvatures.index(max)/L-1/2) #not sure how this will work for if it is curved before midpoint
    delta_c_a = (1/2*(L/2)*max_curv * (L/2 + (1/3)*(L/2)))+ (1/2*L/2 * max_curv * (L/2*(2/3))) #sum(Ai*di) where di is distance from c to centroid of area
    delta_m_a = (1/2 *(L/2) *max_curv*curvature_midpoint_scale) * (1/3*(L/2))
    deflection = delta_c_a/L * (L/2) -delta_m_a
    return deflection
