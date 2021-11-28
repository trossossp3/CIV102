
def V_fail(Tau):
    
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
    V_fail_values = []

    for i in np.arange(0, L, x):

        # I = second_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        # Q = first_moment_of_area(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        # centroid = ybar(height[i], widthtop[i], widthbottom[i], topthickness[i], bottomthickness[i], rightthickness[i], leftthickness[i])
        
        if 0.0 <= y_bar[i] < bottomthickness[i]:
            centroidlocation = "bottom"
        elif bottomthickness[i] <= y_bar[i] < height[i] - topthickness[i]:
            centroidlocation = "middle"
        elif height[i] - topthickness[i] <= y_bar[i] < height[i]:
            centroidlocation = "top"    

        if centroidlocation == "bottom":
            b = widthbottom[i]
        elif centroidlocation == "middle":
            b = rightthickness[i] + leftthickness[i]
        elif centroidlocation == "top":
            b = widthtop[i]
        
        V_fail = I[i]*Tau*b/Q[i]

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

