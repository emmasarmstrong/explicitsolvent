import numpy as np
from math import *

#calculating distance between centres of two spheres
def separation(X1,X2):
    return np.linalg.norm( X1 - X2 )

#function to calculate the range of values for theta and phi for which a new molecule could be placed
def anglelimit(d,r1,r2,R):
    #d distance between two molecules
    #r1 radius of the selected molecule
    #r2 radius of the reference molecule
    #R set distance away the new molecule to be placed
    s = (d**2+(R+r1)**2-(R+r2)**2) / (2.0*d*(R+r1))
    return np.arccos(s)

#conversion from spherical coordinates to cartesian
def cartesian(r,theta,phi):
    x=r*cos(phi)*sin(theta)
    y=r*sin(phi)*sin(theta)
    z=r*cos(theta)
    return x,y,z

#conversion from cartesian to spherical polars
def spherical(x,y,z):
    r=(x**2+y**2+z**2)**(1/2)
    theta=np.arccos(z/r)
    if x<1e-8:
        phi=0
    else:
        phi=np.arctan2(y,x)
    return r,theta,phi

def distance_mat(x):
    n = x.shape[0]
    return [ [separation(x[i], x[j]) for j in range(i, n)] for i in range(n) ]

def periodic_shift(value, lower=0, upper=2*np.pi):
    period = upper - lower
    if value < lower:
        new_value = value + period
    elif value > upper:
        new_value = value - period
    else:
        new_value = value
    return new_value

def rotate_bounds(axis, lower):
    dummyr, theta_ref, phi_ref = spherical(axis[0],axis[1],axis[2])
    theta_lower = periodic_shift(lower - theta_ref)
    theta_upper = periodic_shift(2*np.pi - lower - theta_ref)
    phi_lower = periodic_shift(lower - phi_ref, upper=np.pi)
    phi_upper = periodic_shift(np.pi - lower - phi_ref, upper=np.pi)
    return theta_lower, phi_lower, theta_upper, phi_upper 

def rotate(coords, axis):
    unit_axis = axis / np.linalg.norm(axis)

    bmat_x = np.zeros([3, 3])
    bmat_x[0, 0] = 1.0
    bmat_x[1, 1] = bmat_x[2, 2] = unit_axis[0]
    bmat_x[2, 1]  = np.sqrt(1.0 - unit_axis[0]**2)
    bmat_x[1, 2] = -bmat_x[2, 1]

    bmat_z = np.zeros([3, 3])
    bmat_z[2, 2] = 1.0
    bmat_z[0, 0] = bmat_z[1, 1] = unit_axis[2]
    bmat_z[1, 0] = np.sqrt(1.0 - unit_axis[2]**2)
    bmat_z[0, 1] = -bmat_z[1, 0]

    bmat = np.dot(bmat_x, bmat_z)
    return np.dot(bmat, np.transpose(coords))

def angleselection(x,r,k,molecules):
#radius of the selected molecule                                                                                                            
    r1=r[k]
    #number of reference molecules                                                                                                              
    m=len(molecules)
    if m==0:
        print("M = 0")
        theta=np.random.uniform(0,2*np.pi)
        phi=np.random.uniform(0,np.pi)

    else:
        mol0 = molecules[0]
        r2 = r[mol0]
        d12=separation(x[mol0],x[k])
        #lower limit for angle                                                                                                                  
        lower = anglelimit(d12,r1,r2,s)
        #upper limit for angle                                                                                                                
        lower_theta, lower_phi, upper_theta, upper_phi = rotate_bounds(x[mol0]-x[k], lower) 

        if m==2:
            print("M=2")
            mol1 = molecules[1]
            r3 = r[mol1]
            d13 = separation(x[mol1] , x[k])
            #lower limit                                                                                                                        
            lower2 = anglelimit(d13,r1,r3,s)
            #theta_ref polar angle between reference particles                                                                                  
            #phi_ref azimuth angle between reference particles                                                                                  
            theta_lower2, phi_lower2, theta_upper2, phi_upper2 = rotate_bounds(x[mol1]-x[k], lower2)

            lower_theta = max(lower_theta, theta_lower2)
            lower_phi = max(lower_phi, phi_lower2)
            upper_theta = min(upper_theta, theta_upper2)
            upper_phi = min(upper_phi, phi_upper2)

        theta=np.random.uniform(lower_theta, upper_theta)
        phi=np.random.uniform(lower_phi, upper_phi)

    #convert from spherical polars centred on particle to cartesian centred on the orgin                                                        
    coords = np.array( cartesian(r1+s, theta, phi) )
    coords += x[k, :]
#    coords = rotate(coords, x[k])

    if check_coords(coords, x, r, k, molecules):
        print(theta*180/np.pi, phi*180/np.pi)
        print(lower*180/np.pi, upper_theta*180/np.pi, upper_phi*180/np.pi)

    return coords 

def check_coords(coords, x, r, k, molecules):
    distances = [separation(x[m], coords) for m in molecules]
    distances.append(separation(x[k], coords))
    radii = [r[k], s]
    for m in molecules:
        radii.append(r[m])
    min_r = np.min(radii)
    problem = False
    for d in distances:
        problem = d < 1.95*min_r or problem
    if problem:
        print("Coords", x, coords)
        print("Indices", k, molecules)
        print("Distances", distances, radii)
        print()
    return problem

def placemolecule(x,r,s,l):
    #x list of coordinates of currently placed molecules
    #r list of molecules respective radii
    #s solvent molecule radius
    #l finite box size
    n=len(x)

    #determining molecule for which next placed molcule will be based on
    found=False
    random_flag = False
    indices = [i for i in range(n)]
    ctr = 0
    while not found and ctr < 200:
        distances = []
        molecules = []
        if len(indices) > 0:
            #randomly selecting already placed molecule
            k = np.random.randint( len(indices) )
            k = indices[k]
            #print("Index:", k)
            #calculating distances and storing those for possible reference molecules
    
            for i in [j for j in range(n) if j!=k]:
                d=separation(x[i],x[k])
                max_d = r[i] + r[k] + 2*s
                if d <= max_d:
                    distances.append(d)
                    molecules.append(i)
            indices = [i for i in indices if i != k]
            #print(indices)
        else:
            random_flag = True
            x_random=2*l * np.random.random_sample(3) - l                                                                                     
            #print("Rando:", x_random)
            for i in range(n):                                                                                                                
                d=separation(x[i],x_random)                                                                                                   
                max_d = r[i] + 3*s                                                                                                           
                if d <= max_d:                                                                                                                
                    distances.append(d)                                                                                                       
                    molecules.append(i)                                           

        found = len(distances) < 3
        ctr+=1

    if found:
        if random_flag:
            x = np.append(x, [x_random], axis=0)
        else:
            coords = angleselection(x, r, k, molecules)
            x = np.append(x, [coords], axis=0)
         
        r = np.append(r, s)
    else:
        print("Ctr over 200")
    
    return x, r

if __name__ == '__main__':
    x=np.array( [[-3,0,0], [3,0,0]] )
    r=np.array( [2,2] )
    s=1
    N=100
    l=10

    np.random.seed(42)
    for n in range(N):
        print("Placing molecule ", n)
        x,r = placemolecule(x,r,s,l)

    print(N+2)
    print("Title")
    for c in x:
        print('He', c[0], c[1], c[2])
