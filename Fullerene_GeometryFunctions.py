import numpy as np
from numpy import array, pi, cos, sin, sqrt, linspace
import numpy.linalg as la

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

NA = np.newaxis

from matplotlib import cm, colors

###############################################################################

def split_norm(X, axis=-1):
    R = np.sqrt(np.sum(X*X,axis=axis));
    D = X/R[...,NA]
    return R,D

# X: n x 3, neighbours: n x d
# X[neighbours]: (n x d) x 3
def edge_displacements(X,neighbours):
    Xab = X[neighbours]-X[:,NA]       # Displacement vectors Xv-Xu (n x d x 3)
    return split_norm(Xab)

  
##################################################################################

def corner_cos_angles(Duv): #Duv is a unit vector. 
    ab_hat, ac_hat = Duv, np.roll(Duv,shift=-1,axis=1); #a = np.array(1,2,3), np.roll(a,shift=-1) => (2,3,1)
    return np.sum(ab_hat*ac_hat,axis=-1)

############################################################

def bcd_corner_cos(X,neighbours,on_face):  
    bcd      = neighbours
    bcdplus  = on_face     
    # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    # prev_on_face(a,b)=b- | prev_on_face(a,c)=c- | prev_on_face(a,d)=d-
    
    ##     c-  b+
    ##     |   |  
    ## c+__c   b__b-
    ##      \ /   
    ##       a
    ##       |
    ##       d
    ##      / \
    ##     d- d+
    
    _, ba_hat = split_norm(X[:,NA]    - X[bcd])
    _, bp_hat = split_norm(X[bcdplus] - X[bcd])

    return np.sum(ba_hat*bp_hat,axis=-1)
############################################################

def dihedral_cos_angles(Duv,Ruv):

    Xuv = Duv*Ruv[...,NA]
    # Create the vectors
    ab, ac, ad = Xuv, np.roll(Xuv,shift=-1,axis=1), np.roll(Xuv,shift=-2,axis=1); 
    bc = ac - ab; 
    cd = ad - ac; 
    #Normalize vectors
    cd_hat = cd/ la.norm(cd,axis=-1)[...,NA]
    bc_hat = bc/ la.norm(bc,axis=-1)[...,NA]

    ba_hat = -ab/ la.norm(ab,axis=-1)[...,NA]
    cb_hat = -bc_hat

    # Find angles spaned by vectors
    sin_b = np.sqrt(1 - np.sum(ba_hat*bc_hat,axis=-1)**2)[...,NA]
    sin_c = np.sqrt(1 - np.sum(cb_hat*cd_hat,axis=-1)**2)[...,NA]
                        ###########  cos  ###########
        
    # Plane unit vectors
    nabc = np.cross(ba_hat,bc_hat,axis=-1)/sin_b # N x d x 3
    nbcd = np.cross(cb_hat,cd_hat,axis=-1)/sin_c # N x d x 3

    cosbeta = np.sum(nabc*nbcd,axis=-1) # N x d

    return cosbeta
############################################################
def energy_harmonic(parameter,const,force_const): 
    error = (parameter-const); 
    return np.sum((force_const/2)*(error*error))
############################################################

def energy(X,neighbours,k0,f0):
    Ruv,Duv   = edge_displacements(X,neighbours);
    cos_abc   = corner_cos_angles(Duv)
    dih_cos   = dihedral_cos_angles(Duv,Ruv)

    R0,ang0,dih0,_,_,_,_,_ = k0
    fR,fang,fdih,_,_,_,_,_ = f0
    
    ## the energy for bond length has a factor half 
    ## since the energylength is calculated twice, this does not happen to the angle or dihedral
    
    return  energy_harmonic(Ruv, R0,fR)/2 + energy_harmonic(cos_abc,ang0,fang)  + energy_harmonic(dih_cos, dih0,fdih)



################################################################


def edge_energy_gradient_from_points(X,neighbours,R0,force_const):
    R_ab, ab_hat = edge_displacements(X,neighbours)
    return -2*np.sum((force_const[...,NA]/2)*(R_ab-R0)[...,NA]*ab_hat,axis=1)
############################################################

def corner_cos_gradient(Ruv,Duv,ang0,force_const):
    R_ab = Ruv[...,NA];
    R_ac = np.roll(Ruv,shift=-1,axis=1)[...,NA];
    ab_hat, ac_hat = Duv, np.roll(Duv,shift=-1,axis=1)
    cos_angle = corner_cos_angles(Duv)[...,NA]
    grad = cos_angle*(ab_hat/R_ab + ac_hat/R_ac) - ab_hat/R_ac - ac_hat/R_ab
    return 2*np.sum((force_const[...,NA]/2)*(cos_angle - ang0[...,NA])*grad,axis=1); 

############################################################

def bcd_corner_cos_gradient(X, neighbours, on_face, ang0,force_const):
#     bcd      = neighbours
    bcdplus  = on_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    
    R_ab,  ab_hat  = split_norm(X[:,NA] - X[neighbours])
    R_bbp, bbp_hat = split_norm(X[bcdplus] - X[neighbours])
    cos_angle  = np.sum(ab_hat*bbp_hat,axis=-1)
    grad = (bbp_hat - ab_hat*cos_angle[...,NA])/R_ab[...,NA]
    
    return 2*np.sum(((force_const[...,NA]/2)*(cos_angle[...,NA] - ang0[...,NA]))*grad, axis=1)

############################################################



def dihedral_cos_energy_gradient_from_points(X,neighbours,dih0,force_const):
    Ruv,Duv = edge_displacements(X,neighbours) # points from a - > b,c,d
    
    Rab = Ruv[...,NA]
    Xuv = Duv*Rab
    # Create the vectors
    Xab, Xac, Xad = Xuv, np.roll(Xuv,shift=-1,axis=1), np.roll(Xuv,shift=-2,axis=1); 
    Xbc = Xac - Xab; 
    Xcd = Xad - Xac; 
    
    Rcd = la.norm(Xcd,axis=-1)[...,NA]
    Rbc = la.norm(Xbc,axis=-1)[...,NA]
    #Normalize vectors
    cd_hat = Xcd/ Rcd
    bc_hat = Xbc/ Rbc
    
    cb_hat = -bc_hat
    
    ba_hat = -Xab/ Rab
    # ad_hat = ad/ la.norm(ad,axis=-1)[...,NA]
    # Find angles spaned by vectors
    cos_b = np.sum(ba_hat*bc_hat,axis=-1)[...,NA]
    cos_c = np.sum(cb_hat*cd_hat,axis=-1)[...,NA]
    
    
    sin_b = np.sqrt(1 - cos_b**2)
    sin_c = np.sqrt(1 - cos_c**2)

    # Plane unit vectors
    nabc = np.cross(ba_hat,bc_hat,axis=-1)/sin_b # N x d x 3
    nbcd = np.cross(cb_hat,cd_hat,axis=-1)/sin_c # N x d x 3

    #Angle between planes
    cos_beta = np.sum(nabc*nbcd,axis=-1)[...,NA]
    
    cot_b = cos_b/sin_b
    
    
    d_dih_a = (np.cross(bc_hat,nbcd)/(sin_b*Rab)) - ((ba_hat*cos_beta)/Rab) + ((cot_b*cos_beta)/(sin_b*Rab)) * (bc_hat - ba_hat*cos_b)


    return 2*np.sum((force_const[...,NA]/2)*(cos_beta-dih0[...,NA])*d_dih_a,axis=1)   

############################################################
def outer_dih_a_gradient(X, neighbours, next_face, prev_face, dih,force_const):
    bcd_p = next_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    bcd_m = prev_face
    
    Rba, ba_hat = split_norm(X[:,NA]   - X[neighbours]) # b --> a
    Rbp, bp_hat = split_norm(X[bcd_p]  - X[neighbours]) # b --> p
    Rbm, bm_hat = split_norm(X[bcd_m]  - X[neighbours]) # b --> m
    
    Rba = Rba[...,NA]
    Rbp = Rbp[...,NA]
    Rbm = Rbm[...,NA]
    
    Rab = Rba
    ab_hat = -ba_hat
    
    Xba, Xbp, Xbm = Rba*ba_hat, Rbp*bp_hat, Rbm*bm_hat
    
    Xam = Xbm - Xba
    Ram = np.linalg.norm(Xam,axis=-1)[...,NA]
    am_hat = Xam / Ram
    
    ma_hat, Rma = -am_hat, Ram
    
    Xmp = Xbp - Xbm
    Rmp = np.linalg.norm(Xmp,axis=-1)[...,NA]
    mp_hat = Xmp / Rmp
    
    cos_a = np.sum(ab_hat*am_hat,axis=-1)[...,NA]
    cos_m = np.sum(ma_hat*mp_hat,axis=-1)[...,NA]
    
    sin_a = np.sqrt(1 - cos_a**2)
    sin_m = np.sqrt(1 - cos_m**2)
    
    
    nbam_hat = np.cross(ab_hat,am_hat)/sin_a
    namp_hat = np.cross(ma_hat,mp_hat)/sin_m
    
    cos_beta = np.sum(nbam_hat*namp_hat,axis=-1)[...,NA]
    
    cot_a = cos_a / sin_a
    cot_m = cos_m / sin_m

    d_dih_a = (ab_hat*cos_beta / Rab) - (np.cross(am_hat,namp_hat)/(Rab*sin_a)) + (am_hat*cos_beta/Ram) -  (np.cross(namp_hat,ab_hat)/(Ram*sin_a)) + (cot_a*cos_beta/sin_a)*( (ab_hat*cos_a/Rab) - (am_hat/Rab) + (am_hat*cos_a/Ram) - (ab_hat/Ram) ) + (np.cross(mp_hat,nbam_hat)/(Rma*sin_m)) - (ma_hat*cos_beta/Rma) + (cot_m*cos_beta/sin_m)*( (mp_hat/Rma) - (ma_hat*cos_m/Rma) );       
 
    return 2*np.sum((force_const[...,NA]/2)*(cos_beta-dih[...,NA])*d_dih_a,axis=1)



def outer_dih_m_gradient(X, neighbours, next_face, prev_face, dih,force_const):
    bcd_p = next_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    bcd_m = prev_face
    
    Rba, ba_hat  = split_norm(X[:,NA]    - X[neighbours]) # b --> a
    Rbp, bp_hat   = split_norm(X[bcd_p]  - X[neighbours]) # b --> p
    Rbm, bm_hat   = split_norm(X[bcd_m]  - X[neighbours]) # b --> m
    
    Rba = Rba[...,NA]
    Rbp = Rbp[...,NA]
    Rbm = Rbm[...,NA]
    
    Xba, Xbp, Xbm = Rba*ba_hat, Rbp*bp_hat, Rbm*bm_hat
      
    mb_hat, Rmb = -bm_hat, Rbm
    
    Xmp = Xbp - Xbm
    Rmp = np.linalg.norm(Xmp,axis=-1)[...,NA]
    mp_hat = Xmp / Rmp
    
    pm_hat, Rpm = -mp_hat, Rmp
    
    Xpa = Xba - Xbp
    Rpa = np.linalg.norm(Xpa,axis=-1)[...,NA]
    pa_hat = Xpa / Rpa
    
    cos_m = np.sum(mb_hat*mp_hat,axis=-1)[...,NA]
    cos_p = np.sum(pm_hat*pa_hat,axis=-1)[...,NA]
    
    sin_m = np.sqrt(1 - cos_m**2)
    sin_p = np.sqrt(1 - cos_p**2)
    
    nbmp_hat = np.cross(mb_hat,mp_hat)/sin_m
    nmpa_hat = np.cross(pm_hat,pa_hat)/sin_p
    
    cos_beta = np.sum(nbmp_hat*nmpa_hat,axis=-1)[...,NA]
    
    cot_p = cos_p / sin_p
    
    d_dih_m = (np.cross(nbmp_hat,pm_hat)/(Rpa*sin_p)) - (pa_hat*cos_beta/Rpa) + \
              (cot_p*cos_beta/(sin_p*Rpa))*(pm_hat - pa_hat*cos_p)
    
    return 2*np.sum((force_const[...,NA]/2)*(cos_beta-dih[...,NA])*d_dih_m,axis=1)


def outer_dih_p_gradient(X, neighbours, next_face, prev_face, dih,force_const): 
    bcd_p = next_face     # next_on_face(a,b)=b+ | next_on_face(a,c)=c+ | next_on_face(a,d)=d+
    bcd_m = prev_face
    
    Rba,  ba_hat  = split_norm(X[:,NA] - X[neighbours]) # b --> a
    Rbp, bp_hat = split_norm(X[bcd_p] - X[neighbours]) # b --> p
    Rbm, bm_hat = split_norm(X[bcd_m] - X[neighbours]) # b --> m
    
    Rba = Rba[...,NA]
    Rbp = Rbp[...,NA]
    Rbm = Rbm[...,NA]
    
    Xba, Xbp, Xbm = Rba*ba_hat, Rbp*bp_hat, Rbm*bm_hat
    
    pb_hat  = -bp_hat
    Rpb     = Rbp
    
    Xpa = Xba - Xbp
    Rpa = np.linalg.norm(Xpa,axis=-1)[...,NA]
    pa_hat = Xpa / Rpa
    
    ap_hat = -pa_hat
    Rap    = Rpa
    
    Xam = Xbm - Xba
    Ram = np.linalg.norm(Xam,axis=-1)[...,NA]
    am_hat = Xam / Ram
    
    cos_p = np.sum(pb_hat*pa_hat,axis=-1)[...,NA]
    cos_a = np.sum(ap_hat*am_hat,axis=-1)[...,NA]

    sin_p = np.sqrt(1 - cos_p**2)
    sin_a = np.sqrt(1 - cos_a**2)
    
    nbpa_hat = np.cross(pb_hat,pa_hat)/sin_p
    npam_hat = np.cross(ap_hat,am_hat)/sin_a
    
    cos_beta = np.sum(nbpa_hat*npam_hat,axis=-1)[...,NA]
    
    cot_p = cos_p / sin_p
    cot_a = cos_a / sin_a
    
    d_dih_p = (np.cross(npam_hat,pb_hat)/(Rpa*sin_p)) - (pa_hat*cos_beta/Rpa) + \
              (cot_p*cos_beta/(Rpa*sin_p))*(pb_hat - pa_hat*cos_p) + \
              (ap_hat * cos_beta / Rap) - (np.cross(am_hat,nbpa_hat)/(Rap*sin_a)) + (am_hat*cos_beta/Ram) - \
              (np.cross(nbpa_hat,ap_hat)/(Ram*sin_a)) + (cot_a*cos_beta/sin_a) *  ((ap_hat*cos_a/Rap) - (am_hat/Rap) + (am_hat*cos_a/Ram)-               (ap_hat/Ram))
    
    return 2*np.sum((force_const[...,NA]/2)*(cos_beta-dih[...,NA])*d_dih_p,axis=1)
############################################################

def gradient(X,neighbours,next_on_face, prev_on_face,k0,f0, UseAll):
    
    R0,ang0,dih0,ang0_p,ang0_m,dih0_a,dih0_m,dih0_p = k0
    fR,fang,fdih,fang_p,fang_m,fdih_a,fdih_m,fdih_p = f0
    
    Ruv,Duv = edge_displacements(X,neighbours)
        
    d_R     = edge_energy_gradient_from_points(X,neighbours,R0,fR)
    d_ang   = corner_cos_gradient(Ruv,Duv,ang0,fang)
    d_dih   = dihedral_cos_energy_gradient_from_points(X,neighbours,dih0,fdih)

    
    if UseAll == 0:
        grad    = d_R + d_ang
        
    if UseAll == 1:
        d_bp    = bcd_corner_cos_gradient(X, neighbours, next_on_face, ang0_p,fang_p)
        d_bm    = bcd_corner_cos_gradient(X, neighbours, prev_on_face, ang0_m,fang_m)
        grad    = d_R + d_ang + d_dih + (d_bp + d_bm)
    
    if UseAll == 2:
        d_bp    = bcd_corner_cos_gradient(X, neighbours, next_on_face, ang0_p,fang_p)
        d_bm    = bcd_corner_cos_gradient(X, neighbours, prev_on_face, ang0_m,fang_m)
        diha    = outer_dih_a_gradient(X, neighbours, next_on_face, prev_on_face, dih0_a,fdih_a)
        dihm    = outer_dih_m_gradient(X, neighbours, next_on_face, prev_on_face, dih0_m,fdih_m)
        dihp    = outer_dih_p_gradient(X, neighbours, next_on_face, prev_on_face, dih0_p,fdih_p)
        
        grad = d_R + d_ang + d_dih + (d_bp + d_bm) + (diha + dihm + dihp)
       
    
    return grad

#################################################################

def Bisection_search(d,X,neighbours,next_on_face, prev_on_face,k0,f0,a,b,N_max,eps,N,UseAll): #Bisection
    #d is search dir
    #X is geometry
    # k0,f0 are the parameters
    #(a,b) brackets
    #N_max max iterations
    #eps is threshold.
    
    
    
    count = 0
    nBrack = 0 
    
    dfc = np.sum(-d*d) # Check size on f'(a)

    
    while np.sign(np.sum(gradient(X + b*d,neighbours,next_on_face, prev_on_face,k0,f0,UseAll)*d)) < 0:
        b *= 1.5
        nBrack +=1
            
    while np.abs(dfc) > eps:
        count+= 1
        c =  (a+b)/2 # Take the point between a=0 and b=x 
        Xc = X + c*d # Take a step in the direction of the gradient
        
        fc  =  gradient(Xc,neighbours,next_on_face, prev_on_face,k0,f0,UseAll)  # Calculate new gradint
        
        dfc = np.sum(fc*d) 
        if count > N_max: 
            #raise Exception(f'linesearch did not converge in {N_max} steps')
            #print(f'At itt: {N} the linesearch did not converge in {N_max} steps')
            return c, Xc, -fc

        if dfc < 0: #If gradient is negatice choose c as new a 
            a = c;
        else:
            b = c; # if gradient is positive chose c as new b
    
    return c, Xc, -fc

def linspace_search(d,X,neighbours,k0,f0,a,b,N,next_on_face, prev_on_face,UseAll): 
    nBrack = 0 
       
    while np.sign(np.sum(gradient(X + b*d,neighbours,next_on_face, prev_on_face,k0,f0,UseAll)*d)) < 0:
        a = b
        b *= 1.5
        nBrack +=1

    d /=np.linalg.norm(d)
    brack = np.linspace(a,b,N)

    dc = brack[1] - brack[0]
    E = array([energy(X + k*d,neighbours,k0,f0) for k in brack])
    c = brack[np.argmin(E)]
    Xc = X + c*d
        
    fc = gradient(Xc,neighbours,next_on_face, prev_on_face,k0,f0,UseAll)
    return c, Xc, -fc

##############################################################################################
def PolakRibiere(r0,r1):
    r0 = r0.reshape(-1);
    r1 = r1.reshape(-1)
    return np.dot(r1.T,r1 - r0)/np.dot(r0.T,r0)

def FletcherReeves(r0,r1):
    r0 = r0.reshape(-1)
    r1 = r1.reshape(-1)
    return np.dot(r1.T,r1)/np.dot(r0.T,r0)
##############################################################################################

def lines_for_plotting(X,neighbours,next_on_face,prev_on_face,fig_number):
    ## Plotting function, only connectivity
    
    neigh     = X[neighbours]
    nextplus  = X[next_on_face]
    prevminus = X[prev_on_face]
    
    a = X
    b,c,d = neigh[:,0,:], neigh[:,1,:], neigh[:,2,:]
    
    bp,cp,dp = nextplus[:,0,:], nextplus[:,1,:], nextplus[:,2,:]
    bm,cm,dm = prevminus[:,0,:], prevminus[:,1,:], prevminus[:,2,:]
    
    bb = np.array([a,b,bm,b,bp])
    cc = np.array([a,c,cm,c,cp])
    dd = np.array([a,d,dm,d,dp])
    
    _,n,_ = bb.shape
    fig = plt.figure(fig_number)
    ax = fig.add_subplot(111, aspect='equal', projection='3d')
    ax.scatter(X[:,0],X[:,1],X[:,2], color="r")

    for i in range(n):
        ax.plot(bb[:,i,0], bb[:,i,1], bb[:,i,2], color='k')#,linewidth=1)
        ax.plot(cc[:,i,0], cc[:,i,1], cc[:,i,2], color='k')#,linewidth=1)
        ax.plot(dd[:,i,0], dd[:,i,1], dd[:,i,2], color='k')#,linewidth=1)
    

    ax.set_xlabel('x [Å]',fontsize=15)
    ax.set_ylabel('y [Å]',fontsize=15)
    ax.set_zlabel('z [Å]',fontsize=15);


    ax.get_autoscale_on()
    plt.show()
    return 


def plotting_faces(geom,pentagons,hexagons,ap=0.3,ah=0.3,ax_off=True):
    ## Plots pentagon and hexagon information on geometry
    fig = plt.figure()
    ax = Axes3D(fig)

    a,b = np.min(geom,axis=0),np.max(geom,axis=0)
    a,b = np.min(a),np.max(b)
    ax.set_xlim([a,b])
    ax.set_ylim([a,b])
    ax.set_zlim([a,b])

    pentaA = Poly3DCollection(geom[pentagons], linewidths=1, alpha=ap,edgecolors='k')
    hexaA  = Poly3DCollection(geom[hexagons], linewidths=1, alpha=ah ,edgecolors='k')

    hexaA.set_facecolor('orange')
    pentaA.set_facecolor('blue')

    ax.set_xlabel('x [Å]',fontsize=15)
    ax.set_ylabel('y [Å]',fontsize=15)
    ax.set_zlabel('z [Å]',fontsize=15);

    ax.add_collection3d(hexaA)
    ax.add_collection3d(pentaA)

    if ax_off:
        ax.axis('off')
    plt.show()
    return 

##################################################################################################

def Parameters_from_geometry(face_right,neigh,MyParam = False):
    ## Calculates the force constants and parameters for the fullerene
    ## From the face information and connectivity
    
    
    (R55, R56, R66) = (1.479,1.458, 1.401); # Bond-lengths in Ångstrom
    R0_constants = np.array([[R55, R56],
                             [R56, R66]]);
    R0 = R0_constants[np.roll(face_right,shift=1,axis=1)-5, face_right-5]; # N x d #shift = +1, we need previous face

    (ang5,ang6) = (108*pi/180,120*pi/180); #bond-angles in degrees
    ang_constants = np.array([ang5,ang6])
    ang0 = np.cos(ang_constants[face_right-5])
    ang0_p = ang0
    ang0_m = np.roll(ang0,shift=1,axis=1) # Need previous face so shift (+1) is correct

    dih_constants = cos(np.array([0.652358,0.509674,0.509674,0.345123,0.615841,0.417884,0.417884,0]).reshape(2,2,2))


    f0 = face_right-5
    f1 = np.roll(f0,shift=-1,axis=1)
    f2 = np.roll(f0,shift=-2,axis=1)
    #################################################
    dih0 = dih_constants[[f0[:,0],f1[:,0],f2[:,0]],[f0[:,1],f1[:,1],f2[:,1]],[f0[:,2],f1[:,2],f2[:,2]]].transpose()
    #################################################
    
    if MyParam == False:
        (fpp,fhp,fhh) = (260.0,390.0,450.0);
        (fp,fh) = (100.0,100.0)
        (fppp,fpph,fphh,fhhh) = (35.0,65.0,85.0,270.0)
    
    if MyParam == True:
        (fpp,fhp,fhh) = (260.0,353.377,518.992);
        (fp,fh) = (207.924,216.787)
        (fppp,fpph,fphh,fhhh) = (35.0,65.0,3.772,270.0)


    fR_constants = np.array([[fpp, fhp],
                             [fhp, fhh]]);
    fR = fR_constants[np.roll(face_right,shift=1,axis=1)-5, face_right-5];

    fang_constants = np.array([fp,fh])
    fang = fang_constants[face_right-5]

    fang_p = fang
    fang_m = np.roll(fang,shift=1,axis=1)

    fdih_constants = np.array([[[fppp,fpph],
                               [fpph,fphh]],
                              [[fpph,fphh],
                               [fphh,fhhh]]])
    fdih = fdih_constants[face_right-5, np.roll(face_right,shift=1,axis=1)-5, np.roll(face_right,shift=2,axis=1)-5] 

    bfr,cfr,dfr = sort_neighbour_faces(face_right,neigh)
    fdih_b = fdih_constants[bfr-5, np.roll(bfr,shift=1,axis=1)-5, np.roll(bfr,shift=2,axis=1)-5]
    fdih_c = fdih_constants[cfr-5, np.roll(cfr,shift=1,axis=1)-5, np.roll(cfr,shift=2,axis=1)-5]
    fdih_d = fdih_constants[dfr-5, np.roll(dfr,shift=1,axis=1)-5, np.roll(dfr,shift=2,axis=1)-5]

    f0b = bfr-5
    f1b = np.roll(f0b,shift=-1,axis=1)
    f2b = np.roll(f0b,shift=-2,axis=1)
    dih0_b = dih_constants[[f0b[:,0],f1b[:,0],f2b[:,0]],[f0b[:,1],f1b[:,1],f2b[:,1]],[f0b[:,2],f1b[:,2],f2b[:,2]]].transpose()
    fdih_b = fdih_constants[bfr-5, np.roll(bfr,shift=-1,axis=1)-5, np.roll(bfr,shift=-2,axis=-1)-5]
    #################################################
    f0c = cfr-5
    f1c = np.roll(f0c,shift=-1,axis=1)
    f2c = np.roll(f0c,shift=-2,axis=1)
    dih0_c = dih_constants[[f0c[:,0],f1c[:,0],f2c[:,0]],[f0c[:,1],f1c[:,1],f2c[:,1]],[f0c[:,2],f1c[:,2],f2c[:,2]]].transpose()
    fdih_c = fdih_constants[cfr-5, np.roll(cfr,shift=-1,axis=1)-5, np.roll(cfr,shift=-2,axis=1)-5]
    # #################################################
    f0d = dfr-5
    f1d = np.roll(f0d,shift=-1,axis=1)
    f2d = np.roll(f0d,shift=-2,axis=1)
    dih0_d = dih_constants[[f0d[:,0],f1d[:,0],f2d[:,0]],[f0d[:,1],f1d[:,1],f2d[:,1]],[f0d[:,2],f1d[:,2],f2d[:,2]]].transpose()
    fdih_d = fdih_constants[dfr-5, np.roll(dfr,shift=-1,axis=1)-5, np.roll(dfr,shift=-2,axis=1)-5]
    ###################################################33
    dih_a = array([dih0_b[:,0],dih0_c[:,0],dih0_d[:,0]]).T
    dih_m = array([dih0_b[:,1],dih0_c[:,1],dih0_d[:,1]]).T
    dih_p = array([dih0_b[:,2],dih0_c[:,2],dih0_d[:,2]]).T

    fdih_a = array([fdih_b[:,0],fdih_c[:,0],fdih_d[:,0]]).T
    fdih_m = array([fdih_b[:,1],fdih_c[:,1],fdih_d[:,1]]).T
    fdih_p = array([fdih_b[:,2],fdih_c[:,2],fdih_d[:,2]]).T
    
    #print('R0,ang0,dih0,ang0_p,ang0_m,dih0_b,dih0_c,dih0_d,fR,fang,fdih,fang_p,fang_m,fdih_a,fdih_m,fdih_p')
    return R0,ang0,dih0,ang0_p,ang0_m,dih_a,dih_m,dih_p,fR,fang,fdih,fang_p,fang_m,fdih_a,fdih_m,fdih_p


########################################################################

def fit_plane(points):
    N = len(points)
    CM = (np.sum(points,axis=0)/N)
    M = points - CM[NA,:]
    U, S, V = la.svd(M)
    n = V[np.argmin(S)]
    ev = np.min(S)
    err = np.mean(np.abs((M@n)))
    return n, err, ev


############################################################################
## Functions for reading Gaussian16 ".log" files.

import re
def read_geometry(file):
    
    f = open(file, "r")
        
    lines = f.readlines()
    
    N = None
    for i in range(len(lines)):
        if isinstance(int(re.search(r'\d+', file).group()),int):
            N = int(re.search(r'\d+', file).group())

        elif "NAtoms=" in lines[i]:
            N = re.search(r'\d+', lines[i]).group()

        else:
            print('What is the number of atoms N ?')

            
        if "Standard orientation:" in lines[i]:
            coord_lines = lines[i+5:i+5+N]
            coord = np.array([l.split()[-3:] for l in coord_lines],dtype=float)
            return coord
    return None



def get_Hessian(file,N):
    f = open(file, "r")
    txt = f.read()
    f.close()
    hess = re.split("Hessian after L703:", txt)[1] 
    hess = re.split("Leave Link  703", hess)[0]
    hess = hess.strip()
    column_index_line = re.compile("(\s*[0-9]+\s+[0-9]+\s+[0-9]+\s+[0-9]+\s+[0-9]+\s+)")
 
    Hessian = np.zeros([N,N,3,3])
    
    for line in hess.splitlines():
        if column_index_line.match(line):
            ## line is an column index line ##
            column_id = np.fromstring(re.sub("\s+", ",",  line.strip()), sep=",", dtype=int) - 1
        else:
            ## line is not an column index line ##
            line = re.sub("D", "E", line)
            data_line = np.fromstring(re.sub("\s+", ",",  line.strip()), sep=",")
                  
            if len(data_line) == 0:
                return Hessian
            
            row_id = int(data_line[0]) - 1
            
            data_line = data_line[1:]
            for i in range(len(data_line)):
                atom_row = row_id//3
                atom_col = column_id[i]//3

                Hessian[atom_row,atom_col, row_id%3, column_id[i]%3] = data_line[i]
                Hessian[atom_col,atom_row, row_id%3,column_id[i]%3] = data_line[i]
    return Hessian


