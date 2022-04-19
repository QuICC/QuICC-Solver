"""
This module includes important functions to rotate an hdf5 state 
and perform the integral of vorticity from the spectral velocity field

"""
import numpy as np
import h5py
import sys
sys.path.append('/home/nicolol/workspace/QuICC/Python/')
from quicc.geometry.spherical import shell_radius
from quicc.projection import shell

##Wigner matrix 

def dlmb(L):
    #INPUT:
    #L       Maximum angular degree
    #OUTPUT:
    #D       Lower right quarter of the Wigner D-matrix m>=0
    #d       Masters' concatenated output
    # Computes matrix elements for spherical harmonic polar rotation around
    # the y-axis over 90 degrees only. 
    #Rotation matrix
    #  D_{mm'}(a,b,g) = exp(-ima) d_{mm'}(b) exp(-im'g) 
    # but we factor the rotation itself into:
    #    R(a,b,g)=R(a-pi/2,-pi/2,b)R(0,pi/2,g+pi/2)
    # thus we only need to compute d_{mm'} for b=90.
    # After a code by T. Guy Masters.
    # See also McEwen, 2006.
    # Last modified by fjsimons-at-alum.mit.edu, 08/05/2008
                                  
    d=np.zeros(np.sum( ( np.array(range(L+1)) + 1 )**2 ) )
    # Initialize using D&T C.115.
    # l = 0 
    d[0] = 1
    # l = 1 
    if L >= 1:
        d[1] = 0
        d[2] = 1.0 / np.sqrt(2)
        d[3] = -1.0 / np.sqrt(2)
        d[4]= 0.5 ;
        
    #pointer index
    ind = 5
    #factor
    f1= 0.5
    #Lwait = 100
    #Recursions
    for l in range(2,L+1):
        lp1 = l + 1
        knd = ind + lp1 
        #print(knd)
        fl2p1 = l + lp1
        vect = np.array( range(1,l+1) )
        f = np.sqrt( vect * ( fl2p1 - vect) )
        f1 = f1 * ( 2.0 * l - 1.0 ) / ( 2.0 * l )
        #print f1
        d[knd-1] = -np.sqrt(f1)
        d[knd-2] = 0
        for i in range(2,l+1):
            j = knd-i
            #print('j=',j)
            d[j-1] = -f[i-2] * d[j+1] / f[i-1]
        #print d

        #Positive N (bottom triangle)
        f2 = f1
        g1 = l 
        g2 = lp1
        for N in range(1,l+1):
            knd = knd + lp1
            en2 = N + N
            g1 = g1 + 1
            g2 = g2 - 1
            f2 = f2 * g2 / g1
            #print(f2)
            d[knd - 1] = -np.sqrt(f2)
            d[knd - 2] = d[knd-1]*en2/f[0]
            #print d[knd-2]
            for i in range(2, l-N+1):
                j = knd-i
                d[j-1] = ( en2 * d[j] - f[i-2] * d[j+1] ) / f[i-1]
                #print d[j-1]

        #Fill upper triangle and fix signs
        for j in range(1,l+1):
            for m in range(j,l+1):
                d[ind+m*lp1+j-1]=d[ind+j*lp1+m-l-1]

        isn=1+np.mod(l,2)
        for n in range(0,l+1):
            knd=ind+n*lp1
            for i in range(isn,lp1+1,2):
                d[knd+i-1]=-d[knd+i-1]
        ind=ind+lp1*lp1;

    #Now let's rearrange the coefficients as 1x1, 2x2, 3x3 etc rotation
    #matrices.
    cst=1;
    D=np.empty(L+1,dtype=object)
    #Start of coefficient sequence; need transpose!
    for l in range(1,L+2):
        #Leo: This line doesn't work !!!
        #print l
        #print(len(d[cst-1:cst+l*l-1]))
        #print(np.reshape(d[cst-1:cst+l*l-1],(l,l)))
        D[l-1]=np.reshape(d[cst-1:cst+l*l-1],(l,l))
        cst=cst+l*l
        #print(cst)
    return (D,d)

def printS(string):
    print(string)

def computeAverageVorticity(state):
    f = h5py.File(state, 'r') #read state
    nR=int(np.array(f['/truncation/spectral/dim1D'], dtype='double'))+1

    #Toroidal Velocity
    dataT=np.array(f['/velocity/velocity_tor'].value)

    # obtain the map factor for the Tchebyshev polynomials
    ro = f['/physical/ro'].value
    rratio = f['/physical/rratio'].value
    a, b = shell_radius.linear_r2x(ro,rratio)

    # compute the weight of the operator
    ri = ro* rratio
    delta = (f['/physical/ekman'].value)**.5*10
    volume = ((ro-delta)**5 - (ri+delta)**5 )/ (5. * (3. / (4 *np.pi))**.5)

    # define boundary conditions
    bc = {0:0, 'cr':2}
    R2 = shell_radius.r2(nR+2, a, b, bc)
    bc['cr'] = 1
    R1 = shell_radius.r1(nR+3, a, b, bc)
    I1 = shell_radius.i1(nR+4, a, b, bc)

    Pi = shell.proj(nR+4, a, b, np.array([ro-delta, ri+delta]))

    proj_vec = Pi*I1*R1*R2
    proj_vec = (proj_vec[0,:]-proj_vec[1,:])/volume
    omegax = -(proj_vec*dataT[2,:][:,0])[0]*2**.5
    omegay = (proj_vec* dataT[2,:][:,1])[0]*2**.5
    omegaz = (proj_vec* dataT[1,:][:,0])[0]

    return np.array([omegax, omegay, omegaz])




def rotateState(state, omega, field):
    f = h5py.File(state, 'r+') #read state
    LL=int(np.array(f['/truncation/spectral/dim2D'], dtype='double'))
    MM=int(np.array(f['/truncation/spectral/dim3D'], dtype='double'))
    NN=int(np.array(f['/truncation/spectral/dim1D'], dtype='double'))

    #Toroidal Velocity
    #data=np.array(f[field].value)
    data=f[field].value
    """
    #Poloidal Velocity
    dataP=np.array(f['/velocity/velocity_pol'].value)
    #Codensity
    dataC=np.array(f['/Codensity/Codensity'].value)
    """



    ##Computing rotation axis and allocating memory 
    #f.keys()
    #print('running some stuff')
    DD=dlmb(LL)

    ###Determine rotation axis from spectral coefficients
    ##Leo: Check indices this again
    #divide by normalization of n(1,1) mode 
    #(-1) from condon-shortley phase

    #Data[l+l*m, n, real-imag]
    #1. / norm(W(n=0,l=1)) = sqrt(pi) / 2 = 0.88627
    #n=0
    #Y10=dataT[1,0][0]/0.886227 #OK! 
    #divide by normalization of W(n=0,l=1)=sqrt(pi) / 4 = 0.443113
    #Y11R=-dataT[2,0][0]/0.443113  #OK
    #Y11I=dataT[2,0][1]/0.443113 #OK

    #n=0
    """
    Y10=dataT[1,0][0]*2.0/np.sqrt(np.pi)
    Y11R=-dataT[2,0][0]*4.0/np.sqrt(np.pi)  #OK
    Y11I=dataT[2,0][1]*4.0/np.sqrt(np.pi) #OK
    """

    #n=1
    #no ekman removal
    #d = 0.
    """
    #Ekman layer removal
    d = 10.0*np.sqrt(E)
    ekmanR = (-2+40*(d-2)*d)/(7*np.sqrt(np.pi))
    Y10 = Y10 + dataT[1,1][0]*ekmanR 
    Y11R = Y11R - dataT[2,1][0]*2*ekmanR  #OK
    Y11I = Y11I + dataT[2,1][1]*2*ekmanR #OK

    #n=2
    ekmanR = 10.*(-5. + 4.*(-2. + d)*d*(11.+28.*(-2.+d)*d))/(63.*np.pi)
    Y10 = Y10 + dataT[1,2][0]*ekmanR 
    Y11R = Y11R - dataT[2,2][0]*2*ekmanR  #OK
    Y11I = Y11I + dataT[2,2][1]*2*ekmanR #OK

    #+Omega Leo: Adding Omega_0 to the flow because we are in the mantle frame?
    #Y10=Y10+1.0/E #For viscous time-scale
    Y10=Y10+1.0 #For rotation time-scale
    """

    #How to remove solid body rotation??

    #f.close()
    # replace the



    #[phi0, theta0]: rotation axis of fluid
    if (omega[1]>=0 and omega[0]>=0):
        phi0=np.arctan(omega[1]/omega[0])
    elif (omega[0]<0):
        phi0=np.arctan(omega[1]/omega[0])+np.pi
    else:
        phi0=np.arctan(omega[1]/omega[0])+2*np.pi

    theta0=np.arctan(np.sqrt(omega[0]**2+omega[1]**2)/(omega[2]+1))

    #### determine axis of fluid from uniform vorticity
    # using thetat0 and phi0 determined from mean vorticity
    # theta0=0.837615978648410;
    # phi0=5.569620149823718;
    #Euler Angles
    #Alpha: 1st rotation about Z
    #R = R(alpha-pi/2, -pi/2, beta)* R(0,pi/2, gamma+pi/2)
    alp=-phi0-np.pi/2
    #Beta: 2nd rotation about X
    bta=-theta0
    #Gamma: 3rd rotation about Z
    gam=0


    #e(i*m*alpha) = cos(m*alpha) + i * sin(m*alpha)
    #calp = cos(m*alpha)
    #salp = sin(m*alpha)
    mmVect=np.array(range(0,int(MM)+1))
    llVect=np.array(range(0,int(LL)+1))

    calp=np.cos(mmVect*alp)
    salp=np.sin(mmVect*alp)

    #cbta = cos(l*bta)
    #sbta = sin(l*bta)
    cbta=np.cos(llVect*bta)
    sbta=np.sin(llVect*bta)
    
    #Schmidt normalization factors
    NF=np.array([])
    ind =  0
    for l in range(0, int(LL)+1):
            for m in range(0,min(l,MM)+1):
                if (m==0):
                    #print(l,m)
                    #NF[ind]=1.0
                    #NF=1.0
                    NF=np.append(NF,1.0) #Schmidt normalization for m = 0 
                else:
                    #print(l,m)
                    #NF[ind]=np.sqrt(2)/2;
                    #NF=np.sqrt(2)/2
                    #NF=np.append(NF,np.sqrt(2)/2) #Schmidt normalization for m != 0
                    NF = np.append(NF, 1.)  # Schmidt normalization for m != 0
                
                ind=ind+1;
            
    #ind = 529 (same as matlab)

    #Rotating EPM results

    for n in range(0, NN+1 ): #N worland polynomials
    #for n in range(1,2): #First polynomial
        #Toroidal variables
        alp=np.zeros((int(((LL+2)*(LL+1))/2),2))
        bta=np.zeros((int(((LL+2)*(LL+1))/2),2))
        #Poloidal Variables
        alpP=np.zeros((int(((LL+2)*(LL+1))/2),2))
        btaP=np.zeros((int(((LL+2)*(LL+1))/2),2))
        #Codensity variables
        alpC=np.zeros((int(((LL+2)*(LL+1))/2),2))
        btaC=np.zeros((int(((LL+2)*(LL+1))/2),2))
        
        ind=0
        for l in range(0, LL+1 ):
        #for l in range(0,3+1): #debug
            for m in range(0, min(l,MM)+1 ): #Leo fix: m from 0 to l
                #print ('ind, l, m', ind, l , m )
                  ### Azimuthal rotation (Z) by (ALPHA - PI/2)
                #Dividing by normalisation 
                Cos = data[ind,n][0] / NF[ind] #index, n, real
                Sin = data[ind,n][1] / NF[ind] #index, n, imag
                #Rotate about Z by alpha
                alp[ind,0] = Cos * calp[m] + Sin * salp[m]
                alp[ind,1] = Sin * calp[m] - Cos * salp[m]

                """
                #Dividing poloidal by normalisation 
                CosP = dataP[ind,n][0] / NF[ind] #real
                SinP = dataP[ind,n][1] / NF[ind] #imag
                #Rotate about Z by alpha 
                alpP[ind,0] = CosP * calp[m] + SinP*salp[m]
                alpP[ind,1] = SinP * calp[m] - CosP*salp[m]
                
                #Dividing Codensity by normalisation 
                CosC = dataC[ind,n][0] / NF[ind] #real
                SinC = dataC[ind,n][1] / NF[ind] #imag
                #Rotate about Z by alpha 
                alpC[ind,0] = CosC * calp[m] + SinC*salp[m]
                alpC[ind,1] = SinC * calp[m] - CosC*salp[m]
                """
                
                ind=ind+1
                
            ###### tilt by bta
            li=l+1
            i,j=np.meshgrid(range(1,li+1),range(1,li+1))
            #checkboard * 2 
            IC=(((i+j)+li%2)%2)*2
            IC[:,0]=1 #m=0
            #inverted checkboard *2
            IS=((((i+j)+li%2)%2)<1)*2
            IS[:,0]=1 #m=0
            
            # STEP 1: PASSIVE colatitudinal (Y) rotation over -PI/2
            #Toroidal
            Cp = np.dot( DD[0][l].T * IC, \
                         alp[ int((l+1)*(l+2) /2 -l-1) : int((l+1) * (l+2) / 2), 0]) #X Tor
            Sp = np.dot( DD[0][l].T * IS, \
                         alp[ int((l+1)*(l+2) /2 - l - 1) : int((l+1) * (l+2) / 2), 1]) #Y Tor
            """
            #Poloidal
            CpP = np.dot( DD[0][l].T * IC, \
                         alpP[ (l+1)*(l+2) /2 - l-1 : (l+1)*(l+2) /2, 0]) #X Pol
            SpP = np.dot( DD[0][l].T * IS, \
                         alpP[ (l+1)*(l+2) /2 - l-1 : (l+1)*(l+2) /2, 1]) #Y Pol
            
            #Codensity
            CpC = np.dot( DD[0][l].T * IC, \
                         alpC[ (l+1)*(l+2) /2 - l-1 : (l+1)*(l+2) /2, 0]) #X Pol
            SpC = np.dot( DD[0][l].T * IS, \
                         alpC[ (l+1)*(l+2) /2 - l-1 : (l+1)*(l+2) /2, 1]) #Y Pol
            """
            
            # STEP 2: PASSIVE azimuthal (Z) rotation over BETA 
            Cpp = Cp * cbta[0:l+1].T + Sp * sbta[0:l+1].T
            Spp = Sp * cbta[0:l+1].T - Cp * sbta[0:l+1].T

            """
            CppP = CpP * cbta[0:l+1].T + SpP * sbta[0:l+1].T
            SppP = SpP * cbta[0:l+1].T - CpP * sbta[0:l+1].T
            
            CppC = CpC * cbta[0:l+1].T + SpC * sbta[0:l+1].T
            SppC = SpC * cbta[0:l+1].T - CpC * sbta[0:l+1].T
            """

            #print('CppT2:', CppT) 
            #matlab: 0, 37.1985 
            #python:  -1.00485917e-14,   3.72031861e+01
            #more or less Ok 
            
            # STEP 3: PASSIVE colatitudinal (Y) rotation over PI/2
            #Toroidal
            Cpp = np.dot( DD[0][l] * IC, Cpp)
            Spp = np.dot( DD[0][l] * IS, Spp)

            """
            #Poloidal
            CppP = np.dot( DD[0][l] * IC, CppP)
            SppP = np.dot( DD[0][l] * IS, SppP)
            
            #Codensity
            CppC = np.dot( DD[0][l] * IC, CppC)
            SppC = np.dot( DD[0][l] * IS, SppC)
            """

            
            #print('CppT3:', CppT) 
            #print('DD[0][li-1]', DD[0][li-1] )
            
            ###### STEP4: PASSIVE azimuthal rotation over (GAMMA +PI/2)
            #Toroidal
            bta[int((l+1)*(l+2)/2-l-1):int((l+1)*(l+2)/2),0] = Cpp
            bta[int((l+1)*(l+2)/2-l-1):int((l+1)*(l+2)/2),1] = Spp

            """
            #Poloidal
            btaP[(l+1)*(l+2)/2-l-1:(l+1)*(l+2)/2,0] = CppP
            btaP[(l+1)*(l+2)/2-l-1:(l+1)*(l+2)/2,1] = SppP
            
            #Codensity
            btaC[(l+1)*(l+2)/2-l-1:(l+1)*(l+2)/2,0] = CppC
            btaC[(l+1)*(l+2)/2-l-1:(l+1)*(l+2)/2,1] = SppC
            """

            #print('CppT4:',CppT)
            #matlab: 52.6066, 0
                 
        ind=0
        for l in range(0, LL+1 ):
            for m in range(0, min(l,MM)+1 ): #Leo fix: m from 0 to min(l,M)
                #print("ind, l , m ", ind, l, m)
                #Multiplying by Schmidt normalization 
                #Toroidal
                data[ind,n][0]=bta[ind,0]*NF[ind]
                data[ind,n][1]=bta[ind,1]*NF[ind]
                """
                #Poloidal
                dataP[ind,n][0]=btaP[ind,0]*NF[ind]
                dataP[ind,n][1]=btaP[ind,1]*NF[ind]
                
                #Codensity
                dataC[ind,n][0]=btaC[ind,0]*NF[ind]
                dataC[ind,n][1]=btaC[ind,1]*NF[ind]
                """
                ind=ind+1
                
        #print('dataT:',dataT[6,0])

    rotatedState = data#(dataT, dataP, dataC)
    #f[field].value = rotatedState
    f[field][:] = data
    f.close()
    return rotatedState

def selectModes(state, modes, field):
    f = h5py.File(state, 'r+') #read state
    LL=int(np.array(f['/truncation/spectral/dim2D'], dtype='double'))
    MM=int(np.array(f['/truncation/spectral/dim3D'], dtype='double'))
    NN=int(np.array(f['/truncation/spectral/dim1D'], dtype='double'))

    # impose single modes on field
    data=f[field].value
    """
    #Poloidal Velocity
    dataP=np.array(f['/Velocity/VelocityPol'].value)
    #Codensity
    dataC=np.array(f['/Codensity/Codensity'].value)
    """

    #E=f['PhysicalParameters/E'].value
    f.close()
    
    #Schmidt normalization factors
    NF=np.array([])
    ind =  0
    

    ind=0
    for l in range(0, LL+1 ):
        for m in range(0, min(l,MM)+1 ): #Leo fix: m from 0 to min(l,M)
            if m not in modes:

                #Toroidal
                data[ind,:][0]=0.
                data[ind,:][1]=0.
                """
                #Poloidal
                dataP[ind,n][0]=0.
                dataP[ind,n][1]=0.
            
                #Codensity
                dataC[ind,n][0]=0.
                dataC[ind,n][1]=0.
                """

            ind=ind+1
                
        #print('dataT:',dataT[6,0])

    #singleMode = (dataT, dataP, dataC)
    f[field][:] = data
    f.close()
    return data

def removeModes(state, modes, field):
    f = h5py.File(state, 'r+') #read state
    LL=int(np.array(f['/Truncation/L'], dtype='double'))
    MM=int(np.array(f['/Truncation/M'], dtype='double'))
    NN=int(np.array(f['/Truncation/N'], dtype='double'))

    #Toroidal Velocity
    data=f['/velocity/velocity_tor'].value
    """
    #Poloidal Velocity
    dataP=np.array(f['/velocity/velocity_pol'].value)
    #Codensity
    #dataC=np.array(f['/Codensity/Codensity'].value)
    """

    #Schmidt normalization factors
    NF=np.array([])
    ind =  0

    ind=0
    for l in range(0, LL+1 ):
        for m in range(0, min(l,MM)+1 ): #Leo fix: m from 0 to min(l,M)
            if m in modes:
                data[ind,:][0]=0.
                data[ind,:][1]=0.
                """
                #Poloidal
                dataP[ind,n][0]=0.
                dataP[ind,n][1]=0.
            
                #Codensity
                dataC[ind,n][0]=0.
                dataC[ind,n][1]=0.
                """

            ind=ind+1
                
        #print('dataT:',dataT[6,0])

    f[field][:] = data
    f.close()
    return data

def correctRotation(state, toBeCorrected = ['/velocity/velocity_tor','/velocity/velocity_pol']):
    try:
        omega = computeAverageVorticity(state)
        for field in toBeCorrected:
            rotateState(state,omega, field)



    except Exception as e:
        print(e)
        pass