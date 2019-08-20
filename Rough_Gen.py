import numpy as np
import numpy.ma as ma
import os
import argparse

parser = argparse.ArgumentParser(description='Produce ASCII raster of roughness values a la Casas et al, 2010')
parser.add_argument('ZRef_File', help='ASCII Raster of D values to calculate Epsilon')
parser.add_argument('Depth_File', help='ASCII Raster of h values to caclulate Epsilon')
parser.add_argument('N_File', help='Name for ASCII Raster output of roughness values')
parser.add_argument('Background_N',help='Background Roughness value of channel, where there is no vegetation')
args = parser.parse_args()

# Define function to load Z ref values in
def LoadZRef(ZRef_File):
    # Load array of Zref
    zref = np.genfromtxt(ZRef_File,delimiter = ' ',skip_header = 6)
    # Remove all negative values and NoData
    zref[zref < 0.] = 0.
    return zref

# Define function to load depth values in
def LoadDepth(Depth_File):
    # Load array of depth
    depth = np.genfromtxt(Depth_File,delimiter = ' ',skip_header = 6)
    depth = depth[1:,:-1]
    depth[depth<0] = 0.
    return depth

#Create divide function to avoid dividing by zero
def divide_zero(n, d):
    return n / d if d else 0


# Define function to calculate array of Epsilon values
def Epsilon(zref,depth):
    epsilon = np.zeros_like(zref)
    for i in range(len(zref[:,0])):
        for j in range(len(zref[0,:])):
            epsilon[i,j] = divide_zero(depth[i,j],zref[i,j])
    epsilon[(epsilon < 1.) & (epsilon != 0.)] = 1.
    epsilon[epsilon > 7.] = 0.
    return epsilon

# Define function to calculate array of function of alpha and epsilon
def Func(alpha,epsilon):
    function = np.zeros_like(epsilon)
    for i in range(len(epsilon[:,0])):
        for j in range(len(epsilon[0,:])):
            div_term = divide_zero(1.,epsilon[i,j])
            function[i,j] = 1.+(alpha*(div_term)*np.log((np.cosh((1./alpha)-((1./alpha)*epsilon[i,j])))/(np.cosh(1./alpha))))
    return function

# Define function to calculate array of roughness values
def NCalc(Cu,g,depth,function):
    n = (depth**(1./6.)/((g**(1./2.))*Cu*function))
    return n

# Define function to read in raster ASCII header and create output raster ASCII
def ReadWriteRaster(inraster,outraster,array):
    from numpy import savetxt
    with open(inraster,"r") as file:
        header = {}
        for i in range(6):
            header[i]=file.readline()
    file.close()
    with open(outraster,'w') as file:
        for i in range(6):
            file.write(header[i])
    file.close()
    savetxt('vals.txt',array)
    datafile = open('vals.txt','r')
    data = datafile.read()
    datafile.close()
    with open(outraster,'a') as file:
        file.write(data)
    file.close()

def Main(ZRef_File,Depth_File,N_File,Background_N):
    # Define constants
    Cu = 4.5
    a = 1.
    g = 9.8 # m/s^2
    # Load in values from rasters as arrays
    zref = LoadZRef(ZRef_File)
    depth = LoadDepth(Depth_File)
    print('depth is ' + str(depth.max()) + ' ' + str(depth.mean())+ ' ' + str(depth.min()))
    # Produced masked array of roughness
    eps = Epsilon(zref,depth)
    print(eps.shape)
    print('Zref is '+str(zref.max()) + ' ' + str(zref.mean())+ ' ' + str(zref.min()))
    print('Eps is '+str(eps.max()) + ' ' + str(eps.mean())+ ' ' + str(eps.min()))
    func = Func(a,eps)
    print('func is '+str(func.max()) + ' ' + str(func.mean())+ ' ' + str(func.min()))
    n = NCalc(Cu,g,depth,func)
    print('N is '+str(n.max()) + ' ' + str(n.mean())+ ' ' + str(n.min()))
    # Fill masked values with default of 0.04
    Background_N = float(Background_N)
    n[n == 0.] = Background_N
    print('N is '+str(n.max()) + ' ' + str(n.mean())+ ' ' + str(n.min()))
    n[n < Background_N] = Background_N
    print('N is '+str(n.max()) + ' ' + str(n.mean())+ ' ' + str(n.min()))
    print(str(n[n > 0.04]))
    ReadWriteRaster(ZRef_File,N_File,n)

Main(args.ZRef_File,args.Depth_File,args.N_File,args.Background_N)


