#-----------------------------------Required Modules starts---------------------------------------------->
import numpy as np
from netCDF4 import Dataset  
import matplotlib.pyplot as plt
import pylab as pylab
from scipy.optimize import *
import math
from numpy import pi
import os, sys
from lmfit import Model
import pandas as pd
from IPython.core.display import display
from itertools import product
from sklearn import linear_model
from sklearn.linear_model import LinearRegression, Ridge, BayesianRidge
from matplotlib.mlab import griddata



#---------------------------------Required Modules ends---------------------------------------------------->

#-------------------------- Main Attributes printing functions--------------------------------------------->

def ncdump(nc_fid, verb=True):
    
    print "=============================================="
    """
    ncdump outputs dimensions, variables and their attribute information.
    """
    def print_ncattr(key):
        
        """
        Prints the NetCDF file attributes for a given key
        """
        try:
            
            #print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                #print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    nc_attrs = nc_fid.ncattrs()
    if verb:
        #print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    if verb:
        #print "NetCDF dimension information:"
        for dim in nc_dims:
            #print "\tName:", dim 
            #print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    nc_vars = [var for var in nc_fid.variables]
    if verb:
       # print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                #print '\tName:', var
                #print "\t\tdimensions:", nc_fid.variables[var].dimensions
                #print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

#main  printing functions ends --------------------------------------------------------------> 

#-------------------------------various array declarations starts------------------------------>
ed_final=[]
fs_new=[]
amp=[]
hour=[]
time=[]
minute=[]
second=[]
ed_new=[]
lon_new=[]
alt_new=[]
file_new=[]
lon_final=[]
time_final=[]
sub_lon=[]
sub_time=[]
#--------------------------------array declaration ends------------------------------------->


#----------------------------file management loop starts---------------------------------------->
j=.00123
while j<.047:
	c=2013.000+j
	path = "./projfiles/"+str(c)[:8]+"/"
	day1 = os.listdir( path )
	
    
	for m in range(len(day1)):	
		nc_f = day1[m]
		
	
		nc_fid = Dataset(nc_f, 'r')  
		nc_attrs, nc_dims, nc_vars = ncdump(nc_fid)
		mt = nc_fid.variables['MSL_alt'][:] 
		glt = nc_fid.variables['GEO_lat'][:]
		gln = nc_fid.variables['GEO_lon'][:]
		oa1 = nc_fid.variables['OCC_azi'][:] 
		tc1 = nc_fid.variables['TEC_cal'][:]
		ed1 = nc_fid.variables['ELEC_dens'][:]

		file = day1[m]
		
		
		ncfile =Dataset(file,'r')
		msl_alt = np.array(ncfile.variables['MSL_alt'][:],dtype=np.float32)
		temp = np.array(ncfile.variables['TEC_cal'][:],dtype=np.float32)
		ed = np.array(ncfile.variables['ELEC_dens'][:],dtype=np.float32)
		oa = np.array(ncfile.variables['OCC_azi'][:],dtype=np.float32)
		g_lat = np.array(ncfile.variables['GEO_lat'][:],dtype=np.float32)
		g_lon = np.array(ncfile.variables['GEO_lon'][:],dtype=np.float32)
		ho = nc_fid.getncattr('hour')
		mi = nc_fid.getncattr('minute')
		sec = nc_fid.getncattr('second')
		fs=nc_fid.getncattr('fileStamp')
		
        #----------------------main project loop starts---------------------------------------->
		
		for i in range(len(msl_alt)):
			if 199<(msl_alt[i])<201  and -10<=g_lat[i]<=10 and ed[i]>=0:
				hour.append(ho)
				minute.append(mi)
				second.append(sec)
				lon_new.append(g_lon[i])
				ed_new.append(ed[i])
				fs_new.append(fs)
				
		#-------------------------main project loop ends---------------------------------------->


	#---------------------------------------curve fitting code starts---------------------------->
	"""
	def sinefunction(x, a, b, c):
	    return a + b * np.sin(x*np.pi/180.0 + c)

	smodel = Model(sinefunction)
	result = smodel.fit(y, x=xdeg, a=0, b=70000, c=0)
	amp.append(abs(result.params['b'].value))
	#print "Amplitude",amp
	#print(result.fit_report())
	pylab.xlabel('Longitude------>')
	pylab.ylabel('E_density------>')
	
	plt.plot(xdeg, y, 'o', label='data')
	plt.plot(xdeg, result.best_fit,'*', label='fit')
	plt.legend()
	#plt.show()
	"""
	#-------------------------------------curve fitting code ends---------------------------------->
	j=j+.00100

#----------------------------------------file management loop ends---------------------------------->

for i in range(len(hour)):
	time.append((hour[i])+(float(minute[i])/60)+(second[i]/3600))

data={'Time':time,'Longitude':lon_new,'Electron Density':ed_new,'fileStamp':fs_new}
dataTable=pd.DataFrame(data)
#display(dataTable)	

#---------------------------------------------Floor functioning code starts-------------------------->
count=0
n=0
ed_sum=0
h=0
l=0
inner=outer=0
while outer<8:
	inner=0
	l=0
	while inner<9:
		for i in range(len(lon_new)):
			if (-180+l)<=lon_new[i]<=(-140+l) and (0+h)<=time[i]<=(3+h):
				count=count+1
				if count>1:
					lon_final.append(-180+l)
					time_final.append(0+h)
					break
		if count<=1:
			sub_lon.append(-180+l)
			sub_time.append(0+h)
		l=l+40
		ed_sum=0
		inner=inner+1
		count=0
	h=h+3
	outer=outer+1
#---------------------------------2D array creation starts---------------------------------------------->


#--------------------------------------------Floor Functioning code ends ------------------------------->

#-------------------------------Mean Electron Density Calculation starts-------------------------------->
count=0
n=0
h=0
l=0
inner=outer=0
while outer<8:
	inner=0
	l=0
	while inner<9:
		for i in range(len(lon_new)):
			if (-180+l)<=lon_new[i]<=(-140+l) and (0+h)<=time[i]<=(3+h):
				ed_sum=ed_sum+ed_new[i]
				count=count+1

		if count>1 and (ed_sum/float(count))>0:
			ed_final.append(ed_sum/count)
			count=0

		l=l+40
		ed_sum=0
		inner=inner+1
	h=h+3
	outer=outer+1

# ------------------------------Mean Electron Density Calculation ends---------------------------------->


#-------------------------------Machine Learning starts--------------------------------------------->

X1 = np.vstack((lon_final, time_final)).T
X2=np.vstack((sub_lon, sub_time)).T
Y=ed_final
reg = linear_model.BayesianRidge()
reg.fit(X1, Y)
BayesianRidge(alpha_1=1e-05, alpha_2=1e-05, compute_score=False, copy_X=True,
       fit_intercept=True, lambda_1=1e-06, lambda_2=1e-06, n_iter=300,
       normalize=False, tol=0.001, verbose=False)


for x in X2:
	ed_final.append(reg.predict ([x]))

lon_final=lon_final+sub_lon
time_final=time_final+sub_time

#---------------------------------Machine Learning Ends------------------------------------------------->
	
print lon_final
print time_final

#-----------------------------------Raw graph between longitude and time start--------------------------> 

fig, ax = plt.subplots()
ax.set_axisbelow(True)
ax.xaxis.set_ticks(np.arange(-180,181, 40))
ax.yaxis.set_ticks(np.arange(0,25, 3))
plt.xlabel('Longitude(deg)')
plt.ylabel('time(hours)')
plt.title('Within +/-5 degree latitude')
plt.plot(lon_new,time, 'ro')

ax.grid(linestyle='-', linewidth='2', color='green')
plt.grid(True)
plt.show()


#-------------------------------------Raw graph between longitude and time ends---------------------------> 

#----------------------------------Floor Functioned graph plotting starts--------------------------------->

fig, ax = plt.subplots()
ax.set_axisbelow(True)
ax.xaxis.set_ticks(np.arange(-180,181, 40))
ax.yaxis.set_ticks(np.arange(0,25, 3))
plt.xlabel('Longitude(deg)')
plt.ylabel('time(hours)')
plt.title('Within +/-5 degree latitude')
plt.plot(lon_final,time_final, 'go')
plt.plot(sub_lon,sub_time,'ro')
ax.grid(linestyle='-', linewidth='0.5', color='black')
plt.grid(True)
plt.show()

#------------------------------------Floor Functioned graph plotting starts-------------------------->


#----------------------------------------------Scatter plotting starts------------------------------->

fig, ax = plt.subplots()
ax.set_axisbelow(True)
plt.scatter(lon_final,time_final,40,c=ed_final)
ax.xaxis.set_ticks(np.arange(-180,181, 40))
ax.yaxis.set_ticks(np.arange(0,25, 3))
plt.xlabel('Longitude(deg)')
plt.ylabel('time(hours)')
plt.title('Within +/-5 degree latitude')
ax.grid(linestyle='-', linewidth='.5', color='red')
cbar=plt.colorbar()
cbar.set_label("Electron Density(/cm3)",labelpad=2)
plt.show()

#------------------------------------------------Scatter plotting ends--------------------------------->

#--------------------------------------contour plot starts------------------------------------------->
"""
X, Y = np.meshgrid(lon_final, time_final)
z = ed_final
plt.figure()
Z = griddata(lon_final, time_final, z, X, Y,interp='linear')
levels = [-9.21, -4.61, -2.3, 0.5]
contour = plt.contour(X, Y, Z, levels)
plt.tricontour(lon_final,time_final, z, levels)

#cp = plt.contourf(X, Y, Z)
#plt.colorbar(cp)
plt.title('Filled Contours Plot')
#plt.xlabel('x (cm)')
#plt.ylabel('y (cm)')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

x = [4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10]
y = [-6.95138e-06, -9.09998e-07, 8.24384e-06, 3.4941e-06, 5.08276e-06, 7.82652e-06, -4.7378e-06, -1.40027e-05, -1.62638e-05, -3.97604e-06, 3.19294e-06, 2.50123e-06, -4.13063e-06, -6.60289e-06, -4.02982e-06, -2.32882e-06, -3.86464e-06, -1.09167e-05, -9.42387e-06, -3.48118e-07, 5.22e-06, 8.74445e-06, 1.35842e-05, 2.33632e-05, 2.71328e-05, 1.747e-05, 2.32177e-06, -7.25386e-06, -9.75881e-06, -2.99633e-06, 1.19281e-06, -4.24077e-06, -7.4252e-06, -4.54435e-07, 1.03078e-05, 1.14579e-05, 3.90613e-06, -4.77174e-06, -9.25321e-06, -8.36579e-06, -3.0257e-06, -1.69309e-06, -5.36534e-06, -4.01092e-06, 1.20577e-06, 5.13284e-06, 5.06792e-06, 4.81178e-06, 5.9607e-06, 6.70492e-06, 3.45118e-06, 2.51942e-06, 1.23012e-06, 2.09802e-06, 1.44658e-06, -8.93274e-08, -5.14753e-06, -9.93717e-06, -7.91692e-06, -4.12816e-06, -6.33457e-06]
z = [-0.63, -0.02, -1.05, -0.22, -0.51, -1.26, -0.53, -4.97, -7.32, -0.44, -0.30, -0.19, -0.55, -1.48, -0.58, -0.20, -0.57, -5.00, -3.84, -0.01, -1.17, -3.31, -8.13, -22.59, -30.52, -13.69, -0.27, -2.88, -5.48, -0.49, -0.08, -0.86, -2.94, -0.01, -4.32, -5.49, -0.69, -1.15, -4.70, -3.73, -0.44, -0.14, -1.49, -0.77, -0.07, -1.08, -1.04, -0.93, -1.42, -1.76, -0.51, -0.25, -0.07, -0.18, -0.09, -0.00, -1.08, -5.03, -2.64, -0.65, -1.65]

levels = [-9.21, -4.61, -2.3, 0.5]  # Levels must be increasing.
plt.tricontour(x, y, z, levels)
plt.show()
"""