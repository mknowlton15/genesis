#set up constants matrix
#retrieve water data for each scene desired and create water matrix
#do the matrix math to get the phi values
#find Tsen value using landsat constants and the thermal band radiance
#use Tsen just calculated and thermal band radiance to calculate planck's function factors
#for each scene calculate the temperature using the calculated phi values and the radiance values at each pixel in the mask
import sys, os
import arcpy
import numpy


c = [[.14714, -.15583, 1.1234], [-1.1836, -.37607, -.52894], [-.04554, 1.8719, -.39071]]
w = 5 #give a year and doy as parameters then retrieve w from a csv
b = [[w*w], [w], [1]]
m = numpy.dot(c, b)
ls = 5 #need to actually pull the path radiance of the thermal band here
ts = 1284.3/(numpy.log((671.62/(ls))+1))
by = 1290
y = (ts*ts)/(by*ls)
d = ts-((ts*ts)/by)
s = y*((1/.96)*(((m[0])*ls)+m[1])+m[2])+d
print s