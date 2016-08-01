import numpy
c = [[.14714, -.15583, 1.1234], [-1.1836, -.37607, -.52894], [-.04554, 1.8719, -.39071]]
w=.7085
b = [[w*w], [w], [1]]
m = numpy.dot(c, b)	
ls = (104*0.055158)+1.2378 #numbers for before 4/2/2007 (after that date .055376 and 1.18)
ts = 1260.56/(numpy.log((607.76/(ls))+1))
by = 1256
y = (ts*ts)/(by*ls)
d = ts-((ts*ts)/by)
t = y*(1.052631578947368*((m[0]*ls)+m[1])+m[2])+d