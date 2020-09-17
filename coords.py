import math
import cairo
import numpy as np
delta=(1+math.sqrt(2))
pi=math.pi
def inflate(tiles) :
    newtiles=[]
    for tile in tiles :
        #print tile
        typ,coords=tile
        if(typ==0) :
            a,b,c=coords
            p1= a +(c-a)/delta;
            p2= a + (b-a)*1.0/(1.0+delta)
            p3= a + (b-a)*delta/(1.0+delta)
            p4= c + (b-c)*1.0/delta
            p5= p1 + (b-a)*1.0/(1.0+delta)
            newtiles.append([1,[p5,p1,a,p2]])
            newtiles.append([1,[p3,p4,c,p5]])
            newtiles.append([0,[c,p1,p5]])
            newtiles.append([0,[p3,p2,p5]])
            newtiles.append([0,[b,p4,p3]])
        elif(typ==1) :
            a,b,c,d=coords
            p1= c+(d-c)*1.0/delta
            p2= c+(b-c)*1.0/delta
            p3= a+(b-a)*1.0/delta
            p4= a+(d-a)*1.0/delta
            p5= p1+(a-d)*1.0/delta
            p6= p4+ (b-a)*1.0/delta

            newtiles.append([1,[p5,p1,c,p2]])
            newtiles.append([1,[b,p6,d,p5]])
            newtiles.append([1,[p6,p4,a,p3]])
            newtiles.append([0,[d,p4,p6]])
            newtiles.append([0,[b,p3,p6]])
            newtiles.append([0,[d,p1,p5]])
            newtiles.append([0,[b,p2,p5]])
    return newtiles


#create square to start
#tiles=[]
#r=math.sqrt(2)/2.0
#a=-r*(1+1j)
#b=r*(1+1j)
#c=r*(1-1j)
#d=r*(-1+1j)
#tiles.append([0,[a,b,c]])
#tiles.append([0,[a,b,d]])

#create central octagon
tiles=[]
for i in range(8) :
    a=0
    p1=math.cos(2*pi*i/8) + 1j*math.sin(2*pi*i/8)
    p2=math.cos(2*pi*(i+1)/8) + 1j*math.sin(2*pi*(i+1)/8)
    d=2*math.cos(pi/8.0)*(math.cos(2*pi*(2*i+1)/16) + 1j*math.sin(2*pi*(2*i+1)/16))
    if(i%2==0) :
        b,c=p1,p2
    else :
        b,c=p2,p1
    tiles.append([1,[a,b,d,c]])
    e=b*(1+math.sqrt(2))
    f=c*(1+math.sqrt(2))
    tiles.append([0,[e,b,d]])
    tiles.append([0,[f,c,d]])
    abd=abs(d)
    g=d*(abs(d)+abs(b-c))/abs(d)
    tiles.append([1,[f,d,e,g]])


n_infl=2;
for i in range(n_infl) :
  tiles=inflate(tiles)

typ,coords=tiles[0]
if(typ==0) :
  a,b,c=coords
else :
  a,b,c,d=coords

size=abs(b-c)
sinv=size**-1.0

fname="tiles_inf%d.dat"%(n_infl)
f=open(fname,'w')
print len(tiles)
for tile in tiles :
    typ,coords=tile
    if(typ==0) :
      a,b,c=coords
      print>>f,"%d %.16f %.16f %.16f %.16f %.16f %.16f %f %f"%(typ,sinv*a.real,sinv*a.imag,sinv*b.real,sinv*b.imag,sinv*c.real,sinv*c.imag,-1,-1)
    else :
      a,b,c,d=coords
      print>>f,"%d %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f"%(typ,sinv*a.real,sinv*a.imag,sinv*b.real,sinv*b.imag,sinv*c.real,sinv*c.imag,sinv*d.real, sinv*d.imag)
f.close()
    




