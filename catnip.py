#!/bin/python

import sys
import os
import random
import math
import time

start=time.time()

help=False

if "-h" in sys.argv:
    help=True

try:
    sys.argv[1]
except:
    help=True

needs=['-c','-ci','-p']

if not help:
    for thing in needs:
        if thing not in sys.argv:
            print "missing",thing
            help=True

if help:
    print ""
    print "       catnip is used to put residues from one file into another"
    print ""
    print "-c   : input system coordinate file (.gro) to which molecules are added"
    print "-oc  : [ optional | catnip.gro ] output system coordinate file"
    print "-p   : input topology file (.top) which corresponds to the system coordinate file"
    print "-op  : [ optional | catnip.top ] output topology file (.top)"
    print "-ci  : input molecule file (.top) which is put in the system file"
    print "-dl  : [ optional | some.file ] List of residues that can be removed. 'no RES' prevents a hardcoded residue from being removed."
    print "-n   : [ optional | prompted ] number of new molecules to put into the system"
    print "-g   : [ optional | prompted ] group name, number, or 'interface RES1 (or ceil) RES2 x||y||z'"
    print "-ndx : [ optional | index.ndx ] specify index file"
    print "-h   : helpful info"
    sys.exit()
# useful info ^

deleted=0
i=0
outbox="catnip.gro"
outtop="catnip.top"
# default outs
delgroups=[]
indgroups=[]
groupies=[]
thisgroup=[]
r3=[0,1,2]
g1=6.02314
g2=3.14159
catdel=[]

try:
    sys.path.append('/home/Tony/program-workshop/')
    import f77math
    def sqdist(this,that):
        return f77math.f77sqdist(this[0],this[1],this[2],that[0],that[1],that[2],box[0],box[1],box[2])

    def vectorize(this,that):
        return list(f77math.f77vectorize(this[0],this[1],this[2],that[0],that[1],that[2],box[0],box[1],box[2]))

except:
    print "could not import fortran math"
    print "using pythonic math instead, go grab a snack"
# sooper-dooper optimized (jk)
# roughly 2x slower than f77math
    def sqdist(this,that):
# return squared distance between points
        xyzd=[min([(this[f]-that[f])**2,(this[f]-that[f]-box[f])**2]) for f in r3]
        return sum(xyzd)

    def vectorize(this,that):
# return vector from one point to another
        if abs(that[0]-this[0]) > abs(box[0]-abs(that[0]-this[0])):
            if that[0] > this[0]:
                xdv = box[0]-abs(that[0]-this[0])
            else:
                xdv = (this[0]-that[0])-box[0]
        else:
            xdv = this[0] - that[0]
        if abs(that[1]-this[1]) > abs(box[1]-abs(that[1]-this[1])):
            if that[1] > this[1]:
                ydv = box[1]-abs(that[1]-this[1])
            else:
                ydv = (this[1]-that[1])-box[1]
        else:
            ydv = this[1] - that[1]
        if abs(that[2]-this[2]) > abs(box[2]-abs(that[2]-this[2])):
            if that[2] > this[2]:
                zdv = box[2]-abs(that[2]-this[2])
            else:
                zdv = (this[2]-that[2])-box[2]
        else:
            zdv = this[2] - that[2]
        return [xdv,ydv,zdv]
 
def move(this,that):
# move this molecule away from that molecule
    a=sqdist(xyz[boxnum.index(rngesus)],this)
    b=sqdist(xyz[boxnum.index(rngesus)],that)
    c=sqdist(this,that)
    templist=oldxyz
    try:
        psi=math.acos((c-a-b)/(-2.*math.sqrt(a*b)))
    except:
        psi = 0
    if abs(psi) > 0.1 or attempt <= 4:
        try:
            dxyz=vectorize(this,that)
            n=math.sqrt(sum([f**2 for f in dxyz]))
            dxyz=[f/n for f in dxyz]
        except:
            dxyz=[1.0,0.0,0.0]
        dxyz=[0.25*f for f in dxyz]
        tryxyz=[[templist[i-1][j]+dxyz[j]+xyz[boxnum.index(rngesus)][j] for j in r3] for i in newat]
        return tryxyz
# translate it away
    else:
        a1=that[0]-this[0]
        a2=that[1]-this[1]
        a3=that[2]-this[2]
        b1=that[0]-xyz[boxnum.index(rngesus)][0]
        b2=that[1]-xyz[boxnum.index(rngesus)][1]
        b3=that[2]-xyz[boxnum.index(rngesus)][2]
        vx=(a2*b3)-(a3*b2)
        vy=-(a1*b3)+(a3*b1)
        vz=(a1*b2)-(a2*b1)
        n=math.sqrt(vx**2+vy**2+vz**2)
        try:
            vx=vx/n
            vy=vy/n
            vz=vz/n
        except:
            vx=1.0
            vy=0.0
            vz=0.0
        if psi > 0:
            R=[[math.cos(0.3)+vx*vx*(1-math.cos(0.3)),vx*vy*(1-math.cos(0.3))-vz*math.sin(0.3),vx*vx*(1-math.cos(0.3))+vy*math.sin(0.3)],[vx*vy*(1-math.cos(0.3))+vz*math.sin(0.3),math.cos(0.3)+vy*vy*(1-math.cos(0.3)),vy*vz*(1-math.cos(0.3))-vx*math.sin(0.3)],[vx*vz*(1-math.cos(0.3))-vy*math.sin(0.3),vy*vz*(1-math.cos(0.3))+vx*math.sin(0.3),math.cos(0.3)+vz*vz*(1-math.cos(0.3))]]
        else:
            R=[[math.cos(-0.3)+vx*vx*(1-math.cos(-0.3)),vx*vy*(1-math.cos(-0.3))-vz*math.sin(-0.3),vx*vx*(1-math.cos(-0.3))+vy*math.sin(-0.3)],[vx*vy*(1-math.cos(-0.3))+vz*math.sin(-0.3),math.cos(-0.3)+vy*vy*(1-math.cos(-0.3)),vy*vz*(1-math.cos(-0.3))-vx*math.sin(-0.3)],[vx*vz*(1-math.cos(-0.3))-vy*math.sin(-0.3),vy*vz*(1-math.cos(-0.3))+vx*math.sin(-0.3),math.cos(-0.3)+vz*vz*(1-math.cos(-0.3))]]
        i=0
        while i < len(templist):
            dx=templist[i][0]*R[0][0]+templist[i][1]*R[0][1]+templist[i][2]*R[0][2] 
            dy=templist[i][0]*R[1][0]+templist[i][1]*R[1][1]+templist[i][2]*R[1][2]
            dz=templist[i][0]*R[2][0]+templist[i][1]*R[2][1]+templist[i][2]*R[2][2]
            templist[i]=[dx,dy,dz]
            i+=1 
             
        tryxyz=[[templist[i-1][j]+xyz[boxnum.index(rngesus)][j] for j in r3] for i in newat]
        return tryxyz
# rotate it
         
while i < len(sys.argv):
    if sys.argv[i] == "-c":
        box=sys.argv[i+1]
    if sys.argv[i] == "-oc":
        outbox=sys.argv[i+1]
    if sys.argv[i] == "-p":
        topfile=sys.argv[i+1]
    if sys.argv[i] == "-op":
        outtop=sys.argv[i+1]
    if sys.argv[i] == "-ci":
        molfile=sys.argv[i+1]
    if sys.argv[i] == "-dl":
        dfile=sys.argv[i+1]
    if sys.argv[i] == "-n":
        tries=int(sys.argv[i+1])
    if sys.argv[i] == "-g":
        group=sys.argv[i+1]
        if group == "interface":
            g1=sys.argv[i+2]
            g2=sys.argv[i+3]
            try:
                lim=sys.argv[i+4]
                if lim not in ['x','y','z']: 
                    print lim,'not x, y, or z; setting to z'
                    lim='z'
            except:
                lim='z'
    if sys.argv[i] == "-ndx":
        ndxfile=sys.argv[i+1]
  
    i+=1
# read command line stuff 
  
try:
    tries
except:
    tries=raw_input("How many would you like added? ")

tries=int(tries)
# if not given a command line option, it will ask how many to add

if tries < 1:
    print "must add at least once"
    sys.exit()

try:
    if os.path.isfile(dfile):
        delfile=open(dfile,"r+")
        for line in delfile:
            thing = line.split(',')
            for thang in thing:
                if thang not in delgroups:
                    delgroups.append(thang) 
        delfile.close()
# read file for deleting stuff
    else:
        delgroups=dfile.split(",")
# or read commandline stuff

except:
    pass

if "SOL" not in delgroups:
    delgroups.append("SOL")
if "CHX" not in delgroups:
    delgroups.append("CHX")
if "CHCL" not in delgroups:
    delgroups.append("CHCL")
for thing in delgroups:
    if thing.split()[0] == "no":
        while thing.split()[1] in delgroups:
            delgroups.remove(thing.split()[1])
    while delgroups.count(thing) > 1:
        delgroups.remove(thing)
# add or remove things from the list

print "instructions for removable molecules:",delgroups

try:
    ndxfile
except:
    if os.path.isfile("index.ndx"):
        print "removing old index.ndx"
        os.remove("index.ndx")
        print "specify ndx file with -ndx to use an old ndx file"
    os.system("echo q | gmx make_ndx -f "+box+" &>/dev/null")
    ndxfile="index.ndx"
# make new index if there isn't one

index=open(ndxfile,"rw+")

groups=-1

print "reading index file..."

for line in index:
    try:
        if line.split()[0] == '[':
            groupies.append([])
            groupies[groups]=thisgroup
            indgroups.append(line)
            thisgroup=[]
            groups=groups+1
        else:
            thisgroup=thisgroup+(line.split())
    except:
        pass

groupies[groups]=thisgroup
# read index for group name and the atoms in each group

indgroups=[thing[2:-3] for thing in indgroups]
# remove brackets

index.close()

try:
    group
except:
    for thing in indgroups:
        print indgroups.index(thing),thing
    group = raw_input("Where would you like to place the new molecules? ")
    if group == "interface":
        g1=raw_input("group one: ")
        g2=raw_input("group two: ") 
        lim=raw_input("x||y||z: ")
try:
    group=indgroups[int(group)]
except:
    pass
try: 
    g1=indgroups[int(g1)]
except:
    pass
try:
    g2=indgroups[int(g2)]
except:
    pass
try:
    lim=int(lim)
except:
    try:
        lim
        if lim == "x":
            lim=0
        elif lim == "y":
            lim=1
        elif lim == "z":
            lim=2
        else:
            print "3rd option after -g must be integer between 0 and 2 or x, y, or z"
            sys.exit()
    except:
        pass
# making sure you know what group you want removed 

boxfile=open(box,"r+")

print "reading box file..."

boxlines=boxfile.readlines()

molnum=[]
boxmol=[]
boxele=[]
boxnum=[]
xyz=[]
box=[]

boxatoms=int(boxlines[1])
for j in boxlines[2:-1]:
    molnum.append(int(j[:5]))
    boxmol.append(j[5:9].strip())
    boxele.append(j[12:15].strip())
    boxnum.append(int(j[15:20]))
    xyz.append(j[21:45].split())

xyz=[[float(n) for n in m] for m in xyz]

box=boxlines[-1].split()
box=[float(i) for i in box]
# read system box stuff

boxfile.close()

if group == "interface":
    print "finding interface"
    i=0
    c1=0
    c2=0.
    com1=0
    com2=0
    groups=groups+1
    if g1 not in boxmol:
        try:
            while i < boxatoms:
                if str(boxnum[i]) in groupies[indgroups.index(g1)]:
                    com1=com1+xyz[i][lim]
                    c1+=1
                i+=1
        except:
            pass
    else:
        while i < boxatoms:
            if boxmol[i] == g1:
                com1=com1+xyz[i][lim]
                c1+=1
            i+=1
    i=0
    if g2 not in boxmol:
        try:
            while i < boxatoms:
                if str(boxnum[i]) in groupies[indgroups.index(g2)]:
                    com2=com2+xyz[i][lim]
                    c2+=1
                i+=1
        except:
            pass
    else:
        while i < boxatoms:
            if boxmol[i] == g2:
                com2=com2+xyz[i][lim]
                c2+=1
            i+=1
    try:
        com1=com1/float(c1)
    except:
        if g1 == "ceil":
            com1=box[lim]
        else:
            com1=0.
    try:
        com2=com2/float(c2)
    except:
        if g2 == "ceil":
            com2=box[lim]
        else:
            com2=0.

    if com1 > com2:
        top=com1
        bot=com2
    else:
        top=com2
        bot=com1

    indgroups.append("interface")
    thing=[]
    for i in boxnum:
        if xyz[boxnum.index(i)][lim] < top and xyz[boxnum.index(i)-1][lim] > bot:
            thing.append(i)
    groupies.append([])
    groupies[groups]=thing
# find interface if needed

llptc=open(molfile,"r+")

newnum=[]
newmol=[]
newele=[]
newat=[]
newxyz=[]

print "reading new molecule file..."

newlines=llptc.readlines()

newatoms=int(newlines[1])

for s in newlines[2:-1]:
    newnum.append(int(s[:5].strip()))
    newmol.append(s[5:9].strip())
    newele.append(s[12:15].strip())
    newat.append(int(s[15:20].strip()))
    newxyz.append(s[21:45].split())

llptc.close()

added=0

newxyz=[[float(j) for j in i] for i in newxyz]

cx=0.
cy=0.
cz=0.

for thing in newxyz:
    cx=cx+thing[0]/len(newxyz)
    cy=cy+thing[1]/len(newxyz)
    cz=cz+thing[2]/len(newxyz)

oldxyz=[[p[0]-cx, p[1]-cy, p[2]-cz] for p in newxyz]
    
while added < tries:

    added=added+1

    rngesus=int(groupies[indgroups.index(group)][random.randrange(0,len(groupies[indgroups.index(group)]))])
# pick the chosen one

    print "placing new molecule",added,"by",rngesus
    tryxyz=[[xyz[boxnum.index(rngesus)][j]+oldxyz[i-1][j] for j in r3] for i in newat]

    attempt=0
    deleted=0

    print "adjusting orientation..."

    i=0
    k=len(xyz)
    kk=len(tryxyz)

    while i < k:
        j=0
        if boxmol[i] not in delgroups:
            while j < kk:
                if sqdist(xyz[i],tryxyz[j]) < 0.06:
                    attempt=attempt+1
                    tryxyz=move(tryxyz[j],xyz[i])
                    i=0
                    j=0
                if attempt > 30:
                    break
                j+=1
        i+=1
        if attempt > 30:
            break
# it's got some tries to move the molecule away if it's too close to something it can't delete
      
    print "removing molecules..."

    j=0
  
    while j < k:
        i=0
        if boxmol[j] in delgroups and molnum[j] not in catdel:
            while i < kk:
                if sqdist(tryxyz[i],xyz[j]) < 0.012:
                    catdel.append(molnum[j])
                    deleted=deleted+molnum.count(molnum[j])
                    break
                i+=1
        j+=1
# delete atoms that are too close

    mols=molnum[len(molnum)-1]
    for i in newat:
         ik=newat.index(i)
         molnum.append(newnum[ik]+mols)
         boxmol.append(newmol[ik])
         boxele.append(newele[ik])
         boxnum.append(i+boxatoms)
         xyz.append(tryxyz[ik])

    boxatoms=boxatoms+newatoms-deleted
# update number of atoms in the box

print "writing output files..."

catcoord=open(outbox,"w")
thetop=open(topfile,"r")
cattop=open(outtop,"w")

templist=[]

for line in thetop:
    templist.extend(line)
    try:
        if line.split()[0]+line.split()[1]+line.split()[2] == "[molecules]":
            break
# copy everything from the old topology up until [ molecules ]
    except:
        pass

thetop.close()

cattop.write(''.join(templist))


catcoord.write("catnip conffile\n")
catcoord.write(str(boxatoms)+"\n")
# write new number of atoms to system

nextmol=1
nextatom=0
molsofthat=0

i=0
k=len(boxnum)
templist=[]

while i < k:
    if molnum[i] not in catdel:
        if i > 1:
            if molnum[i] != molnum[i-1]:
                nextmol+=1
                molsofthat+=1
# renumber atoms and molecules so they are all consecutive
            if boxmol[i] != boxmol[i-1]:
                cattop.write(boxmol[i-1]+"    "+str(molsofthat)+"\n")
                molsofthat=0
# write new number of things to topology
        nextatom+=1
        templist.extend(("%5i%-4s  %4s%5i  %6.3f  %6.3f  %6.3f\n" % (nextmol,boxmol[i],boxele[i],nextatom,xyz[i][0],xyz[i][1],xyz[i][2])))
    i+=1
    if i == k:
        molsofthat+=1
        cattop.write(boxmol[i-1]+"    "+str(molsofthat)+"\n")
# write all atoms to the new file

catcoord.write(''.join([j for j in templist]))

catcoord.close()
cattop.close()

os.system("echo '    '"+str(box[0])+"'    '"+str(box[1])+"'    '"+str(box[2])+" >> "+outbox)
# echo the last line so vmd can open the new file?

stop=time.time()

print "took",str(stop-start)+"s"
