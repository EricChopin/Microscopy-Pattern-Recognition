#################################################################################################
#Authors: Eric and Yvan Chopin , August 2021
#Zones separation identification program (Immunofluorescence microscopy)
# 
#This program performs some data analysis on two microscopy images of the same sample
#(comparison of zone-separations intensities between two 16 bits 1200x1200 tiff images 
#of the same sample under two different laser lights, triggering different fluorescence emissions). 
#Zone separation patterns are searched on images *_C0.tiff and their intensities compared with *_C1.tiff 
#
#This program was realized as an exercice to help a student after a short term traineeship 
#in biochemistry (focusing on microtubules in the cytoskeleton). 
#The images were taken in immunofluorescence microscopy by a PhD student  
#################################################################################################
from PIL import Image, ImageFilter, ImageDraw, ImageFont
import tifffile
import numpy as np  
import math
from matplotlib import pyplot as plt 
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import datetime

###################################
# Basic functions 
###################################

# function used to log what the program is doing 
def addlog(s,t=0):
    if t==1:
        with open(fld2+logfile,"a") as tfile:
            tfile.write(str(datetime.datetime.now()) + ":"+s+"\n")
    else:
        with open(fld2+logfile,"a") as tfile:
            tfile.write(s+"\n")

# function used to log modifications in a specific zone of the image for debug purposes
def scanpx(x,y):
#    if x>= 326 and x<=349 and y>=604 and y<=627:
#        addlog("zone="+str(zoneMaxNb)+" ajout x="+str(x)+" y="+str(y)+" cur subzone="+str(subzones[x+dimx*y]))
    return 0
	
# draw a square around identified zones
def zoneSquare_indentification(xmin,xmax,ymin,ymax,n):
    def drawsquare(xmin,xmax,ymin,ymax):	
        for x in range(xmin, xmax+1):
            img1.putpixel((x,ymin),(0,255,0))
            img1.putpixel((x,ymax),(0,255,0))       
        for y in range(ymin, ymax+1):
            img1.putpixel((xmin,y),(0,255,0))
            img1.putpixel((xmax,y),(0,255,0))

    drawsquare(xmin,xmax,ymin,ymax)
    draw = ImageDraw.Draw(img1)
# write zone number on the picture in hexadecimal format using font "Courrier New"
    font = ImageFont.truetype('cour.ttf', 12)
    draw.text((xmax+7,ymin),hex(n)[2:],(255,0,0),font=font)

# compute x,y position of a pixel from its position within the tiff file
def OneToTwo(p):
    x = p%dimx
    return x, int((p-x)/dimx)

#At position p in the image, look at pixel on a circle of radius r pixels around p 
#return mean and sigma of intensity of pixels on this circle
#(not used so far... may be used to identify round zones like phase separations)
#in practise we used a simpler characterization of zone separations in small squares 
def circletest(p,r):
    x,y=OneTwoTwo(p)
    if x>=r and x< dimx-r and y>=r and y<dimy-r :
        intensity=[]
        for itheta in range(8*r):
            xp = x+int(r*math.cos(np.pi*itheta/4/r))
            yp = y+int(r*math.sin(np.pi*itheta/4/r))
            intensity.append(tr[xp+dimx*yp])
        return mean(intensity), var(intensity)**0.5
    else:
        return 0,0

#not used: for debug purposes to draw specific zones on the image
def drawsq(xmin,xmax,ymin,ymax):
    global img2
    for x in range(xmin, xmax+1):
        img2.putpixel((x,ymin),(0,255,0))
        img2.putpixel((x,ymax),(0,255,0))       
    for y in range(ymin, ymax+1):
        img2.putpixel((xmin,y),(0,255,0))
        img2.putpixel((xmax,y),(0,255,0))
 
#gaussian function to smooth the image to get a more reliable identification of zones
def GaussianFilter(P):
    epsilon = 0.05
    sigma = P*1.0/math.sqrt(-2*math.log(epsilon))
    h = np.zeros((2*P+1,2*P+1))
    som = 0
    for m in range(-P,P+1):
        for n in range(-P,P+1):
            h[m+P][n+P] = math.exp(-(n*n+m*m)/(2*sigma*sigma))
            som += h[m+P][n+P]
    h = h/som
    return h

# intensities in the original tiff file are expanded linearly 
# to get a visible image in a png file. NB; computations are performed 
# on the original intensities from the tiff file, not from the "linearly mapped" image
def linearmap(ar):
    tmin = mini
    tmax = maxi
    k = 65535.0/(tmax-tmin)
    for i in range(dimx*dimy):
        ar[i] = tmin+int(math.floor(k*(ar[i]-tmin)))
    return ar

#From an array describing an histogram, find the position of the maximum.
def maxhist(t):
    tmin = min(t)
    tmax = max(t)
    h=[]
    for i in range(10):
        h.append(0)
    for iter in range(5):
        d = (tmax-tmin)/10
        for i in range(10):
            h[i]=0
        nb=0
        imax=10
        hmax=0 
        for item in t:
            if item>=tmin and item<tmax:
                i = int(math.floor((item-tmin)/d))
                nb+=1
                h[i]+=1
                if h[i]>hmax:
                    hmax=h[i]
                    imax=i
        tmin = tmin + imax*d
        tmax = tmin + d
    return tmin, tmax

# compute average intensity on the border of a square of width = 2*r+1
def average(x,y,r):
    sum=0.0
    np=0
    for i1 in range(-r,r+1):
        if x+i1>-1 and x+i1<dimx and y-r>-1:
            np+=1
            sum +=tr[x+i1+dimx*(y-r)]
        if x+i1>-1 and x+i1<dimx and y+r<dimy:
            np+=1
            sum +=tr[x+i1+dimx*(y+r)]
    
    for j1 in range(-(r-1),r):
        if y+j1>-1 and y+j1<dimy and x-r>-1:
            np+=1
            sum +=tr[x-r+dimx*(y+j1)]
        if y+j1>-1 and y+j1<dimy and x+r<dimx:
            np+=1
            sum +=tr[x+r+dimx*(y+j1)]

    return sum/np

# compute average intensity on the border of a square of width = 2*r+1
def average2(ar,zon):
    sum=0.0
    np=0
    for x1 in range(zoneinfo[zon][0]-3,zoneinfo[zon][1]+4):
        np+=1
        sum +=ar[x1+dimx*(zoneinfo[zon][2]-3)]
        np+=1
        sum +=ar[x1+dimx*(zoneinfo[zon][3]+3)]
    
    for y1 in range(zoneinfo[zon][2]-3,zoneinfo[zon][3]+4):
        np+=1
        sum +=ar[zoneinfo[zon][0]-3+dimx*y1]
        np+=1
        sum +=ar[zoneinfo[zon][1]+3+dimx*y1]

    return sum/np



def testpix(ima,x,y,pxinit):
    global px 
    if x<1 or x>dimx-1 or y<1 or y>dimy-1:
        return 1
    pix1=ima.getpixel((x,y))
    px1=float(tr[x+dimx*y])
    meanval=average(x,y,7)
    if px1 < cellBg or px < cellBg or px1 < meanval or px < meanval:
        return 1
    if (pix1==(255,0,0) or pix1==(255,255,0) or pix1==(0,255,0) or abs(px-px1) < 0.2*(px-meanval)) and px1-meanval > 0.3*(pxinit-meanval):
       px=px1
       if not(pix1==(0,255,0)):
           ima.putpixel((x,y),(255,255,0))

       scanpx(x,y)

       if subzones[x+dimx*y] == 0:
           subzones[x+dimx*y]= zoneMaxNb
           zones[x+dimx*y]= zoneMaxNb
       else: 
           # updating the zone number of all the subzones mapping to the zone of the pixel we have touched
           subzo = subzones[x+dimx*y]
           init_zone = zonemap[subzo]
           if init_zone != zoneMaxNb:
               for zo in range(len(zonemap)):
                   if zonemap[zo]==init_zone :
                       zonemap[zo]= zoneMaxNb
                       #addlog("remap at "+str(x)+","+str(y)+" subzo="+str(subzo)+" => "+str(zonemap[subzo]) +" new mapping "+str(zo)+" => "+str(zonemap[zo]))
       return 0
    else:
       return 1

#test if the pixel is either yellow (inner part of the zone), or red (maximum of intensity), or green (zone boundary) 
def pixOK(p):
    if p==(255,0,0) or p==(0,255,0) or p==(255,255,0):
        return 1
    else:
        return 0

#interpolation function used by the Sobel algorithm that detects the zone boundaries
def interpolation(ar,x,y):
    i = math.floor(x)
    j = math.floor(y)
    t = x-i
    u = y-j
    u1 = 1.0-u
    t1 = 1.0-t
    if j==dimy-1:
        if i==dimx-1:
            return ar[i+dimx*j]
        return t*ar[i+dimx*j]+t1*ar[i+dimx*(j+1)]
    if i==dimx-1:
        return u*ar[i+dimx*j]+u1*ar[i+1+dimx*j]
    return t1*u1*ar[i+dimx*j]+t*u1*ar[i+1+dimx*j]+t*u*ar[i+1+dimx*(j+1)]+t1*u*ar[i+dimx*(j+1)]

########################################################################################################
#
#          MAIN PROGRAM STARTS HERE
#
########################################################################################################
#Some file path definition

# path to the tiff images to be processed  
fld="C:/tiffimages/"
#path where the log file is stored
fld2="C:/tiffimages/python/"
#name of the log file
logfile="zone_separation_recognition.log"


####################################################################
#image to be processed (we set here the file name of the *_C0.tiff 
#image, then the file name of the corresponding *_C1.tiff image 
#is deduced from the first filename
####################################################################
filename="rpe1wt bfp clip noc - 21_XY1625146426_Z0_T0_C0.tiff"
#####################################################################

#Opening the first tiff image
img = tifffile.imread(fld+filename)

#name of intermediate and final result files 
filename2 ="zlin_" + filename
filename_bmp = filename2.replace(".tiff",".bmp")
filename_txt = filename2.replace(".tiff",".txt")



##########################################################
 
#dimx=img2.size[0]
#dimy=img2.size[1]
dimx=1200
dimy=1200
print(" Start for "+filename+" : " + str(datetime.datetime.now()))
addlog(" Start for "+filename+" : " + str(datetime.datetime.now()))
addlog("Image size: " + str(dimx)+","+str(dimy))

#check that the tiff file is indeed a 16 bits grey levels tiff file
if img.dtype==np.uint16:
    # r= Gaussian filter radius
    r=3
    h=GaussianFilter(r)
    #print("Gaussian filter matrix:")
    #print(h)
# raw data from original tiff image
    tr1=[]
# array of zones parts numbers
    subzones=[] 
    zoneMaxNb=1
# array of zones numbers
    zones=[]

#  mapping subzone=> zone: zonemap[subzone]=zone 
# initialization with a zero because 0 is a non valid value for a zone number, 
# so that in the subzone array a pixel with a zero value means that it is not associated 
# to a zone
    zonemap=[0]

#copying raw data of image in array tr1
    for i in range(img.size):
        val = float(img.flat[i])
        tr1.append(val)
        subzones.append(0)
        zones.append(0)

    trmin = min(tr1)
    trmax = max(tr1)
    k = 65535.0/(trmax-trmin)

###################### .bmp (bitmap) image generation (coloured version of the first tiff image)
    img2 = Image.new(mode="RGB", size=(dimx, dimy))

    for i in range(img.size):
        j=i%1200
        z=int(math.floor(k*(img.flat[i]-trmin)))

        R = ((z & 0b1<<15)>>15)*(0b1<<7)+((z & 0b1<<14)>>14)*(0b1<<6)+((z & 0b1<<13)>>13)*(0b1<<5)+((z & 0b1<<12)>>12)*(0b1<<4)+((z & 0b1<<11)>>11)*(0b1<<3)+((z & 0b1<<8)>>8)*(0b1<<2)+((z & 0b1<<5)>>5)*(0b1<<1)+((z & 0b1<<2)>>2)  
        G = ((z & 0b1<<15)>>15)*(0b1<<7)+((z & 0b1<<14)>>14)*(0b1<<6)+((z & 0b1<<13)>>13)*(0b1<<5)+((z & 0b1<<12)>>12)*(0b1<<4)+((z & 0b1<<10)>>10)*(0b1<<3)+((z & 0b1<<7)>>7)*(0b1<<2)+((z & 0b1<<4)>>4)*(0b1<<1)+((z & 0b1<<1)>>1)
        B = ((z & 0b1<<9)>>9)*(0b1<<3)+((z & 0b1<<6)>>6)*(0b1<<2)+((z & 0b1<<3)>>3)*(0b1<<1)+(z & 0b1)    

        img2.putpixel((j,int((i-j)/1200)),(R,G,B))

#Saving coloured image ( .bmp )
    img2.save(fld+filename_bmp)
    img1 = img2.copy()


# data smoothed by the gaussian filter
    tr=[]
# gradient data (Sobel algorithm)
    grd=[]


#Performing a gaussian filtering
    print("Applying the gaussian filter of radius="+str(r)+" : "+ str(datetime.datetime.now()))
    for i in range(img.size):
        if i >dimx*r and i < dimx*(dimy-r) and i%dimx > r-1 and i%dimx < dimx-r:
            val=0.0
            for i1 in range(-r,r+1):
                for j1 in range(-r,r+1):
                    val += h[j1+r][i1+r] * tr1[i+i1+dimx*j1]
        else:
            val=tr1[i]
        tr.append(val)

    mini = np.min(tr)
    maxi = np.max(tr)
    
    tr2=[0.0 for elt in range(120*120)]

    print("Computing averages on 10x10 pads")
    for i in range(img.size):
         col,line = OneToTwo(i)
         col2 = int(math.floor(col/10))
         line2 = int(math.floor(line/10))
         tr2[col2+120*line2] += img.flat[i]

    for i in range(len(tr2)):
        tr2[i] = tr2[i]/100.0
    mini2 = np.min(tr2)
    maxi2 = np.max(tr2)
    addlog("min / max of averaged pads : "+str(mini2)+" => "+str(maxi2))
    a,b = maxhist(tr2)
    addlog("Interval of maximum intensity for the overall background:")
    addlog(str(a)+" => "+str(b))

    d=b-mini2

    tr3=[]
    for i in range(len(tr2)):
        col2 = i%120
        line2 = int((i-col2)/120)
        if tr2[i] > b+3*d:
            tr3.append(tr2[i])
    a2,b2 = maxhist(tr3)
    cellBg = b2
    addlog("Cell background=" + str(cellBg))
    # for debugging, histogram of cell background
    #plt.plot(x3, tr3, color='red', linewidth=1)
    #plt.axvline(x=cellBg,color="yellow",linestyle="--")
    
    #basic gradient 
    #for i in range(img.size):
    #    if i >dimx-1 and i < dimx*(dimy-1) and i%dimx > 0 and i%dimx < dimx-1 and tr[i]> cellBg:
    #        grd.append( ((tr[i-1]-tr[i+1])**2 +(tr[i+dimx]-tr[i-dimx])**2)**(0.5) )
    #plt.hist(grd,1000)

    print("Applying Sobel filter " +str(datetime.datetime.now()))
    addlog("Applying Sobel filter " +str(datetime.datetime.now()))
    #Sobel:
    h2x=[[-0.25,0,0.25],[-0.5,0,0.5],[-0.25,0,0.25]]
    h2y=[[-0.25,-0.5,-0.25],[0.0,0,0.0],[0.25,0.5,0.25]]
    tr4=[]
    tr4b=[]
    tr4c=[]
    tr4s=[]
#applying a Sobel filter to detect the zones boundaries
    for i in range(img.size):
        sobelx=0.0
        sobely=0.0
        G=0.0
        theta=0.0
        cos=0.0
        sin=0.0
        if i >dimx and i < dimx*(dimy-1) and i%dimx > 0 and i%dimx < dimx-1:
            for i1 in range(-1,2):
                for j1 in range(-1,2):
                    sobelx+= h2x[j1+1][i1+1] * tr[i+i1+dimx*j1]
                    sobely+= h2y[j1+1][i1+1] * tr[i+i1+dimx*j1]
            gradient = sobelx+sobely*1j
            G = np.absolute(gradient)
            theta=np.angle(gradient)
            cos = math.cos(theta)
            sin = math.sin(theta)
        tr4.append(G)
        tr4c.append(cos)
        tr4s.append(sin)
    grdmean = np.mean(tr4)
    addlog("Gradient mean value:"+str(grdmean))
    
    print("End Sobel part 1 (gradient modulus and angle computation)")
    for j in range(2,dimy-2):
        for i in range(2,dimx-2):
#            if tr4[i+dimx*j]>2000:
            if tr4[i+dimx*j]>10*grdmean :
                cos = tr4c[i+dimx*j]
                sin = tr4s[i+dimx*j]
                g1 = interpolation(tr4,i+cos,j+sin)
                g2 = interpolation(tr4,i-cos,j-sin)
                if (tr4[i+dimx*j]>g1) and (tr4[i+dimx*j]>g2):
                    img2.putpixel((i,j),(0,255,0))
    print("End Sobel part 2 (maximas of gradients)")
#fin detection des bords de zones 


# the "Alpin" algorithm (I found that this name suits well for our peak detection algorithm... ;-) )
    n2=0
    m2=0
    for j in range(1,dimy-1):
        for i in range(1,dimx-1):
            px=float(tr[i+dimx*j])
            n2=0
            m2=0
#     we only consider points > 2x cell Banckground 
            if px > cellBg*2:
                for j2 in range(-1,2):
                    for i2 in range(-1,2):
                        if i2!=0 or j2!=0:
                            px2 = float(tr[i+i2+dimx*(j+j2)])
                            if px2 < px:
                                n2+=1
                            else:
                                if px2 < px*1.05 :
                                    m2+=1
            if n2+m2==8 and m2<4:
                img2.putpixel((i,j),(255,0,0))

# for debugging purposes => 3D display of the intensity of a small part of the image to study details
#    xsurf=[]
#    ysurf=[]
#    zsurf=[]
#    for j in range(1,dimy-1):
#        for i in range(1,dimx-1):
#            #if i>750 and i<800 and j>630 and j<650:
#            #if i>385 and i<397 and j>515 and j<535:
#            if i>250 and i<300 and j>350 and j<440:
#                px=float(tr[i+dimx*j])
#                xsurf.append(float(i))
#                ysurf.append(float(j))
#                zsurf.append(px)
#    fig = plt.figure()
#    ax = Axes3D(fig)
#    surf = ax.plot_trisurf(xsurf, ysurf, zsurf, cmap=cm.jet, linewidth=0.1, vmin=0, vmax=10000)
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    plt.show()
#    exit()

    print("Zones colouring "+str(datetime.datetime.now()))
# Zones colouring
    nbscanned=0
    for j in range(1,dimy-1):
        for i in range(1,dimx-1):
            pix=img2.getpixel((i,j))
            pxinit=0.0
            if pix==(255,0,0) or pix==(0,255,0): 
                zonemap.append(zoneMaxNb)
                subzones[i+dimx*j]=zoneMaxNb
                scanpx(i,j)
                px=float(tr[i+dimx*j])
                pxinit=px
                updown=1
                endscanup = 0
                nbfound=1
                
                i1=i
                j1=j
                while endscanup==0 and nbscanned< 100000000 : 
                    nbscanned+=1
                    j1+=updown
                    if testpix(img2,i1,j1,pxinit)==1:
                        i1+=1
                        if testpix(img2,i1,j1,pxinit)==1:
                            updown=-updown
                            j1+=updown
                            if testpix(img2,i1,j1,pxinit)==1:
                                j1+=updown
                                if testpix(img2,i1,j1,pxinit)==1:
                                    endscanup=1
                                else:
                                    nbfound+=1
                            else:
                                nbfound+=1
                        else:
                            nbfound+=1
                    else:
                        nbfound+=1

                px=float(tr[i+dimx*j])
                updown=-1
                endscanup = 0
                
                i1=i
                j1=j
                while endscanup==0 and nbscanned< 100000000 : 
                    nbscanned+=1
                    j1+=updown
                    if testpix(img2,i1,j1,pxinit)==1:
                        i1-=1
                        if testpix(img2,i1,j1,pxinit)==1:
                            updown=-updown
                            j1+=updown
                            if testpix(img2,i1,j1,pxinit)==1:
                                j1+=updown
                                if testpix(img2,i1,j1,pxinit)==1:
                                    endscanup=1
                                else:
                                    nbfound+=1
                            else:
                                nbfound+=1
                        else:
                            nbfound+=1
                    else:
                        nbfound+=1

                #print("subzone: " + str(zoneMaxNb) + " nb pixels="+str(nbfound)+" end time="+str(datetime.datetime.now()))
                if zoneMaxNb%50==0:
                    print("subzone: " + str(zoneMaxNb) + " ("+str(i)+","+str(j)+ ") nb pixels="+str(nbfound))
                zoneMaxNb+=1

#    print("Zone map:")                
#    print(zonemap[0:100])

# updating zone Nb for the current subzone and the other subzones we have touched, which acquire the same zone number as the current subzone:
    addlog("updating merged zones "+str(datetime.datetime.now()))
    for z in range(img.size):
        zo = subzones[z]
        if subzones[z] !=0 :
           zones[z]= zonemap[zo]  

# filling holes (5 iterations are performed)
    addlog("filling holes surrounded by selected pixels")
    for iter in range(4):
        npxAdded=0
        addlog("Iteration: "+str(iter+1)+" : "+str(datetime.datetime.now()))
        for j in range(1,dimy-1):
            for i in range(1,dimx-1):
                pix=img2.getpixel((i,j))
                if pixOK(pix)==0:
                    pE=pixOK(img2.getpixel((i+1,j)))
                    pO=pixOK(img2.getpixel((i-1,j)))
                    pN=pixOK(img2.getpixel((i,j-1)))
                    pS=pixOK(img2.getpixel((i,j+1)))
                    pNE=pixOK(img2.getpixel((i+1,j-1)))
                    pNO=pixOK(img2.getpixel((i-1,j-1)))
                    pSE=pixOK(img2.getpixel((i+1,j+1)))
                    pSO=pixOK(img2.getpixel((i-1,j+1)))
                    cross1 = pE+pO+pN+pS
                    cross2 = pNE+pNO+pSE+pSO
                    if cross1==4 or cross1+cross2>6 or (cross1==3 and cross2==3):
                        img2.putpixel((i,j),(255,255,0))
                        npxAdded+=1
                        #attributing zone number: taking the max zone number of surrounding 
                        #pixels, but normally they should all have the same zone number
                        neighbour_zones=[]
                        for i1 in range(-1,2):
                            for j1 in range(-1,2):
                                if (i1!=0 or j1 !=0) and subzones[i+i1+dimx*(j+j1)]>0:
                                    neighbour_zones.append(zonemap[subzones[i+i1+dimx*(j+j1)]])
                        if len(neighbour_zones)>0:
                            maxzone = max(neighbour_zones)
                            minzone = max(neighbour_zones)

                            if minzone != maxzone:
                                addlog("zone anomaly in "+str(i) + ","+str(j)+" minzone="+str(minzone)+" max="+str(maxzone))
                            subzones[i+dimx*j]= maxzone

                        
        addlog("Pixels added:"+str(npxAdded))

# merge of apparently isolated points 
    print("Merging isolated points "+str(datetime.datetime.now()))
    for j in range(1,dimy-1):
        for i in range(1,dimx-1):
            subzo=subzones[i+dimx*j]
            if subzo !=0 :
                center_zone=zonemap[subzo]
                neighbour_zones=[]
                same_zone_found=0
                for i1 in range(-1,2):
                    for j1 in range(-1,2):
                        if (i1!=0 or j1 !=0) and subzones[i+i1+dimx*(j+j1)]>0:
                            curzone=zonemap[subzones[i+i1+dimx*(j+j1)]]
                            if curzone==center_zone: 
                                same_zone_found=1
                            else:
                                if curzone!=0:
                                    neighbour_zones.append(curzone)

                if len(neighbour_zones)>0 and same_zone_found==0:
                    maxzone = max(neighbour_zones)
                    subzones[i+dimx*j]= maxzone
   
    addlog("finding zone squares "+str(datetime.datetime.now()))
    zoneinfo=[]
    for zo in range(len(zonemap)):
        #zoneinfo format:
        #zoneinfo[zone_number][0] = xmin of the rectangle
        #zoneinfo[zone_number][1] = xmax of the rectangle
        #zoneinfo[zone_number][2] = ymin of the rectangle
        #zoneinfo[zone_number][3] = ymax of the rectangle
        #zoneinfo[zone_number][4] = nb of pixels in this zone (not the area of the rectangle)
        #zoneinfo[zone_number][5] = total intensity in this zone (for the "*_C0.tiff" image)
        #zoneinfo[zone_number][6] = zone number as displayed in the final resulting image
        #zoneinfo[zone_number][7] = 1 if the zone is apparently a phase separation, 0 otherwise
        #zoneinfo[zone_number][8] = total intensity in this zone (for the "*_C1.tiff" image)
        #zoneinfo[zone_number][9] = average background around the zone (for the "*_C0.tiff" image)
        #zoneinfo[zone_number][10] = average background around the zone (for the "*_C1.tiff" image)
        zoneinfo.append([dimx,0,dimy,0,0,0.0,0,0,0.0,0.0,0.0])

    #computing zoneinfos data => rectangle, nb of concerned pixels, total intensity
    for z in range(img.size):
        zo = zones[z]
        if zo == zonemap[zo] and zo!=0:
            x,y = OneToTwo(z)
            if x < zoneinfo[zo][0]:
                zoneinfo[zo][0] = x
            if x > zoneinfo[zo][1]:
                zoneinfo[zo][1] = x
            if y < zoneinfo[zo][2]:
                zoneinfo[zo][2] = y
            if y > zoneinfo[zo][3]:
                zoneinfo[zo][3] = y
            zoneinfo[zo][4]+=1
            zoneinfo[zo][5]+=float(tr1[z])

# merging small zones
    addlog("Merging zones < 15 pixels "+str(datetime.datetime.now()))
    for zo in range(len(zonemap)):
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][4]>0 and zoneinfo[zo][4]<15:
            min_nbpixels_outerzone=999999
            dist_to_center=float(dimx)
            matching_zone=0
            for zi in range(len(zonemap)):
                if zi == zonemap[zi] and zi!=0 and zoneinfo[zi][4]>0 and zoneinfo[zi][0]<=zoneinfo[zo][0] and zoneinfo[zi][1] >= zoneinfo[zo][1] and zoneinfo[zi][2]<=zoneinfo[zo][2] and zoneinfo[zi][3]>=zoneinfo[zo][3]:
                    if zoneinfo[zi][4] < min_nbpixels_outerzone:
                        min_nbpixels_outerzone=zoneinfo[zi][4]
                        matching_zone=zi
            if matching_zone !=0:
                dist_to_center=(( (zoneinfo[matching_zone][0]+zoneinfo[matching_zone][1])/2.0 - (zoneinfo[zo][0]+zoneinfo[zo][1])/2.0)**2+( (zoneinfo[matching_zone][2]+zoneinfo[matching_zone][3])/2.0 - (zoneinfo[zo][2]+zoneinfo[zo][3])/2.0)**2)**(0.5)
                if dist_to_center < ((zoneinfo[matching_zone][0]-zoneinfo[matching_zone][1])**2+(zoneinfo[matching_zone][2]-zoneinfo[matching_zone][3])**2)**(0.5)/3.0:
                    #print("merging zone "+str(zo)+" into zone "+str(matching_zone)+" pixels="+str(zoneinfo[zo][4]))
                    zonemap[zo]=matching_zone
                    zoneinfo[matching_zone][4]+=zoneinfo[zo][4]
                    zoneinfo[matching_zone][5]+=zoneinfo[zo][5]
                    zoneinfo[zo][4]=0
                    #udpating zones array 
                    for pxu in range(img.size):
                        if subzones[pxu]==zo:
                            zones[pxu] = matching_zone


# Identification of zones separation (criterias are: width <110 pixels, and the ratio length/width must be close to 1)
    addlog("Indentifying phase separations "+str(datetime.datetime.now()))
    for zo in range(len(zonemap)):
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][4]>0 and zoneinfo[zo][4]<110:
            zo_w = float(abs(zoneinfo[zo][1]-zoneinfo[zo][0]))
            zo_h = float(abs(zoneinfo[zo][3]-zoneinfo[zo][2]))
            if zo_h*0.7 < zo_w < zo_h*1.2:
                zoneinfo[zo][7]=1
                zoneinfo[zo][9]=average2(tr,zo)


    ##################################################################################################################
    #computing corresponding intensities in the 2nd channel image (needs to recompute its own cell background level)
    ##################################################################################################################
    img3 = tifffile.imread(fld+filename.replace("_C0.tiff","_C1.tiff"))
#copying raw data of 2nd image in array tr5
    tr5=[]
    for i in range(img3.size):
        val = float(img3.flat[i])
        tr5.append(val)

    tr6=[0.0 for elt in range(120*120)]

    print("Computing averages on 10x10 pads of 2nd image"+str(datetime.datetime.now()))
    for i in range(img3.size):
         col,line = OneToTwo(i)
         col2 = int(math.floor(col/10))
         line2 = int(math.floor(line/10))
         tr6[col2+120*line2] += img3.flat[i]

    for i in range(len(tr6)):
        tr6[i] = tr6[i]/100.0
    mini2_C1 = np.min(tr6)
    maxi2_C1 = np.max(tr6)
    addlog("min / max of averaged pads on 2nd image: "+str(mini2_C1)+" => "+str(maxi2_C1))
    a_C1,b_C1 = maxhist(tr6)
    addlog("Interval of maximum for the overall background of 2nd image:")
    addlog(str(a_C1)+" => "+str(b_C1))

    d_C1=b_C1-mini2_C1

    tr7=[]
    for i in range(len(tr6)):
        col2 = i%120
        line2 = int((i-col2)/120)
        if tr6[i] > b_C1+3*d_C1:
            tr7.append(tr6[i])
    a2_C1,b2_C1 = maxhist(tr7)
    cellBg2 = b2_C1
    addlog("Cell background (2nd image)=" + str(cellBg2))

    print("Collecting phase separation intensities on 2nd image : "+str(datetime.datetime.now()))
    for z in range(img3.size):
        zo = zones[z]
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][7]==1:
            zoneinfo[zo][8]+=float(tr5[z])

    for zo in range(len(zonemap)):
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][4]>0 and zoneinfo[zo][7]==1:
            zoneinfo[zo][10]=average2(tr5,zo)


    ##################################################################################################################
    #drawing squares surrounding zones and write the zones infos in a text file
    ##################################################################################################################
    zoneNumber=0
    print("Printing zone squares and result data in a text file : "+str(datetime.datetime.now()))

    for zo in range(len(zonemap)):
        #addlog(str(zo)+","+str(zonemap[zo]))
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][4]>14:
            zoneNumber+=1
            zoneinfo[zo][6]=zoneNumber
            with open(fld+filename_txt,"a") as text_file:
                text_file.write(hex(zoneinfo[zo][6])[2:]+","+str(zoneinfo[zo][0])+","+str(zoneinfo[zo][1])+","+str(zoneinfo[zo][2])+","+str(zoneinfo[zo][3])+","+str(zoneinfo[zo][4])+","+str(zoneinfo[zo][5])+","+str(zoneinfo[zo][8])+","+str(zoneinfo[zo][9])+","+str(zoneinfo[zo][10])+","+str(zoneinfo[zo][7])+","+str(zo)+"\n")
               
            zoneSquare_indentification(zoneinfo[zo][0]-1,zoneinfo[zo][1]+1,zoneinfo[zo][2]-1,zoneinfo[zo][3]+1,zoneNumber)

    # writing rsult data in a .txt file 
    for zo in range(len(zonemap)):
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][4]>0 and zoneinfo[zo][4]<15:
            zoneNumber+=1
            zoneinfo[zo][6]=zoneNumber
            with open(fld+filename_txt,"a") as text_file:
                text_file.write(hex(zoneinfo[zo][6])[2:]+","+str(zoneinfo[zo][0])+","+str(zoneinfo[zo][1])+","+str(zoneinfo[zo][2])+","+str(zoneinfo[zo][3])+","+str(zoneinfo[zo][4])+","+str(zoneinfo[zo][5])+","+str(zoneinfo[zo][8])+","+str(zoneinfo[zo][9])+","+str(zoneinfo[zo][10])+","+str(zoneinfo[zo][7])+","+str(zo)+"\n")
               

    #Saving result images:
    #img2 = coloured image with yellow/green/red pixels
    #img1 = coloured image with zones surrounded by a green rectangle and a zone number in red expressed in hexadecimal format to save space on the image

    #img2.save(fld+filename_bmp.replace(".bmp","_interm.bmp")) 
    img1.save(fld+filename_bmp.replace(".bmp","_result.bmp"))   


    ##################################################################################################################
    #scatterplot
    ##################################################################################################################
     
    scx=[]
    scy=[]
    #  b= overall background on 1st image, b_C1 = idem for the 2nd image
    for zo in range(len(zonemap)):
        if zo == zonemap[zo] and zo!=0 and zoneinfo[zo][7]==1:
            scx.append((zoneinfo[zo][5]-zoneinfo[zo][4]*zoneinfo[zo][9])/b)
            scy.append((zoneinfo[zo][8]-zoneinfo[zo][4]*zoneinfo[zo][10])/b_C1)
    plt.scatter(scx,scy)
    plt.show()


    ##################################################################################################################
    #end
    ##################################################################################################################
     
    print(" End computing at " + str(datetime.datetime.now()))
    addlog(" End computing at " + str(datetime.datetime.now()))
    # display final image on screen:
    #img1.show()

