#coding:utf-8

import numpy as np

class maxSimRegionMerging:
    """docstring for maxSimRegionMerging"""
    def __init__(self, figHandle,orgImage,bwSegmentResultImage,bwMarkerImage,Label):
        super(maxSimRegionMerging, self).__init__()
        self.figHandle                 # the handle of the interface
        self.orgImage                  # the original image
        self.bwSegmentResultImage
        self.bwMarkerImage               # the binary image with object and background marker
        self.Label                     # 
    
    def maxSimRegionMerging(self):

        Label2 = Label
        #hdl = get(figHandle,'userData')    # 
        #set(figHandle,'currentAxes',hdl.figHandle_Axes2)               # set hdl.figHandle_Axes2 as the current axes
        mergeLabel = Label                                               # the label of region merging
        mergingTimes = 1                                                 # the number of the region merging
    
        regionNum = Label.max()                                      # the number of the region
        redbins = 16
        greenbins = 16
        bluebins = 16                            # the quantization of image
        indImage = self.quanimage(orgImage,redbins,greenbins,bluebins)       
    
        [height,width] = Label.shape
        Region = self.InitRegion(Label,indImage,bwMarkerImage,height,width)     # initilize the region. Region is a struct variable
        SimTable = self.CompSim(Label,height,width,Region)                    # SimTable is a table which saves the similarity between the regions
    
        MergeTable = -1*ones(1,regionNum)                                # initilize the MergeTable
    
        k = 0                       # k is the times of region merging. the initial value is 0.
        while True:
        k = k + 1              
            flag = 0                # 0£ºno region merging  1:region merging
            k1 = 0                  # k1 is the times of region merging in the first stage
            while True:               # region merging in the first stage        
                k1 = k1 + 1           # 
                flag1 = 0           # 0£ºno region merging  1:region merging
                for i in range(regionNum):  #         
                    # markerType   0 non-marker region  1 background  2 foreground
                    # background marker regions merge with non-marker regions as
                    # possible in the first stage
                    if Region[i].markerType==1:                        # region i is background
                        # 
                        for j in range(regionNum):                             # 
                            if i!=j and SimTable[i,j]>0 and Region[j].markerType!=2:        # region i and region j is adjacent and region i belongs to non-marker region
                                index = self.MaxSimIndex(j,SimTable,regionNum)               # 
                                if i==index:                                            # 
                                    if i<j:                         
                                        MergeTable = self.MergeRecord(MergeTable,i,j,regionNum)                                
                                        flag1 = 1                                       # 
                                        flag = 1                                        # 
                                    elif i>j:
                                        MergeTable = self.MergeRecord(MergeTable,j,i,regionNum)
                                        flag = 1                                        # 
                                        flag1 = 1                                       # 
                if flag1==0:
                    break
                Label2,regionNum,Region = self.MergePostProc2(Label2,MergeTable,regionNum,Region)
                SimTable = CompSim(Label2,height,width,Region)
    
               # 
                MergeTable=-1*ones(1,regionNum)
                ImageE=drawEdge(orgImage,Label2)
    
                #set(hdl.segImage, 'XData',[1 width],'YData',[1 height],'Cdata',ImageE)drawnow
    
                if k==1:
                    str2='st round'
                elif k==2:
                    str2='nd round'        
                #title([num2str(k1),'th ','merging in 1th stage ','( ',num2str(k),str2,')'])drawnow
                #pause(1)
            # -------------------------- region merging of the first stage merges ends  --------------------------
    
            if flag==0: 
                break
    
            # -------------------------- region merging of the second stage merge starts --------------------------
            flag2 = 1
            k1 = 0
            while True:   #   non-marker regions are merged each other in the second stage
                k1 = k1 + 1 
                flag2 = 0
    
                for i in range(regionNum):
                    if Region[i].markerType==0:     # for non-marker regions
                        # 
                        for j in range(regionNum):
                            # region i is not marked and region i and region is is adjacent
                            if i!=j and Region(j).markerType==0 and SimTable(i,j)>0:
                                index = self.MaxSimIndex(j,SimTable,regionNum)
                                if i==index:
                                    if i<j:
                                        MergeTable = self.MergeRecord(MergeTable,i,j,regionNum)
                                        flag = 1          
                                        flag2 = 1         
                                    else:
                                        MergeTable = self.MergeRecord(MergeTable,j,i,regionNum)
                                        flag = 1
                                        flag2 = 1
    
                if flag2==0
                    break
    
               #
                Label2,regionNum,Region = self.MergePostProc2(Label2,MergeTable,regionNum,Region)
                SimTable = CompSim(Label2,height,width,Region)
    
                MergeTable = -1*ones([1,regionNum])
    
                ImageE = drawEdge(orgImage,Label2)
                #set(hdl.segImage, 'XData',[1 width],'YData',[1 height],'Cdata',ImageE)drawnow
    
                if k==1
                    str2='st round'
                elif k==2
                    str2='nd round'                
                #title([num2str(k1),'th ','merging in 2th stage ','(',num2str(k),str2,')'])drawnow
                #pause(1)
            if flag==0:
                break
    
        # extract the object
        bwSegmentResultImage = 0.0*Label
        for i in range(regionNum):
            if Region[i].markerType!=1:
                bwSegmentResultImage(find(Label2==i))=1
    
        ImageE=drawEdge(orgImage,bwSegmentResultImage)
        #set(hdl.segImage, 'XData',[1 width],'YData',[1 height],'Cdata',ImageE)drawnow
        #title(['Segmentation by MSRM'])
        return D,ImageE

    def CompSim(self,I,height,width,Region):

        RegionNum = I.max()
        SimTable = np.zeros([RegionNum,RegionNum],int)
    
        for i in range(height):
            for j in range(width):
                index = I[i,j]
                SimTable[index,index] = 1
                for k in range(-1,2,1)
                    for l in range(-1,2,1)
                        row = i + k
                        col = j + l
                        if row>0 and row<=height and col>0 and col<=width:
                            if index != I[row,col]:
                                SimTable[index,I[row,col]] = 1
    
        for i in range(RegionNum):
            for j in range(RegionNum):
                if i!=j and SimTable[i,j]>0:
                    H1 = sqrt(Region[i].rgbHistogram/Region[i].area)
                    H2 = sqrt(Region[j].rgbHistogram/Region[j].area)           
                    SimTable[i,j] = H1*H2.conj().T + 0.000001
        return SimTable

    def drawEdge(self,RGB_I,L):
    #draw the edge between two different region
    
        H = L.shape[0]
        W = L.shape[1]
    
        DirX = [1,1,0,-1,-1,-1, 0, 1]
        DirY = [0,1,1, 1, 0,-1,-1,-1]
    
        ImageE = RGB_I
    
        for i in range(H):
            for j in range(W):
                for k in range(8):
                    y = i + DirY[k]
                    x = j + DirX[k]
                    if x>=1 and x<=W and y>=1 and y<=H:
                        if L[i,j] < L[y,x]:        
                            ImageE[i,j,:] = 255
        return ImageE

    def fillAllZero(self,I,height,width):
    #fill the zeor in the image I and ready to region merging
    
        dirX = [1,1,0,-1,-1,-1, 0, 1]
        dirY = [0,1,1, 1, 0,-1,-1,-1]
    
        I2 = I
        for i in range(width):
            for j in range(height):
                Flag = 0
                if I[j,i] == 0:
                    for k in range(8):
                        curX = i + dirX[k]
                        curY = j + dirY[k]
                        if curX>=1 and curX<=width and curY>=1 and curY<=height:
                            if I[curY,curX] != 0:
                                I2[j,i] = I[curY,curX]
                                flag = 1
                                break
        return I2

    def InitRegion(self,I,indImage,Mask,height,width):
    # initialize region information and calcluate the area, markerType and
    # rgbHistogram
    
        RegionNum = max(for i in I) 
        for i in range(RegionNum):             
            Region[i].area = 0
            Region[i].markerType = 0    
            Region[i].rgbHistogram = np.zeros(1,max(for i in indImage),int)
    
        for i in range(height):
            for j in range(width):
                index=I[i,j]        
                Region[index].area=Region[index].area+1 
                Region[index].rgbHistogram[1,indImage[i,j]] = Region[index].rgbHistogram[1,indImage[i,j]] + 1
                Region[index].markerType=max(Region[index].markerType,Mask[i,j])
        return Region

    def MaxSimIndex(self,index,SimTable,regionNum):
    
        maxSim = -1
        for i in range(regionNum):
            if i != index:
                sim = SimTable[i,index]
                if maxSim<sim:
                    maxSim = sim
                    maxIndex = i
        return maxIndex

    def MergePostProc2(self,I,record,N,Region):
    
        I2 = np.zeros(I.shape[0],I.shape[1],int)   
        k = 0
        for i in range(N):
            # 1. 
            if record[i] == -1:
                k = k + 1
                ind = (I==i)         
                I2[ind] = k               
            
                Region2[k].rgbHistogram=Region[i].rgbHistogram
                Region2[k].markerType=Region[i].markerType
                Region2[k].area=Region[i].area
    
            elif record[i] == -2:
                k = k + 1
                ind = (I==i)        
                I2[ind] = k
            
                Region2[k].rgbHistogram=Region[i].rgbHistogram
                Region2[k].markerType=Region[i].markerType
                Region2[k].area=Region[i].area
            
                for j in range(N):              
                    if record[j] == i:     
                        ind = (I==j)
                        I2[ind] = k
    
                        Region2[k].rgbHistogram = Region2[k].rgbHistogram+Region[j].rgbHistogram
                        Region2[k].markerType = max(Region2[k].markerType,Region[j].markerType)      %
                        Region2[k].area = Region2[k].area+Region[j].area               
                
        newRegionNum = k
        return I2,newRegionNum,Region2

    def MergeRecord(self,record,index1,index2,regionNum):
    # merging index1 and index2 by updating record
    
        newRecord = record                 
        # 1.
        if record[index1]==-1 and record[index2]==-1:
            newRecord[index1] = -2         
            newRecord[index2] = index1    
        # 2.
        elif record[index1]==-1 and record[index2]==-2:
            newRecord[index1] = -2         
            newRecord[index2] = index1     
            for i in range(regionNum):           
                if record[i]==index2:
                    newRecord[i] = index1  
        # 3.
        elif record[index1]==-1 and record[index2]>0:
            rootNode = record[index2]      
            if index1>rootNode:            
                newRecord[index1] = rootNode
            elif index1<rootNode:
                newRecord[index1] = -2     
                newRecord[rootNode] = index1 
                for i in range(regionNum):
                    if record[i]==rootNode:
                        newRecord[i] = index1
        # 4.
        elif record[index1]==-2 and record[index2]==-1:
            newRecord[index2] = index1          
        # 5.
        elif record[index1]==-2 and record[index2]==-2:
            newRecord[index2] = index1          
            for i in range(regionNum):
                if record[i]==index2:
                    newRecord[i] = index1
        # 6.
        elif record[index1]==-2 and record[index2]>0:
            rootNode = record[index2]
            if index1<rootNode:
                newRecord[rootNode] = index1    
                for i in range(regionNum):
                    if record[i]==rootNode:    
                        newRecord[i] = index1
            elif index1>rootNode:       
                newRecord[index1] = rootNode    
                for i in range(regionNum):            
                    if record[i]==index1:      
                        newRecord[i] = rootNode
        #7. 
        elif record[index1]>0 and record[index2]==-2:
            rootNode = record[index1]           
            newRecord[index2] = rootNode        
            for i in range(regionNum):
                if record[i]==index2:           
                    newRecord[i] = rootNode     
        #8. 
        elif record[index1]>0 and record[index2]==-1:
            rootNode = record[index1]           
            newRecord[index2] = rootNode        
        #9. 
        elif record[index1]>0 and record[index2]>0:
            rootNode1 = record[index1]
            rootNode2 = record[index2]
            if rootNode1>rootNode2:             
                newRecord[rootNode1] = rootNode2
                for i in range(regionNum):       
                    if record[i]==rootNode1    
                        newRecord[i] = rootNode2
            elif rootNode1<rootNode2:
                newRecord[rootNode2] = rootNode1
                for i in range(regionNum):    
                    if record[i]==rootNode2:    
                        newRecord[i] = rootNode1
        return newRecord

    def quanimage(self,image,redBins,greenBins,blueBins):
    
        image = double(image)
    
        M = image.shape[0]
        N = image.shape[1]
        indImage = np.zeros([M,N],int)
    
        red = 256/redBins        
        green = 256/greenBins
        blue = 256/blueBins
    
        for y in range(M):
            for x in range(N):
                i = math.floor(image[y,x,1]/red) + 1
                j = math.floor(image[y,x,2]/green) + 1
                k = math.floor(image[y,x,3]/blue) + 1        
                index = red*green*(i-1) + green*(j-1) + k 
                indImage[y,x] = index
        return indImage

    def bwlabels(self,img,temp,x,y,n):#种子点(x,y)，递归找出与该点为同一目标的点
        w=len(img[0])
        h=len(img)
        
        num=0
        pt=[]
        #图像的九种情况
        if (x!=0 and x!=w-1) and (y!=0 and y!=h-1):
            for i in range(3):
                for j in range(3):          
                    if img[y-1+i][x-1+j]==1 and temp[y-1+i][x-1+j]==0:
                        num+=1
                        temp[y-1+i][x-1+j]=n
                        t=[y-1+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==0 and (y!=0 and y!=h-1):#第一列
            for i in range(3):
                for j in range(2):          
                    if img[y-1+i][x+j]==1 and temp[y-1+i][x+j]==0:
                        num+=1
                        temp[y-1+i][x+j]=n
                        t=[y-1+i,x+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==0 and y==0:#第一个点
            for i in range(2):
                for j in range(2):          
                    if img[y+i][x+j]==1 and temp[y+i][x+j]==0:
                        num+=1
                        temp[y+i][x+j]=n
                        t=[y+i,x+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif (x!=0 and x!=w-1) and y==0:#第一行
            for i in range(2):
                for j in range(3):          
                    if img[y+i][x-1+j]==1 and temp[y+i][x-1+j]==0:
                        num+=1
                        temp[y+i][x-1+j]=n
                        t=[y+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==w-1 and y==0:#第一行最后一个点
            for i in range(2):
                for j in range(2):          
                    if img[y+i][x-1+j]==1 and temp[y+i][x-1+j]==0:
                        num+=1
                        temp[y+i][x-1+j]=n
                        t=[y+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==w-1 and (y!=h-1 and y!=0):#最后一列
            for i in range(3):
                for j in range(2):          
                    if img[y-1+i][x-1+j]==1 and temp[y-1+i][x-1+j]==0:
                        num+=1
                        temp[y-1+i][x-1+j]=n
                        t=[y-1+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==w-1 and y==h-1:#最后一点
            for i in range(2):
                for j in range(2):          
                    if img[y-1+i][x-1+j]==1 and temp[y-1+i][x-1+j]==0:
                        num+=1
                        temp[y-1+i][x-1+j]=n
                        t=[y-1+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif (x!=0 and x!=w-1) and y==h-1:#最后一行
            for i in range(2):
                for j in range(3):          
                    if img[y-1+i][x-1+j]==1 and temp[y-1+i][x-1+j]==0:
                        num+=1
                        temp[y-1+i][x-1+j]=n
                        t=[y-1+i,x-1+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)
        elif x==0 and y==h-1:#最后一行第一个点
            for i in range(2):
                for j in range(2):          
                    if img[y-1+i][x+j]==1 and temp[y-1+i][x+j]==0:
                        num+=1
                        temp[y-1+i][x+j]=n
                        t=[y-1+i,x+j]
                        pt.append(t)
            if num==0:
                return temp
            else:
                for k in range(num):
                    self.bwlabels(img,pt[k][1],pt[k][0],n)

    def mybwLabels(self,img):#对图像中的点进行标记
        w=len(img[0])
        h=len(img)
        temp=np.zeros(a.shape[1])
        temp=np.tile(temp,[a.shape[0],1]) 
        num=0
        for i in range(h):
            for j in range(w):
                if img[i][j]==1 and temp[i][j]==0:#未标记的非0点
                    num += 1
                    temp = self.bwlabels(img,temp,j,i,num)
        return temp

    def setLabel(self,segImage):
        # number each region
        #imBW = im2bw(segImage)
        ret,imBW = cv2.threshold(segImage,127,255,cv2.THRESH_BINARY)  
        #Label = bwlabel(imBW) 
        Label = self.mybwLabels(imBW)
        print(temp)                                 
        [H,W] = Label.shape
        while True:                                       
            Label = self.fillAllZero(Label,H,W)
            if (Label==0) is None:
                break
        return Label

    def setMarkerImage(self,Image,imMask):
    # Image: original image
    # imMask: Mask
    
        [indy,indx] = nonzero(imMask)              
        N = indy.shape[0]                      
    
        for i in range(N):
            row = indy[i]
            col = indx[i]
            if imMask[row,col]==1: 
                Image[row,col,:] = 0
                Image[row,col,3] = 255
            elif imMask[row,col]==2:
                Image[row,col,:] = 0
                Image[row,col,2] = 255
        return Image

    def setMaskImage(self,markerImage,imMask,isObjectMarker,isBackgroundMarker):
        index = nonzero(imMask==1)
        if isObjectMarker==1:
            markerImage[index] = 2         
        elif isBackgroundMarker==1:
            markerImage[index] = 1         
        return markerImage