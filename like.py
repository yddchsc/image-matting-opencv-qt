import cv2
 
def createHist(img): 
    #cv2.CvtColor(img,img,cv2.CV_BGR2HSV) 
    b_plane = cv2.CreateImage((img.width,img.height), 8, 1)
    g_plane = cv2.CreateImage((img.width,img.height), 8, 1)
    r_plane = cv2.CreateImage((img.width,img.height), 8, 1)
 
    cv2.split(img,b_plane,g_plane,r_plane,None)
    planes = [b_plane, g_plane, r_plane]
     
    bins = 4
    b_bins = bins
    g_bins = bins
    r_bins = bins
 
    hist_size = [b_bins,g_bins,r_bins]
    b_range = [0,255]
    g_range = [0,255]
    r_range = [0,255]
 
    ranges = [b_range,g_range,r_range]
    hist = cv2.CreateHist(hist_size, cv2.CV_HIST_ARRAY, ranges, 1)
    cv2.CalcHist([cv2.GetImage(i) for i in planes], hist)
    cv2.NormalizeHist(hist,1)
    return hist
 
def imgcompare(image1,image2):
    img1 = cv2.LoadImage(image1)
    hist1 = createHist(img1)
    img2 = cv2.LoadImage(image2)
    hist2 = createHist(img2)
    return cv2.CompareHist(hist1,hist2,cv2.CV_COMP_CORREL)
     
print(imgcompare("2.jpg","0010.jpg"))
