# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 15:42:48 2020

@author: 87198
"""

import csv
import numpy as np
import copy
import sys


epsl0=8.8542e-12
epslRPdms=2.8
structureName='chapter5-otft-'
cDict={}
for n in range(0,21):
    fileName=structureName+'%02d.csv'%n
    #设置求解坐标范围
    stepTime=float(n)/20
    x=300
    y=300
    z=300-200*stepTime
    C=0
    
    print '**********'+fileName+'**********'
    print 'Step time=%.3f' % stepTime
    
    #根据csv文件中每个单元对应节点的evf void出现次数最多的值作为该单元的evf void，然后初始化epsilon_r
    #xNum yNum zNum为各方向上单元数
    meshSize=5
    xNum=int(x/meshSize)
    yNum=int(y/meshSize)
    zNum=int(z/meshSize)
    epslR=np.ones((xNum+2,yNum+2,zNum+2)) #四周预留以供设置边界条件
    
    with open(fileName,'r')as csvFile:
        next(csvFile)
        reader=csv.reader(csvFile)
        for row in reader:
            if float(row[6])<x and float(row[7])<y and float(row[8])<z and float(row[6])>0 and float(row[7])>0 and float(row[8])>0: #判断单元是否在电场求解范围内
                epslR[int((float(row[6])+meshSize/2)/meshSize),int((float(row[7])+meshSize/2)/meshSize),int((float(row[8])+meshSize/2)/meshSize)]=1+(epslRPdms-1)*(1-float(row[12])) #根据EVF_VOID赋给单元epsilon值
    
    """
    有限差分法求电势分布，进而由电场能求电容
    """
    #在z方向上均匀初始化电势场
    phi=np.zeros_like(epslR)
    for n in range(1,zNum+1):
        phi[:,:,n].fill(float(n)/zNum)
    #设置第一类边界条件
    phi[:,:,0]=0
    phi[:,:,zNum+1]=1
    epslR[:,:,0]=epslR[:,:,1]
    epslR[:,:,zNum+1]=epslR[:,:,zNum]
    #初始化误差和迭代数
    error=1
    iteration=0
    #迭代求解
    print "Iteration bigins!"
    while iteration<1000 and error>1e-4: #迭代直至1000代或误差小于1e-4
        iteration+=1
        phiLast=copy.deepcopy(phi)
        #遍历每个单元，更新其电势值
        for i in range(1,xNum+1):
            for j in range(1,yNum+1):
                for k in range(1,zNum+1):
                    phi[i,j,k]=(phi[i+1,j,k]+phi[i-1,j,k]+phi[i,j+1,k]+phi[i,j-1,k]+phi[i,j,k+1]+phi[i,j,k-1]+((epslR[i+1,j,k]-epslR[i-1,j,k])*(phi[i+1,j,k]-phi[i-1,j,k])+(epslR[i,j+1,k]-epslR[i,j-1,k])*(phi[i,j+1,k]-phi[i,j-1,k])+(epslR[i,j,k+1]-epslR[i,j,k-1])*(phi[i,j,k+1]-phi[i,j,k-1]))/(4*epslR[i,j,k]))/6
        #更新第二类边界条件
        phi[0,:,:]=phi[2,:,:]
        phi[xNum+1,:,:]=phi[xNum-1,:,:]
        phi[:,0,:]=phi[:,2,:]
        phi[:,yNum+1,:]=phi[:,yNum-1,:]
        #求当前误差
        error=abs(phi-phiLast).max()
        sys.stdout.write('\r'+'iteration = %3d' % iteration + '   ' +'error = %.5f' % error)
        sys.stdout.flush()
        
    print ''
    print "Iteration finished!"
    
    W=0 #初始化电场能
    #对电场能遍历求和
    for x in range(1,xNum+1):
        for y in range(1,yNum+1):
            for z in range(1,zNum+1):
                W+=(epsl0*epslR[x,y,z]/2)*(((phi[x+1,y,z]-phi[x-1,y,z])/2)**2+((phi[x,y+1,z]-phi[x,y-1,z])/2)**2+((phi[x,y,z+1]-phi[x,y,z-1])/2)**2)*(meshSize*1e-6)
    C=2*W #求解电容
    
    """
    串并联方法求电容
    
    cXY=np.zeros((xNum,yNum))
    for xx in range(1,xNum+1):
        for yy in range(1,yNum+1):
            cXY[xx-1,yy-1]=1/sum(1/(epsl0*epslR[xx,yy,1:-1]*meshSize*1e-6)) #z方向上串联
    print zNum
    C=sum(sum(cXY)) #xy方向上并联
    """
    
    cDict[stepTime]=C*10**12
    print "C = %.8f pF" %(C*10**12)