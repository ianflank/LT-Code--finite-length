import math
import numpy as np
from release import release

k=int(input("k 為 encoding distribution degree 最大值:"))
u=k
n=int(input("收到封包數:"))


#建立ideal soliton distribution
distribution=np.zeros(k+1)
distribution[1]=np.float32(1/k)
for i in range(2,k+1):
    distribution[i]=np.float32(1/(i*(i-1)))


#建立num_arrays個數的子矩陣，為檢驗某一state其上一層所有可能轉移state
num_arrays =n*((2*n)*(k+1))+(n)*(k+1)+k
arrays = [np.zeros([n+1,n+1,u+1], dtype=np.float32) for _ in range(num_arrays+1)]

#建立此一矩陣，為顯示經計算後各state發生機率
probability=np.zeros([n+1,n+1,u+1], dtype=np.float32)
for Q in range(0,n+1):
    probability[n-Q][Q][u]=math.pow(1/k,Q)*math.pow(1-1/k,n-Q)*math.comb(n,Q)

#在probability矩陣中標示初始state其下一層的所有可能state的值標記為0.5
def under_layer(c,r,u):
    s=i=0
    r1=r-1
    c1=c
    while(c1 >=0):
        while((s <=2) & (r1 >=0)):
            probability[c1][r1][u-1]=0.5
            while(i==0):
                r2=r1
                i=i+1
            s=s+1
            r1=r1-1
        i=s=0
        c1=c1-1
        r1=r2+1
    print("end layer")
    print(probability)

#為計算state的發生機率，並將其上一層所有可能轉變state標記
#x=j
#y=i
def value_probability(c1,r1,u1):
    h1=0
    if(u1>=1):
        for x in range (0,n-c1+1):
            y=0
            while(r1+1+y-x<=n):
                if((r1+y-x>=0)&(r1-x>=0)):
                    h=probability[c1+x][r1+1+y-x][u1+1]*state_transfer(x,y,c1,r1,u1+1)
                    arrays[c1*((2*n)*(u+1))+r1*(u+1)+u1][c1+x,r1+1+y-x,u1+1]=100
                else:
                    h=0
                h1=h1+h
                y=y+1
  #print(f"({c1},{r1},{u1}),{c1*((2*n)*(u+1))+r1*(u+1)+u1}")
        probability[c1][r1][u1]=h1



#為計算轉移機率
def state_transfer (j,i,c1,r1,u1):
    if u1>1:
        value=np.float32(math.comb(c1+j,j)* math.pow(Pu,j)* math.pow(1-Pu,c1)*math.comb(r1+i-j,i)* math.pow(1/u1,i)* math.pow(1-1/u1,r1-j))
        return value
    else:
        value=1
        return value




#將所有state利用迴圈計算其發生機率
x=n
u2=u
L=1
D=p=0
while(u2>2):
    Pu=release (distribution,k,u2)
    for j in range (0,n+1):
        for i in range (0,n+1):
            if(j+i<=x):
                value_probability(j,i,u2-1)
    x=x-1
    u2=u2-1
    print(u2)

while( u2==2):
   Pu=release (distribution,k,u2)
   for i in range (0,n+1):
        print("loading...")
        value_probability(0,i,u2-1)
   u2=u2-1

if(u2==1):
   while(L<=n):
        p=p+probability[2][L][2]
        D=D+probability[0][L][1]
        L=L+1

   print(f"D為{D}為解碼成功率")

#將值以以四捨五入方式取到小數點第4位
rounded_array = np.round(probability, decimals=4)
print(rounded_array)


output_file_path = "C:/Users/xxxx.txt"




#将三维矩阵写入文本文件
with open(output_file_path, 'w') as file:
    for matrix_2d in rounded_array:
        np.savetxt(file, matrix_2d, fmt='%f', delimiter='\t')
        file.write('\n')  # 在每个2D矩阵之间添加空行

while True:
      O=int(input("1.問機率  /  2.所有轉移情形:"))
      if(O==1):
            A,B,C=map(int,input("欲知state的probability:(c,r,u)").split())
            print(probability[A][B][C])
      else:
            c1,r1,u1=map(int,input("確認狀態的u+1的所有轉移可能性").split())
            print(arrays[c1*((2*n)*(u+1))+r1*(u+1)+u1])

