import math
import numpy as np
from release import release

k=int(input("k 為 encoding distribution degree 最大值:"))
u=k



#c,r=map(int,input("c、r:").split())
# print("c=",c)
# print("r=",r)

# # 建立ideal soliton distribution
# distribution=np.zeros(k+1)
# distribution[1]=float(1/k)
# for i in range(2,k+1):
#  distribution[i]=float(1/(i*(i-1)))


# 建立Robust soliton distribution
phi=np.zeros(k+1)
phi[1]=float(1/k)
for i in range(2,k+1):
 phi[i]=float(1/(i*(i-1)))


delta=float(input("delta 為可接受失敗率:"))
G=float(input("G 為何"))
R=G*np.log(k/delta)*math.pow(k,1/2)
t=np.zeros(k+1)

i=1
while(i<=k):
 if(i<=round((k/R))-1):
  t[i]=float(R/(i*k))
  print(f"t[i]={float(R/(i*k))}")
 elif(i==round(k/R)):
  t[i]=float(R*np.log(R/delta)/k)
  print(f"t[i]={float(R*np.log(R/delta)/k)}")
  print(f"解{R*np.log(R/delta)}")
 else:
  t[i]=0
  print(f"t[i]={0}")
 i=i+1

i=1
beta1=0
f=0
while(i<=k):
 f=f+t[i]
 beta=phi[i]+t[i]
 beta1=beta1+beta
#  print(beta1)
#  print(f)
 i=i+1


distribution=np.zeros(k+1)
for i in range(1,k+1):
 distribution[i]=(phi[i]+t[i])/beta1
 print(f"distribution[{i}]={distribution[i]}")
 i=i+1
    # print(f"i={i},down_first_item={down_first_item}")

print(f"n={(k*beta1)},beta1={beta1},f={f},k/R={k/R},R={R}")
print(f"N={k+2*np.log(R/delta)*R}")
n=int(input("收到封包數:"))

# 建立此一矩陣，為顯示經計算後各state發生機率
probability=np.zeros([n+1,n+1,u+1], dtype=float)
for Q in range(0,n+1):
 probability[n-Q][Q][u]=math.pow(distribution[1],Q)*math.pow(1-distribution[1],n-Q)*math.comb(n,Q)

# 在probability矩陣中標示初始state其下一層的所有可能state的值標記為0.5
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

# 為計算state的發生機率，並將其上一層所有可能轉變state標記
# x=j
# y=i
def value_probability(c1,r1,u1):
  h1=0
  if(u1>=1):
   for x in range (0,n-c1+1):
    y=0
    while(r1+1+y-x<=n):
       if((r1+y-x>=0)&(r1-x>=0)):
        h=probability[c1+x][r1+1+y-x][u1+1]*state_transfer(x,y,c1,r1,u1+1) #
       else:
        h=0
       h1=h1+h
       y=y+1
  #  print(f"({c1},{r1},{u1}),{c1*((2*n)*(u+1))+r1*(u+1)+u1}")
   probability[c1][r1][u1]=h1


# def value_probability(c1, r1, u1):
#     h1 = 0
#     if u1 >= 1:
#         x_values = np.arange(0, n-c1+1)
#         y_values = np.arange(0, n)

#         x, y = np.meshgrid(x_values, y_values)
#         valid_indices = (r1 + 1 + y - x <= n) & (r1 + y - x >= 0) & (r1 - x >= 0)

#         h_values = probability[c1 + x, r1 + 1 + y - x, u1 + 1] * state_transfer(x, y, c1, r1, u1 + 1)
#         h_values[~valid_indices] = 0

#         h1 = np.sum(h_values)

#     probability[c1, r1, u1] = h1



  # elif(u1==1):
  #  for x in range (0,n-c1+1):
  #   y=0
  #   while(r1+1+y-x<=n):
  #      if((r1+y-x>=0)&(r1-x>=0)):
  #       h=probability[c1+x][r1+1+y-x][u1+1]*state_transfer_1(x,y,c1,r1,u1+1) #
  #      else:
  #       h=0
  #      h1=h1+h
  #      y=y+1
  # #  print(f"({c1},{r1},{u1}),{c1*((2*n)*(u+1))+r1*(u+1)+u1}")
  #  probability[c1][r1][u1]=h1


  # else:
  #  z=1
  #  while(r1+z<=n):
  #   h=probability[0][r1+z][u1+1]
  #   h1=h1+h
  #   z=z+1
  # #  print(f"({c1},{r1},{u1}),{c1*((2*n)*(u+1))+r1*(u+1)+u1}")
  #  probability[0][r1][u1]=h1



# 為計算轉移機率
def state_transfer (j,i,c1,r1,u1):
    if u1>1:
     value=np.float32(math.comb(c1+j,j)* math.pow(Pu,j)* math.pow(1-Pu,c1)*math.comb(r1+i-j,i)* math.pow(1/u1,i)* math.pow(1-1/u1,r1-j))
     return value
    else:
     value=1
     return value

# 為計算Pu的機率
# def release (distribution,k,u):
#   numerator=0
#   denominator=0
#   for i in range(1,k+1):
#     a=float(math .comb(k,i))
#     down_first_item=float(distribution[i]*math .comb(k-u,i)/a)
#     down_second_item=float(distribution[i]*math .comb(k-u,i-1)*u/a)
#     if(i==1):
#         numerator=numerator
#     else:
#         upper=float(distribution[i]*math .comb(u,2)*math .comb(k-u,i-2)*(2/u)/a)
#         temp1=upper
#         numerator=numerator+temp1
#     temp2=down_first_item+down_second_item
#     denominator=denominator+temp2

#   denominator=1-denominator
#   Pu=float(numerator/denominator)
#   return Pu




# 將所有state利用迴圈計算其發生機率
x=n
# under_layer(c,r,u)
u2=u
D=p=0
L=1
while(u2>2):
  Pu=release (distribution,k,u2)
  for j in range (0,n+1):
   for i in range (0,n+1):
      if(j+i<=x):
      #  print("loading...")
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
# 將值以以四捨五入方式取到小數點第4位
# rounded_array = np.round(probability, decimals=4)
print(probability)


output_file_path = "C:/Users/88697/Desktop/LTcode.txt"
# np.savetxt(output_file_path, rounded_array, fmt='%d', delimiter='\t')



# 将三维矩阵写入文本文件
with open(output_file_path, 'w') as file:
    for matrix_2d in probability:
        np.savetxt(file, matrix_2d, fmt='%f', delimiter='\t')
        file.write('\n')  # 在每个2D矩阵之间添加空行

print(f"n={n}")
print(f"1+R={1+R}")
while True:
        A,B,C=map(int,input("欲知state的probability:(c,r,u)").split())
        print(probability[A][B][C])


