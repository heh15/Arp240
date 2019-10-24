import matplotlib.pyplot as plt
import numpy as np

x=np.array(range(5))
y=np.array(range(5))
z=y+2


fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(x,y)


fig1=plt.figure()
plt.plot(x,z)

X,Y=[],[]
for lines in plt.gca.getlines():
    for x,y in lines.get_xydata():
        X.append(x)
        Y.append(y)

# plt.gca.getlines()

plt.plot(X,Y,'g')
