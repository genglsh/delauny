import os  
import os.path  
import gzip  
import math
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D

def dot(a, b):

    return a.x*b.x + a.y*b.y + a.z*b.z

def cross(a, b):

    return Point( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x )

class Point:

    def __init__(s, x, y, z=0):
        s.x = x
        s.y = y
        s.z = z

    def __repr__(s):

        return "( " + str(s.x) + ", " + str(s.y) + " )"

    def __add__(s, b):

        return Point(s.x+b.x, s.y+b.y)

    def __sub__(s, b):

        return Point(s.x-b.x, s.y-b.y)

    def __mul__(s, b):

        return Point(b*s.x, b*s.y)

    __rmul__ = __mul__

    def IsIn(self, t):
        ''' Checks to see if p is in t using the barycenter method'''
        b = t.v[1] - t.v[0]
        c = t.v[2] - t.v[0]
        d = self   - t.v[0]

        det = c.x*b.y-c.y*b.x
        u = (d.x*b.y-d.y*b.x)//float(det)
        v = (c.x*d.y-c.y*d.x)//float(det)
        return u >= 0 and v >= 0 and u + v < 1

    def IsInCircumcircleOf(self, T):

        a = T.v[0] - T.v[2]
        b = T.v[1] - T.v[2]

        # Ref: https://en.wikipedia.org/wiki/Circumscribed_circle#Circumcircle_equations
        z = cross(a,b)
        p0 = cross(dot(a,a)*b-dot(b,b)*a, z)*(0.5/dot(z,z)) + T.v[2]

        r2 = 0.25*dot(a, a)*dot(b,b)*dot(a-b, a-b)/dot(z, z)

        #print "IsInC"
        #print self, p0
        #print sqrt(r2), "\n"
        #print 

        return dot(self-p0, self-p0) <= r2

class Triangle:


    def __init__(self, a, b, c):
        
        self.v = [None]*3
        self.v[0] = a
        self.v[1] = b
        self.v[2] = c

        self.neighbour = [None]*3    # Adjacent triangles


    def __repr__(s):

        '''
        return '<%s, [%s, %s, %s]>' % (
                hex(id(s)), 
                hex(id(s.neighbour[0])), 
                hex(id(s.neighbour[1])), 
                hex(id(s.neighbour[2])))
        '''
        return '< ' + str(s.v) + ' >'

    #def vOppositeOf(self, T):
    #
    #   return self.v[(self.neighbour.index(T))] 

    #def ReplaceNeighbour(self, A, B):
    #   
    #    self.neighbour[self.neighbour.index(A)] = B
        
    def SetEdge(self, edge, T):
        'Set the edge neighbour that matches "edge" to T'

        temp_v = self.v + self.v[0:1]
        for i in range(3):
            if edge[0] == temp_v[i] and edge[1] == temp_v[i+1]:
                self.neighbour[(i+2)%3] = T
                return
        print ('This function should never get this far')
        print (edge)
        print (temp_v)
        print (T)

class Delaunay_Triangulation:
    """Bowyer Watson"""

    def __init__(self):

        # Create a two triangle 'frame'
        
#         a = Point(0, 0)
#         b = Point(99, 0)
#         c = Point(99, 99)
#         d = Point(0, 99)
#         a = Point(-80, -115)
#         b = Point(100, -115)
#         c = Point(100, 95)
#         d= Point(-80, 95)
        a = Point(110, 400)
        b = Point(1040, 400)
        c = Point(1040, 1000)
        d= Point(110, 1000)

        T1 = Triangle(a, d, b)
        T2 = Triangle(c, b, d)

        T1.neighbour[0] = T2
        T2.neighbour[0] = T1

        self.triangles = [T1, T2]


    def AddPoint(self, p):
       
        bad_triangles = []

        # Search for the triangle where the point is.
        ''' For now I am just doing a naive search,
        but I hope to replace this with an initial guess
        and a BFS'''
        for T in self.triangles:
            
            if p.IsInCircumcircleOf(T):
                bad_triangles.append(T)
        
        # print("查看当前坏三角形的个数， 如果坏三角形始终为2或者1，说明一般不会出现三个德洛内三角形外接圆相遮的情况")
        # print(len(bad_triangles))

        # Find the convex hull of the bad triangles.
        # Expressed a list of edges (point pairs) in ccw order
        boundary = self.Boundary(bad_triangles)
		
		#edge 是3个元素， 公共边的两个定点还有与之对应的三角形。


        for T in bad_triangles:
            self.triangles.remove(T)

        # Retriangle to hole
        new_triangles = []
        for edge in boundary:
            T = Triangle(p, edge[0], edge[1])

            T.neighbour[0] = edge[2]                   # To neighbour
            if T.neighbour[0]:
                T.neighbour[0].SetEdge(edge[1::-1], T)     # from neighbour

            new_triangles.append(T)

        # Link the new triangles
        N = len(new_triangles)
        for i, T in enumerate(new_triangles):
            T.neighbour[2] = new_triangles[(i-1) % N]   # back
            T.neighbour[1] = new_triangles[(i+1) % N]   # forward
   
        self.triangles.extend(new_triangles)

      
    def Boundary(self, bad_triangles):

        # Start with a triangle at random
        T = bad_triangles[0]
        edge = 0

        boundary = []

        while True:
            
            if len(boundary) > 1:
                if boundary[0] == boundary[-1]:
                    break

            if T.neighbour[edge] in bad_triangles:

                last = T
                T = T.neighbour[edge]

                edge = (T.neighbour.index(last) + 1) % 3 

            else:   # Found an edge that is on the boundary
                # Add to list
                boundary.append((T.v[(edge+1)%3], T.v[(edge+2)%3], T.neighbour[edge]))
                edge = (edge + 1) % 3

        return boundary[:-1]

    def export(self):

        ps = [p for t in self.triangles for p in t.v ]

        xs = [p.x for p in ps]
        ys = [p.y for p in ps]
        zs = [p.z for p in ps]
        #xs = list(set(xs))
        #ys = list(set(ys))

        ts = [(ps.index(t.v[0]), ps.index(t.v[1]), ps.index(t.v[2])  ) for t in self.triangles]

        return xs, ys, zs, ts


  
  
# def read_gz_file(path):
#     if os.path.exists(path):
#         with gzip.open(path, 'r') as pf:
#             for line in pf:
#                 yield line
#     else:
#         print('the path [{}] is not exist!'.format(path))
#
# con = read_gz_file(os.path.join('3DTEC', 'rangeimages/90333d51.abs.gz'))
# row  = 0
# col = 0
# for key, data in enumerate(con):
#     if key > 2:
#         data_str = data.decode('utf-8').strip().split(' ')
#         if key == 3:
#             flag_array = np.array(data_str).reshape(480, 640)
#             row, col = flag_array.shape
#         elif key == 4:
#             x_array = np.array(data_str).reshape(480, 640)
#         elif key == 5:
#             y_array = np.array(data_str).reshape(480, 640)
#         elif key == 6:
#             z_array = np.array(data_str).reshape(480, 640)
# select_x = list()
# select_y = list()
# select_z = list()
# dir_x = 0
# dir_y = 0
# step = 8
# while dir_x < row:
#     while dir_y < col:
#         if flag_array[dir_x][dir_y] == '1':
#             select_x.append(float(x_array[dir_x][dir_y]))
#             select_y.append(float(y_array[dir_x][dir_y]))
#             select_z.append(float(z_array[dir_x][dir_y]))
#         dir_y += step
#     dir_x += step
#     dir_y = 0
# row  = 0
# col = 0
# select_x = list()
# select_y = list()
# select_z = list()
all_points = np.loadtxt(os.path.join('data', 'points.txt'))
select_x = list(all_points[:, 0])
select_z = list(all_points[:, 1])
select_y = list(all_points[:, 2])
# print(all_points)
# print(select_z)
print(len(select_x))
print(max(select_x), min(select_x))
print(max(select_y), min(select_y))
print(max(select_z), min(select_z))
DT = Delaunay_Triangulation()
# 经过该步骤构建了 一个Delaunay_Triangulation 对象
for x, y, z in zip(select_x, select_y, select_z):
    DT.AddPoint(Point(x, y, z))

XS, YS, ZS, TS = DT.export()

ax = plt.figure().add_subplot(111, projection='3d')
# ax.set_xlim(-100000, 100000)
# ax.set_ylim(-100000, 100000)
# ax.set_zlim(-1600, -1500)
ax.set_xlim(100, 2000)
ax.set_zlim(0, 1000)
ax.set_ylim(400, 1000)
ax.set_xlabel('X label')
ax.set_ylabel('Y label')
ax.set_zlabel('Z label')
len_c = len(TS)
for i in range(len_c):
    if i % 100 == 0 :
        print(i, len_c)
    i_0 = TS[i][0]
    i_1 = TS[i][1]
    i_2 = TS[i][2]
    try:
        if(ZS[i_0] == 0 or ZS[i_1] == 0 or ZS[i_2] == 0):
            continue
        x1 = [XS[i_0], XS[i_1]]
        y1 = [YS[i_0], YS[i_1]]
        z1 = [ZS[i_0], ZS[i_1]]
        x2 = [XS[i_0], XS[i_2]]
        y2 = [YS[i_0], YS[i_2]]
        z2 = [ZS[i_0], ZS[i_2]]
        x3 = [XS[i_2], XS[i_1]]
        y3 = [YS[i_2], YS[i_1]]
        z3 = [ZS[i_2], ZS[i_1]]


        ax.plot(x1, y1, z1, c='r')
        ax.plot(x2, y2, z2, c='r')
        ax.plot(x3, y3, z3, c='r')
    except:
        print(len(ZS), i_0, i_1)
plt.show()