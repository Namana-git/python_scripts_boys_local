import numpy as np
import math

# read a file with filename ud.dat

# Now the tt state has 36 eigenvectors on which the localized states have non negligible projection
# so tt is an array of size 36
# ss is also an array of size 36
# group_3 is now an array of size 16
# group_5 is an array of size 24
# group_4 is an array of size 24    
# sing is an array of size 16

# and the corresponding loc_arrays are also of same size

nv = 4
nc = 4

n_tt = nc*(nc-1)*nv*(nv-1)//4
n_ss = nc*(nc-1)*nv*(nv-1)//4
n_grp3 = nc*nv 
n_grp5 = nc*(nc-1)*nv//2
n_grp4 = nc*nv*(nv-1)//2
nsing = nc*nv

tt = np.zeros((n_tt), dtype=float)
tt_i = np.zeros((n_tt), dtype=int)
tt_j = np.zeros((n_tt), dtype=int)
tt_a = np.zeros((n_tt), dtype=int)
tt_b = np.zeros((n_tt), dtype=int)

ss = np.zeros((n_ss), dtype=float)
ss_i = np.zeros((n_ss), dtype=int)
ss_j = np.zeros((n_ss), dtype=int)
ss_a = np.zeros((n_ss), dtype=int)
ss_b = np.zeros((n_ss), dtype=int)

grp_3 = np.zeros((n_grp3), dtype=float)
grp_3_i = np.zeros((n_grp3), dtype=int)
grp_3_a = np.zeros((n_grp3), dtype=int)

grp_4 = np.zeros((n_grp4), dtype=float)
grp_4_i = np.zeros((n_grp4), dtype=int)
grp_4_a = np.zeros((n_grp4), dtype=int)
grp_4_b = np.zeros((n_grp4), dtype=int)

grp_5 = np.zeros((n_grp5), dtype=float)
grp_5_i = np.zeros((n_grp5), dtype=int)
grp_5_j = np.zeros((n_grp5), dtype=int)
grp_5_a = np.zeros((n_grp5), dtype=int)

sing = np.zeros((nsing), dtype=float)
sing_i = np.zeros((nsing), dtype=int)
sing_a = np.zeros((nsing), dtype=int)
# tt.dat file contains ijab indices and the corresponding tt values
# read all the above values from tt.dat

with open('tt.dat', 'r') as f:
    for x in range(n_tt):
        line = f.readline()
        values = line.split()
        tt_i[x] = int(values[0])
        tt_j[x] = int(values[1])
        tt_a[x] = int(values[2])
        tt_b[x] = int(values[3])
        tt[x] = float(values[4])
        print(tt_i[x], tt_j[x], tt_a[x], tt_b[x], tt[x])
with open('ss.dat', 'r') as f:
    for x in range(n_ss):
        line = f.readline()
        values = line.split()
        ss_i[x] = int(values[0])
        ss_j[x] = int(values[1])
        ss_a[x] = int(values[2])
        ss_b[x] = int(values[3])
        ss[x] = float(values[4])
        print(ss_i[x], ss_j[x], ss_a[x], ss_b[x], ss[x])

with open('grp_3.dat', 'r') as f:
    for x in range(n_grp3):
        line = f.readline()
        values = line.split()
        grp_3_i[x] = int(values[0])
        grp_3_a[x] = int(values[1])
        grp_3[x] = float(values[2])
        print(grp_3_i[x], grp_3_a[x], grp_3[x])

with open('grp_4.dat', 'r') as f:
    for x in range(n_grp4):
        line = f.readline()
        values = line.split()
        grp_4_i[x] = int(values[0])
        grp_4_a[x] = int(values[1])
        grp_4_b[x] = int(values[2])
        grp_4[x] = float(values[3])
        print(grp_4_i[x], grp_4_a[x], grp_4_b[x], grp_4[x])
with open('grp_5.dat', 'r') as f:
    for x in range(n_grp5):
        line = f.readline()
        values = line.split()
        grp_5_i[x] = int(values[0])
        grp_5_j[x] = int(values[1])
        grp_5_a[x] = int(values[2])
        grp_5[x] = float(values[3])
        print(grp_5_i[x], grp_5_j[x], grp_5_a[x], grp_5[x])

with open('sing.dat', 'r') as f:
    for x in range(nsing):
        line = f.readline()
        values = line.split()
        sing_i[x] = int(values[0])
        sing_a[x] = int(values[1])
        sing[x] = float(values[2])
        print(sing_i[x], sing_a[x], sing[x])






# map tt indices with ijab and ss indices with ijab and grp_3 indices with ia and grp_4 indices with iab and grp_5 indices with ija and sing indices with ia


loc_tt = np.zeros((n_tt), dtype=float)
loc_ss = np.zeros((n_ss), dtype=float)
loc_grp_3 = np.zeros((n_grp3), dtype=float)
loc_grp_4 = np.zeros((n_grp4), dtype=float)
loc_grp_5 = np.zeros((n_grp5), dtype=float)
loc_sing = np.zeros((nsing), dtype=float)

    

Uh= np.zeros((nv,nv), dtype=float)
Ul = np.zeros((nc,nc), dtype=float)
with open('Uh.dat', 'r') as f:
    for i in range(nv):
        #line = f.readline()
        #values = line.split()
        for j in range(nv):
            line = f.readline()
            values = line.split()
            Uh[i,j] = float(values[2])

with open('Ul.dat', 'r') as f:
    for i in range(nc):
        for j in range(nc):
            line = f.readline()
            values = line.split()
            Ul[i,j] = float(values[2])
            
        #print()


# read the values from the file

#Uh[1,0] =   0.707166092116449  # HOMO-1  contribution to left
#Uh[0,0] =  0.7070474652813279  # HOMO contribution to up orbital
#Uh[1,1] =  0.7070474652813277  # HOMO-1 contribution to right
#Uh[0,1] =  -0.7071660921164493 # HOMO contribution to right
#
#Ul[0,0] =  0.7070616462336817 # LUMO contribution to left
#Ul[1,0] =  -0.7071519132586113
#Ul[0,1] =  0.7071519132586109
#Ul[1,1] = 0.7070616462336821

# HOMO
# 0 (1) --> right up 
# 1 (2) --> left up
# 2 (3) --> right down
# 3 (4) --> left down


# LUMO
# 0 (1) --> right down
# 1 (2) --> left up
# 2 (3) --> right up
# 3 (4) --> left down



loc_tt_i = np.zeros((n_tt),dtype=int)
loc_tt_j = np.zeros((n_tt),dtype=int)
loc_tt_a = np.zeros((n_tt),dtype=int)
loc_tt_b = np.zeros((n_tt),dtype = int)

loc_ss_i = np.zeros((n_ss),dtype=int)
loc_ss_j = np.zeros((n_ss),dtype=int)
loc_ss_a = np.zeros((n_ss),dtype=int)
loc_ss_b = np.zeros((n_ss),dtype = int)

loc_grp_3_i = np.zeros((n_grp3),dtype=int)
loc_grp_3_a = np.zeros((n_grp3),dtype = int)


loc_grp_4_i  = np.zeros((n_grp4),dtype=int)
loc_grp_4_a = np.zeros((n_grp4),dtype=int)
loc_grp_4_b =  np.zeros((n_grp4),dtype=int)

loc_grp_5_i = np.zeros((n_grp5),dtype=int)
loc_grp_5_j = np.zeros((n_grp5),dtype=int)
loc_grp_5_a = np.zeros((n_grp5),dtype=int)

loc_sing_i = np.zeros((nsing),dtype=int)
loc_sing_a = np.zeros((nsing),dtype=int)

c = 0
for j in range(nc):                 # Fortran: j = 1,nc
    for i in range(j+1, nc):        # Fortran: i = j+1,nc
        for b in range(nv):      # Fortran: beta = 1,nv
            for a in range(b+1, nv):
                loc_tt_i[c] = i+1
                loc_tt_j[c] = j+1
                loc_tt_a[c] = a+1
                loc_tt_b[c] = b+1
                loc_ss_i[c] = i+1
                loc_ss_j[c] = j+1
                loc_ss_a[c] = a+1
                loc_ss_b[c] = b+1
                c = c+1
c1 =0
for i in range(nc): 
    for a in range(nv):
        loc_grp_3_i[c1] = i+1
        loc_grp_3_a[c1] = a+1
        loc_sing_i[c1] = i+1
        loc_sing_a[c1] = a+1
        c1 = c1+1


c2 = 0

for i in range(nc):
    for b in range(nv):
        for a in range(b+1, nv):
            loc_grp_4_i[c2] = i+1
            loc_grp_4_a[c2] = a+1
            loc_grp_4_b[c2] = b+1
            c2 = c2+1

c3 = 0
for j in range(nc):
    for i in range(j+1, nc):
        for a in range(nv):
            loc_grp_5_i[c3] = i+1
            loc_grp_5_j[c3] = j+1
            loc_grp_5_a[c3] = a+1
            c3 = c3+1




            





for x in range(n_tt):  # loc orbitals
     loc_i = loc_tt_i[x] -1
     loc_j = loc_tt_j[x] -1
     loc_a = loc_tt_a[x] -1
     loc_b = loc_tt_b[x] -1
     for y in range(n_tt): # front orbitals
         i = tt_i[y] -1
         j = tt_j[y] -1
         a = tt_a[y] -1
         b = tt_b[y] -1
         loc_tt[x] = loc_tt[x] +  ((Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_a]*Uh[b,loc_b]) - (Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[b,loc_b]) - (Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_b]*Uh[b,loc_a]) + (Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_b]*Uh[b,loc_a]))*tt[y] 
#
#loc_tt =  ((Ul[1,1]*Ul[0,0]*Uh[1,1]*Uh[0,0]) - (Ul[1,0]*Ul[0,1]*Uh[1,1]*Uh[0,0]) - (Ul[1,1]*Ul[0,0]*Uh[1,0]*Uh[0,1]) + (Ul[1,0]*Ul[0,1]*Uh[1,0]*Uh[0,1]))*tt 
#
#
## has contribution from ss,grp3_0000, grp3_1111,grp3_0011,grp3_1100, grp5_1000, grp5_1011, grp4_0010, grp4_1110             

for x in range(n_ss):
      loc_i = loc_tt_i[x] -1
      loc_j = loc_tt_j[x] -1
      loc_a = loc_tt_a[x] -1
      loc_b = loc_tt_b[x] -1 
      for y in range(n_ss):
          i = ss_i[y] -1
          j = ss_j[y] -1
          a = ss_a[y] -1
          b = ss_b[y] -1
          loc_ss[x] = loc_ss[x] +  ((Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_a]*Uh[b,loc_b]) + (Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[b,loc_b]) + (Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_b]*Uh[b,loc_a]) + (Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_b]*Uh[b,loc_a]))*ss[y]
      for z in range(n_grp3):
          i = grp_3_i[z] -1
          a = grp_3_a[z] -1
          loc_ss[x] = loc_ss[x] + 2*(Ul[i,loc_i]*Ul[i,loc_j]*Uh[a,loc_a]*Uh[a,loc_b])*grp_3[z]
      for r in range(n_grp5):
          i = grp_5_i[r] -1
          j = grp_5_j[r] -1
          a = grp_5_a[r] -1
          loc_ss[r] = loc_ss[r] + 2/np.sqrt(2)*((Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_a]*Uh[a,loc_b])+(Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[a,loc_b]))*grp_5[r]
      for s in range(n_grp4):
          i = grp_4_i[s] -1
          a = grp_4_a[s] -1
          b = grp_4_b[s] -1
          loc_ss[s] = loc_ss[s] + 2/np.sqrt(2)*((Ul[i,loc_i]*Ul[i,loc_j]*Uh[a,loc_a]*Uh[b,loc_b])+(Ul[i,loc_i]*Ul[i,loc_j]*Uh[a,loc_b]*Uh[b,loc_a]))*grp_4[s]

#loc_ss = ((Ul[1,1]*Ul[0,0]*Uh[1,1]*Uh[0,0]) + (Ul[1,0]*Ul[0,1]*Uh[1,1]*Uh[0,0]) + (Ul[1,1]*Ul[0,0]*Uh[1,0]*Uh[0,1]) + (Ul[1,0]*Ul[0,1]*Uh[1,0]*Uh[0,1]))*ss + \
#         2*(Ul[0,1]*Ul[0,0]*Uh[0,1]*Uh[0,0])*grp_3_0000 + \
#         2*(Ul[1,1]*Ul[1,0]*Uh[1,1]*Uh[1,0])*grp_3_1111 + \
#         2*(Ul[1,1])*Ul[1,0]*Uh[0,1]*Uh[0,0]*grp_3_1100 + \
#         2*(Ul[0,1]*Ul[0,0]*Uh[1,1]*Uh[1,0])*grp_3_0011 + \
#         2/np.sqrt(2)*((Ul[1,1]*Ul[0,0]*Uh[0,1]*Uh[0,0])+(Ul[1,0]*Ul[0,1]*Uh[0,0]*Uh[0,1]))*grp_5_1000 + \
#         2/np.sqrt(2)*((Ul[1,1]*Ul[0,0]*Uh[1,1]*Uh[1,0])+(Ul[1,0]*Ul[0,1]*Uh[1,0]*Uh[1,1]))*grp_5_1011 + \
#         2/np.sqrt(2)*((Ul[0,1]*Ul[0,0]*Uh[1,1]*Uh[0,0])+(Ul[0,0]*Ul[0,1]*Uh[1,0]*Uh[0,1]))*grp_4_0010 + \
#         2/np.sqrt(2)*((Ul[1,1]*Ul[1,0]*Uh[1,1]*Uh[0,0])+(Ul[1,0]*Ul[1,1]*Uh[1,0]*Uh[0,1]))*grp_4_1110

for x in range(n_grp3):
    loc_i = loc_grp_3_i[x] -1
    loc_a = loc_grp_3_a[x] -1
    for y in range(n_ss):
        i = ss_i[y] -1
        j = ss_j[y] -1
        a = ss_a[y] -1
        b = ss_b[y] -1
        loc_grp_3[x] = loc_grp_3[x] + 2*(Ul[i,loc_i]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[b,loc_a])*ss[y]
    for z in range(n_grp3):
        i = grp_3_i[z] -1
        a = grp_3_a[z] -1
        loc_grp_3[x] = loc_grp_3[x] + (Ul[i,loc_i]*Ul[i,loc_i]*Uh[a,loc_a]*Uh[a,loc_a])*grp_3[z]
    for r in range(n_grp5):
        i = grp_5_i[r] -1
        j = grp_5_j[r] -1
        a = grp_5_a[r] -1
        loc_grp_3[x] = loc_grp_3[x] + 2/np.sqrt(2)*((Ul[i,loc_i]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[a,loc_a]))*grp_5[r]
    for s in range(n_grp4):
        i = grp_4_i[s] -1
        a = grp_4_a[s] -1
        b = grp_4_b[s] -1
        loc_grp_3[x] = loc_grp_3[x] + 2/np.sqrt(2)*((Ul[i,loc_i]*Ul[i,loc_i]*Uh[a,loc_a]*Uh[b,loc_a]))*grp_4[s]





    


#
## loc_grp_3 has contributions from ss,grp_3_0000, grp_3_1111,grp_3_0011,grp_3_1100, grp5_1000, grp5_1011, grp4_0010, grp4_1110
#loc_grp_3_0000 = 2*(Ul[1,0]*Ul[0,0]*Uh[1,0]*Uh[0,0])*ss + \
#                (Ul[0,0]*Ul[0,0]*Uh[0,0]*Uh[0,0])*grp_3_0000 + \
#                (Ul[1,0]*Ul[1,0]*Uh[1,0]*Uh[1,0])*grp_3_1111 + \
#                (Ul[0,0]*Ul[0,0]*Uh[1,0]*Uh[1,0])*grp_3_0011 + \
#                (Ul[1,0]*Ul[1,0]*Uh[0,0]*Uh[0,0])*grp_3_1100 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[0,0]*Uh[0,0]*Uh[0,0]))*grp_5_1000 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[0,0]*Uh[1,0]*Uh[1,0]))*grp_5_1011 + \
#                2/np.sqrt(2)*((Ul[0,0]*Ul[0,0]*Uh[1,0]*Uh[0,0]))*grp_4_0010 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[1,0]*Uh[1,0]*Uh[0,0]))*grp_4_1110
#
#loc_grp_3_1111 = 2*(Ul[1,1]*Ul[0,1]*Uh[1,1]*Uh[0,1])*ss + \
#                (Ul[0,1]*Ul[0,1]*Uh[0,1]*Uh[0,1])*grp_3_0000 + \
#                (Ul[1,1]*Ul[1,1]*Uh[1,1]*Uh[1,1])*grp_3_1111 + \
#                (Ul[0,1]*Ul[0,1]*Uh[1,1]*Uh[1,1])*grp_3_0011 + \
#                (Ul[1,1]*Ul[1,1]*Uh[0,1]*Uh[0,1])*grp_3_1100 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[0,1]*Uh[0,1]*Uh[0,1]))*grp_5_1000 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[0,1]*Uh[1,1]*Uh[1,1]))*grp_5_1011 + \
#                2/np.sqrt(2)*((Ul[0,1]*Ul[0,1]*Uh[1,1]*Uh[0,1]))*grp_4_0010 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[1,1]*Uh[1,1]*Uh[0,1]))*grp_4_1110 
#
#loc_grp_3_0011 = 2*(Ul[1,0]*Ul[0,0]*Uh[1,1]*Uh[0,1])*ss + \
#                (Ul[0,0]*Ul[0,0]*Uh[0,1]*Uh[0,1])*grp_3_0000 + \
#                (Ul[1,0]*Ul[1,0]*Uh[1,1]*Uh[1,1])*grp_3_1111 + \
#                (Ul[0,0]*Ul[0,0]*Uh[1,1]*Uh[1,1])*grp_3_0011 + \
#                (Ul[1,0]*Ul[1,0]*Uh[0,1]*Uh[0,1])*grp_3_1100 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[0,0]*Uh[0,1]*Uh[0,1]))*grp_5_1000 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[0,0]*Uh[1,1]*Uh[1,1]))*grp_5_1011 + \
#                2/np.sqrt(2)*((Ul[0,0]*Ul[0,0]*Uh[1,1]*Uh[0,1]))*grp_4_0010 + \
#                2/np.sqrt(2)*((Ul[1,0]*Ul[1,0]*Uh[1,1]*Uh[0,1]))*grp_4_1110
#
#loc_grp_3_1100 = 2*(Ul[1,1]*Ul[0,1]*Uh[1,0]*Uh[0,0])*ss + \
#                (Ul[0,1]*Ul[0,1]*Uh[0,0]*Uh[0,0])*grp_3_0000 + \
#                (Ul[1,1]*Ul[1,1]*Uh[1,0]*Uh[1,0])*grp_3_1111 + \
#                (Ul[0,1]*Ul[0,1]*Uh[1,0]*Uh[1,0])*grp_3_0011 + \
#                (Ul[1,1]*Ul[1,1]*Uh[0,0]*Uh[0,0])*grp_3_1100 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[0,1]*Uh[0,0]*Uh[0,0]))*grp_5_1000 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[0,1]*Uh[1,0]*Uh[1,0]))*grp_5_1011 + \
#                2/np.sqrt(2)*((Ul[0,1]*Ul[0,1]*Uh[1,0]*Uh[0,0]))*grp_4_0010 + \
#                2/np.sqrt(2)*((Ul[1,1]*Ul[1,1]*Uh[1,0]*Uh[0,0]))*grp_4_1110
#
## loc_grp_5 has contributions from ss,grp_3_0000, grp_3_1111,grp_3_0011,grp_3_1100, grp5_1000, grp5_1011, grp4_0010, grp4_1110
#


for x in range(n_grp5):
        loc_i = loc_grp_5_i[x] -1
        loc_j = loc_grp_5_j[x] -1
        loc_a = loc_grp_5_a[x] -1
        for y in range(n_ss):
            i = ss_i[y] -1
            j = ss_j[y] -1
            a = ss_a[y] -1
            b = ss_b[y] -1
            loc_grp_5[x] = loc_grp_5[x] + np.sqrt(2)*((Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_a]*Uh[b,loc_a])+(Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[b,loc_a]))*ss[y]
        for z in range(n_grp3):
            i = grp_3_i[z] -1
            a = grp_3_a[z] -1
            loc_grp_5[x] = loc_grp_5[x] + np.sqrt(2)*((Ul[i,loc_i]*Ul[i,loc_j]*Uh[a,loc_a]*Uh[a,loc_a]))*grp_3[z]
        for r in range(n_grp5):
            i = grp_5_i[r] -1
            j = grp_5_j[r] -1
            a = grp_5_a[r] -1
            loc_grp_5[x] = loc_grp_5[x] + ((Ul[i,loc_i]*Ul[j,loc_j]*Uh[a,loc_a]*Uh[a,loc_a])+ (Ul[i,loc_j]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[a,loc_a]))*grp_5[r]
        for s in range(n_grp4):
            i = grp_4_i[s] -1
            a = grp_4_a[s] -1
            b = grp_4_b[s] -1
            loc_grp_5[x] = loc_grp_5[x] + 2*((Ul[i,loc_i]*Ul[i,loc_j]*Uh[a,loc_a]*Uh[b,loc_a]))*grp_4[s]



#loc_grp_5_1000 = np.sqrt(2)*((Ul[1,1]*Ul[0,0]*Uh[1,0]*Uh[0,0]) + (Ul[1,0]*Ul[0,1]*Uh[1,0]*Uh[0,0]))*ss + \
#                 np.sqrt(2)*(Ul[0,1]*Ul[0,0]*Uh[0,0]*Uh[0,0])*grp_3_0000 + \
#                 np.sqrt(2)*(Ul[1,1]*Ul[1,0]*Uh[1,0]*Uh[1,0])*grp_3_1111 + \
#                 np.sqrt(2)*(Ul[0,1]*Ul[0,0]*Uh[1,0]*Uh[1,0])*grp_3_0011 + \
#                 np.sqrt(2)*(Ul[1,1]*Ul[1,0]*Uh[0,0]*Uh[0,0])*grp_3_1100 + \
#                 (Ul[1,1]*Ul[0,0]*Uh[0,0]*Uh[0,0] + Ul[1,0]*Ul[0,1]*Uh[0,0]*Uh[0,0])*grp_5_1000 + \
#                 (Ul[1,1]*Ul[0,0]*Uh[1,0]*Uh[1,0] + Ul[1,0]*Ul[0,1]*Uh[1,0]*Uh[1,0])*grp_5_1011 + \
#                 2*(Ul[0,1]*Ul[0,0]*Uh[1,0]*Uh[0,0])*grp_4_0010 + \
#                 2*(Ul[1,1]*Ul[1,0]*Uh[1,0]*Uh[0,0])*grp_4_1110
#
#
#
#loc_grp_5_1011 =  np.sqrt(2)*((Ul[1,1]*Ul[0,0]*Uh[1,1]*Uh[0,1])+(Ul[1,0]*Ul[0,1]*Uh[1,1]*Uh[0,1]))*ss + \
#                  np.sqrt(2)*(Ul[0,1]*Ul[0,0]*Uh[0,1]*Uh[0,1])*grp_3_0000 + \
#                  np.sqrt(2)*(Ul[1,1]*Ul[1,0]*Uh[1,1]*Uh[1,1])*grp_3_1111 + \
#                  np.sqrt(2)*(Ul[0,1]*Ul[0,0]*Uh[1,1]*Uh[1,1])*grp_3_0011 + \
#                  np.sqrt(2)*(Ul[1,1]*Ul[1,0]*Uh[0,1]*Uh[0,1])*grp_3_1100 + \
#                  (Ul[1,1]*Ul[0,0]*Uh[0,1]*Uh[0,1] + Ul[1,0]*Ul[0,1]*Uh[0,1]*Uh[0,1])*grp_5_1000 + \
#                  (Ul[1,1]*Ul[0,0]*Uh[1,1]*Uh[1,1] + Ul[1,0]*Ul[0,1]*Uh[1,1]*Uh[1,1])*grp_5_1011 + \
#                  2*(Ul[0,1]*Ul[0,0]*Uh[1,1]*Uh[0,1])*grp_4_0010 + \
#                  2*(Ul[1,1]*Ul[1,0]*Uh[1,1]*Uh[0,1])*grp_4_1110

for x in range(n_grp4):
        loc_i = loc_grp_4_i[x] -1
        loc_a = loc_grp_4_a[x] -1
        loc_b = loc_grp_4_b[x] -1
        for y in range(n_ss):
            i = ss_i[y] -1
            j = ss_j[y] -1
            a = ss_a[y] -1
            b = ss_b[y] -1
            loc_grp_4[x] = loc_grp_4[x] + np.sqrt(2)*((Ul[i,loc_i]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[b,loc_b])+(Ul[i,loc_i]*Ul[j,loc_i]*Uh[a,loc_b]*Uh[b,loc_a]))*ss[y]
        for z in range(n_grp3):
            i = grp_3_i[z] -1
            a = grp_3_a[z] -1
            loc_grp_4[x] = loc_grp_4[x] + np.sqrt(2)*((Ul[i,loc_i]*Ul[i,loc_i]*Uh[a,loc_a]*Uh[a,loc_b]))*grp_3[z]
        for r in range(n_grp5):
            i = grp_5_i[r] -1
            j = grp_5_j[r] -1
            a = grp_5_a[r] -1
            loc_grp_4[x] = loc_grp_4[x] + 2*((Ul[i,loc_i]*Ul[j,loc_i]*Uh[a,loc_a]*Uh[a,loc_b]))*grp_5[r]
        for s in range(n_grp4):
            i = grp_4_i[s] -1
            a = grp_4_a[s] -1
            b = grp_4_b[s] -1
            loc_grp_4[x] = loc_grp_4[x] + ((Ul[i,loc_i]*Ul[i,loc_i]*Uh[a,loc_a]*Uh[b,loc_b]) + (Ul[i,loc_i]*Ul[i,loc_i]*Uh[a,loc_b]*Uh[b,loc_a]))*grp_4[s]
#loc_grp_4_0010 =  np.sqrt(2)*((Ul[1,0]*Ul[0,0]*Uh[1,1]*Uh[0,0])+(Ul[1,0]*Ul[0,0]*Uh[1,0]*Uh[0,1]))*ss +  \
#                  np.sqrt(2)*(Ul[0,0]*Ul[0,0]*Uh[0,1]*Uh[0,0])*grp_3_0000 + \
#                  np.sqrt(2)*(Ul[1,0]*Ul[1,0]*Uh[1,1]*Uh[1,0])*grp_3_1111 + \
#                  np.sqrt(2)*(Ul[0,0]*Ul[0,0]*Uh[1,1]*Uh[1,0])*grp_3_0011 + \
#                  np.sqrt(2)*(Ul[1,0]*Ul[1,0]*Uh[0,1]*Uh[0,0])*grp_3_1100 + \
#                  2*(Ul[1,0]*Ul[0,0]*Uh[0,1]*Uh[0,0])*grp_5_1000 + \
#                  2*(Ul[1,0]*Ul[0,0]*Uh[1,1]*Uh[1,0])*grp_5_1011 + \
#                  (Ul[0,0]*Ul[0,0]*Uh[1,1]*Uh[0,0] + Ul[0,0]*Ul[0,0]*Uh[1,0]*Uh[0,1])*grp_4_0010 + \
#                  (Ul[1,0]*Ul[1,0]*Uh[1,1]*Uh[0,0] + Ul[1,0]*Ul[1,0]*Uh[1,0]*Uh[0,1])*grp_4_1110
#
#loc_grp_4_1110 =  np.sqrt(2)*((Ul[1,1]*Ul[0,1]*Uh[1,1]*Uh[0,0])+(Ul[1,1]*Ul[0,1]*Uh[1,0]*Uh[0,1]))*ss +  \
#                  np.sqrt(2)*(Ul[0,1]*Ul[0,1]*Uh[0,1]*Uh[0,0])*grp_3_0000 + \
#                  np.sqrt(2)*(Ul[1,1]*Ul[1,1]*Uh[1,1]*Uh[1,0])*grp_3_1111 + \
#                  np.sqrt(2)*(Ul[0,1]*Ul[0,1]*Uh[1,1]*Uh[1,0])*grp_3_0011 + \
#                  np.sqrt(2)*(Ul[1,1]*Ul[1,1]*Uh[0,1]*Uh[0,0])*grp_3_1100 + \
#                  2*(Ul[1,1]*Ul[0,1]*Uh[0,1]*Uh[0,0])*grp_5_1000 + \
#                  2*(Ul[1,1]*Ul[0,1]*Uh[1,1]*Uh[1,0])*grp_5_1011 + \
#                  (Ul[0,1]*Ul[0,1]*Uh[1,1]*Uh[0,0] + Ul[0,1]*Ul[0,1]*Uh[1,0]*Uh[0,1])*grp_4_0010 + \
#                  (Ul[1,1]*Ul[1,1]*Uh[1,1]*Uh[0,0] + Ul[1,1]*Ul[1,1]*Uh[1,0]*Uh[0,1])*grp_4_1110
loc_sing = np.zeros((nsing),dtype=float)
for x in range(nsing):
        loc_i = loc_sing_i[x] -1
        loc_a = loc_sing_a[x] -1
        for y in range(nsing):
            i = sing_i[y] -1
            a = sing_a[y] -1
            loc_sing[x] = loc_sing[x] + (Ul[i,loc_i]*Uh[a,loc_a])*sing[y]  
            print("x,y,z,loc_sing[x]",x,y,((Ul[i,loc_i]*Uh[a,loc_a])*sing[y]),loc_sing[x]) 
             
#loc_00 =    (sing_00)*(Ul[0,0]*Uh[0,0]) + \
#            (sing_11)*(Ul[1,0]*Uh[1,0]) + \
#            (sing_10)*(Ul[1,0]*Uh[0,0]) + \
#            (sing_01)*(Ul[0,0]*Uh[1,0])
#
#loc_11 =    (sing_00)*(Ul[0,1]*Uh[0,1]) + \
#            (sing_11)*(Ul[1,1]*Uh[1,1]) + \
#            (sing_10)*(Ul[1,1]*Uh[0,1]) + \
#            (sing_01)*(Ul[0,1]*Uh[1,1])
#
#loc_10 =    (sing_00)*(Ul[0,1]*Uh[0,0]) + \
#            (sing_11)*(Ul[1,1]*Uh[1,0]) + \
#            (sing_10)*(Ul[1,1]*Uh[0,0]) + \
#            (sing_01)*(Ul[0,1]*Uh[1,0])
#
#loc_01 =   (sing_00)*(Ul[0,0]*Uh[0,1]) + (sing_11)*(Ul[1,0]*Uh[1,1]) + \
#           (sing_10)*(Ul[1,0]*Uh[0,1]) + \
#           (sing_01)*(Ul[0,0]*Uh[1,1])
#
# 
#print(loc_tt)
print("percent of loc_tt = ", np.sum((np.square(loc_tt))))  
print("percent of loc_ss = ", np.sum((np.square(loc_ss))))
sum1 = 0.0
for x in range(n_tt):
    loc_i = loc_tt_i[x] -1
    loc_j = loc_tt_j[x] -1
    loc_a = loc_tt_a[x] -1
    loc_b = loc_tt_b[x] -1

    # only add up contributions from electron  belonging to different dimer
    #  similarly add up contributions from holes belonging to different dimers

    if (((loc_i == 1) and (loc_j == 0)) or ((loc_i == 3) and (loc_j ==0)) or ((loc_i == 2) and (loc_j ==1)) or ((loc_i == 3) and (loc_j ==2))):
         if (((loc_a == 1) and (loc_b == 0)) or ((loc_a == 3) and (loc_b ==0)) or ((loc_a == 2) and (loc_b ==1)) or ((loc_a == 3) and (loc_b ==2))):
                #print("loc_tt = ", loc_tt[x])
                #print("loc_ss = ", loc_ss[x])
                sum1 = sum1 + ((loc_tt[x]*loc_tt[x]) + (loc_ss[x]*loc_ss[x]))

print("ME percent = ", sum1)

print("percent of loc_grp_3 = ", np.sum((np.square(loc_grp_3))))
print("percent of loc_grp_5 = ", np.sum((np.square(loc_grp_5))))
print("percent of loc_grp_4 = ", np.sum((np.square(loc_grp_4))))
print(" DE percent =", np.sum((np.square(loc_grp_3))) + np.sum((np.square(loc_grp_5))) + np.sum((np.square(loc_grp_4))) + (np.sum(np.square(loc_tt)) + np.sum(np.square(loc_ss)) - sum1))


# LE like states 
sum2 = 0.0
for x in range(nsing):
    loc_i = loc_sing_i[x] -1
    loc_a = loc_sing_a[x] -1
    if (((loc_i == 0) and (loc_a == 0)) or ((loc_i == 2) and (loc_a ==2)) or ((loc_i == 0) and (loc_a ==2)) or ((loc_i == 2) and (loc_a ==0)) or ((loc_i == 1) and (loc_a == 1)) or ((loc_i == 3) and (loc_a ==3)) or ((loc_i == 1) and (loc_a ==3)) or ((loc_i == 3) and (loc_a ==1))):
        #print("loc_sing = ", loc_sing[x])
        sum2 = sum2 + (loc_sing[x]*loc_sing[x])


print("LE percent = ", sum2)
print("CT percent = ", np.sum(np.square(loc_sing)) - sum2)









# The question is what combinations add to ME


#print("loc_tt = ", loc_tt)
#print("loc_ss = ", loc_ss)      
#print("loc_grp_3_0000 = ", loc_grp_3_0000)
#print("loc_grp_3_1111 = ", loc_grp_3_1111)
#print("loc_grp_3_0011 = ", loc_grp_3_0011)      
#print("loc_grp_3_1100 = ", loc_grp_3_1100)
#print("loc_grp_5_1000 = ", loc_grp_5_1000)  
#print("loc_grp_5_1011 = ", loc_grp_5_1011)
#print("loc_grp_4_0010 = ", loc_grp_4_0010)
#print("loc_grp_4_1110 = ", loc_grp_4_1110)
#print("loc_00 = ", loc_00)
#print("loc_11 = ", loc_11)
#print("loc_10 = ", loc_10)
#print("loc_01 = ", loc_01)
#print(" LE:",(loc_00)**2 + (loc_11)**2)
#print(" CT:",(loc_01)**2 + (loc_10)**2)
#print(" ME:",(loc_tt*loc_tt) + (loc_ss*loc_ss))
#print(" DE:",(loc_grp_3_0000*loc_grp_3_0000) + (loc_grp_3_1111*loc_grp_3_1111) + (loc_grp_3_0011*loc_grp_3_0011) + (loc_grp_3_1100*loc_grp_3_1100) + (loc_grp_5_1000*loc_grp_5_1000) + (loc_grp_5_1011*loc_grp_5_1011) + (loc_grp_4_0010*loc_grp_4_0010) + (loc_grp_4_1110*loc_grp_4_1110))
