### See end of the code for the plot settings

#Module importation
import numpy as np
import csv
import sys
import math as math
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import lines
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import scipy as scipy
from matplotlib.colors import LogNorm
from pylab import figure, cm
import string
import os
linestyles = ['-', '--', '-.', ':']

script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in

#Pointer initialization

table1_Om01 = []

#Data initialization

OMEGA1_Om01 = []
ARPES_add1_Om01 = []
ARPES_red1_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP00.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3

with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table1_Om01.append(row)


#Size of the imported array

N = len(table1_Om01)

M = len(table1_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA1_Om01.append(float(table1_Om01[i][0].replace(',','.')))

    ARPES_add1_Om01.append(float(table1_Om01[i][1].replace(',','.')))
    
    ARPES_red1_Om01.append(float(table1_Om01[i][2].replace(',','.')))
    

#Pointer initialization

table2_Om01 = []

#Data initialization

OMEGA2_Om01= []
ARPES_add2_Om01 = []
ARPES_red2_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP01.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table2_Om01.append(row)


#Size of the imported array

N = len(table2_Om01)

M = len(table2_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA2_Om01.append(float(table2_Om01[i][0].replace(',','.')))

    ARPES_add2_Om01.append(float(table2_Om01[i][1].replace(',','.')))
    
    ARPES_red2_Om01.append(float(table2_Om01[i][2].replace(',','.')))
    
#Pointer initialization

table3_Om01 = []

#Data initialization

OMEGA3_Om01 = []
ARPES_add3_Om01 = []
ARPES_red3_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP02.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table3_Om01.append(row)


#Size of the imported array

N = len(table3_Om01)

M = len(table3_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA3_Om01.append(float(table3_Om01[i][0].replace(',','.')))

    ARPES_add3_Om01.append(float(table3_Om01[i][1].replace(',','.')))
    
    ARPES_red3_Om01.append(float(table3_Om01[i][2].replace(',','.')))
 
#Pointer initialization

table4_Om01 = []

#Data initialization

OMEGA4_Om01 = []
ARPES_add4_Om01 = []
ARPES_red4_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP03.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table4_Om01.append(row)


#Size of the imported array

N = len(table4_Om01)

M = len(table4_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA4_Om01.append(float(table4_Om01[i][0].replace(',','.')))

    ARPES_add4_Om01.append(float(table4_Om01[i][1].replace(',','.')))
    
    ARPES_red4_Om01.append(float(table4_Om01[i][2].replace(',','.')))
 

#Pointer initialization

table5_Om01 = []

#Data initialization

OMEGA5_Om01 = []

ARPES_add5_Om01 = []
ARPES_red5_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP04.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table5_Om01.append(row)


#Size of the imported array

N = len(table5_Om01)

M = len(table5_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA5_Om01.append(float(table5_Om01[i][0].replace(',','.')))

    ARPES_add5_Om01.append(float(table5_Om01[i][1].replace(',','.')))
    
    ARPES_red5_Om01.append(float(table5_Om01[i][2].replace(',','.')))
    
#Pointer initialization

table6_Om01 = []

#Data initialization

OMEGA6_Om01 = []

ARPES_add6_Om01 = []
ARPES_red6_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP05.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table6_Om01.append(row)


#Size of the imported array

N = len(table6_Om01)

M = len(table6_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA6_Om01.append(float(table6_Om01[i][0].replace(',','.')))

    ARPES_add6_Om01.append(float(table6_Om01[i][1].replace(',','.')))
    
    ARPES_red6_Om01.append(float(table6_Om01[i][2].replace(',','.')))
    

 #Pointer initialization

table7_Om01 = []

#Data initialization

OMEGA7_Om01 = []
ARPES_add7_Om01 = []
ARPES_red7_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP06.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table7_Om01.append(row)


#Size of the imported array

N = len(table7_Om01)

M = len(table7_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA7_Om01.append(float(table7_Om01[i][0].replace(',','.')))

    ARPES_add7_Om01.append(float(table7_Om01[i][1].replace(',','.')))
    
    ARPES_red7_Om01.append(float(table7_Om01[i][2].replace(',','.')))
 
#Pointer initialization

table8_Om01 = []

#Data initialization

OMEGA8_Om01 = []
ARPES_add8_Om01 = []
ARPES_red8_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP07.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table8_Om01.append(row)


#Size of the imported array

N = len(table8_Om01)

M = len(table8_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA8_Om01.append(float(table8_Om01[i][0].replace(',','.')))

    ARPES_add8_Om01.append(float(table8_Om01[i][1].replace(',','.')))
    
    ARPES_red8_Om01.append(float(table8_Om01[i][2].replace(',','.')))
    
#Pointer initialization

table9_Om01 = []

#Data initialization

ARPES_red9_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP08.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table9_Om01.append(row)


#Size of the imported array

N = len(table9_Om01)

M = len(table9_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red9_Om01.append(float(table9_Om01[i][2].replace(',','.')))
    

#Pointer initialization

table10_Om01 = []

#Data initialization


ARPES_red10_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP09.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table10_Om01.append(row)


#Size of the imported array

N = len(table10_Om01)

M = len(table10_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red10_Om01.append(float(table10_Om01[i][2].replace(',','.')))
    

#Pointer initialization

table11_Om01 = []

#Data initialization


ARPES_red11_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP10.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table11_Om01.append(row)


#Size of the imported array

N = len(table11_Om01)

M = len(table11_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red11_Om01.append(float(table11_Om01[i][2].replace(',','.')))
    
#Pointer initialization

table12_Om01 = []

#Data initialization


ARPES_red12_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP11.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table12_Om01.append(row)


#Size of the imported array

N = len(table12_Om01)

M = len(table12_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red12_Om01.append(float(table12_Om01[i][2].replace(',','.')))
    
#Pointer initialization

table13_Om01 = []

#Data initialization

ARPES_red13_Om01 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP12.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table13_Om01.append(row)


#Size of the imported array

N = len(table13_Om01)

M = len(table13_Om01[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red13_Om01.append(float(table13_Om01[i][2].replace(',','.')))
    
###########
###########

#Pointer initialization

table1_Om01_05 = []

#Data initialization

OMEGA1_Om01_05 = []
ARPES_add1_Om01_05 = []
ARPES_red1_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP005.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3

with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table1_Om01_05.append(row)


#Size of the imported array

N = len(table1_Om01_05)

M = len(table1_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA1_Om01_05.append(float(table1_Om01_05[i][0].replace(',','.')))
    
    ARPES_red1_Om01_05.append(float(table1_Om01_05[i][2].replace(',','.')))
    

#Pointer initialization

table2_Om01_05 = []

#Data initialization

OMEGA2_Om01_05= []
ARPES_add2_Om01_05 = []
ARPES_red2_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP015.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table2_Om01_05.append(row)


#Size of the imported array

N = len(table2_Om01_05)

M = len(table2_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA2_Om01_05.append(float(table2_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red2_Om01_05.append(float(table2_Om01_05[i][2].replace(',','.')))
    
#Pointer initialization

table3_Om01_05 = []

#Data initialization

OMEGA3_Om01_05 = []
ARPES_add3_Om01_05 = []
ARPES_red3_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP025.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table3_Om01_05.append(row)


#Size of the imported array

N = len(table3_Om01_05)

M = len(table3_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA3_Om01_05.append(float(table3_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red3_Om01_05.append(float(table3_Om01_05[i][2].replace(',','.')))
 
#Pointer initialization

table4_Om01_05 = []

#Data initialization

OMEGA4_Om01_05 = []
ARPES_add4_Om01_05 = []
ARPES_red4_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP035.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table4_Om01_05.append(row)


#Size of the imported array

N = len(table4_Om01_05)

M = len(table4_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA4_Om01_05.append(float(table4_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red4_Om01_05.append(float(table4_Om01_05[i][2].replace(',','.')))
 

#Pointer initialization

table5_Om01_05 = []

#Data initialization

OMEGA5_Om01_05 = []

ARPES_add5_Om01_05 = []
ARPES_red5_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP045.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table5_Om01_05.append(row)


#Size of the imported array

N = len(table5_Om01_05)

M = len(table5_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA5_Om01_05.append(float(table5_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red5_Om01_05.append(float(table5_Om01_05[i][2].replace(',','.')))
    
#Pointer initialization

table6_Om01_05 = []

#Data initialization

OMEGA6_Om01_05 = []

ARPES_add6_Om01_05 = []
ARPES_red6_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP055.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table6_Om01_05.append(row)


#Size of the imported array

N = len(table6_Om01_05)

M = len(table6_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA6_Om01_05.append(float(table6_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red6_Om01_05.append(float(table6_Om01_05[i][2].replace(',','.')))
    

 #Pointer initialization

table7_Om01_05 = []

#Data initialization

OMEGA7_Om01_05 = []
ARPES_add7_Om01_05 = []
ARPES_red7_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP065.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table7_Om01_05.append(row)


#Size of the imported array

N = len(table7_Om01_05)

M = len(table7_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA7_Om01_05.append(float(table7_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red7_Om01_05.append(float(table7_Om01_05[i][2].replace(',','.')))
 
#Pointer initialization

table8_Om01_05 = []

#Data initialization

OMEGA8_Om01_05 = []
ARPES_add8_Om01_05 = []
ARPES_red8_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP075.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table8_Om01_05.append(row)


#Size of the imported array

N = len(table8_Om01_05)

M = len(table8_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA8_Om01_05.append(float(table8_Om01_05[i][0].replace(',','.')))

    
    
    ARPES_red8_Om01_05.append(float(table8_Om01_05[i][2].replace(',','.')))
    
#Pointer initialization

table9_Om01_05 = []

#Data initialization

ARPES_red9_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP085.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table9_Om01_05.append(row)


#Size of the imported array

N = len(table9_Om01_05)

M = len(table9_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red9_Om01_05.append(float(table9_Om01_05[i][2].replace(',','.')))
    

#Pointer initialization

table10_Om01_05 = []

#Data initialization


ARPES_red10_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP095.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table10_Om01_05.append(row)


#Size of the imported array

N = len(table10_Om01_05)

M = len(table10_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red10_Om01_05.append(float(table10_Om01_05[i][2].replace(',','.')))
    

#Pointer initialization

table11_Om01_05 = []

#Data initialization


ARPES_red11_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP105.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table11_Om01_05.append(row)


#Size of the imported array

N = len(table11_Om01_05)

M = len(table11_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red11_Om01_05.append(float(table11_Om01_05[i][2].replace(',','.')))
    
#Pointer initialization

table12_Om01_05 = []

#Data initialization


ARPES_red12_Om01_05 = []

rel_path = "DATA/J02_samefreq_Om04/ARPES_PHS_J02_samefreq_WP115.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table12_Om01_05.append(row)


#Size of the imported array

N = len(table12_Om01_05)

M = len(table12_Om01_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red12_Om01_05.append(float(table12_Om01_05[i][2].replace(',','.')))
    
    
    
##########
#########
    

#Pointer initialization

table16_Om02 = []

#Data initialization

Wfreq_Om02 = []
WPLASMA_PM_Om02 = []
WPLASMA_2M_Om02 = []
WPLASMA_2P_Om02 = []
row_number=9

rel_path = "DATA/Contrib/Contrib_Om02_V2.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:

    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table16_Om02.append(row)


#Size of the imported array

N = len(table16_Om02)

M = len(table16_Om02[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02.append(float(table16_Om02[i][0].replace(',','.')))
    WPLASMA_PM_Om02.append(float(table16_Om02[i][6].replace(',','.')))
    WPLASMA_2P_Om02.append(float(table16_Om02[i][7].replace(',','.')))
    WPLASMA_2M_Om02.append(float(table16_Om02[i][8].replace(',','.')))
    

    
############

ARPES_tot_Om01 = np.zeros((1601,25))

for i in range(0,1601):
    
    ARPES_tot_Om01[i][0]=ARPES_red1_Om01[i]
    ARPES_tot_Om01[i][1]=ARPES_red1_Om01_05[i]
    ARPES_tot_Om01[i][2]=ARPES_red2_Om01[i]
    ARPES_tot_Om01[i][3]=ARPES_red2_Om01_05[i]
    ARPES_tot_Om01[i][4]=ARPES_red3_Om01[i]
    ARPES_tot_Om01[i][5]=ARPES_red3_Om01_05[i]
    ARPES_tot_Om01[i][6]=ARPES_red4_Om01[i]
    ARPES_tot_Om01[i][7]=ARPES_red4_Om01_05[i]
    ARPES_tot_Om01[i][8]=ARPES_red5_Om01[i]
    ARPES_tot_Om01[i][9]=ARPES_red5_Om01_05[i]
    ARPES_tot_Om01[i][10]=ARPES_red6_Om01[i]
    ARPES_tot_Om01[i][11]=ARPES_red6_Om01_05[i]
    ARPES_tot_Om01[i][12]=ARPES_red7_Om01[i]
    ARPES_tot_Om01[i][13]=ARPES_red7_Om01_05[i]
    ARPES_tot_Om01[i][14]=ARPES_red8_Om01[i]
    ARPES_tot_Om01[i][15]=ARPES_red8_Om01_05[i]
    ARPES_tot_Om01[i][16]=ARPES_red9_Om01[i]
    ARPES_tot_Om01[i][17]=ARPES_red9_Om01_05[i]
    ARPES_tot_Om01[i][18]=ARPES_red10_Om01[i]
    ARPES_tot_Om01[i][19]=ARPES_red10_Om01_05[i]
    ARPES_tot_Om01[i][20]=ARPES_red11_Om01[i]
    ARPES_tot_Om01[i][21]=ARPES_red11_Om01_05[i]
    ARPES_tot_Om01[i][22]=ARPES_red12_Om01[i]
    ARPES_tot_Om01[i][23]=ARPES_red12_Om01_05[i]
    ARPES_tot_Om01[i][24]=ARPES_red13_Om01[i]
        
    
    
#Pointer initialization

table17_Om02 = []

#Data initialization

Wfreq_Om02_1peak = []
WPLASMA_PM_Om02_1peak = []
WPLASMA_2M_Om02_1peak = []
WPLASMA_2P_Om02_1peak = []

row_number=9

rel_path = "DATA/Contrib/Contrib_Om02_V2.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table17_Om02.append(row)


#Size of the imported array

N = len(table17_Om02)

M = len(table17_Om02[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_1peak.append(float(table17_Om02[i][0].replace(',','.')))
    WPLASMA_PM_Om02_1peak.append(float(table17_Om02[i][6].replace(',','.')))
    WPLASMA_2P_Om02_1peak.append(float(table17_Om02[i][7].replace(',','.')))
    WPLASMA_2M_Om02_1peak.append(float(table17_Om02[i][8].replace(',','.')))





########################
########################

    
    
    
#Pointer initialization

table1_Om01_U4J = []

#Data initialization

OMEGA1_Om01_U4J = []
ARPES_add1_Om01_U4J = []
ARPES_red1_Om01_U4J = []



row_number=3

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP00.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table1_Om01_U4J.append(row)


#Size of the imported array

N = len(table1_Om01_U4J)

M = len(table1_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red1_Om01_U4J.append(float(table1_Om01_U4J[i][2].replace(',','.')))
    

#Pointer initialization

table2_Om01_U4J = []

#Data initialization

OMEGA2_Om01_U4J= []
ARPES_add2_Om01_U4J = []
ARPES_red2_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP01.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table2_Om01_U4J.append(row)


#Size of the imported array

N = len(table2_Om01_U4J)

M = len(table2_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red2_Om01_U4J.append(float(table2_Om01_U4J[i][2].replace(',','.')))
    
#Pointer initialization

table3_Om01_U4J = []

#Data initialization

ARPES_red3_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP02.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table3_Om01_U4J.append(row)


#Size of the imported array

N = len(table3_Om01_U4J)

M = len(table3_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red3_Om01_U4J.append(float(table3_Om01_U4J[i][2].replace(',','.')))
 
#Pointer initialization

table4_Om01_U4J = []

#Data initialization


ARPES_red4_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP03.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table4_Om01_U4J.append(row)


#Size of the imported array

N = len(table4_Om01_U4J)

M = len(table4_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red4_Om01_U4J.append(float(table4_Om01_U4J[i][2].replace(',','.')))
 

#Pointer initialization

table5_Om01_U4J = []

#Data initialization


ARPES_red5_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP04.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table5_Om01_U4J.append(row)


#Size of the imported array

N = len(table5_Om01_U4J)

M = len(table5_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red5_Om01_U4J.append(float(table5_Om01_U4J[i][2].replace(',','.')))
    
#Pointer initialization

table6_Om01_U4J = []

#Data initialization


ARPES_red6_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP05.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table6_Om01_U4J.append(row)


#Size of the imported array

N = len(table6_Om01_U4J)

M = len(table6_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red6_Om01_U4J.append(float(table6_Om01_U4J[i][2].replace(',','.')))
    

 #Pointer initialization

table7_Om01_U4J = []

#Data initialization

ARPES_red7_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP06.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table7_Om01_U4J.append(row)


#Size of the imported array

N = len(table7_Om01_U4J)

M = len(table7_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red7_Om01_U4J.append(float(table7_Om01_U4J[i][2].replace(',','.')))
 
#Pointer initialization

table8_Om01_U4J = []

#Data initialization

ARPES_red8_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP07.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:
    
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table8_Om01_U4J.append(row)


#Size of the imported array

N = len(table8_Om01_U4J)

M = len(table8_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red8_Om01_U4J.append(float(table8_Om01_U4J[i][2].replace(',','.')))
    
#Pointer initialization

table9_Om01_U4J = []

#Data initialization

ARPES_red9_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP08.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table9_Om01_U4J.append(row)


#Size of the imported array

N = len(table9_Om01_U4J)

M = len(table9_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red9_Om01_U4J.append(float(table9_Om01_U4J[i][2].replace(',','.')))
    

#Pointer initialization

table10_Om01_U4J = []

#Data initialization


ARPES_red10_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP09.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table10_Om01_U4J.append(row)


#Size of the imported array

N = len(table10_Om01_U4J)

M = len(table10_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red10_Om01_U4J.append(float(table10_Om01_U4J[i][2].replace(',','.')))
    

#Pointer initialization

table11_Om01_U4J = []

#Data initialization


ARPES_red11_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP10.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table11_Om01_U4J.append(row)


#Size of the imported array

N = len(table11_Om01_U4J)

M = len(table11_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red11_Om01_U4J.append(float(table11_Om01_U4J[i][2].replace(',','.')))
    
#Pointer initialization

table12_Om01_U4J = []

#Data initialization


ARPES_red12_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP11.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table12_Om01_U4J.append(row)


#Size of the imported array

N = len(table12_Om01_U4J)

M = len(table12_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red12_Om01_U4J.append(float(table12_Om01_U4J[i][2].replace(',','.')))
    
#Pointer initialization

table13_Om01_U4J = []

#Data initialization

ARPES_red13_Om01_U4J = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP12.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table13_Om01_U4J.append(row)


#Size of the imported array

N = len(table13_Om01_U4J)

M = len(table13_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red13_Om01_U4J.append(float(table13_Om01_U4J[i][2].replace(',','.')))
    


###########
##########

#Pointer initialization

table1_Om01_U4J_05 = []

#Data initialization

OMEGA1_Om01_U4J_05 = []
ARPES_add1_Om01_U4J_05 = []
ARPES_red1_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP005.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3

with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table1_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table1_Om01_U4J_05)

M = len(table1_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red1_Om01_U4J_05.append(float(table1_Om01_U4J_05[i][2].replace(',','.')))
    

#Pointer initialization

table2_Om01_U4J_05 = []

#Data initialization

OMEGA2_Om01_U4J_05= []
ARPES_add2_Om01_U4J_05 = []
ARPES_red2_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP015.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table2_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table2_Om01_U4J_05)

M = len(table2_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red2_Om01_U4J_05.append(float(table2_Om01_U4J_05[i][2].replace(',','.')))
    
#Pointer initialization

table3_Om01_U4J_05 = []

#Data initialization

ARPES_red3_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP025.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table3_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table3_Om01_U4J_05)

M = len(table3_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red3_Om01_U4J_05.append(float(table3_Om01_U4J_05[i][2].replace(',','.')))
 
#Pointer initialization

table4_Om01_U4J_05 = []

#Data initialization


ARPES_red4_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP035.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:




    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table4_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table4_Om01_U4J_05)

M = len(table4_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red4_Om01_U4J_05.append(float(table4_Om01_U4J_05[i][2].replace(',','.')))
 

#Pointer initialization

table5_Om01_U4J_05 = []

#Data initialization


ARPES_red5_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP045.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table5_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table5_Om01_U4J_05)

M = len(table5_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red5_Om01_U4J_05.append(float(table5_Om01_U4J_05[i][2].replace(',','.')))
    
#Pointer initialization

table6_Om01_U4J_05 = []

#Data initialization


ARPES_red6_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP055.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:















    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table6_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table6_Om01_U4J)

M = len(table6_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red6_Om01_U4J_05.append(float(table6_Om01_U4J_05[i][2].replace(',','.')))
    

 #Pointer initialization

table7_Om01_U4J_05 = []

#Data initialization

ARPES_red7_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP065.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:











    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table7_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table7_Om01_U4J)

M = len(table7_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red7_Om01_U4J_05.append(float(table7_Om01_U4J_05[i][2].replace(',','.')))
 
#Pointer initialization

table8_Om01_U4J_05 = []

#Data initialization

ARPES_red8_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP075.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:










    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table8_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table8_Om01_U4J)

M = len(table8_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red8_Om01_U4J_05.append(float(table8_Om01_U4J_05[i][2].replace(',','.')))
    
#Pointer initialization

table9_Om01_U4J_05 = []

#Data initialization

ARPES_red9_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP085.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:










    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table9_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table9_Om01_U4J)

M = len(table9_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):



    
    ARPES_red9_Om01_U4J_05.append(float(table9_Om01_U4J_05[i][2].replace(',','.')))
    

#Pointer initialization

table10_Om01_U4J_05 = []

#Data initialization


ARPES_red10_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP095.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:










    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table10_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table10_Om01_U4J)

M = len(table10_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    
    ARPES_red10_Om01_U4J_05.append(float(table10_Om01_U4J_05[i][2].replace(',','.')))
    

#Pointer initialization

table11_Om01_U4J_05 = []

#Data initialization


ARPES_red11_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP105.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:










    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table11_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table11_Om01_U4J_05)

M = len(table11_Om01_U4J_05[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    
    ARPES_red11_Om01_U4J_05.append(float(table11_Om01_U4J_05[i][2].replace(',','.')))
    
#Pointer initialization

table12_Om01_U4J_05 = []

#Data initialization


ARPES_red12_Om01_U4J_05 = []

rel_path = "DATA/J0_samefreq_Om04/ARPES_PHS_J00_samefreq_WP115.csv"
abs_file_path = os.path.join(script_dir, rel_path)


row_number=3
with open (abs_file_path, newline='') as csvfile:










    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table12_Om01_U4J_05.append(row)


#Size of the imported array

N = len(table12_Om01_U4J)

M = len(table12_Om01_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):
    
    ARPES_red12_Om01_U4J_05.append(float(table12_Om01_U4J_05[i][2].replace(',','.')))
    
    
    
    
ARPES_tot_Om01_U4J = np.zeros((1601,25))

for i in range(0,1601):
    
    ARPES_tot_Om01_U4J[i][0]=ARPES_red1_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][1]=ARPES_red1_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][2]=ARPES_red2_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][3]=ARPES_red2_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][4]=ARPES_red3_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][5]=ARPES_red3_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][6]=ARPES_red4_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][7]=ARPES_red4_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][8]=ARPES_red5_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][9]=ARPES_red5_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][10]=ARPES_red6_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][11]=ARPES_red6_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][12]=ARPES_red7_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][13]=ARPES_red7_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][14]=ARPES_red8_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][15]=ARPES_red8_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][16]=ARPES_red9_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][17]=ARPES_red9_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][18]=ARPES_red10_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][19]=ARPES_red10_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][20]=ARPES_red11_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][21]=ARPES_red11_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][22]=ARPES_red12_Om01_U4J[i]
    ARPES_tot_Om01_U4J[i][23]=ARPES_red12_Om01_U4J_05[i]
    ARPES_tot_Om01_U4J[i][24]=ARPES_red13_Om01_U4J[i]
    

#Pointer initialization

table16_Om02 = []

#Data initialization


Wfreq_Om02 = []
WPLASMA_PM_Om02 = []
WPLASMA_2P_Om02 = []
WPLASMA_2M_Om02 = []
row_number=7

rel_path = "DATA/Contrib/Contrib_U5J_2peak_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)


with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table16_Om02.append(row)


#Size of the imported array

N = len(table16_Om02)

M = len(table16_Om02[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02.append(float(table16_Om02[i][6].replace(',','.')))
    WPLASMA_PM_Om02.append(float(table16_Om02[i][3].replace(',','.')))
    WPLASMA_2P_Om02.append(float(table16_Om02[i][4].replace(',','.')))
    WPLASMA_2M_Om02.append(float(table16_Om02[i][5].replace(',','.')))
    
#Pointer initialization

table16_Om02_WP = []

#Data initialization


Wfreq_Om02_WP = []
WPLASMA_PM_Om02_WP = []
WPLASMA_2P_Om02_WP = []
row_number=7

rel_path = "DATA/Contrib/Contrib_U5J_2peak_WP_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table16_Om02_WP.append(row)


#Size of the imported array

N = len(table16_Om02_WP)

M = len(table16_Om02_WP[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_WP.append(float(table16_Om02_WP[i][6].replace(',','.')))
    WPLASMA_PM_Om02_WP.append(float(table16_Om02_WP[i][3].replace(',','.')))
    WPLASMA_2P_Om02_WP.append(float(table16_Om02_WP[i][4].replace(',','.')))
    
    
    
    
#Pointer initialization

table17_Om02 = []

#Data initialization

Wfreq_Om02_1peak = []
WPLASMA_PM_Om02_1peak = []
WPLASMA_2P_Om02_1peak = []
WPLASMA_2M_Om02_1peak = []

row_number=7

rel_path = "DATA/Contrib/Contrib_J0_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table17_Om02.append(row)


#Size of the imported array

N = len(table17_Om02)

M = len(table17_Om02[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_1peak.append(float(table17_Om02[i][6].replace(',','.')))
    WPLASMA_PM_Om02_1peak.append(float(table17_Om02[i][3].replace(',','.')))
    WPLASMA_2P_Om02_1peak.append(float(table17_Om02[i][4].replace(',','.')))
    WPLASMA_2M_Om02_1peak.append(float(table17_Om02[i][5].replace(',','.')))
    

#Pointer initialization

table17_Om02_WP = []

#Data initialization

Wfreq_Om02_1peak_WP = []
WPLASMA_PM_Om02_1peak_WP = []
WPLASMA_2P_Om02_1peak_WP = []

row_number=7

rel_path = "DATA/Contrib/Contrib_J0_WP_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table17_Om02_WP.append(row)


#Size of the imported array

N = len(table17_Om02_WP)

M = len(table17_Om02_WP[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_1peak_WP.append(float(table17_Om02_WP[i][6].replace(',','.')))
    WPLASMA_2P_Om02_1peak_WP.append(float(table17_Om02_WP[i][4].replace(',','.')))
    
    
    
#Pointer initialization

table17_Om02_U4J = []

#Data initialization

Wfreq_Om02_1peak_U4J = []
WPLASMA_PM_Om02_1peak_U4J = []
WPLASMA_2P_Om02_1peak_U4J = []
WPLASMA_2M_Om02_1peak_U4J = []

row_number=7

rel_path = "DATA/Contrib/Contrib_U5J_1peak_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table17_Om02_U4J.append(row)


#Size of the imported array

N = len(table17_Om02_U4J)

M = len(table17_Om02_U4J[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_1peak_U4J.append(float(table17_Om02_U4J[i][6].replace(',','.')))
    WPLASMA_PM_Om02_1peak_U4J.append(float(table17_Om02_U4J[i][3].replace(',','.')))
    WPLASMA_2P_Om02_1peak_U4J.append(float(table17_Om02_U4J[i][4].replace(',','.')))
    WPLASMA_2M_Om02_1peak_U4J.append(float(table17_Om02_U4J[i][5].replace(',','.')))

    
    #Pointer initialization

table17_Om02_U4J_WP = []

#Data initialization

Wfreq_Om02_1peak_U4J_WP = []
WPLASMA_PM_Om02_1peak_U4J_WP = []
WPLASMA_2P_Om02_1peak_U4J_WP = []


row_number=7

rel_path = "DATA/Contrib/Contrib_U5J_1peak_WP_samefreq.csv"
abs_file_path = os.path.join(script_dir, rel_path)

with open (abs_file_path, newline='') as csvfile:


    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table17_Om02_U4J_WP.append(row)


#Size of the imported array

N = len(table17_Om02_U4J_WP)

M = len(table17_Om02_U4J_WP[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()


#DATA importation

for i in range(1,N):
    
    Wfreq_Om02_1peak_U4J_WP.append(float(table17_Om02_U4J_WP[i][6].replace(',','.')))
    WPLASMA_PM_Om02_1peak_U4J_WP.append(float(table17_Om02_U4J_WP[i][3].replace(',','.')))
    WPLASMA_2P_Om02_1peak_U4J_WP.append(float(table17_Om02_U4J_WP[i][4].replace(',','.')))
    
############

#Pointer initialization

table_Hub4 = []

#Data initialization

OMEGA_Hub4 = []
ARPES_red_Hub4 = []

rel_path = "DATA/Hubbard/ARPES_PHS_hubbard_J02_WP00_SIG8_NT24_V2.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3
with open (abs_file_path, newline='') as csvfile:





    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table_Hub4.append(row)


#Size of the imported array

N = len(table_Hub4)

M = len(table_Hub4[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):

    OMEGA_Hub4.append(float(table_Hub4[i][0].replace(',','.')))
    
    ARPES_red_Hub4.append(float(table_Hub4[i][2].replace(',','.')))
    
    

#Pointer initialization

table_Hub0 = []

#Data initialization

OMEGA_Hub0 = []
ARPES_red_Hub0 = []

rel_path = "DATA/Hubbard/ARPES_PHS_hubbard_J00_WP00_SIG8_NT24_V2.csv"
abs_file_path = os.path.join(script_dir, rel_path)

row_number=3

with open (abs_file_path, newline='') as csvfile:



    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    for row in spamreader:

        table_Hub0.append(row)


#Size of the imported array

N = len(table_Hub0)

M = len(table_Hub0[0]) 


#Check of the array structur

if M != row_number :

    print("Problem number of colum")

    sys.exit()



#DATA importation

for i in range(1,N):


    OMEGA_Hub0.append(float(table_Hub0[i][0].replace(',','.')))
    
    ARPES_red_Hub0.append(float(table_Hub0[i][2].replace(',','.')))
    


###############
fontsize = 14

mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['font.size'] = 8  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 5.5
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['xtick.major.width'] = .7
mpl.rcParams['ytick.major.width'] = .7
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['figure.titlesize'] = 8


mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage[helvet]{sfmath}',
   
]

    
fig = plt.figure()

gs = mpl.gridspec.GridSpec(2, 5)

ax1 = fig.add_subplot(gs[0, 2:6])
ax2 = fig.add_subplot(gs[1, 2:6])
ax3 = fig.add_subplot(gs[0, 0:2])
ax4 = fig.add_subplot(gs[1, 0:2])
axs = [ax1, ax2, ax3, ax4]


fig.set_size_inches(6., 5.)

plot_vals= np.arange(0,6.0001,1)
plot_wp = []
omega_values = []

for i in range(0,1601):
    omega_values = np.append(omega_values,[(-15/1600)*i])
    
omega_values = np.flip(omega_values)
for i in range (0,25):
        plot_wp = np.append(plot_wp, [(3/24)*i])


plot_omega = [-15, -10.5, -6.5, -2.5, 0]

boxProps = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.15)
boxProps2 = dict(boxstyle='square', facecolor='white', alpha=1., linewidth=0., fill=True, pad=0.3)


sm = plt.cm.ScalarMappable(cmap='terrain', norm=LogNorm(vmin=0.00001, vmax=2.5))
CS = axs[0].pcolormesh(plot_wp, omega_values , ARPES_tot_Om01_U4J, cmap='terrain', norm=LogNorm(), shading = 'terrain')

axs[0].plot(Wfreq_Om02_1peak, WPLASMA_PM_Om02_1peak,'red')
axs[0].plot(Wfreq_Om02, WPLASMA_2M_Om02_1peak,'red')
axs[0].plot(Wfreq_Om02_1peak_WP, WPLASMA_2P_Om02_1peak_WP,'red')

axs[0].set_xlabel(r"light-matter coupling $[\omega_P \, / \, \omega_{\rm{phot}}]$", fontsize=10)

axs[0].set_yticks(plot_omega)
axs[0].xaxis.set_visible(False)
axs[0].yaxis.set_visible(False)


CS = axs[1].pcolormesh(plot_wp, omega_values , ARPES_tot_Om01, cmap='terrain', norm=LogNorm(), shading = 'terrain')


axs[1].set_xlabel(r"light-matter coupling $[\omega_{\rm P} \, / \, \omega_{\rm{phot}}]$", fontsize=10)
axs[1].set_xticks(plot_vals)
axs[1].set_ylabel(r"$\omega\left[J\right]$", fontsize=10)
axs[1].set_yticks(plot_omega)
axs[1].yaxis.set_visible(False)

axs[1].tick_params(axis = 'x', labelsize = 10)



axs[1].plot(Wfreq_Om02, WPLASMA_PM_Om02,'red',linestyle='--')
axs[1].plot(Wfreq_Om02_WP, WPLASMA_2P_Om02_WP,'red',linestyle='--')
axs[1].plot(Wfreq_Om02, WPLASMA_2M_Om02,'red',linestyle='--')

axs[1].plot(Wfreq_Om02_1peak_U4J, WPLASMA_PM_Om02_1peak_U4J,'red')
axs[1].plot(Wfreq_Om02, WPLASMA_2M_Om02_1peak_U4J,'red')
axs[1].plot(Wfreq_Om02_1peak_U4J_WP, WPLASMA_2P_Om02_1peak_U4J_WP,'red')





## Set dashed lines

axs[0].axhline(y=-2.409375, color='black',linewidth='1.5')


axs[1].axhline(y=-2.128125, color='black',linewidth='1.5')


axs[1].axhline(y=-4.125,linestyle='dashed', color='black',linewidth='1.5')


####################
####################


axs[2].plot(ARPES_red1_Om01_U4J,omega_values,'red', color='indigo')
axs[2].fill_betweenx(OMEGA_Hub0,ARPES_red_Hub0,1e-10,color='silver', label="uncoupled \n Hubbard")
axs[2].set_xscale('log')
axs[2].xaxis.set_visible(False)
axs[2].set_ylabel(r"$\omega \, \left[J\right]$", fontsize=10)
axs[2].set_ylim(-15,0)
axs[2].set_xlim(1e1,9.99e-8)
axs[2].set_yticks(plot_omega)

legend = axs[2].legend(fontsize = 10, loc = 'lower left', bbox_to_anchor=(.1, 0.), edgecolor = 'black', ncol = 1)
legend.get_frame().set_alpha(0.)
legend.get_frame().set_boxstyle('Square', pad=0.1)
legend.get_frame().set_linewidth(0.0)


axs[2].text(5000, -0.5, r"$\rm{a.)}$",fontweight='bold', fontsize=10, color='black', alpha=1., bbox=boxProps)
axs[3].text(5000, -0.5, r"$\rm{b.)}$",fontweight='bold', fontsize=10, color='black', alpha=1., bbox=boxProps)

axs[0].text(3.1, -0.5, r"$\rm{c.)}$",fontweight='bold', fontsize=10, color='black', alpha=1., bbox=boxProps)
axs[1].text(3.1, -0.5, r"$\rm{d.)}$",fontweight='bold', fontsize=10, color='black', alpha=1., bbox=boxProps)


axs[2].tick_params(axis = 'y', labelsize = 10)
axs[2].tick_params(axis = 'x', labelsize = 10)
axs[0].text(3.2, -7.5,r'$J=0$',rotation=90, fontsize = 10)
axs[1].text(3.2, -7.5,r'$U=5 \, J$',rotation=90, fontsize = 10)


##############


rect1 = mpl.patches.Rectangle((0.365, -3.45),
                                     0.35, 0.6,
                                     color ='white')

rect2 = mpl.patches.Rectangle((0.17, -3.2),
                                     0.35, 0.6,
                                     color ='white')
axs[0].add_patch(rect1)
axs[1].add_patch(rect2)

###############


axs[0].annotate('', xy=(0.25,-2.4), xytext=(0.25, -5.9),arrowprops=dict(arrowstyle="<->", color='black'))
axs[0].annotate('', xy=(1.5,-2.4), xytext=(1.5, -7.5),arrowprops=dict(arrowstyle="<->", color='black'))
axs[0].annotate('', xy=(2.4,-2.4), xytext=(2.4, -13.3),arrowprops=dict(arrowstyle="<->", color='black'))

axs[0].text(1.6, -5.8,r'$ \omega_+ + \omega_-$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[0].text(2.5, -12,r'$ 2 \, \omega_+$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[0].text(0.35, -3.5,r'$ 2 \, \omega$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[0].text(0.63, -3.63,r'$ -$', fontsize = 6, va='center',alpha=1., bbox=boxProps2)


axs[1].annotate('', xy=(0.25,-4), xytext=(0.25, -7.7),arrowprops=dict(arrowstyle="<->", color='black',linestyle= 'dashed'))
axs[1].annotate('', xy=(1.7,-4), xytext=(1.7, -9.1),arrowprops=dict(arrowstyle="<->", color='black',linestyle= 'dashed'))
axs[1].annotate('', xy=(2.4,-4), xytext=(2.4, -15),arrowprops=dict(arrowstyle="<->", color='black',linestyle= 'dashed'))


axs[1].text(1.35, -3.2,r'$ \omega_+ + \omega_-$', fontsize = 10,  va='center',alpha=1., bbox=boxProps)
axs[1].text(2.5, -12.3,r'$ 2 \, \omega_+$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[1].text(0.15, -3.2,r'$ 2 \, \omega $', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[1].text(0.43, -3.33,r'$ -$', fontsize = 6, va='center',alpha=1., bbox=boxProps2)


axs[2].text(2.5, -5.2,r'$ 2 \, \omega_{\rm phon}$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[2].text(3e-3, -9.5,r'$ 4 \, \omega_{\rm phon}$', fontsize = 10, va='center',alpha=1., bbox=boxProps)

axs[3].text(9e-4, -4.5,r'$ 2 \, \omega_{\rm phon}$', fontsize = 10, va='center',alpha=1., bbox=boxProps)
axs[3].text(3, -7.,r'$ 2 \, \omega_{\rm phon}$', fontsize = 10, va='center',alpha=1., bbox=boxProps)


#############
############
 

axs[3].set_xlabel(r"$A(\omega)$", fontsize=10)
axs[3].set_xscale('log')
axs[3].set_ylim(-15,-0)
axs[3].set_xlim(1e1,9.99e-8)
axs[3].plot(ARPES_red1_Om01,omega_values,'blue',label=r"$U = 5 \, J$", color='indigo')
axs[3].fill_betweenx(OMEGA_Hub4,ARPES_red_Hub4,1e-10,color='silver')
axs[3].set_ylabel(r"$\omega \, \left[J\right]$", fontsize=10)
axs[3].tick_params(axis = 'both', labelsize = 10)
axs[3].set_yticks(plot_omega)

axs[2].annotate('', xy=(1e-2,-2.4), xytext=(1e-2, -6.5),arrowprops=dict(arrowstyle="<->", color='black'))
axs[2].annotate('', xy=(1e-5,-2.4), xytext=(1e-5, -10.5),arrowprops=dict(arrowstyle="<->", color='black'))

axs[3].annotate('', xy=(3e-3,-2.2), xytext=(3e-3, -6.8),arrowprops=dict(arrowstyle="<->", color='black'))
axs[3].annotate('', xy=(0.01,-4), xytext=(0.01, -8.5),arrowprops=dict(arrowstyle="<->", color='black'))


###################
###################

plt.draw()
p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()



ax_cbar = fig.add_axes([ 1.02, p0[0]-0.3, 0.03, p1[2]-p0[0]+0.3 ])

fig.colorbar(sm, cax=ax_cbar, orientation='vertical')
ax_cbar.text(0.001, 4,r'$A(\omega)$', fontsize = 10)


##########
##########


fig.tight_layout()
fig.subplots_adjust(wspace=0., hspace=0.1)
#plt.savefig('False_color.png', format='png', bbox_inches='tight', dpi = 600)
plt.show()


    

    
    



