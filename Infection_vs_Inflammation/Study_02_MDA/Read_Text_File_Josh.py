import numpy as np

file_name='Goldenberg_M5_S5_2016_12_27_T2exDCE.txt'

f1 = open(file_name,'r') 
lines=f1.readlines()

for original_lines in lines:
    new_lines=original_lines.split()
    print(new_lines[1::3])


myarray=np.random.rand(1,10)
num_el=10
A=np.zeros_like(num_el,)

for rows in A:
    print(rows)

#B=np.zeros(10,)
