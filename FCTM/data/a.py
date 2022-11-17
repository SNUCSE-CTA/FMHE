import os

#for year in range(1907,2019):
path="seoul/"
num_element=0
temp_list=[]
start_year=1971
for fil in sorted(os.listdir(path)):
    if fil[0:7]=='SURFACE' and int(fil[20:24])>=start_year:
        f=open(path+fil,"r", encoding='utf8', errors="ignore")
        f.readline()
        for line in f:
        #print(line[1])
        #print(int(10*float(line.split(',')[2])),' ')
            k=line.split(',')[2]
            if k is not '':
                num_element=num_element+1
                temp_list.append(str(float(k)))
        #print(fil)
        f.close()
print(num_element)
print(" ".join(temp_list))
