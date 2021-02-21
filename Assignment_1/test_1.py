import numpy as np 
import pandas as pd
import math
import seaborn as sns
import matplotlib.pyplot as plt
import csv

# %matplotlib inline

df = pd.read_csv("out.csv")
print(df['p'])

method = []

for i in range(len(df['p'])):
	method.append((i)%3 + 1)

df['method'] = method

print(df['method'])

PT = [16, 25, 36, 49, 64]
NT = [256, 1024, 4096, 16384, 65536, 262144, 1048576]

iteration = []

for i in range(len(PT)):
	for j in range(len(NT)):
		for k in range(3):
			for l in range(5):
				iteration.append(l)

df['iteration'] = iteration

lg = []

for i in range(len(df['p'])):
	lg.append(math.log(df['dp'][i],2))
df['lg'] = lg

print(df['iteration'])

# for i in range(len(df['time'])):
# 	df['time'][i] = round(df['time'][i],2)

l1 = df[(df['p']==16) & (df['method']==1)]
print(l1)

sns.set(style="whitegrid")
# p = sns.boxplot(x = l1['time'], y = l1['lg'], data = l1,width=0.1, linewidth=1, notch=False, color = "red", hue='method')
# plt.set(yscale="log")
data = l1['time']
sns.boxplot(y = data, x = df['lg'])
# p.get_figure()
plt.yscale('log')
plt.show()

# for p in PT:
# 	for n in NT:
# 		for k in range(3):
# 			s = []
# 			for i in range(len(df['p'])):
# 				if df['p'][i]==p and df['dp'][i]==n and (k+1)==df['iteration'][i]:
# 					s.append(df['time'][i])
# 			seaborn.boxplot(x = s, y = n)