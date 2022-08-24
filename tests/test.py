from itertools import combinations
import pandas as pd

A = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
temp = combinations(A, 2)
temp = list(temp)
count = 0
for i in list(temp):
    # print(type(i))
    count = count+1
    # print(i)
# print(count)


df = pd.DataFrame(temp, columns=["atom1", "atom2"])
df = df.assign
df["new"][0] = 1
print(df)
a1 = (df["atom1"][44],df["atom2"][44])
print(a1)

print("-------------------")