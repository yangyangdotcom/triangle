from itertools import combinations
A = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
temp = combinations(A, 2)
count = 0
for i in list(temp):
    print(type(i))
    count = count+1
    print(i)
print(count)