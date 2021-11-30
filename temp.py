import numpy as np
arr =[1,2,3]
arr1 = [4,5,6]
arr2 = [7,8,9]
arr2 = np.concatenate((arr, arr1, arr2),axis=None)

print(arr2)