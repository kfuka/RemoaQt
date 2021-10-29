import numpy as np
import glob

npz_files = glob.glob("./data_base/0000000/*.npz")
# print(npz_files)

np_array = np.load(npz_files[0])
print(np_array.files)

roi1 = np_array["roi_center_1"]
print(roi1)

img1 = np_array["img1"]
print(img1)
