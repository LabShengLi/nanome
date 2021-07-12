import os
import matplotlib.font_manager as font_manager
font_manager._rebuild()

import matplotlib.pyplot as plt

from matplotlib import font_manager
for font in font_manager.fontManager.ttflist:
    print(font)

plt.rcParams['font.sans-serif'] = ['Arial']

# from nanocompare.global_config import pic_base_dir
# print(plt.rcParams["font.family"])
# print(plt.rcParams['font.sans-serif'])

# plt.rcParams.update({'font.sans-serif': ['Arial']})
#
# print(plt.rcParams["font.family"])
# print(plt.rcParams['font.sans-serif'])

# number of employees of A
emp_count = [3, 20, 50, 200, 350, 400]
year = [2014, 2015, 2016, 2017, 2018, 2019]
# plot a line chart
fig, ax = plt.subplots()
ax.plot(year, emp_count)
# set axis titles
ax.set_xlabel("Year")
ax.set_ylabel("Employees")
# set chart title
ax.set_title("Employee Growth at A")
# plt.show()

# outfn = os.path.join(pic_base_dir, "test.jpg")
plt.savefig("test.jpg", format='jpg', bbox_inches='tight', dpi=600)
