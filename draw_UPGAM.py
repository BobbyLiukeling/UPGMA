# -*- coding: utf-8 -*-
# @Time : 2021/5/21 0021 19:29
# @Author : Bobby_Liukeling
# @File : draw_UPGAM.py

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb

df = pd.read_csv("runtimes.csv")
df = df.sort_values(by = "Unnamed: 0")
# pdb.set_trace()
#画布
fig, ax = plt.subplots()
#使用Axes.plot作图
x = df["Unnamed: 0"] # 数据X轴取值范围
y1 = df["UPGMA"]
ax.plot(x, y1, label='UPGMA') # 作y1 = x 图，并标记此线名为linear
y2 = df["UPGMA_Double"]
ax.plot(x, y2, label='UPGMA_Double') #作y2 = x^2 图，并标记此线名为quadratic
y3 = df["UPGMA_Triple"]
ax.plot(x, y3, label='UPGMA_Triple') # 作y3 = x^3 图，并标记此线名为cubic
ax.set_xlabel('dataset') #设置x轴名称 x label
ax.set_ylabel('runtime') #设置y轴名称 y label
ax.set_title('the running time of the code') #设置图名为Simple Plot
ax.legend() #自动检测要在图例中显示的元素，并且显示

plt.show()



