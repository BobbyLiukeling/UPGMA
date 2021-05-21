# # -*- coding: utf-8 -*-
# # @Time : 2021/4/3 0003 20:54
# # @Author : Bobby_Liukeling
# # @File : design_data.py
#
# '''
# 设计一组随机数据，数据大小为0-100，小数点后保留一位，数据量为10,000
# 数据代表两组生物序列之间的距离，数值越小则距离越小
# '''
import numpy as np
import pdb
import pandas as pd

def get_data(shape):
    np.random.seed(0)
     #矩阵大小
    variance = 12.5 #方差
    mean = 50 #均值
    a = np.random.normal(mean,variance,(shape,shape)) #产生随机数
    # pdb.set_trace()

    #对数据进行归一化处理，将数据数值集中到1（1,100）
    a = ((a-a.min())/(a.max()-a.min()))*100

    # 保留小数点后两位
    prices = 2
    asign = lambda t: np.around(t, prices)
    a = asign(a)

    # 保留下三角区域的0
    a =  pd.DataFrame(a,columns=np.arange(0,shape,1))
    a = a.replace(0,'0')

    b = np.tril(a,-1) #下三角矩阵
    data  = pd.DataFrame(b,columns=np.arange(0,shape,1))

    data = data.replace(0, np.nan)#将上三角值置为空
    data = data.replace('0',0)#还原下三角0 值

    # data.to_csv("randn_data{}.csv".format(shape),index=None)
    data.to_csv("randn_data{}.csv".format(shape),index=None)
    return data

# class data_position:
#     def __init__(self,data,row,loc):
#         self.data = data #数据
#         self.row = row #横坐标
#         self.loc = loc #纵坐标



if __name__ == "__main__":
    shape = 100
    data = get_data(shape)
    pdb.set_trace()