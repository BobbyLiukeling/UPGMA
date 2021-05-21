# -*- coding: utf-8 -*-
# @Time : 2021/5/21 0021 16:49
# @Author : Bobby_Liukeling
# @File : UPGMA_runtime.py


'''
结论：
由于研究人员的水平限制，在大数据集上的运行时间会缩短一半以上，理论上应该是减少一个数量级的
但是在小数据集上的表现却恰恰相反
将修改的部分单独比较，重新编写
'''
import pandas as pd
import time
import pdb
import numpy as np
import math



def UPGMA(data,shape):
    # print(data)

    for i in range(0,shape-2,1):
        data = UpMatrix(data)
        # print(data.columns)


    #最后一项需要单独合并
    # tree_list = "({},{})".format(data.columns[0],data.columns[1])
    #
    # print(tree_list)
    # Draw_tree(tree_list,shape)

def UpMatrix(data):
    #对矩阵的一次更新
    #data的index和column是都是string类型数据

    # a = time.time()
    values = pd.DataFrame(columns=['min_value','min_index',]) #记录最小值级位置的值


    values['min_value'] = data.min(axis=1) #每行最小值
    values['min_index'] = data.idxmin(axis=1) #每行最小值列下标

    min = 101 #全局最小值
    index_row = '' #全局最小值行坐标
    index_clo = ''#全局最小值列坐标
    for i,index_c,index_r in zip(values['min_value'],values['min_index'],values.index):
        if min>i:
            min = i
            index_clo = index_c #列坐标
            index_row = index_r #行坐标
    # pdb.set_trace()
    #更新矩阵距离s
    #此合并方法会比之前的快20倍
    data = combined(data=data, index_row=index_row, index_clo=index_clo)
    # d = time.time()

    # print("rate:",b-a," ",d-c," ", (d-c)/(b-a))
    return data

def Draw_tree(tree_list,shape):
    from Bio import Phylo
    from io import StringIO
    import pylab

    handle = StringIO(tree_list)
    tree = Phylo.read(handle, "newick")
    # Phylo.draw_graphviz(tree)
    # print(tree)
    # tree.rooted = True
    # Phylo.draw_ascii(tree)
    # pylab.savefig
    # pdb.set_trace()
    # pdb.set_trace()
    Phylo.draw(tree,size=shape)
    # print("XXx")
    # Phylo.write([tree,], "example-both.xml")

def combined(data,index_row,index_clo):

    list_renew = (((data.loc[index_row] + data[index_clo]).fillna(value=0) + (data.loc[index_row] + data.loc[index_clo]).fillna(value=0) + (data[index_row] + data[index_clo]).fillna(value=0)) / 2)
    list_renew = list_renew.drop([index_clo, index_row], axis=0)

    # 将被合并的两行内容从list_renew和data中删除
    data = data.drop([index_clo, index_row], axis=1)
    data = data.drop([index_clo, index_row], axis=0)

    # 将更新的数据添加到data中
    data.loc["({},{})".format(index_row, index_clo)] = list_renew
    data["({},{})".format(index_row, index_clo)] = np.nan
    return data

#将全局更新改为局部更新
def UPMatrix_PLUS(data,values):
    # 对矩阵的一次更新
    # data的index和column是都是string类型数据

    min = 101  # 全局最小值
    index_row = ''  # 全局最小值行坐标
    index_clo = ''  # 全局最小值列坐标
    for i, index_c, index_r in zip(values['min_value'], values['min_index'], values.index):
        if min > i:
            min = i
            index_clo = index_c  # 列坐标
            index_row = index_r  # 行坐标
    # 更新矩阵距离
    list_renew = (((data.loc[index_row] + data[index_clo]).fillna(value=0) + (
                data.loc[index_row] + data.loc[index_clo]).fillna(value=0) + (data[index_row] + data[index_clo]).fillna(
        value=0)) / 2).drop([index_clo, index_row], axis=0)
    data = data.drop([index_clo, index_row], axis=1)
    data = data.drop([index_clo, index_row], axis=0)

    # 将更新的数据添加到data中
    columns_name = "({},{})".format(index_row, index_clo)
    data.loc[columns_name] = list_renew
    data[columns_name] = np.nan
    # 删除被选中合并的两行

    # pdb.set_trace()
    #将values中涉及index_row，index_clo的数据都进行更新，重新查找最小值和最小值下标


    values = values.drop([index_clo, index_row], axis=0)

    for i, index in zip(values['min_index'], values.index):
        if i == index_row or i == index_clo:  # 待更新数据
            temp_loc = data.loc[index]
            if index == data.index[0] and temp_loc.isnull().all(): #修改之后可以避免一个for循环，比之前快近7%
            # if temp_loc.isnull().all():
                temp_index = np.nan
                temp_min = np.nan
            else:
                temp_min = temp_loc.min()  # 更新的每行最小值
                temp_index = temp_loc[temp_loc == temp_min].index.values[0]  # 找下标新的每行最小值下标
            values.loc[index] = [temp_min, temp_index]

    #增加合并的一行
    row_min = np.min(list_renew) # 每行最小值
    row_index_min = data.loc[columns_name][data.loc[columns_name]==row_min].index.values[0] # 每行最小值列下标
    values.loc[columns_name] = [row_min,row_index_min] #添加最后一行

    # pdb.set_trace()

    return data,values

def UPGMA_PLUS (data,shape):
    values = pd.DataFrame()  # 记录最小值级位置的值
    values['min_value'] = data.min(axis=1)  # 每行最小值
    values['min_index'] = data.idxmin(axis=1)  # 每行最小值列下标
    for i in range(0,shape-2,1):
        data, values = UPMatrix_PLUS(data, values)
    # tree_list = "({},{})".format(data.columns[0],data.columns[1])
    # print(tree_list)
    # Draw_tree(tree_list,shape)

#每次合并两项
def UPMatrix_Double(data):
    # a = time.time()
    values = pd.DataFrame(columns=['min_value', 'min_index'])  # 记录最小值级位置的值
    values['min_value'] = data.min(axis=1)  # 每行最小值
    values['min_index'] = data.idxmin(axis=1)  # 每行最小值列下标

    # b = time.time()

    # c = time.time()
    min_1 = 101  # 全局最小值
    min_2 = 102 # 全局第二最小值
    index_row = ''  # 全局最小值行坐标
    index_clo = ''  # 全局最小值列坐标
    index_row_2 = ''  # 全局第二最小值行坐标
    index_clo_2 = ''  # 全局第二最小值列坐标
    for i, index_c, index_r in zip(values['min_value'], values['min_index'], values.index):
        if min(min_1,min_2) > i:#比二者都小，赋值给最小的
            min_1 = i
            index_clo = index_c  # 列坐标
            index_row = index_r  # 行坐标
        elif i>min_1: #i 介于min_1he min_2之间
            min_2 = i
            index_clo_2 = index_c  # 列坐标
            index_row_2 = index_r  # 行坐标

    #判断交集

    if (
        len(set([index_clo,index_row,index_clo_2,index_row_2]))<4
    ): #有交集,不进行加速合并
        # pdb.set_trace()
        data = combined(data=data, index_row=index_row, index_clo=index_clo)
    else: #无交集
        # 更新矩阵距离s
        data = combined(data=data, index_row=index_row, index_clo=index_clo)
        data = combined(data=data, index_row=index_row_2, index_clo=index_clo_2)
    return data

def UPGMA_Double(data,shape):

    # while len(data)>2:
    #     data = UPMatrix_Double(data)
        # print(data.columns)

    while True:
        try:
            data = UPMatrix_Double(data)
        except:
            break


    # #最后一项需要单独合并
    # tree_list = "({},{})".format(data.columns[0],data.columns[1])
    # print(tree_list)
    # Draw_tree(tree_list,shape)

#每次合并三项
def UPMatrix_DoubleX(data):
    # a = time.time()
    values = pd.DataFrame(columns=['min_value', 'min_index'])  # 记录最小值级位置的值
    values['min_value'] = data.min(axis=1)  # 每行最小值
    values['min_index'] = data.idxmin(axis=1)  # 每行最小值列下标

    # b = time.time()

    # c = time.time()
    min_1 = 101  # 全局最小值
    min_2 = 102 # 全局第二最小值
    min_3 = 103
    index_row = ''  # 全局最小值行坐标
    index_clo = ''  # 全局最小值列坐标
    index_row_2 = ''  # 全局第二最小值行坐标
    index_clo_2 = ''  # 全局第二最小值列坐标
    index_row_3 = ''  # 全局第三最小值行坐标
    index_clo_3 = ''  # 全局第三最小值列坐标
    for i, index_c, index_r in zip(values['min_value'], values['min_index'], values.index):
        if min(min_1,min_2,min_3) > i:#比二者都小，赋值给最小的
            min_1 = i
            index_clo = index_c  # 列坐标
            index_row = index_r  # 行坐标
        elif min(min_2,min_3)>i:#i 介于min_1he min_2之间
            min_2 = i
            index_clo_2 = index_c  # 列坐标
            index_row_2 = index_r  # 行坐标

        elif i<min_3: #i 介于min_3he min_2之间
            min_3 = i
            index_clo_3 = index_c  # 列坐标
            index_row_3 = index_r  # 行坐标

    #判断交集




    if (
        len(set([index_clo,index_row,index_clo_2,index_row_2,index_clo_3,index_row_3]))<6
    ): #有交集,不进行加速合并
        # pdb.set_trace()
        data = combined(data = data,index_row = index_row,index_clo = index_clo)

    else: #无交集
        # 更新矩阵距离s
        data = combined(data=data, index_row=index_row, index_clo=index_clo)
        data = combined(data=data, index_row=index_row_2, index_clo=index_clo_2)
        data = combined(data=data, index_row=index_row_3, index_clo=index_clo_3)

    return data

def UPGMA_DoubleX(data,shape):

    # while len(data)>2:
    #     data = UPMatrix_Double(data)
        # print(data.columns)

    while True:
        try:
            data = UPMatrix_DoubleX(data)
        except:
            break


    #最后一项需要单独合并
    # tree_list = "({},{})".format(data.columns[0],data.columns[1])
    # print(tree_list)
    # Draw_tree(tree_list,shape)


def draw_contrast(df):
    import matplotlib.pyplot as plt
    import numpy as np

    labels = ['10', '100', '1000']
    # pdb.set_trace()
    UPGMA = df['UPGMA']
    UPGMA_Double = df['UPGMA_Double']
    UPGMA_Triple = df['UPGMA_Triple']
    #
    # men_means = [20, 34, 30, 35, 27]
    # women_means = [25, 32, 34, 20, 25]

    # pdb.set_trace()

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width, UPGMA, width, label='UPGMA')
    rects2 = ax.bar(x , UPGMA_Double, width, label='UPGMA_Double')
    rects3 = ax.bar(x + width, UPGMA_Triple, width, label='UPGMA_Triple')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('runtime')
    ax.set_title('runtime of the UPGMA')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    ax.bar_label(rects1, padding=2)
    ax.bar_label(rects2, padding=2)
    ax.bar_label(rects3, padding=2)

    fig.tight_layout()

    plt.show()


def draw():
    import matplotlib.pyplot as plt
    import numpy as np

    labels = ['G1', 'G2', 'G3', 'G4', 'G5']
    men_means = [20, 34, 30, 35, 27]
    women_means = [25, 32, 34, 20, 25]

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width / 2, men_means, width, label='Men')
    rects2 = ax.bar(x + width / 2, women_means, width, label='Women')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Scores')
    ax.set_title('Scores by group and gender')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    fig.tight_layout()

    plt.show()


def run():
    shape = 1000  # 指定读取数据的规模
    data = pd.read_csv("randn_data{}.csv".format(shape))

    # data = pd.DataFrame(data,index=map(lambda x:str(x) ,np.arange(0,shape,1))) #设置index,且index为字符串
    data = pd.DataFrame(data, index=np.arange(0, shape, 1))  # 设置index,且index为字符串
    data.index = map(lambda x: str(x), np.arange(0, shape, 1))

    df = pd.DataFrame(columns=["UPGMA", "UPGMA_Double", "UPGMA_Triple"])


    batch = 100
    while(len(data)>=batch):
        shape = len(data)
        a = time.time()
        UPGMA(data, shape)
        b = time.time()
        time1 = b - a
        print("UPGMA:", time1)

        e = time.time()
        UPGMA_Double(data, shape)
        f = time.time()
        time2 = f - e
        print("UPGMA_Double: ", time2)

        g = time.time()
        UPGMA_DoubleX(data, shape)
        h = time.time()
        time3 = h - g
        print("UPGMA_DoubleX : ", time3)

        #将运行时间添加到数据中
        times = [time1,time2,time3]
        df.loc[len(data)] =times

        # 减少数据量
        data = data.iloc[:shape - batch, :shape - batch]





    df.to_csv('runtimes.csv')
    print(df)


if __name__ == '__main__':
    run()


