# -*- coding: utf-8 -*-
# @Time : 2021/4/3 0003 20:54
# @Author : Bobby_Liukeling
# @File : UPGMA+.py

'''
结论：
由于研究人员的水平限制，在大数据集上的运行时间会缩短一半以上，理论上应该是减少一个数量级的
但是在小数据集上的表现却恰恰相反
'''
import pandas as pd
import time
import pdb
import numpy as np
import math



def UPGMA(data,shape):

    for i in range(0,shape-2,1):
        # a = time.time()
        data = UpMatrix(data)
        # b = time.time()
        # print("UPGMA",b-a)
        #
        # c = time.time()
        # data, values = UPMatrix_PLUS(data, values)
        # d = time.time()
        # print("UPGMA+:", d - c)
        # pdb.set_trace()

    #最后一项需要单独合并
    tree_list = "({},{})".format(data.columns[0],data.columns[1])
    # print(tree_list)

    # Draw_tree(tree_list,shape)

    # return  tree_list


def UpMatrix(data):
    #对矩阵的一次更新
    #data的index和column是都是string类型数据
    values = pd.DataFrame() #记录最小值级位置的值

    values['min_value'] = data.min(axis=1) #每行最小值
    values['min_index'] = data.idxmin(axis=1) #每行最小值列下标



    # c = time.time()
    min = values['min_value'].min() #全局最小值
    index_row = values['min_value'][values['min_value'] == min].index.values[0] #全局最小值行下标，str类型
    index_clo = values['min_index'][index_row] #全局最小值列下标,str类型数据
    # d = time.time()
    # print("d-c:",d-c)


    # pdb.set_trace()
    #更新矩阵距离
    list_renew = []
    # pdb.set_trace()

    for index,value in data.iteritems(): #逐列遍历
        try:
            if(index == index_clo or index == index_row): #被合并行不更新
                pass
            else:
                # pdb.set_trace()

                if (pd.notna(data[index][index_row]) and pd.notna(data[index][index_clo])):
                    list_renew.append(np.around((data[index][index_row] + data[index][index_clo]) / 2, 2))
                elif(pd.notna(data[index][index_row]) and pd.isna(data[index][index_clo])):
                    list_renew.append(np.around((data[index][index_row] + data[index_clo][index]) / 2,2))
                elif(pd.isna(data[index][index_row]) and pd.notna(data[index][index_clo])):
                    list_renew.append(np.around((data[index_row][index] + data[index][index_clo]) / 2,2))
                elif (pd.isna(data[index][index_row]) and pd.isna(data[index][index_clo])):
                    list_renew.append(np.around((data[index_row][index] + data[index_clo][index]) / 2, 2))
                else:
                    pdb.set_trace()
        except Exception as e:
            print("Exception:",e)
            pdb.set_trace()

        # 将被合并的两行内容从list_renew和data中删除


    # print("2:", data.index)
    data = data.drop(index_clo, axis=1)
    data = data.drop(index_row, axis=1)
    data = data.drop(index_clo, axis=0)
    data = data.drop(index_row, axis=0)
    # print("3:", data.index)

    #将更新的数据添加到data中
    # pdb.set_trace()

    data.loc["({},{})".format(index_row,index_clo)] =list_renew
    data["({},{})".format(index_row,index_clo)] = np.nan

    # pdb.set_trace()



    return data
    # print("4:", data.index)

    # pdb.set_trace()


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

def UPMatrix_PLUS(data,values):
    # 对矩阵的一次更新
    # data的index和column是都是string类型数据

    # min_row = np.arange(0,shape,1).tolist() #记录行位置，便于删除行

    min = values['min_value'].min()  # 全局最小值
    index_row = values['min_value'][values['min_value'] == min].index.values[0]  # 全局最小值行下标，str类型
    index_clo = values['min_index'][index_row]  # 全局最小值列下标,str类型数据

    # 更新矩阵距离
    list_renew = []
    # pdb.set_trace()

    for index, value in data.iteritems():  # 逐列遍历
        try:
            if (index == index_clo or index == index_row):  # 被合并行不更新
                pass
            else:
                # pdb.set_trace()

                if (pd.notna(data[index][index_row]) and pd.notna(data[index][index_clo])):
                    list_renew.append(np.around((data[index][index_row] + data[index][index_clo]) / 2, 2))
                elif (pd.notna(data[index][index_row]) and pd.isna(data[index][index_clo])):
                    list_renew.append(np.around((data[index][index_row] + data[index_clo][index]) / 2, 2))
                elif (pd.isna(data[index][index_row]) and pd.notna(data[index][index_clo])):
                    list_renew.append(np.around((data[index_row][index] + data[index][index_clo]) / 2, 2))
                elif (pd.isna(data[index][index_row]) and pd.isna(data[index][index_clo])):
                    list_renew.append(np.around((data[index_row][index] + data[index_clo][index]) / 2, 2))
                else:
                    pdb.set_trace()
        except Exception as e:
            print("Exception:", e)
            pdb.set_trace()

        # 将被合并的两行内容从list_renew和data中删除

    # print("2:", data.index)
    # pdb.set_trace()
    data = data.drop(index_clo, axis=1)
    data = data.drop(index_row, axis=1)
    data = data.drop(index_clo, axis=0)
    data = data.drop(index_row, axis=0)
    # print("3:", data.index)

    # 将更新的数据添加到data中
    # pdb.set_trace()

    columns_name = "({},{})".format(index_row, index_clo)
    data.loc[columns_name] = list_renew
    data[columns_name] = np.nan

    # 删除被选中合并的两行
    values = values.drop(index_clo, axis=0)
    values = values.drop(index_row, axis=0)

    # pdb.set_trace()
    #将values中涉及index_row，index_clo的数据都进行更新，重新查找最小值和最小值下标
    for i,index in zip(values['min_index'],values.index):
        if i == index_row or i == index_clo: #待更新数据
            # pdb.set_trace()
            temp_min = data.loc[index].min() #更新的每行最小值
            temp_index =  data.loc[index][data.loc[index]==temp_min].index.values
            if len(temp_index) == 0: #整行为空
                temp_index = np.nan
                temp_min = np.nan
            else:
                temp_index = temp_index[0]#更新的每行最小值下标
            values.loc[index] = [temp_min,temp_index]


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
        # a = time.time()
        # data = UpMatrix(data)
        # b = time.time()
        # print("UPGMA", b - a)

        # c = time.time()
        data, values = UPMatrix_PLUS(data, values)
        # d = time.time()
        # print("UPGMA+:", d - c)
        pdb.set_trace()

    #最后一项需要单独合并
    tree_list = "({},{})".format(data.columns[0],data.columns[1])
    # print(tree_list)

    # Draw_tree(tree_list,shape)
    pass


if __name__ == '__main__':

    shape = 100 #指定读取数据的规模
    data = pd.read_csv("randn_data{}.csv".format(shape) )

    # data = pd.DataFrame(data,index=map(lambda x:str(x) ,np.arange(0,shape,1))) #设置index,且index为字符串
    data = pd.DataFrame(data,index=np.arange(0,shape,1)) #设置index,且index为字符串
    data.index = map(lambda x:str(x) ,np.arange(0,shape,1))

    start = time.time()
    UPGMA(data,shape)
    end = time.time()
    print("UPGMA:",end - start)

    # start = time.time()
    # UPGMA_PLUS(data, shape)
    # end = time.time()
    # print("UPGMA_PLUS:", end - start)


    # pdb.set_trace()
