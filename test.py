# -*- coding: utf-8 -*-
# @Time : 2021/4/7 0007 14:16
# @Author : Bobby_Liukeling
# @File : test.py

from numba import cuda

#设置模拟器，本主机上没有英伟达显卡
# SET NUMBA_ENABLE_CUDASIM = 1

def cpu_print():
    print("print by cpu.")

@cuda.jit
def gpu_print():
    # GPU核函数
    print("print by gpu.")

def main():
    gpu_print[1, 2]()
    cuda.synchronize()
    cpu_print()

if __name__ == "__main__":
    main()
