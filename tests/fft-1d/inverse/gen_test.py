#!/usr/bin/env python3

from random import uniform
import os,sys

ntest = int(0)

def test_1():
    global ntest
    ntest += 1
    f=open('%.2d.in'%(ntest),'w')
    n=1
    print(n,file=f)
    print("1",file=f)
    f.close()

def test_2():
    global ntest
    ntest += 1
    f=open('%.2d.in'%(ntest),'w')
    n=2
    print(n,file=f)
    print("1 1",file=f)
    f.close()

def test_3():
    global ntest
    ntest += 1
    f=open('%.2d.in'%(ntest),'w')
    n=2
    print(n,file=f)
    print("0 0",file=f)
    f.close()

def test_4(n):
    global ntest
    ntest += 1
    f=open('%.2d.in'%(ntest),'w')
    print(n,file=f)
    for i in range(n):
        print(uniform(-10,10),end=' ',file=f)
    f.close()



if __name__ == "__main__":
    path=sys.argv[1]
    os.chdir(path)
    test_1()
    test_2()
    test_3()
    test_4(128)
    test_4(1024)
    test_4(100)
    test_4(2**13)
    test_4(8191)
