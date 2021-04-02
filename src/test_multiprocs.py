from multiprocessing import Manager
manager = Manager()
l = manager.list([i*i for i in range(10)])

d1={1:3, 2:[1,2,3]}
dm1=manager.dict(d1)

print(l)

print(dm1)