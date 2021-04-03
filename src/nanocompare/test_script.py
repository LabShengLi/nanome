import os, psutil
process = psutil.Process(os.getpid())
print(process.memory_info().rss/1000/1000)  # in bytes GB

