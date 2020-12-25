import subprocess

if __name__ == '__main__':
    outDir = '/projects/li-lab/yang'
    ret = subprocess.Popen("ls {}".format(outDir), shell=True, stdout=subprocess.PIPE).stdout.read()
    print(ret.decode('utf-8'))

    # subprocess.Popen("echo nihao wohao dajiahao".format(outDir), shell=True, stdout=subprocess.PIPE).communicate()
    # print(ret)



