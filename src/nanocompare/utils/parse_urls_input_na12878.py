import urllib  # the lib that handles the url stuff
import urllib.request
import re

in_html = "https://raw.githubusercontent.com/nanopore-wgs-consortium/NA12878/master/nanopore-human-genome/rel_3_4.md"

if __name__ == '__main__':
    with urllib.request.urlopen(in_html) as response:
        html = response.read()
    html_str = str(html)

    ## Cope with: http://s3.amazonaws.com/nanopore-human-wgs/rel3-fast5-chr1.part01.tar
    flist = re.findall(r"http://s3\.amazonaws\.com/nanopore-human-wgs/rel\d-fast5-chr[0-9XYM]+\.part\d+\.tar", html_str)

    print(len(flist))

    with open("na12878.filelist.txt", 'w') as outf:
        for fn in flist:
            outf.write(f'{fn}\n')








