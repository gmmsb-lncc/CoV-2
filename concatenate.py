import shutil

with open('B.1.1.7.afa','wb') as wfd:
    for f in [
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt1_MSA.afa',
    'B.1.1.7-UK-pt1_MSA.afa'

        		]:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)
