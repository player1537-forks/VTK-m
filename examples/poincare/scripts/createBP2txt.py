import sys
import os.path

dir = sys.argv[0][0: sys.argv[0].rfind('/')]
posFileNm = sys.argv[1]
negFileNm = sys.argv[2]
outFileNm = sys.argv[3]

num = 0
while True :

    pFile = posFileNm + '.%03d.bp' % num
    nFile = negFileNm + '.%03d.bp' % num

    if not os.path.exists(pFile) or not os.path.exists(nFile) :
        break

    pOutFile = 'pos.' + outFileNm + '.%03d' % num
    nOutFile = 'neg.' + outFileNm + '.%03d' % num

    cmd = 'python3 %s/bp2txt.py %s %s 1' % (dir, pFile, pOutFile)
    print cmd

    cmd = 'python3 %s/bp2txt.py %s %s 1' % (dir, nFile, nOutFile)
    print cmd


    ##combine them.
    cmd = 'tail -n +2 %s.txt > .tmp.file' % nOutFile
    print cmd
    cmd = 'cat %s.txt .tmp.file > both.%s.%03d.txt' % (pOutFile, outFileNm, num)
    print cmd
    cmd = 'rm .tmp.file'
    print cmd

    num = num+1
