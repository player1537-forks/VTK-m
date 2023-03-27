import sys


outNm = sys.argv[1]
isPos = int(sys.argv[2])
outFile = open(outNm,'w')


h = 0.01


runCmd = 'jsrun -n1 -a1 -c2 -g1 ./examples/poincare/Poincare'
runCmd = runCmd + ' --gpuParams 256 128'
runCmd = runCmd + ' --dir ../ITER-data'
runCmd = runCmd + ' --numPunc 1000'
runCmd = runCmd + ' --psiRange 0.99 1.0 3'
runCmd = runCmd + ' --thetaRange 0 360 3'
#runCmd = runCmd + ' --thetaRange 3'
#runCmd = runCmd + ' --psiVals 0.99 0.995 0.999'
#runCmd = runCmd + ' --thetaVals 252.5 253.5'

whichType = ''
if isPos : whichType = 'pos'
else : whichType = 'neg'



header ='''#!/bin/bash
#BSUB -env "all"
#BSUB -P PHY122
#BSUB -J Poincare
#BSUB -W 2:00
#BSUB -nnodes 1
#BSUB -o %s-run-%%J.log
#BSUB -q debug ''' % whichType

outFile.write(header)
outFile.write('\n\n')

xgc3DFileRange = range(406, 470, 2)

cnt = 0
for i in xgc3DFileRange :
    xgcFile = 'xgc.3d.%05d.bp' % i

    if isPos == 1 :
        output = 'pos.out.%03d' % cnt
        step = ' --stepSize %f' % h
        outStr = runCmd + ' --output ' + output + step + ' --3DFile %s' % xgcFile
        outFile.write(outStr + ' & \n')
    else:
        output = 'neg.out.%03d' % cnt
        step = ' --stepSize %f' % -h
        outStr = runCmd + ' --output ' + output + step + ' --3DFile %s' % xgcFile
        outFile.write(outStr + ' & \n')

    cnt = cnt+1
    if cnt % 6 == 0 or i == xgc3DFileRange[-1]:
        outFile.write('wait \n\n')
