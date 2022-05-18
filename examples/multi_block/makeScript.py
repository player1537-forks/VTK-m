header = '''#!/bin/bash
#BSUB -env "all"
#BSUB -P PHY122
#BSUB -J Poincare
#BSUB -W 0:30
#BSUB -nnodes 1
#BSUB -o run-%J.log
#BSUB -e run-%J.log
#BSUB -q debug

cd /gpfs/alpine/proj-shared/csc143/pugmire/vtkm/multiblock'''

print(header)


dims = [8,16,32,64,128]
numBlocks = [8, 16, 32, 64, 128]
numTask = [1, 2, 4, 8, 16, 32, 64]

basecmd = 'jsrun -n1 -a1 -c2 -g1 ./build/examples/multi_block/MultiBlock'
basecmd = basecmd + ' --field tangle --isolevels 20 --gpu'

def dumpCmd(runcmd, cnt) :
  #print(runcmd + '&')
  print(runcmd)
  if cnt % 6 == 0 :
     print('wait\n\n')
     

cnt = 0
for d in dims :
    for n in numBlocks :
        cmd = basecmd + ' --tangle %d %d %d %d' % (n, d, d, d)
        print('echo \n\n\n')
        print('echo "dims= %d n= %d"\n' % (d, n))
        
        #do a serial task
        print('echo serial')
        runcmd = cmd + ' --threading serial'
        runcmd = runcmd +  ' --output ./output/out.%d.txt' % 0 # (cnt%6)
#        cnt = cnt+1
#        dumpCmd(runcmd, cnt)
        
        for t in numTask :
            print('echo task= ', t)
            runcmd = cmd + ' --threading task %d'% t
            runcmd = runcmd +  ' --output ./output/out.%d.txt' % 0 # (cnt%6)
#            cnt = cnt+1
#            dumpCmd(runcmd, cnt)

        for t in numTask :
            print('echo taskOld= ', t)
            runcmd = cmd + ' --threading taskOld %d'% t
            runcmd = runcmd +  ' --output ./output/out.%d.txt' % 0 # (cnt%6)
            cnt = cnt+1
            dumpCmd(runcmd, cnt)
           

