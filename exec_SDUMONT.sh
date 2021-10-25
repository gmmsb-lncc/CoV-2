module load python/3.8.2  
echo 'Python LOADED'


srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.427_429/B.1.427_B.1.429-USA-California_MSA.afa  ./B.1.427_429/B.1.427_429_cosensus.csv ./B.1.427_429/B.1.427_429_canonical_cosensus.fasta 
echo 'B.1.427_429 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.525/B.1.525-UK-Nigeria_MSA.afa  ./B.1.525/B.1.525_cosensus.csv ./B.1.525/B.1.525_canonical_cosensus.fasta 
echo 'B.1.525 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.526/BB.1.526-NewYork_MSA.afa  ./B.1.526/B.1.526_cosensus.csv ./B.1.526/B.1.526_canonical_cosensus.fasta 
echo 'B.1.526 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.617.1/B.1.617.1-India_MSA.afa  ./B.1.617.1/B.1.617.1_cosensus.csv ./B.1.617.1/B.1.617.1_canonical_cosensus.fasta 
echo 'B.1.617.1 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./C.37/C.37-Peru_MSA.afa  ./C.37/C.37_cosensus.csv ./C.37/C.37_canonical_cosensus.fasta 
echo 'C.37 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./P.1/P.1-Brazil-Japan_MSA.afa  ./P.1/P.1_cosensus.csv ./P.1/P.1_canonical_cosensus.fasta 
echo 'P.1 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./P.2/P.2-Brazil_MSA.afa  ./P.2/P.2_cosensus.csv ./P.2/P.2_canonical_cosensus.fasta 
echo 'P.2 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./P.3/P.3-Philippines_MSA.afa  ./P.3/P.3_cosensus.csv ./P.3/P.3_canonical_cosensus.fasta 
echo 'P.3 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/1/B.1.1.7-UK-pt1_MSA.afa  ./B.1.1.7/1/B.1.1.7-UK-pt1_cosensus.csv ./B.1.1.7/1/B.1.1.7-UK-pt1_canonical_cosensus.fasta
echo 'B.1.1.7-UK-pt1 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/2/B.1.1.7-UK-pt2_MSA.afa  ./B.1.1.7/2/B.1.1.7-UK-pt2_cosensus.csv ./B.1.1.7/2/B.1.1.7-UK-pt2_canonical_cosensus.fasta 
echo 'B.1.1.7-UK-pt2 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/3/B.1.1.7-UK-pt3_MSA.afa  ./B.1.1.7/3/B.1.1.7-UK-pt3_cosensus.csv ./B.1.1.7/3/B.1.1.7-UK-pt3_canonical_cosensus.fasta 
echo 'B.1.1.7-UK-pt3 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/4/B.1.1.7-UK-pt4_MSA.afa  ./B.1.1.7/4/B.1.1.7-UK-pt4_cosensus.csv ./B.1.1.7/4/B.1.1.7-UK-pt4_canonical_cosensus.fasta 
echo 'B.1.1.7-UK-pt4 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/5/B.1.1.7-UK-pt5_MSA.afa  ./B.1.1.7/5/B.1.1.7-UK-pt5_cosensus.csv ./B.1.1.7/5/B.1.1.7-UK-pt5_canonical_cosensus.fasta 
echo 'B.1.1.7-UK-pt5 FINISH'

srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py ./B.1.1.7/6/B.1.1.7-UK-pt6_MSA.afa  ./B.1.1.7/6/B.1.1.7-UK-pt6_cosensus.csv ./B.1.1.7/6/B.1.1.7-UK-pt6_canonical_cosensus.fasta 
echo 'B.1.1.7-UK-pt6 FINISH'
