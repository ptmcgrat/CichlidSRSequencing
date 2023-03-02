import subprocess as sp
import shlex

sp.run('source activate pipeline && echo $PWD && source deactivate', shell=True)

# sp.run(shlex.split('conda init bash'), shell=True)
# sp.run(shlex.split('conda deactivate'), shell=True)