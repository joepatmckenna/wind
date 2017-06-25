from datetime import datetime, timedelta

HOUR = timedelta(hours=1)
write_dir = 'submit/'

t = datetime(2005,8,14,0)
idx = 0
while t < datetime(2005,9,15,0):

    filename = write_dir + '%03i.sh' % (idx,)

    with open(filename, 'w+') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -J "wind%03i"\n' % (idx,))
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH --ntasks-per-node=32\n')
        if t < datetime(2005,9,9,0):
            f.write('#SBATCH -p backfill\n')
        else:
            f.write('#SBATCH -p backfill2\n')
        f.write('#SBATCH --mail-type="ALL"\n')
        f.write('#SBATCH -t 02:00:00\n')
        f.write('#SBATCH -o out.txt\n')
        f.write('python write_flow.py %s' % (script, t.strftime('%Y-%m-%d-%H'),))

    t += 3*HOUR
    idx += 1




