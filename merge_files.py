# submit to run over all files of a generator
import sys
import os 
import glob


input_dir = "OUTPUT"

files = glob.glob(input_dir + "/output_*.root")

files = set([''.join([i for i in f.split("/")[-1] if not i.isdigit()]) for f in files])
print(files)

for generator in files:
  if not "dpm" in generator: continue
  if os.path.exists(os.path.join(input_dir,generator)):
    os.system('rm {}'.format(os.path.join(input_dir,generator)) )
  os.system('hadd {} {}'.format(os.path.join(input_dir,generator), os.path.join(input_dir,generator).replace(".root","*")))
