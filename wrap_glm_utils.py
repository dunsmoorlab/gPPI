import argparse
from fearmem_glm_utils import *

parser = argparse.ArgumentParser(description='Func args')

parser.add_argument('-s','--subj',default=None,type=int)

args = parser.parse_args()

os.system('export OMP_NUM_THREADS=48')
afni_fwhmx(args.subj)
