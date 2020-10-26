import argparse
# from fearmem_glm_utils import *
from afni_gPPI import *

parser = argparse.ArgumentParser(description='Func args')

parser.add_argument('-s','--subj',default=None,type=int)

args = parser.parse_args()

os.system('export OMP_NUM_THREADS=48')
# afni_fwhmx(args.subj)
# basic_model_reg_smooth(args.subj)
reg_smooth_gPPI(args.subj)