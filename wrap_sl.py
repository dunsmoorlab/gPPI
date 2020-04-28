import argparse
from sl_er_rsa import *


parser = argparse.ArgumentParser(description='Func args')

parser.add_argument('-s','--subj',default=None,type=int)

args = parser.parse_args()

ER_searchlight(args.subj,run=True)
