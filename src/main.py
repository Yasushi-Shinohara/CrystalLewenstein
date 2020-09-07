#!/usr/bin/python
# coding: UTF-8
# This is created 2020/06/08 by Y. Shinohara
# This is lastly modified 2020/06/08 by Y. Shinohara #This part is highly doubtable because of my lazyness
import time
ts = time.time()
from modules.print_funcs import print_header, print_footer, print_midtime, print_endtime
from modules.functions import *
from modules.parameters import parameter_class
print_header()
import sys
import numpy as np
import math
import ctypes as ct
from modules.constants import *


tt = time.time()
print_midtime(ts,tt)

te = time.time()
#print_endtime(ts,tt,te,param.Nt)


print_footer() 
sys.exit()

