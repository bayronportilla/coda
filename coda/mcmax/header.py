from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
import os,sys,fnmatch
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from mcmax3dpy import read as mread, plot as mplot
from scipy.integrate import romb,simps,quad,simpson
from astropy.stats import median_absolute_deviation
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from prodimopy.interface1D.infile import write
from prodimopy.read import read_prodimo
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
from scipy.stats import norm
import pandas as pd
from gofish import imagecube
from astropy.convolution import convolve,Gaussian2DKernel
from matplotlib.colors import LogNorm
#from astropy.modeling import models
from scipy import optimize
