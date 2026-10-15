"""Power/load cases deliberately different for training and hold-out verification."""
import numpy as np


def make_times(duration=6.,dt=.15):
    n=int(round(duration/dt))
    return np.linspace(0.,n*dt,n+1)


def training_power(times):
    x=np.empty((len(times),4))
    for i,t in enumerate(times):
        if t < .9:
            x[i]=[100,100,100,100]
        elif t < 2.7:
            x[i]=[250,100,150,100]
        elif t < 4.2:
            x[i]=[100,250,100,200]
        else:
            x[i]=[200,200,200,200]
    return x


def holdout_power(times):
    x=np.empty((len(times),4))
    for i,t in enumerate(times):
        if t < .6:
            x[i]=[120,110,90,105]
        elif t < 2.1:
            x[i]=[180,280,140,180]
        elif t < 3.9:
            x[i]=[350,150,100,250]
        else:
            x[i]=[160,220,270,130]
    return x


def dashboard_power(times,tdp=800,scenario="step"):
    p=np.full((len(times),4),tdp/8)
    for j,t in enumerate(times):
        if t >= 5:
            p[j]=tdp/4
        if scenario=="hotspot" and 15<=t<35:
            p[j]=tdp/4*np.array([1.75,.75,.75,.75])
        elif scenario=="surge" and 15<=t<35:
            p[j]=tdp/4*1.5
    return p
