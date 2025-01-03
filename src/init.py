import numpy as np
from scipy import constants as pc
global const
global opt

const = {}
const['icm2ifs'] = pc.c*1e2*2.0*np.pi*1e-15
const['hbar'] = pc.h/2.0/np.pi
const['kb'] = pc.k

opt = {}

def setpara(f):
    data = f.readlines()
    for line in data:
        # parse input, assign values to variables
        if "=" in line:
            key, value = line.split(" = ")
            opt[key.strip()] = value.strip()
        f.close()

    # Parse standard options
    opt['temperature'] = float(opt['temperature'])
    opt['Tc'] = float(opt['Tc'])
    opt['N_t'] = int(opt['N_t'])
    opt['wmax_quad'] = float(opt['wmax_quad'])

    # Parse method options
    if opt['method'] == "ID":
        opt['Omega_min'] = float(opt['Omega_min']) 
        opt['Omega_max'] = float(opt['Omega_max'])
        opt['N_w'] = int(opt['N_w'])
        opt['eps'] = float(opt['eps'])
        opt['frank'] = int(opt['frank'])
        #opt['rand'] = bool(int(opt['rand']))
    elif opt['method'] == "BSDO":
        opt['Omega_min'] = float(opt['Omega_min']) 
        opt['Omega_max'] = float(opt['Omega_max'])
        opt['Msp'] = int(opt['Msp'])
        opt['N_w'] = int(opt['N_w'])
    elif opt['method'] == "LOG":
        opt['Msp'] = int(opt['Msp'])
        opt['Omega_max'] = float(opt['Omega_max'])
    elif opt['method'] == "MDM":
        opt['Msp'] = int(opt['Msp'])
        opt['Omega_max'] = float(opt['Omega_max'])
    else:
        print("Invalid method: ", opt['method'])
        print("Please choose from ID, BSDO, LOG, MDM.")
        exit()

    # Parse spectral density options
    if opt['stype'] == "PWR":
        opt['s'] = float(opt['s'])
        opt['alpha'] = float(opt['alpha'])/float(opt['gamc'])
        opt['gamc'] = float(opt['gamc'])
    elif opt['stype'] == "TM" or opt['stype'] == "BO":
        opt['Omg'] = np.array([float(num.strip()) for num in opt['Omg'].split(',')])
        opt['Gam'] = np.array([float(num.strip()) for num in opt['Gam'].split(',')])
        opt['Lam'] = np.array([float(num.strip()) for num in opt['Lam'].split(',')])
    else:
        print("Invalid spectral density type: ", opt['stype'])
        print("Please choose from PWR, TM, BO.")
        exit()