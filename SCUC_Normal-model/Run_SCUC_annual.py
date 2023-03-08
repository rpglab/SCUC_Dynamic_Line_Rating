### Run SCUC
## Author: Jin Lu, University of Houston.
## link: https://rpglab.github.io/resources/

### Import
import os
from UC_function import *
import numpy as np

# read test case
case_name = '123bus_case_final.pkl'  # after the generator cost adjustment
with open(case_name, 'rb') as input:
    case_inst = pickle.load(input)


### form annual test case profiles using annual line rating, wind, solar, load profiles
def run_annual_UC(case_inst,dnum_start,dnum_end):
    for d in range(dnum_start-1,dnum_end):
        day_num = d + 1
        # read annual line rating
        line_annual = np.loadtxt('Line_annual_Dmin.txt')    # line rating file after the line adjustment [255][365]
        for l in range(case_inst.branchtotnum):
            case_inst.branch.rateA[l] = line_annual[l][d]
        # read solar, wind annual
        solar_folder = os.getcwd() + '//' + 'solar_annual'
        solar_flnm = 'solar_annual_D' + str(day_num) + '.txt'
        solar_path = solar_folder + '//' + solar_flnm
        solar_d = np.loadtxt(solar_path)  # [sgen][hr]
        wind_folder = os.getcwd() + '//' + 'wind_annual'
        wind_flnm = 'wind_annual_D' + str(day_num) + '.txt'
        wind_path = wind_folder + '//' + wind_flnm
        wind_d = np.loadtxt(wind_path)  # [wgen][hr]
        # read load profile
        load_folder = os.getcwd() + '//' + 'load_annual'
        load_flnm = 'load_annual_D' + str(day_num) + '.txt'
        load_path = load_folder + '//' + load_flnm
        load_d = np.loadtxt(load_path)  # [hr][bus]
        ## calculate the load for this day
        # list of non wind or solar generators
        reg_gen_list = []
        nonreg_gen_list = []
        for i in range(case_inst.gentotnum):
            if case_inst.gen.fuel_type[i] != 'Wind':
                if case_inst.gen.fuel_type[i] != 'Solar':
                    gen_num = i + 1
                    reg_gen_list.append(gen_num)
            if case_inst.gen.fuel_type[i] == 'Wind':
                gen_num = i + 1
                nonreg_gen_list.append(gen_num)
            if case_inst.gen.fuel_type[i] == 'Solar':
                gen_num = i + 1
                nonreg_gen_list.append(gen_num)
        ## wind&solar output combined with load data
        solar_genidx_list, solar_busnum_list = find_solar_bus(case_inst)
        wind_genidx_list, wind_busnum_list = find_wind_bus(case_inst)
        load_b_t = np.transpose(load_d)
        # count wind power
        for h in range(24):
            for i in range(len(wind_genidx_list)):
                g_idx = wind_genidx_list[i]
                b_idx = case_inst.gen.bus[g_idx] - 1
                load_b_t[b_idx, h] -= wind_d[i, h]
        # count solar power
        for h in range(24):
            for i in range(len(solar_genidx_list)):
                g_idx = solar_genidx_list[i]
                b_idx = case_inst.gen.bus[g_idx] - 1
                load_b_t[b_idx, h] -= solar_d[i, h]
        ## run full SCUC
        UC_case_run1 = build_UC_full(case_inst, load_b_t)
        solve_UC(UC_case_run1, "UCresults.pickle", 'UCcase', day_num)
        # fixed the u_g_t, solve again and obtain dual variable (electricity prices)
        UC_case_run2 = build_UC_full_Run2(case_inst, UC_case_run1, load_b_t)
        solve_UC(UC_case_run2, "UCresults_run2.pickle", "UCcase", day_num)
        print('finish UC for Day ' + str(day_num))

# run SCUC for the selected days
run_annual_UC(case_inst,365,365)


