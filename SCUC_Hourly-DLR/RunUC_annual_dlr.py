### run SCUC with dynamic line rating profiles

### Import
import os
import numpy as np
import shutil
from UC_function_DLR import *

# read test case
case_name = '123bus_case_final.pkl'  # after the line adjustment
with open(case_name, 'rb') as input:
    case_inst = pickle.load(input)

### function to read the DLR result
def read_dlr(case_inst):
    dlr_result = np.zeros((365,24,case_inst.branchtotnum))
    dlr_folder = os.getcwd() + '//' + 'dynamic_rating_result'
    for l in range(case_inst.branchtotnum):
        dlr_flnm = 'dynamic_rating_L' + str(l + 1) + '.txt'
        dlr_path = dlr_folder + '//' + dlr_flnm
        dlr_l = np.loadtxt(dlr_path)
        for d in range(365):
            for h in range(24):
                dlr_result[d,h,l] = dlr_l[d,h]
    return dlr_result

### form annual test case profiles using hourly DLR, wind, solar, load profiles
def run_annual_UC_dlr(case_inst,dnum_start,dnum_end):
    for d in range(dnum_start-1,dnum_end):
        day_num = d + 1
        # read annual line rating
        line_dlr = read_dlr(case_inst)    # line rating file after the line adjustment and dlr [d,h,l]
        line_d = line_dlr[d,:,:]    # [hour24][line255]
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
        line_l_t = np.transpose(line_d)
        UC_case_run1 = build_UC_full(case_inst, load_b_t,line_l_t)
        solve_UC(UC_case_run1, "UCresults.pickle", 'UCcase', day_num)
        # fixed the u_g_t, solve again and obtain dual variable (electricity price)
        UC_case_run2 = build_UC_full_Run2(case_inst, UC_case_run1, load_b_t,line_l_t)
        solve_UC(UC_case_run2, "UCresults_run2.pickle", "UCcase", day_num)
        print('finish UC for Day ' + str(day_num))

# run SCUC with DLR for the selected days
run_annual_UC_dlr(case_inst,365,365)

