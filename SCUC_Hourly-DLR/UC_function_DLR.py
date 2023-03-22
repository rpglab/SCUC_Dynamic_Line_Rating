### Unit commitment function: full SCUC version
### run UC twice: second run will obtain electricity price
### consider dynamic line rating

### Author: Jin Lu, University of Houston.

### Import
from pyomo.environ import *
#import pyomo.environ as pyo
import pickle
import numpy as np
import os
from formpyomo_UC import *

## function to find the solar gen&bus list
def find_solar_bus(case_inst):
    solar_genidx_list = []
    solar_busnum_list = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] == 'Solar':  # the solar power plants
            solar_genidx_list.append(i)  # solar gen index in the Gen
            bus_num = case_inst.gen.bus[i]
            solar_busnum_list.append(bus_num)
    return solar_genidx_list, solar_busnum_list

## function to find the wind gen&bus list
def find_wind_bus(case_inst):
    wind_genidx_list = []
    wind_busnum_list = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] == 'Wind':  # the wind power plants
            wind_genidx_list.append(i)  # windgen index in the Gen
            bus_num = case_inst.gen.bus[i]
            wind_busnum_list.append(bus_num)
    return wind_genidx_list, wind_busnum_list

# function to build the SCUC with DLR (return pyomo case instance)
def build_UC_full(case_instance,load_b_t,line_l_t):
    ### Prepare and form UC pyomo data file
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    for i in range(case_instance.gentotnum):
        if case_instance.gen.fuel_type[i] != 'Wind':
            if case_instance.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_instance.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
        if case_instance.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC(case_instance, time_totnum, load_b_t, reg_gen_list)    # full UC file 'formpyomo_UC.dat'
    ### Python set
    # hour time set (except first hour)
    Time_23 = []
    for i in range(23):
        Time_23.append(i + 2)  # range is 2-24

    # BaseMVA
    BaseMVA = 100

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    model.u_g_t = Var(model.GEN, model.TIME, domain=Binary) # u_g_t
    model.v_g_t = Var(model.GEN, model.TIME, domain=Binary)  # v_g_t
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Renewable Curtailment Var
    model.rnwcur_b_t = Var(model.BUS, model.TIME,domain=NonNegativeReals)

    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t]*BaseMVA + model.gen_cost_NL[g] * model.u_g_t[g, t] + model.gen_cost_SU[g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj

    model.object = Objective(rule=objfunction, sense=minimize)

    # Renewable curtailment constraint
    def rnwcur_f_1(model, b, t):
        if load_b_t[b - 1][t - 1] >= 0:
            return model.rnwcur_b_t[b, t] == 0
        else:
            return model.rnwcur_b_t[b, t] <= -load_b_t[b - 1][t - 1]/BaseMVA
    model.rnwcur_cons_1 = Constraint(model.BUS, model.TIME, rule=rnwcur_f_1)

    ## Generator power and reserve constraints
    # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g]/BaseMVA * model.u_g_t[g, t] <= model.p_g_t[g, t]
        #return 0 <= model.p_g_t[g, t]
    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g]/BaseMVA * model.u_g_t[g, t]
    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t reserve constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g]/BaseMVA * model.u_g_t[g, t]
    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t reserve constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0
    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right
    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]
    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]
    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * line_l_t[k-1,t-1]/BaseMVA

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * line_l_t[k-1,t-1]/BaseMVA

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]/BaseMVA + model.rnwcur_b_t[b,t]
        return nodal_balance_left == nodal_balance_right
    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]/BaseMVA

    model.gen_rr1_cons = Constraint(model.GEN, Time_23, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]/BaseMVA
    model.gen_rr2_cons = Constraint(model.GEN, Time_23, rule=gen_rr2_f)

    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 1:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=1
        else:
            return model.v_g_t[g, t] >= model.u_g_t[g, t] - model.u_g_t[g, t - 1]
    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo

# function to build the SCUC case with fixed u_g_t (return pyomo case instance)
def build_UC_full_Run2(case_instance,UC_case,load_b_t,line_l_t):
    ### Prepare and form UC pyomo data file
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    for i in range(case_instance.gentotnum):
        if case_instance.gen.fuel_type[i] != 'Wind':
            if case_instance.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_instance.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
        if case_instance.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)

    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC(case_instance, time_totnum, load_b_t, reg_gen_list)
    ### Python set
    # hour time set (except first hour)
    Time_23 = []
    for i in range(23):
        Time_23.append(i + 2)  # range is 2-24

    # BaseMVA
    BaseMVA = 100

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    # model.u_g_t = Var(model.GEN, model.TIME, domain=Binary)
    model.v_g_t = Var(model.GEN, model.TIME,domain=Binary)  # v_g_t is simplified to non-binary variable
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Renwable Curtailment Var
    model.rnwcur_b_t = Var(model.BUS, model.TIME,domain=NonNegativeReals)

    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t]*BaseMVA + model.gen_cost_NL[g] * UC_case.u_g_t[g, t]() + model.gen_cost_SU[g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj
    model.object = Objective(rule=objfunction, sense=minimize)

    # Renewable curtailment constraint
    def rnwcur_f_1(model, b, t):
        if load_b_t[b - 1][t - 1] >= 0:
            return model.rnwcur_b_t[b, t] == 0
        else:
            return model.rnwcur_b_t[b, t] <= -load_b_t[b - 1][t - 1]/BaseMVA
    model.rnwcur_cons_1 = Constraint(model.BUS, model.TIME, rule=rnwcur_f_1)

    ## Generator power and reserve constraints
    # # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g]/BaseMVA * UC_case.u_g_t[g, t]() <= model.p_g_t[g, t]
    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g]/BaseMVA * UC_case.u_g_t[g, t]()
    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t ramping constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g]/BaseMVA * UC_case.u_g_t[g, t]()
    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t ramping constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0
    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right
    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]
    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]
    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * line_l_t[k-1,t-1]/BaseMVA

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * line_l_t[k-1,t-1]/BaseMVA

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]/BaseMVA + model.rnwcur_b_t[b,t]
        return nodal_balance_left == nodal_balance_right
    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]/BaseMVA
    model.gen_rr1_cons = Constraint(model.GEN, Time_23, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]/BaseMVA
    model.gen_rr2_cons = Constraint(model.GEN, Time_23, rule=gen_rr2_f)
    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 1:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=1
        else:
            return model.v_g_t[g, t] >= UC_case.u_g_t[g, t]() - UC_case.u_g_t[g, t - 1]()
    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo

# function to do pyomo solving
def solve_UC(UC_case,pickle_filename,case_nm,day_num):
    # set the solver
    #UC_solver = SolverFactory('glpk', executable='C:\\Users\\jlu27\\Desktop\\glpk-4.65\\w64\\glpsol.exe')
    UC_solver = SolverFactory('conopt',
                                executable='C:\\Users\\lujin\\OneDrive - University Of Houston\\work folder\\ampl new\\ampl_mswin64\\conopt.exe')
    UC_solver.options.mipgap = 0.01
    results = UC_solver.solve(UC_case)
    print('the solution is found')
    # display solution
    print("\nresults.Solution.Status: " + str(results.Solution.Status))
    print("\nresults.solver.status: " + str(results.solver.status))
    print("\nresults.solver.termination_condition: " + str(results.solver.termination_condition))
    print("\nresults.solver.termination_message: " + str(results.solver.termination_message) + '\n')

    # save solution to pickle file
    #pickle.dump(results, open(pickle_filename, "wb"))
    # write result
    write_UCresult_day(UC_case, case_nm,day_num)

### function to write the SCUC results
def write_UCresult_day(UC_case, case_nm, day_num):
    BaseMVA = 100
    ### print power flow result
    flnm_pf = case_nm + '_pf.txt'
    fdpath = os.getcwd() + '\\UC_results' + '\\Day' + str(day_num)
    flpath = fdpath + '\\' + flnm_pf
    # Check whether the specified path exists or not
    isExist = os.path.exists(fdpath)
    if not isExist:
        os.makedirs(fdpath)  # Create a new directory because it does not exist
    f = open(flpath, 'w')
    # print p_k_t
    for k in UC_case.LINE:
        p_k_t_str = ''
        for t in UC_case.TIME:
            if t>=1:    # delete the initial hour
                p_k_t_str = p_k_t_str + str(int(UC_case.p_k_t[k, t]()*BaseMVA)) + ' '
        f.write(p_k_t_str)
        f.write('\n')
    f.write('\n')
    f.close()
    # print p_g_t
    flnm_pgt = case_nm + '_pgt.txt'
    flpath = fdpath + '\\' + flnm_pgt
    f = open(flpath, 'w')
    for g in UC_case.GEN:
        p_g_t_str = ''
        for t in UC_case.TIME:
            if t >= 1:  # delete the initial hour
                p_g_t_str = p_g_t_str + str(int(UC_case.p_g_t[g, t]()*BaseMVA)) + ' '
        f.write(p_g_t_str)
        f.write('\n')
    f.write('\n')
    f.close()
    # print power flow percentage
    flnm_pf = case_nm + '_pfpct.txt'
    flpath = fdpath + '\\' + flnm_pf
    f = open(flpath, 'w')
    for k in UC_case.LINE:
        p_k_t_pct_str = ''
        for t in UC_case.TIME:
            if t >= 1:
                p_k_t = UC_case.p_k_t[k, t]()
                p_k_max = UC_case.line_Pmax[k]
                # print(str(case_pyomo.line_Pmax[k])+' ')
                p_k_t_pct_str = p_k_t_pct_str + str(p_k_t*BaseMVA / p_k_max) + ' '
        f.write(p_k_t_pct_str)
        f.write('\n')
    f.write('\n')
    f.close()
    # print LMP
    flnm_lmp = case_nm + '_lmp.txt'
    flpath = fdpath + '\\' + flnm_lmp
    f = open(flpath, 'w')
    for b in UC_case.BUS:
        lmp_str = ''
        for t in UC_case.TIME:
            if t >= 1:
                nodal_balance_cons = getattr(UC_case, 'nodal_balance_cons')
                lmp = UC_case.dual.get(nodal_balance_cons[b,t])
                lmp_str += str(lmp/BaseMVA) + ' '
        f.write(lmp_str)
        f.write('\n')
    f.write('\n')
    f.close()
    # print operational cost
    flnm_opcost = case_nm + '_Opcost.txt'
    flpath = fdpath + '\\' + flnm_opcost
    f = open(flpath, 'w')
    Opcost_str = str(UC_case.object())
    f.write(Opcost_str)
    f.close()