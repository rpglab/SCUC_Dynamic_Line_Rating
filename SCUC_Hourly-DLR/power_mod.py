### Create a Class for matpower/pypower case
'''
GenericModel(power module), Bus, Gen, Branch, Gencost
input is index(start from 0), output is num(start from 1)
GenericModel_add add the power flow for branch
'''
from pypower import loadcase
import numpy as np
import pickle

### classic power system case class
class GenericModel:
    def __init__(self, py_case):
        self.case_dict = loadcase.loadcase(py_case)
        self.baseMVA = self.case_dict['baseMVA']
        self.bus_array = self.case_dict['bus']
        self.gen_array = self.case_dict['gen']
        self.branch_array = self.case_dict['branch']
        if 'gencost' in self.case_dict:
            self.gencost_array = self.case_dict['gencost']
            self.gencost = self.loadgencost(self.gencost_array)
        self.bus = self.loadbus(self.bus_array)
        self.gen = self.loadgen(self.gen_array)
        self.branch = self.loadbranch(self.branch_array)
        self.bustotnum = self.bus_array.shape[0]             #total bus number
        self.gentotnum = self.gen_array.shape[0]  # total bus number
        self.branchtotnum = self.branch_array.shape[0]  # total bus number

    # loadbus function can load the bus array data to the class Bus
    def loadbus(self,bus_data):
        bus_num = bus_data[:, 0].astype(int)
        bus_type = bus_data[:, 1].astype(int)
        bus_Pd = bus_data[:, 2]
        bus_Qd = bus_data[:, 3]
        bus_Gs = bus_data[:, 4]
        bus_Bs = bus_data[:, 5]
        bus_area = bus_data[:, 6]
        bus_Vm = bus_data[:, 7]
        bus_Va = bus_data[:, 8]
        bus_baseKV = bus_data[:, 9]
        bus_zone = bus_data[:, 10]
        bus_Vmax = bus_data[:, 11]
        bus_Vmin = bus_data[:, 12]
        bus_loaded = Bus(bus_num, bus_type, bus_Pd, bus_Qd, bus_Gs, bus_Bs, bus_area, bus_Vm, bus_Va, bus_baseKV,
                         bus_zone, bus_Vmax, bus_Vmin)
        return bus_loaded

    def loadgen(self,gen_data):
        gen_bus = gen_data[:, 0].astype(int)
        gen_Pg = gen_data[:, 1]
        gen_Qg = gen_data[:, 2]
        gen_Qmax = gen_data[:, 3]
        gen_Qmin = gen_data[:, 4]
        gen_Vg = gen_data[:, 5]
        gen_mBase = gen_data[:, 6]
        gen_status = gen_data[:, 7]
        gen_Pmax = gen_data[:, 8]
        gen_Pmin = gen_data[:, 9]
        gen_Pc1 = gen_data[:, 10]
        gen_Pc2 = gen_data[:, 11]
        gen_Qc1min = gen_data[:, 12]
        gen_Qc1max = gen_data[:, 13]
        gen_Qc2min = gen_data[:, 14]
        gen_Qc2max = gen_data[:, 15]
        gen_ramp_agc = gen_data[:, 16]
        gen_ramp_10 = gen_data[:, 17]
        gen_ramp_30 = gen_data[:, 18]
        gen_ramp_q = gen_data[:, 19]
        gen_apf = gen_data[:, 20]
        gen_loaded = Gen(gen_bus, gen_Pg, gen_Qg, gen_Qmax, gen_Qmin, gen_Vg, gen_mBase, gen_status, gen_Pmax, gen_Pmin,
                         gen_Pc1, gen_Pc2, gen_Qc1min, gen_Qc1max, gen_Qc2min, gen_Qc2max, gen_ramp_agc, gen_ramp_10,
                         gen_ramp_30, gen_ramp_q, gen_apf)
        return gen_loaded

    def loadbranch(self,branch_data):
        branch_fbus = branch_data[:, 0].astype(int)
        branch_tbus = branch_data[:, 1].astype(int)
        branch_r = branch_data[:, 2]
        branch_x = branch_data[:, 3]
        branch_b = branch_data[:, 4]
        branch_rateA = branch_data[:, 5]
        branch_rateB = branch_data[:, 6]
        branch_rateC = branch_data[:, 7]
        branch_ratio = branch_data[:, 8]
        branch_angle = branch_data[:, 9]
        branch_status = branch_data[:, 10]
        branch_angmin = branch_data[:, 11]
        branch_angmax = branch_data[:, 12]
        branch_loaded = Branch(branch_fbus, branch_tbus, branch_r, branch_x, branch_b, branch_rateA, branch_rateB,
                               branch_rateC, branch_ratio, branch_angle, branch_status, branch_angmin, branch_angmax)
        return branch_loaded

    def loadgencost(self,gencost_data):
        gencost_type = gencost_data[:, 0]
        gencost_startup = gencost_data[:, 1]
        gencost_shutdown = gencost_data[:, 2]
        gencost_loaded = Gencost(gencost_type, gencost_startup, gencost_shutdown)
        return gencost_loaded

class Bus:
    def __init__(self, num, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, Vmax, Vmin):
        self.num = num
        self.type = type
        self.Pd = Pd
        self.Qd = Qd
        self.Gs = Gs
        self.Bs = Bs
        self.area = area
        self.Vm = Vm
        self.Va = Va
        self.baseKV = baseKV
        self.zone = zone
        self.Vmax = Vmax
        self.Vmin = Vmin

class Gen:
    def __init__(self,bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2, Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf):
        self.bus = bus
        self.Pg = Pg
        self.Qg = Qg
        self.Qmax = Qmax
        self.Qmin = Qmin
        self.Vg = Vg
        self.mBase = mBase
        self.status = status
        self.Pmax = Pmax
        self.Pmin = Pmin
        self.Pc1 = Pc1
        self.Pc2 = Pc2
        self.Qc1min = Qc1min
        self.Qc1max = Qc1max
        self.Qc2min = Qc2min
        self.Qc2max = Qc2max
        self.ramp_agc = ramp_agc
        self.ramp_10 = ramp_10
        self.ramp_30 = ramp_30
        self.ramp_q = ramp_q
        self.apf = apf

class Branch:
    def __init__(self, fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax):
        self.fbus = fbus
        self.tbus = tbus
        self.r = r
        self.x = x
        self.b = b
        self.rateA = rateA
        self.rateB = rateB
        self.rateC = rateC
        self.ratio = ratio
        self.angle = angle
        self.status = status
        self.angmin = angmin
        self.angmax = angmax

class Gencost:
    def __init__(self, type, startup, shutdown):
        self.type = type
        self.startup = startup
        self.shutdown = shutdown

### class add the power flow data for branch
class GenericModel_add(GenericModel):
    def loadbranch(self,branch_data):
        branch_fbus = branch_data[:, 0].astype(int)
        branch_tbus = branch_data[:, 1].astype(int)
        branch_r = branch_data[:, 2]
        branch_x = branch_data[:, 3]
        branch_b = branch_data[:, 4]
        branch_rateA = branch_data[:, 5]
        branch_rateB = branch_data[:, 6]
        branch_rateC = branch_data[:, 7]
        branch_ratio = branch_data[:, 8]
        branch_angle = branch_data[:, 9]
        branch_status = branch_data[:, 10]
        branch_angmin = branch_data[:, 11]
        branch_angmax = branch_data[:, 12]
        branch_Pf = branch_data[:, 13]
        branch_Qf = branch_data[:, 14]
        branch_Pt = branch_data[:, 15]
        branch_Qt = branch_data[:, 16]
        branch_loaded = Branch_add(branch_fbus, branch_tbus, branch_r, branch_x, branch_b, branch_rateA, branch_rateB,
                               branch_rateC, branch_ratio, branch_angle, branch_status, branch_angmin, branch_angmax,branch_Pf,branch_Qf,branch_Pt,branch_Qt)
        return branch_loaded

class Branch_add(Branch):
    def __init__(self, fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax,Pf,Qf,Pt,Qt):
        self.fbus = fbus
        self.tbus = tbus
        self.r = r
        self.x = x
        self.b = b
        self.rateA = rateA
        self.rateB = rateB
        self.rateC = rateC
        self.ratio = ratio
        self.angle = angle
        self.status = status
        self.angmin = angmin
        self.angmax = angmax
        self.Pf = Pf
        self.Qf = Qf
        self.Pt = Pt
        self.Qt = Qt

### the case class used for backbone case power flow verification
### parameters for pf: PV/PQ bus is initialized, load is the aggregated loads from the 2000-bus case
class power_pf_mod(GenericModel):
    def __init__(self, bus_data, line_data, gen_data, load_data):
        self.bus = power_pf_mod.form_bus_data(self, bus_data, load_data)
        self.branch = power_pf_mod.form_branch_data(self, line_data)
        self.gen = power_pf_mod.form_gen_data(self, gen_data)
        self.bustotnum = len(self.bus.num)
        self.branchtotnum = len(self.branch.fbus)
        self.gentotnum = len(self.gen.bus)

    def form_bus_data(self,bus_data,load_data):
        num = bus_data.bus_num              ### bus_num for bckbn_case is the num for full_case (backbone_argis_data line 21)
        type = np.ones(bus_data.bus_totnum) ### all bus is set to PQ bus, the first bus is set to PV bus and will convert to reference bus (pypower.bustypes line 47)
        Pd = np.zeros(bus_data.bus_totnum)
        Qd = np.zeros(bus_data.bus_totnum)
        for i in range(load_data.load_totnum):
            bus_num = load_data.bus_num[i]
            for j in range(bus_data.bus_totnum):
                b_num = num[j]
                if b_num==bus_num:
                    bus_index = j
            Pd[bus_index] += load_data.MW_value[i]      # aggregate the Active load on a bus
            Qd[bus_index] += load_data.Mvar_value[i]    # aggregate the Reactive load on a bus
        Gs = np.zeros(bus_data.bus_totnum)
        Bs = np.zeros(bus_data.bus_totnum)
        area = np.zeros(bus_data.bus_totnum)
        Vm = np.ones(bus_data.bus_totnum)       # Initial voltage is 1 pu
        Va = np.full( (bus_data.bus_totnum), 0)       # Initial angle is 0
        baseKV = bus_data.nominal_KV
        zone = np.zeros(bus_data.bus_totnum)
        Vmax = np.full( (bus_data.bus_totnum), 1.06)    # max voltage set to 1.06
        Vmin = np.full( (bus_data.bus_totnum), 0.94)    # min voltage set to 0.94
        bus = Bus(num, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, Vmax, Vmin)
        return bus

    def form_branch_data(self,line_data):
        fbus = line_data.fbus_num
        tbus = line_data.tbus_num
        r = line_data.r
        x = line_data.x
        b = line_data.b
        rateA = line_data.MVA_limit
        rateB = np.zeros(line_data.line_totnum)
        rateC = np.zeros(line_data.line_totnum)
        ratio = np.zeros(line_data.line_totnum)
        angle = np.zeros(line_data.line_totnum)
        status = np.ones(line_data.line_totnum)
        angmin = np.full((line_data.line_totnum),-360)      # minimum angle difference set to -360
        angmax = np.full((line_data.line_totnum),360)      # maximum angle difference set to 360
        branch = Branch(fbus,tbus,r,x,b,rateA,rateB,rateC,ratio,angle,status,angmin,angmax)
        return branch

    def form_gen_data(self,gen_data):
        bus = gen_data.bus_num
        Pg = gen_data.P_default
        Qg = np.zeros(gen_data.gen_totnum)      ### Reactive output is set to 0
        Qmax = gen_data.Q_max
        Qmin = gen_data.Q_min
        Pmax = gen_data.P_max
        Pmin = gen_data.P_min
        #Vg = gen_data.V_set         ### voltage set point use origin data
        Vg = np.full((gen_data.gen_totnum), 1.02)  ### voltage set point use origin data
        mBase = np.full( (gen_data.gen_totnum), 100)          # MVA base is 100
        status = np.ones(gen_data.gen_totnum)
        Pc1 = np.zeros(gen_data.gen_totnum)
        Pc2 = np.zeros(gen_data.gen_totnum)
        Qc1min = np.zeros(gen_data.gen_totnum)
        Qc1max = np.zeros(gen_data.gen_totnum)
        Qc2min = np.zeros(gen_data.gen_totnum)
        Qc2max = np.zeros(gen_data.gen_totnum)
        ramp_agc = np.zeros(gen_data.gen_totnum)
        ramp_10 = np.zeros(gen_data.gen_totnum)
        ramp_30 = np.zeros(gen_data.gen_totnum)
        ramp_q = np.zeros(gen_data.gen_totnum)
        apf = np.zeros(gen_data.gen_totnum)
        gen = Gen(bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2, Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf)
        return gen

### Texas case class with more data(latitude,longitude,bus_num_convert,bus_other_data)
## origin bus/branch numbers are converted
class Texas_bckbn_class(GenericModel_add):
    def __init__(self,bckbn_bus,bckbn_branch_index,GenericModel_add_case):
        # same data in GenericModel_add
        self.bus = GenericModel_add_case.bus
        self.gen = GenericModel_add_case.gen
        self.branch = GenericModel_add_case.branch
        # total num
        self.bustotnum = len(self.bus.num)
        self.branchtotnum = len(self.branch.fbus)
        self.gentotnum = len(self.gen.Pg)
        # add data through bckbn_bus data class
        self.bus.latitude = bckbn_bus.latitude
        self.bus.longitude = bckbn_bus.longitude
        self.bus.name = bckbn_bus.bus_name
        self.bus.area_name = bckbn_bus.area_name
        self.bus.sub_num = bckbn_bus.sub_num
        self.bus.nominal_KV = bckbn_bus.nominal_KV
        # save the origin bus_num to num_convert
        self.bus.num_convert = bckbn_bus.bus_num
        # change the bus_num to num start from 1
        bus_num_list = []
        for i in range(self.bustotnum):
            bus_num = i + 1
            bus_num_list.append(bus_num)
        self.bus.num = bus_num_list
        # save the origin branch index to index_convert
        self.branch.index_convert = bckbn_branch_index
        # change fbus tbus gen_bus to num start from 1
        branch_fbus_list = []
        branch_tbus_list = []
        for i in range(self.branchtotnum):
            fbus_origin_num = self.branch.fbus[i]
            tbus_origin_num = self.branch.tbus[i]
            for j in range(self.bustotnum):
                if self.bus.num_convert[j]==fbus_origin_num:
                    fbus_num = j + 1
                if self.bus.num_convert[j]==tbus_origin_num:
                    tbus_num = j + 1
            branch_fbus_list.append(fbus_num)
            branch_tbus_list.append(tbus_num)
        self.branch.fbus = branch_fbus_list
        self.branch.tbus = branch_tbus_list
        # convert the origin gen bus num to new gen bus num
        gen_bus_list = []
        for i in range(self.gentotnum):
            genbus_origin_num = self.gen.bus[i]
            for j in range(self.bustotnum):
                if self.bus.num_convert[j]==genbus_origin_num:
                    genbus_num = j + 1
            gen_bus_list.append(genbus_num)
        self.gen.bus = gen_bus_list

### Texas case class with more data(latitude,longitude,bus_num_convert,bus_other_data)
class Texas_aggregate_class(GenericModel_add):
    def __init__(self,GenericModel_add_case,aggregate_bus):
        # same data in GenericModel_add
        self.bus = GenericModel_add_case.bus
        self.gen = GenericModel_add_case.gen
        self.branch = GenericModel_add_case.branch
        # total num
        self.bustotnum = len(self.bus.num)
        self.branchtotnum = len(self.branch.fbus)
        self.gentotnum = len(self.gen.Pg)
        # add data through aggregate_bus data class
        self.bus.latitude = aggregate_bus.latitude
        self.bus.longitude = aggregate_bus.longitude
        self.bus.name = aggregate_bus.bus_name
        self.bus.area_name = aggregate_bus.area_name
        self.bus.sub_num = aggregate_bus.sub_num
        self.bus.nominal_KV = aggregate_bus.nominal_KV
        # save the origin bus_num to num_convert
        self.bus.num_convert = self.bus.num
        # change the bus_num to num start from 1
        bus_num_list = []
        for i in range(self.bustotnum):
            bus_num = i + 1
            bus_num_list.append(bus_num)
        self.bus.num = bus_num_list
        # change fbus tbus gen_bus to num start from 1
        branch_fbus_list = []
        branch_tbus_list = []
        for i in range(self.branchtotnum):
            fbus_origin_num = self.branch.fbus[i]
            tbus_origin_num = self.branch.tbus[i]
            for j in range(self.bustotnum):
                if self.bus.num_convert[j]==fbus_origin_num:
                    fbus_num = j + 1
                if self.bus.num_convert[j]==tbus_origin_num:
                    tbus_num = j + 1
            branch_fbus_list.append(fbus_num)
            branch_tbus_list.append(tbus_num)
        self.branch.fbus = branch_fbus_list
        self.branch.tbus = branch_tbus_list
        gen_bus_list = []
        for i in range(self.gentotnum):
            genbus_origin_num = self.gen.bus[i]
            for j in range(self.bustotnum):
                if self.bus.num_convert[j]==genbus_origin_num:
                    genbus_num = j + 1
            gen_bus_list.append(genbus_num)
        self.gen.bus = gen_bus_list


### Data Structure for Texas case
## define the substation class
class Substation:
    def __init__(self, sub_num,sub_name,area_name,latitude,longitude,max_KV):
        self.sub_totnum = len(sub_num)
        self.sub_num = sub_num
        self.sub_name = sub_name
        self.area_name = area_name
        self.latitude = latitude
        self.longitude = longitude
        self.max_KV = max_KV

## define the bus data class
class Bus_data:
    def __init__(self, bus_num,bus_name,area_name,sub_num,nominal_KV):
        self.bus_totnum = len(bus_num)
        self.bus_num = bus_num
        self.bus_name = bus_name
        self.area_name = area_name
        self.sub_num = sub_num
        self.nominal_KV = nominal_KV

## define the bus data class for backbone case
class Bckbn_Bus_data:
    def __init__(self, bus_num,bus_name,area_name,sub_num,nominal_KV,bus_latitude,bus_longitude):
        self.bus_totnum = len(bus_num)
        self.bus_num = bus_num
        self.bus_name = bus_name
        self.area_name = area_name
        self.sub_num = sub_num
        self.nominal_KV = nominal_KV
        self.latitude = bus_latitude
        self.longitude = bus_longitude

## define the line data class
class Line_data:
    def __init__(self, fbus_num,tbus_num,circuit_num,r,x,b,MVA_limit):
        self.line_totnum = len(fbus_num)
        self.fbus_num = fbus_num
        self.tbus_num = tbus_num
        self.circuit_num = circuit_num
        self.r = r
        self.x = x
        self.b = b
        self.MVA_limit = MVA_limit

## define the Load data class
class Load_data:
    def __init__(self, bus_num, MW_value, Mvar_value):
        self.load_totnum = len(bus_num)
        self.bus_num = bus_num
        self.MW_value = MW_value
        self.Mvar_value = Mvar_value

### define the Gen data class
class Gen_data:
    def __init__(self, bus_num, P_min, P_max, Q_min, Q_max, P_default, V_set):
        self.gen_totnum = len(bus_num)
        self.bus_num = bus_num
        self.P_min = P_min
        self.P_max = P_max
        self.Q_min = Q_min
        self.Q_max = Q_max
        self.P_default = P_default
        self.V_set = V_set

# function to load the pickle file
# remember to import power_mod when load case
def load_object(filename):
    with open(filename, 'rb') as input:
        case_object = pickle.load(input)
    return case_object