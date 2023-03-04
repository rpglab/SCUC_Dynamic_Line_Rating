# Power System Day-Ahead Generation Scheduling through SCUC or SCUC-DLR
This program implements day-ahead scheduling models: (i) SCUC with daily/constant line rating, and (ii) SCUC with hourly/dynamic line rating.

This set of program and data is a part of <a class="" target="_blank" href="https://rpglab.github.io/resources/TX-123BT/">our synthetic system (TX-123BT) dataset release</a>.

### Environment Setting:
* Recommended Python Version: Python 3.11
* Required packages: Numpy, pyomo, pypower, pickle
* A third-party solver is required, which can be called by pyomo to solve the SCUC problem.

### Code Folders/Files:
* 'Sample_Codes_SCUC' folder: A standard SCUC model.
	* The load, solar generation, wind generation profiles are provided by 'load_annual', 'solar_annual', 'wind_annual' folders.
	* The daily line rating profiles are provided by 'Line_annual_Dmin.txt'.
	* 'power_mod.py': define the python class for the power system.
	* 'UC_function.py': define functions to build, solve, and save results for pyomo SCUC model.
	* 'formpyomo_UC': define the function to create the input file for pyomo model.
	* 'Run_SCUC_annual': run this file to perform SCUC simulation on the selected days of the TX-123BT profiles.
	* Steps to run SCUC simulation:
		* 1) Set up the python environment. 
		* 2) Set the solver location: 'UC_function.py' => 'solve_UC' function => UC_solver=SolverFactory('solver_name',executable='solver_location')
		* 3) Set the days you want to run SCUC: 'Run_SCUC_annual.py' => last row:  run_annual_UC(case_inst,start_day,end_day)
			* For example: to run SCUC simulations for 125th-146th days in 2019, the last row of the file is 'run_annual_UC(case_inst,125,146)'
			* You can also run a single day's SCUC simulation by using: 'run_annual_UC(case_inst,single_day,single_day)'

* 'Sample_Codes_SCUC_HourlyDLR' folder: The SCUC model consider hourly dynamic line rating (DLR) profiles.
	* The load, solar generation, wind generation profiles are provided by 'load_annual','solar_annual', 'wind_annual' folders.
	* The hourly line rating profiles in 2019 are provided by 'dynamic_rating_result' folder.
	* 'power_mod.py': define the python class for the power system.
	* 'UC_function_DLR.py': define functions to build, solve, and save results for pyomo SCUC model (with hourly DLR).
	* 'formpyomo_UC': define the function to create the input file for pyomo model.
	* 'RunUC_annual_dlr': run this file to perform SCUC simulation (with hourly DLR) on the selected days of the TX-123BT profiles.
	* Steps to run SCUC simulation (with hourly DLR):
		* 1) Set up the python environment. 
		* 2) Set the solver location: 'UC_function_DLR.py' => 'solve_UC' function => UC_solver=SolverFactory('solver_name',executable='solver_location')
		* 3) Set the daily profiles for SCUC simulation: 'RunUC_annual_dlr.py' => last row:  run_annual_UC_dlr(case_inst,start_day,end_day)
			* For example: to run SCUC simulations (with hourly DLR) for 125th-146th days in 2019, the last row of the file is 'run_annual_UC_dlr(case_inst,125,146)'
			* You can also run a single day's SCUC simulation (with hourly DLR) by using: 'run_annual_UC_dlr(case_inst,single_day,single_day)'


### Simulation Results:

The SCUC/SCUC with DLR simulation results are saved in the 'UC_results' folders under the corresponding folder. Files 
under 'UC_results' folder are explained below:
* 'UCcase_Opcost.txt': total operational cost ($)
* 'UCcase_pf.txt': the power flow results (MW). Rows represent lines, columns represent hours.
* 'UCcase_pfpct.txt': the percentage of the power flow to the line capacity (%). Rows represent lines, columns represent hours.
* 'UCcase_pgt.txt': the generators output power (MW). Rows represent conventional generators, columns represent hours.
* 'UCcase_lmp.txt': the locational marginal price ($/MWh). Rows represent buses, columns represent hours.


## Citation:
If you use this dataset and/or the attached sample codes in your work, please cite the following paper:

Jin Lu, Xingpeng Li, Hongyi Li, Taher Chegini, Carlos Gamarra, Y. C. Ethan Yang, Margaret Cook, and Gavin Dillingham, “A Synthetic Texas Backbone Power System with Climate-Dependent Spatio-Temporal Correlated Profiles”, *arXiv*, Feb. 2023.

Paper website: <a class="off" href="/papers/JinLu-TX-123BT/"  target="_blank">https://rpglab.github.io/papers/JinLu-TX-123BT/</a>


## Contributions:
Jin Lu created this program and dataset. Xingpeng Li supervised this work.


## Contact:
If you need any techinical support, please feel free to reach out to <a class="" target="_blank" href="https://rpglab.github.io/people/Jin-Lu/">Jin Lu</a> at jlu27@CougarNet.UH.EDU.

For collaboration, please contact <a class="" target="_blank" href="https://rpglab.github.io/people/Xingpeng-Li/">Dr. Xingpeng Li</a> at xli83@central.uh.edu.

Website: https://rpglab.github.io/


## License:
This work is licensed under the terms of the <a class="off" href="https://creativecommons.org/licenses/by/4.0/"  target="_blank">Creative Commons Attribution 4.0 (CC BY 4.0) license.</a>


## Disclaimer:
The author doesn’t make any warranty for the accuracy, completeness, or usefulness of any information disclosed and the author assumes no liability or responsibility for any errors or omissions for the information (data/code/results etc) disclosed.
