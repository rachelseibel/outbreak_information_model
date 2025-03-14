# Main file running SEIRV_model_heterogeneous over range of parameters from parameters_simulations.csv

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from SEIRV_model_heterogeneous import SEIRV_model_heterogeneous, event, new_day_event, update_vaccine_uptake_rate, alive, foi
from memory_window import memory_window_heterogeneous

# Import simulation parameters
parameters = pd.read_csv('results/heterogeneous/inputs/parameters_scenario1.csv')

# Import initial conditions
initial_conditions = pd.read_csv('results/heterogeneous/inputs/parameters_IC.csv')
# By field name, convert str to int array and append to y0
y0 = []
for i in range(len(initial_conditions.columns)):
    # Split array by comma and remove brackets
    y0_temp = initial_conditions.iloc[0,i].strip(']').strip('[').split(',')
    # Convert str to int
    y0_temp = [float(k) for k in y0_temp]
    # Append to y0
    y0.append(y0_temp)
# Flatten y0
y0 = [item for sublist in y0 for item in sublist]

# Consider splitting the population into groups
num_groups = int(len(y0)/12)

# Time vector
t_start = 0
t_end = 10*365
t_step = 1
t = np.arange(t_start, t_end, t_step)

# Import disease parameters
sarscov2 = pd.read_csv('results/heterogeneous/inputs/parameters_sarscov2.csv')
ebola = pd.read_csv('results/heterogeneous/inputs/parameters_ebola.csv')
influenza = pd.read_csv('results/heterogeneous/inputs/parameters_influenza.csv')

# Initialise results dataframe
df_results = pd.DataFrame(columns=['pathogen','split','memory_window','behaviour_function','r','alpha','vaccine_efficacy','final_cases','final_deaths','final_vaccinated','epidemic_duration'])
df_temporal = pd.DataFrame(columns=['pathogen','split','memory_window','behaviour_function','r','alpha','vaccine_efficacy','t','S','E','I','R','PD','Sv','Ev','Iv','Rv','PDv','C','Cv'])

# Run SEIRV_model over range of parameters
for i in range(len(parameters)):
    # Print progress
    print('Simulation ' + str(i+1) + ' of ' + str(len(parameters)))
    if parameters['pathogen'][i] == 'sarscov2':
        df = sarscov2
        for j in range(len(parameters.columns)):
            df[parameters.columns[j]] = parameters.iloc[i,j]
    elif parameters['pathogen'][i] == 'ebola':
        df = ebola
        for j in range(len(parameters.columns)):
            df[parameters.columns[j]] = parameters.iloc[i,j]
    elif parameters['pathogen'][i] == 'influenza':
        df = influenza
        for j in range(len(parameters.columns)):
            df[parameters.columns[j]] = parameters.iloc[i,j]
    else:
        print('Error: pathogen not found')
    # By field name, convert str to int array and append to y0
    y0 = []
    for x in range(len(initial_conditions.columns)):
        # Split array by comma and remove brackets
        y0_temp = initial_conditions.iloc[0,x].strip(']').strip('[').split(',')
        # Convert str to int
        y0_temp = [float(k) for k in y0_temp]
        # Append to y0
        y0.append(y0_temp)
    # Flatten y0
    y0 = [item for sublist in y0 for item in sublist]
    # Get population split as list
    split = parameters['split'][i].strip(']').strip('[').split(', ')
    # Adjust initial conditions to account for splitting
    for x in range(num_groups):
        y0[x] = float(y0[x]) * float(split[x])
        y0[x+num_groups] = float(y0[x+num_groups]) * float(split[x])
        y0[x+2*num_groups] = float(y0[x+2*num_groups]) * float(split[x])
        y0[x+3*num_groups] = float(y0[x+3*num_groups]) * float(split[x])
        y0[x+4*num_groups] = float(y0[x+4*num_groups]) * float(split[x])
        y0[x+5*num_groups] = float(y0[x+5*num_groups]) * float(split[x])
        y0[x+6*num_groups] = float(y0[x+6*num_groups]) * float(split[x])
        y0[x+7*num_groups] = float(y0[x+7*num_groups]) * float(split[x])
        y0[x+8*num_groups] = float(y0[x+8*num_groups]) * float(split[x])
        y0[x+9*num_groups] = float(y0[x+9*num_groups]) * float(split[x])
        y0[x+10*num_groups] = float(y0[x+10*num_groups]) * float(split[x])
        y0[x+11*num_groups] = float(y0[x+11*num_groups]) * float(split[x])
        if float(split[x]) == 0:
            y0[x+0*num_groups] = 0
            y0[x+num_groups] = 0
            y0[x+2*num_groups] = 0
            y0[x+3*num_groups] = 0
            y0[x+4*num_groups] = 0
            y0[x+5*num_groups] = 0
            y0[x+6*num_groups] = 0
            y0[x+7*num_groups] = 0
            y0[x+8*num_groups] = 0
            y0[x+9*num_groups] = 0
            y0[x+10*num_groups] = 0
            y0[x+11*num_groups] = 0
    # Include four additional columns for memory window
    df['ucases_within_memory'] = [[[1]*num_groups]]
    df['udeaths_within_memory'] = [[[0]*num_groups]]
    df['ucases_memory'] = [[1]]
    df['udeaths_memory'] = [[0]]
    df['cases_within_memory'] = [[[1]*num_groups]]
    df['deaths_within_memory'] = [[[0]*num_groups]]
    df['cases_memory'] = [[1]]
    df['deaths_memory'] = [[0]]
    # Convert df to list by column
    df = df.values.tolist()
    # Solve the SEIRV model
    sol = solve_ivp(SEIRV_model_heterogeneous, [t_start, t_end], y0, t_eval=t, events=[event, new_day_event], args=df, rtol=1e-4, atol=1e-7, max_step=1)
    # Calculate simulation results and save to dataframe by pathogen and behaviour function
    # Calculate final cases and deaths
    final_cases = 0
    final_deaths = 0
    for x in range(num_groups):
        tmp = sol.y[10*num_groups+x][-1] + sol.y[11*num_groups+x][-1]
        final_cases += tmp
        tmp = sol.y[4*num_groups+x][-1] + sol.y[9*num_groups+x][-1]
        final_deaths += tmp

    # Calculate peak cases across all groups
    peak_cases = np.array([0]*num_groups)
    for x in range(num_groups):
        peak_cases[x] = (sol.y[10*num_groups+x] + sol.y[11*num_groups+x]).max()
    peak_cases = peak_cases.sum()

    # Calculate total vaccines administered across all groups
    total_vaccines = 0
    for x in range(num_groups):
        tmp = sol.y[5*num_groups+x][-1] + sol.y[6*num_groups+x][-1] + sol.y[7*num_groups+x][-1] + sol.y[8*num_groups+x][-1] + sol.y[9*num_groups+x][-1]
        total_vaccines += tmp
    
    df_newsim = pd.DataFrame({
        'pathogen': parameters['pathogen'][i],
        'split': parameters['split'][i],
        'memory_window': parameters['memory_window'][i],
        'behaviour_function': parameters['behaviour_function'][i],
        'r': parameters['r'][i],
        'alpha': parameters['alpha'][i],
        'vaccine_efficacy': parameters['vaccine_efficacy'][i],
        'final_cases': final_cases,
        'final_deaths': final_deaths,
        'final_vaccinated': total_vaccines,
        'epidemic_duration': sol.t[-1]
    }, index=[i])
    df_results = pd.concat([df_results, df_newsim], axis=0)

    # Save temporal results to dataframe
    file_size = 0
    for j in range(len(sol.t)):
        df_temporal_new = pd.DataFrame({
            'pathogen': parameters['pathogen'][i],
            'split': parameters['split'][i],
            'memory_window': parameters['memory_window'][i],
            'behaviour_function': parameters['behaviour_function'][i],
            'r': parameters['r'][i],
            'alpha': parameters['alpha'][i],
            'vaccine_efficacy': parameters['vaccine_efficacy'][i],
            't': sol.t[j],
            'S': [[sol.y[0][j], sol.y[1][j], sol.y[2][j]]],
            'E': [[sol.y[3][j], sol.y[4][j], sol.y[5][j]]],
            'I': [[sol.y[6][j], sol.y[7][j], sol.y[8][j]]],
            'R': [[sol.y[9][j], sol.y[10][j], sol.y[11][j]]],
            'PD': [[sol.y[12][j], sol.y[13][j], sol.y[14][j]]],
            'Sv': [[sol.y[15][j], sol.y[16][j], sol.y[17][j]]],
            'Ev': [[sol.y[18][j], sol.y[19][j], sol.y[20][j]]],
            'Iv': [[sol.y[21][j], sol.y[22][j], sol.y[23][j]]],   
            'Rv': [[sol.y[24][j], sol.y[25][j], sol.y[26][j]]],
            'PDv': [[sol.y[27][j], sol.y[28][j], sol.y[29][j]]],
            'C': [[sol.y[30][j], sol.y[31][j], sol.y[32][j]]],
            'Cv': [[sol.y[33][j], sol.y[34][j], sol.y[35][j]]]
        }, index=[file_size])
        df_temporal = pd.concat([df_temporal, df_temporal_new], axis=0)