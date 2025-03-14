# Main file running SEIRV_model over range of parameters from parameters_simulations.csv

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from SEIRV_model_homogeneous import SEIRV_model_homogeneous, event, new_day_event, update_vaccine_uptake_rate, alive, foi

parameters = pd.read_csv('results/homogeneous/inputs/parameters_scenario1a.csv')

# Import initial conditions
initial_conditions = pd.read_csv('results/homogeneous/inputs/parameters_IC.csv')
y0 = initial_conditions.iloc[:,:].values[0]

# Time vector
t_start = 0
t_end = 10*365
t_step = 1
t = np.arange(t_start, t_end, t_step)

# Import disease parameters
sarscov2 = pd.read_csv('results/homogeneous/inputs/parameters_sarscov2.csv')
ebola = pd.read_csv('results/homogeneous/inputs/parameters_ebola.csv')
influenza = pd.read_csv('results/homogeneous/inputs/parameters_influenza.csv')

# Initialise results dataframe
df_results = pd.DataFrame(columns=['pathogen','memory_window','behaviour_function','r','alpha','vaccine_efficacy','final_cases','final_deaths','final_vaccinated','epidemic_duration'])
df_temporal = pd.DataFrame(columns=['pathogen','memory_window','behaviour_function','r','alpha','vaccine_efficacy','t','S','E','I','R','PD','Sv','Ev','Iv','Rv','PDv','C','Cv'])

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
    # Include additional columns for memory window
    df['ucases_memory'] = [[1]]
    df['udeaths_memory'] = [[0]]
    df['cases_memory'] = [[0]]
    df['deaths_memory'] = [[0]]
    # Convert df to list by column
    df = df.values.tolist()
    # Solve the SEIRV model
    sol = solve_ivp(SEIRV_model_homogeneous, [t_start, t_end], y0, t_eval=t, events=[event,new_day_event], args=df, rtol=1e-4, atol=1e-7, max_step=1)
    # Calculate simulation results and save to dataframe by pathogen and behaviour function
    df_newsim = pd.DataFrame({
        'pathogen': parameters['pathogen'][i],
        'memory_window': parameters['memory_window'][i],
        'behaviour_function': parameters['behaviour_function'][i],
        'r': parameters['r'][i],
        'alpha': parameters['alpha'][i],
        'vaccine_efficacy': parameters['vaccine_efficacy'][i],
        'final_cases': (sol.y[10][-1]+sol.y[11][-1]),
        'final_deaths': (sol.y[4][-1]+sol.y[9][-1]),
        'final_vaccinated': (sol.y[5][-1] + sol.y[6][-1] + sol.y[7][-1] + sol.y[8][-1] + sol.y[9][-1]),
        'epidemic_duration': (sol.t[-1])
    }, index=[i])
    df_results = pd.concat([df_results, df_newsim], axis=0)

    # Save temporal results to dataframe
    file_size = 0
    for j in range(len(sol.t)):
        df_temporal_new = pd.DataFrame({
            'pathogen': parameters['pathogen'][i],
            'memory_window': parameters['memory_window'][i],
            'behaviour_function': parameters['behaviour_function'][i],
            'r': parameters['r'][i],
            'alpha': parameters['alpha'][i],
            'vaccine_efficacy': parameters['vaccine_efficacy'][i],
            't': sol.t[j],
            'S': sol.y[0][j],
            'E': sol.y[1][j],
            'I': sol.y[2][j],
            'R': sol.y[3][j],
            'PD': sol.y[4][j],
            'Sv': sol.y[5][j],
            'Ev': sol.y[6][j],
            'Iv': sol.y[7][j],
            'Rv': sol.y[8][j],
            'PDv': sol.y[9][j],
            'C': sol.y[10][j],
            'Cv': sol.y[11][j]
        }, index=[file_size])
        df_temporal = pd.concat([df_temporal, df_temporal_new], axis=0)
        file_size += 1
