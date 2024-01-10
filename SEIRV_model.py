from memory_window import memory_window
import math

def SEIRV_model(t, y, df):
    R0, beta, sigma, gamma, d, td, s1, s2, s3, v, pathogen, m, behaviour_function, r, alpha, vaccine_efficacy, ucases_memory, udeaths_memory, cases_memory, deaths_memory = df

    if (r+alpha) == 0:
        v = 0

    S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, C, Cv = y
    
    # Total population
    N = S + E + I + R + H + Sv + Ev + Iv + Rv + Hv

    if t==0:
        behaviour, ucases_memory, udeaths_memory, cases_memory, deaths_memory = update_vaccine_uptake_rate(t, td, ucases_memory, udeaths_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N, C, Cv)
        if (S+E+R) < 10:
            v = 0
        elif S < 10:
            v = 0
        elif (r + alpha) == 0:
            v = 0
        else:
            v = 0.005*behaviour
        # Update parameters
        df[9] = v
        df[16] = ucases_memory
        df[17] = udeaths_memory
        df[18] = cases_memory
        df[19] = deaths_memory

    # Check if elig/aliv_elig is zero or NaN
    if (S < 10) or ((S + E + R) < 10) or (r+alpha == 0):
        # Equations for the unvaccinated population
        dSdt = -foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S
        dEdt = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S - sigma*E
        dIdt = sigma*E - gamma*I
        dRdt = (1-d)*gamma*I
        dHdt = d*gamma*I
        # Equation for cumulative infectious
        dCdt = sigma*E
        # dCdt = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S

        # Equations for the vaccinated population
        dSvdt = -(1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv
        dEvdt = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv - sigma*Ev
        dIvdt = sigma*Ev - gamma*Iv
        dRvdt = (1-(1-vaccine_efficacy)*d)*gamma*Iv
        dHvdt = (1-vaccine_efficacy)*d*gamma*Iv
        # Equation for cumulative infectious
        dCvdt = sigma*Ev
        # dCvdt = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv
    else:
         # Equations for the unvaccinated population
        dSdt = -foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S - v*S
        dEdt = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S - sigma*E - v*E
        dIdt = sigma*E - gamma*I
        dRdt = (1-d)*gamma*I - v*R
        dHdt = d*gamma*I
        # Equation for cumulative infectious
        dCdt = sigma*E
        # dCdt = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*S

        # Equations for the vaccinated population
        dSvdt = -(1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv + v*S
        dEvdt = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv - sigma*Ev + v*E
        dIvdt = sigma*Ev - gamma*Iv
        dRvdt = (1-(1-vaccine_efficacy)*d)*gamma*Iv + v*R
        dHvdt = (1-vaccine_efficacy)*d*gamma*Iv
        # Equation for cumulative infectious
        dCvdt = sigma*Ev
        # dCvdt = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N)*Sv
    
    return [dSdt, dEdt, dIdt, dRdt, dHdt, dSvdt, dEvdt, dIvdt, dRvdt, dHvdt, dCdt, dCvdt]

# Define the event function
def event(t, y, df):
    return (y[1]+y[2]+y[6]+y[7]) - 1

event.terminal = True
event.direction = -1

def new_day_event(t, y, df):
    # Unpack the dataframe
    R0, beta, sigma, gamma, d, td, s1, s2, s3, v, pathogen, m, behaviour_function, r, alpha, vaccine_efficacy, ucases_memory, udeaths_memory, cases_memory, deaths_memory = df
    # Unpack y
    S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, C, Cv = y
    # Total population
    N = S + E + I + R + H + Sv + Ev + Iv + Rv + Hv
    if new_day_event.previous_time==0:
        new_day_event.previous_time = t
    if (math.floor(t) > math.floor(new_day_event.previous_time) | (t == 0)):
        # Update the vaccine uptake rate at beginning of new day
        if behaviour_function == 'no_vaccine':
            v = 0
        else:
            behaviour, ucases_memory, udeaths_memory, cases_memory, deaths_memory = update_vaccine_uptake_rate(t, td, ucases_memory, udeaths_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N, C, Cv)
            if (S+E+R) < 10:
                v = 0
            elif S < 10:
                v = 0
            elif (r + alpha) == 0:
                v = 0
            else:
                v = 0.005*behaviour
            # Update parameters
            df[9] = v
            df[16] = ucases_memory
            df[17] = udeaths_memory
            df[18] = cases_memory
            df[19] = deaths_memory
    new_day_event.previous_time = t
    return 1

new_day_event.previous_time = 0

def update_vaccine_uptake_rate(t, td, ucases_memory, udeaths_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N, C, Cv):   
    # Total population
    N = S + E + I + R + H + Sv + Ev + Iv + Rv + Hv
    
    # Calculate the number of cases and hospitalisations
    cases = C + Cv
    ucases = C
    deaths = H + Hv
    udeaths = H

    # if (m != 0) | (alpha != 0):
    # Memory window implementation
    ucases, udeaths, cases, deaths, ucases_memory, udeaths_memory, cases_memory, deaths_memory = memory_window(ucases, udeaths, cases, deaths, ucases_memory, udeaths_memory, cases_memory, deaths_memory, m, t, td)

    # Check behaviour function
    if behaviour_function == 'cases':
        theta = ucases/(0.00001+cases)
        b = r + alpha * theta
    elif behaviour_function == 'deaths':
        theta = udeaths/(0.00001+deaths)
        b = r + alpha * theta
    elif behaviour_function == 'fixed':
        b = 1
    elif behaviour_function == 'no_vaccine':
        b = 0
    else:
        print('Error: behaviour function not found')
    return b, ucases_memory, udeaths_memory, cases_memory, deaths_memory

def alive(S, E, I, R, H, Sv, Ev, Iv, Rv, Hv):
    # Calculate the proportion of the population alive
    aliv = S + E + I + R + Sv + Ev + Iv + Rv
    return aliv

def foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, N):
    FOI = beta*(I+Iv)/alive(S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)
    return FOI