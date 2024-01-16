from memory_window import memory_window_heterogeneous
import numpy as np
import math

def SEIRV_model_heterogeneous(t, y, df):

    R0, beta, sigma, gamma, d, td, s1, s2, s3, v, pathogen, split, m, behaviour_function, r, alpha, vaccine_efficacy, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory = df

    # Calculate num_groups
    num_groups = int(len(y)/12)
    # Assign y0 elements to their respective variables as arrays
    S = y[0:(num_groups)]
    E = y[num_groups:(2*num_groups)]
    I = y[2*num_groups:(3*num_groups)]
    R = y[3*num_groups:(4*num_groups)]
    H = y[4*num_groups:(5*num_groups)]
    Sv = y[5*num_groups:(6*num_groups)]
    Ev = y[6*num_groups:(7*num_groups)]
    Iv = y[7*num_groups:(8*num_groups)]
    Rv = y[8*num_groups:(9*num_groups)]
    Hv = y[9*num_groups:(10*num_groups)]
    C = y[10*num_groups:(11*num_groups)]
    Cv = y[11*num_groups:(12*num_groups)]
    
    # Convert r from str to list
    r = np.array(r.strip(']').strip('[').split(', ')).astype(float)

    # Total population
    N = S.sum() + E.sum() + I.sum() + R.sum() + H.sum() + Sv.sum() + Ev.sum() + Iv.sum() + Rv.sum() + Hv.sum()

    # Initialise the differential equations
    dSdt = [0]*num_groups
    dEdt = [0]*num_groups
    dIdt = [0]*num_groups
    dRdt = [0]*num_groups
    dHdt = [0]*num_groups
    dSvdt = [0]*num_groups
    dEvdt = [0]*num_groups
    dIvdt = [0]*num_groups
    dRvdt = [0]*num_groups
    dHvdt = [0]*num_groups
    dCdt = [0]*num_groups
    dCvdt = [0]*num_groups

    if t==0:
        v = [v]*3
        check_split = r
        for i in range(3):
            if check_split[i] == 0:
                v[i] = 0
            else:
                pass
        behaviour, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory = update_vaccine_uptake_rate(t, td, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, C, Cv)

        for i in range(num_groups):
            if (S[i] + E[i] + R[i]) < 10:
                v[i] = 0
            elif S[i] < 10:
                v[i] = 0
            elif (r[i] + alpha) == 0:
                v[i] = 0
            else:
                v[i] = behaviour[i] * 0.005

        # Update parameters
        df[9] = v
        df[17] = ucases_within_memory
        df[18] = udeaths_within_memory
        df[19] = ucases_memory
        df[20] = udeaths_memory
        df[21] = cases_within_memory
        df[22] = deaths_within_memory
        df[23] = cases_memory
        df[24] = deaths_memory

    # Equations for the unvaccinated population
    for i in range(num_groups):
        # Check if elig/aliv_elig is zero or NaN
        if (S[i] < 10) or ((S[i] + E[i] + R[i]) < 10) or (r[i]+alpha == 0):
            dSdt[i] = -foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i]
            dEdt[i] = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i] - sigma*E[i]
            dIdt[i] = sigma*E[i] - gamma*I[i]
            dRdt[i] = (1-d)*gamma*I[i]
            dHdt[i] = d*gamma*I[i]
            # Equation for cumulative infectious
            dCdt[i] = sigma*E[i]
            # dCdt[i] = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i]
            # Equations for the vaccinated population
            dSvdt[i] = -(1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i]
            dEvdt[i] = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i] - sigma*Ev[i]
            dIvdt[i] = sigma*Ev[i] - gamma*Iv[i]
            dRvdt[i] = (1-(1-vaccine_efficacy)*d)*gamma*Iv[i]
            dHvdt[i] = (1-vaccine_efficacy)*d*gamma*Iv[i]
            # Equation for cumulative infectious
            dCvdt[i] = sigma*Ev[i]
            # dCvdt[i] = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i]
        else:
            dSdt[i] = -foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i] - v[i]*S[i]
            dEdt[i] = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i] - sigma*E[i] - v[i]*E[i]
            dIdt[i] = sigma*E[i] - gamma*I[i]
            dRdt[i] = (1-d)*gamma*I[i] - v[i]*R[i]
            dHdt[i] = d*gamma*I[i]
            # Equation for cumulative infectious
            dCdt[i] = sigma*E[i]
            # dCdt[i] = foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*S[i]
            # Equations for the vaccinated population
            dSvdt[i] = -(1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i] + v[i]*S[i]
            dEvdt[i] = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i] - sigma*Ev[i] + v[i]*E[i]
            dIvdt[i] = sigma*Ev[i] - gamma*Iv[i]
            dRvdt[i] = (1-(1-vaccine_efficacy)*d)*gamma*Iv[i] + v[i]*R[i]
            dHvdt[i] = (1-vaccine_efficacy)*d*gamma*Iv[i]
            # Equation for cumulative infectious
            dCvdt[i] = sigma*Ev[i]
            # dCvdt[i] = (1-vaccine_efficacy)*foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)*Sv[i]

    # Create array of differential equations
    diffeqs = [0]*(12*num_groups)
    diffeqs[0:(num_groups)] = dSdt
    diffeqs[num_groups:(2*num_groups)] = dEdt
    diffeqs[2*num_groups:(3*num_groups)] = dIdt
    diffeqs[3*num_groups:(4*num_groups)] = dRdt
    diffeqs[4*num_groups:(5*num_groups)] = dHdt
    diffeqs[5*num_groups:(6*num_groups)] = dSvdt
    diffeqs[6*num_groups:(7*num_groups)] = dEvdt
    diffeqs[7*num_groups:(8*num_groups)] = dIvdt
    diffeqs[8*num_groups:(9*num_groups)] = dRvdt
    diffeqs[9*num_groups:(10*num_groups)] = dHvdt
    diffeqs[10*num_groups:(11*num_groups)] = dCdt
    diffeqs[11*num_groups:(12*num_groups)] = dCvdt

    return diffeqs

# Define the event function
def event(t, y, df):
    # Calculate num_groups
    num_groups = int(len(y)/12)
    # Hetermine the number of active infections (E + Ev + I + Iv)
    active_infs = 0
    for i in range(num_groups):
        active_infs += y[num_groups+i] + y[6*num_groups+i] + y[2*num_groups+i] + y[7*num_groups+i]
    return active_infs - 1

event.terminal = True
event.direction = -1

def new_day_event(t, y, df):
    # Unpack the dataframe
    R0, beta, sigma, gamma, d, td, s1, s2, s3, v, pathogen, split, m, behaviour_function, r, alpha, vaccine_efficacy, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory = df

    # Calculate num_groups
    num_groups = int(len(y)/12)

    # Assign y0 elements to their respective variables as arrays
    S = y[0:(num_groups)]
    E = y[num_groups:(2*num_groups)]
    I = y[2*num_groups:(3*num_groups)]
    R = y[3*num_groups:(4*num_groups)]
    H = y[4*num_groups:(5*num_groups)]
    Sv = y[5*num_groups:(6*num_groups)]
    Ev = y[6*num_groups:(7*num_groups)]
    Iv = y[7*num_groups:(8*num_groups)]
    Rv = y[8*num_groups:(9*num_groups)]
    Hv = y[9*num_groups:(10*num_groups)]
    C = y[10*num_groups:(11*num_groups)]
    Cv = y[11*num_groups:(12*num_groups)]
    
    # Convert r from str to list
    r = np.array(r.strip(']').strip('[').split(', ')).astype(float)

    if (math.floor(t) > math.floor(new_day_event.previous_time) | (t == 0)):
        # Update the vaccine uptake rate at beginning of new day

        num_groups = 3

        # Update the vaccine uptake rate
        if behaviour_function == 'no_vaccine':
            v = [0]*num_groups
        else:
            behaviour, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory = update_vaccine_uptake_rate(t, td, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, C, Cv)

            for i in range(num_groups):
                if (S[i] + E[i] + R[i]) < 10:
                    v[i] = 0
                elif S[i] < 10:
                    v[i] = 0
                elif (r[i] + alpha) == 0:
                    v[i] = 0
                else:
                    v[i] = behaviour[i] * 0.005

    # Update parameters
    df[9] = v
    df[17] = ucases_within_memory
    df[18] = udeaths_within_memory
    df[19] = ucases_memory
    df[20] = udeaths_memory
    df[21] = cases_within_memory
    df[22] = deaths_within_memory
    df[23] = cases_memory
    df[24] = deaths_memory
    new_day_event.previous_time = t
    return 1

new_day_event.previous_time = 0

def update_vaccine_uptake_rate(t, td, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory, m, behaviour_function, r, alpha, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv, C, Cv):

    # Initialise the cases and deaths within each group
    ucases_within = [0]*len(S)
    udeaths_within = [0]*len(S)
    cases_within = [0]*len(S)
    deaths_within = [0]*len(S)

    # Calculate the number of cases and deaths for today within each group
    for i in range(len(S)):
        cases_within[i] = C[i] + Cv[i]
        ucases_within[i] = C[i]
        deaths_within[i] = H[i] + Hv[i]
        udeaths_within[i] = H[i]
    
    # Calculate the number of cases and deaths for today
    cases = C.sum() + Cv.sum()
    ucases = C.sum()
    deaths = H.sum() + Hv.sum()
    udeaths = H.sum()

    # Memory window implementation
    if (m != 0) | (alpha != 0):
        ucases_within, udeaths_within, ucases, udeaths, cases_within, deaths_within, cases, deaths, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory = memory_window_heterogeneous(ucases_within, udeaths_within, ucases, udeaths, cases_within, deaths_within, cases, deaths, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory, m, t, td)

    b = [0]*len(S)
    theta = [0]*len(S)

    for i in range(len(S)):
        # Check behaviour function
        if behaviour_function == 'cases':
            theta[i] = ucases/(0.00001+cases)
            b[i] = r[i] + alpha * theta[i]
        elif behaviour_function == 'deaths':
            theta[i] = udeaths/(0.00001+(deaths))
            b[i] = r[i] + alpha * theta[i]
        elif behaviour_function == 'cases_within':
            theta[i] = ucases_within[i]/(0.00001+cases_within[i])
            b[i] = r[i] + alpha * theta[i]
        elif behaviour_function == 'deaths_within':
            theta[i] = udeaths_within[i]/(0.00001+deaths_within[i])
            b[i] = r[i] + alpha * theta[i]
        elif behaviour_function == 'fixed':
            b[i] = 1
        elif behaviour_function == 'no_vaccine':
            b[i] = 0
        else:
            print('Error: behaviour function not found')
    return b, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory

def alive(S, E, I, R, H, Sv, Ev, Iv, Rv, Hv):
    # Calculate the proportion of the population alive
    aliv = S.sum() + E.sum() + I.sum() + R.sum() + Sv.sum() + Ev.sum() + Iv.sum() + Rv.sum()
    return aliv

def foi(beta, vaccine_efficacy, S, E, I, R, H, Sv, Ev, Iv, Rv, Hv):
    FOI = beta*(I.sum()+Iv.sum())/alive(S, E, I, R, H, Sv, Ev, Iv, Rv, Hv)
    return FOI
