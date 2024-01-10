import numpy as np
from operator import add

def memory_window(ucases, udeaths, cases, deaths, ucases_memory, udeaths_memory, cases_memory, deaths_memory, m, t, td):
    # Append the cases and deaths to the cases_memory and deaths_memory
    ucases_memory.append(ucases)
    cases_memory.append(cases)
    udeaths_memory.append(udeaths)
    deaths_memory.append(deaths)

    if ((t > (m)) and (m > 0)):
        # Check the cases and deaths m days ago
        ucases_m = ucases_memory[-1]-ucases_memory[-(m+1)]
        cases_m = cases_memory[-1]-cases_memory[-(m+1)]
        if t > (m+td+1):
            udeaths_m = udeaths_memory[-1-td]-udeaths_memory[-(m+1)-td]
            deaths_m = deaths_memory[-1-td]-deaths_memory[-(m+1)-td]
        elif (t <= (m+td+1)) & (t > (td+1)):
            udeaths_m = udeaths_memory[-1-td]
            deaths_m = deaths_memory[-1-td]
        else:
            udeaths_m = 0
            deaths_m = 0
    else:
        ucases_m = ucases
        cases_m = cases
        if t > (td+1):
            udeaths_m = udeaths_memory[-1-td]
            deaths_m = deaths_memory[-1-td]
        else:
            udeaths_m = 0
            deaths_m = 0

    return ucases_m, udeaths_m, cases_m, deaths_m, ucases_memory, udeaths_memory, cases_memory, deaths_memory

def memory_window_heterogeneous(ucases_within, udeaths_within, ucases, udeaths, cases_within, deaths_within, cases, deaths, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory, m, t, td):
    
    # Append the cases to the cases_memory
    ucases_within_memory.append(ucases_within)
    cases_within_memory.append(cases_within)
    ucases_memory.append(ucases)
    cases_memory.append(cases)
    # Append the deaths to the deaths_memory
    udeaths_within_memory.append(udeaths_within)
    deaths_within_memory.append(deaths_within)
    udeaths_memory.append(udeaths)
    deaths_memory.append(deaths)

    groups = 3

    if ((t > (m)) and (m > 0)):

        ucases_within_m = [0]*groups
        cases_within_m = [0]*groups
        for y in range(groups):
            ucases_within_m[y] = ucases_within_memory[-1][y]-ucases_within_memory[-m-1][y]
            cases_within_m[y] = cases_within_memory[-1][y]-cases_within_memory[-m-1][y]
        ucases_m = ucases_memory[-1]-ucases_memory[-m-1]
        cases_m = cases_memory[-1]-cases_memory[-m-1]

        if t > (m+td+1):
            udeaths_within_m = [0]*groups
            deaths_within_m = [0]*groups
            for y in range(groups):
                udeaths_within_m[y] = udeaths_within_memory[-1-td][y]-udeaths_within_memory[-m-1-td][y]
                deaths_within_m[y] = deaths_within_memory[-1-td][y]-deaths_within_memory[-m-1-td][y]
            udeaths_m = udeaths_memory[-1-td]-udeaths_memory[-m-1-td]
            deaths_m = deaths_memory[-1-td]-deaths_memory[-m-1-td]
        elif (t <= (m+td+1)) & (t > (td+1)):
            udeaths_within_m = [0]*groups
            deaths_within_m = [0]*groups
            for y in range(groups):
                udeaths_within_m[y] = udeaths_within_memory[-1-td][y]
                deaths_within_m[y] = deaths_within_memory[-1-td][y]
            udeaths_m = udeaths_memory[-1-td]
            deaths_m = deaths_memory[-1-td]
        else:
            udeaths_within_m = [0]*groups
            deaths_within_m = [0]*groups
            udeaths_m = 0
            deaths_m = 0

    else:
        ucases_within_m = ucases_within
        ucases_m = ucases
        cases_within_m = cases_within
        cases_m = cases
        udeaths_within_m = [0]*groups
        deaths_within_m = [0]*groups
        if t > (td+1):
            for y in range(groups):
                udeaths_within_m[y] = udeaths_within_memory[-1-td][y]
                deaths_within_m[y] = deaths_within_memory[-1-td][y]
            udeaths_m = udeaths_memory[-1-td]
            deaths_m = deaths_memory[-1-td]
        else:
            udeaths_within_m = [0]*groups
            deaths_within_m = [0]*groups
            udeaths_m = 0
            deaths_m = 0

    return ucases_within_m, udeaths_within_m, ucases_m, udeaths_m, cases_within_m, deaths_within_m, cases_m, deaths_m, ucases_within_memory, udeaths_within_memory, ucases_memory, udeaths_memory, cases_within_memory, deaths_within_memory, cases_memory, deaths_memory