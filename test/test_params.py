params = dict(

    # Population size
    N=1000,

    # Parameters relating to duration of infection period, including ID period
    av_I_duration=2,
    av_ID_duration=200 / 7,
    inf_red=0.45,
    min_ID=11,

    # Parameters relating to duration of disease period
    av_D_duration=300 / 7,
    dis_red=0.3,
    min_D=1,

    # Parameters relating to lambda function - calculating force of infection
    v_1=1,
    v_2=2.6,
    phi=1.4,
    epsilon=0.5,

    n_inf_sev = 30,
)

# Demography parameters
demog = dict(
    tau=1 / (40 * 52),  # death rate in weeks^-1
    max_age=60 * 52,    # maximum age in population
    mean_age=20 * 52    # mean age in population
)

bet = 0.165582669212017
