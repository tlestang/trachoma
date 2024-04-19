import numpy as np
import matplotlib.pyplot as plt


betavals = [0.2, 0.35, 0.4]
popsize = 1024
nrecords = 1144
nbytes = (nrecords * popsize) // 8

for ibeta, beta in enumerate(betavals):
    inf_records_packed = np.fromfile(
        "infected_state.bin",
        dtype=np.uint8,
        count=nbytes,
        offset=ibeta * nbytes,
    )
    ages_records = np.fromfile(
        "ages.bin",
        dtype=np.uint16,
        count=popsize * nrecords,
        offset=ibeta * popsize * nrecords,
    ).reshape(nrecords, popsize)

    rec_size = popsize // 8
    inf_records = [
        np.unpackbits(
            inf_records_packed[i * rec_size:(i + 1) * rec_size]
        )
        for i in range(nrecords)
    ]
    ninf = [
        np.count_nonzero(
            inf[(ages >= 9 * 52) & (ages <= 15 * 52)]
        )
        for inf, ages in zip(inf_records, ages_records)
    ]

    plt.plot(ninf[::2], label=f"beta = {beta}")

plt.title("Number of infected children aged 9 to 15")
plt.xlabel("weeks")
plt.legend()
plt.show()
