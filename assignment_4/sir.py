#######################################################
# Oleksandra Golub (856706)
#
# Modello SIR (Susceptible-Infected-Recovered) in Python
# Integrazione numerica con metodo di Eulero (for-loop)
#
# dS/dt = -beta * S * I / N
# dI/dt =  beta * S * I / N - gamma * I
# dR/dt =  gamma * I
#######################################################

import os
import numpy as np
import matplotlib.pyplot as plt

# -----------
# Parametri
# -----------

N = 1000        # popolazione totale
beta = 0.5      # tasso di trasmissione
gamma = 0.1     # tasso di guarigione

I0 = 1          # infetti iniziali
R0 = 0          # guariti iniziali
S0 = N - I0 - R0  # suscettibili iniziali

t_max = 160     # tempo totale simulazione (giorni)
dt = 0.1        # passo temporale (Eulero)
times = np.arange(0, t_max + dt, dt)  # 0..t_max con passo dt

R0_basic = beta / gamma
print("R0 =", R0_basic)

# ------------------
# Funzione derivate
# ------------------

def sir_deriv(S, I, R, beta, gamma, N):
    dS = -beta * S * I / N
    dI =  beta * S * I / N - gamma * I
    dR =  gamma * I
    return dS, dI, dR

# ---------------------------
# Eulero (integrazione)
# ---------------------------

S = np.zeros(len(times))
I = np.zeros(len(times))
R = np.zeros(len(times))

S[0] = S0
I[0] = I0
R[0] = R0

for k in range(1, len(times)):
    dS, dI, dR = sir_deriv(S[k-1], I[k-1], R[k-1], beta, gamma, N)

    S[k] = S[k-1] + dS * dt
    I[k] = I[k-1] + dI * dt
    R[k] = R[k-1] + dR * dt

# ----------
# Risultati
# ----------

peak_idx = int(np.argmax(I))
print("Picco infetti =", round(I[peak_idx]), "a t =", round(times[peak_idx], 2), "giorni")
print("Recovered finali =", round(R[-1]))
print("Suscettibili finali =", round(S[-1]))

# ----------------------
# Grafico + salvataggio
# ----------------------

os.makedirs("outputs", exist_ok=True)

plt.figure(figsize=(10, 5))
plt.plot(times, S, color="blue",  lw=2, label="Suscettibili (S)")
plt.plot(times, I, color="red",   lw=2, label="Infetti (I)")
plt.plot(times, R, color="green", lw=2, label="Guariti (R)")

plt.title("Modello SIR (metodo di Eulero)")
plt.xlabel("Tempo (giorni)")
plt.ylabel("Numero di individui")
plt.legend(loc="right")
plt.tight_layout()

out_path = "outputs/sir_python_minimal.png"
plt.savefig(out_path, dpi=150)
plt.close()

print("Grafico salvato in:", out_path)
