"""
Golub Oleksandra 856706

Implementazione dell'algoritmo di Gillespie (SSA - Direct Method)
per simulare sistemi di reazioni chimiche stocastiche
"""

import math
import random


def gillespie_direct(initial_state, stoichiometry, propensity_func, max_time):
    """
    Parametri
    ----------
    initial_state: lista di interi
        Stato iniziale del sistema (numero di molecole per specie)
    stoichiometry: lista dei vettori stechiometrici
        Ogni elemento è una lista che rappresenta il cambiamento 
        per ogni specie quando avviene una certa reazione
        dimensione in n_reazioni per n_specie
    propensity_func: 
        Funzione che, dato lo stato corrente (lista di interi),
        restituisce una lista con le propensità di ciascuna reazione
    max_time : float, tempo massimo di simulazione

    Retorna in uscita
    -------
    times : lista di float
        Tempi degli eventi simulati
    states : lista di liste di interi
        Stati del sistema corrispondenti ai tempi in times
    """
    t = 0.0                        # tempo iniziale
    x = initial_state[:]           # copia dello stato
    times = [t]
    states = [x[:]]

    while t < max_time:
        # si calcola la propensità per ogni reazione
        a = propensity_func(x)
        a0 = sum(a)

        # se nessuna reazione è possibile, tutto si termina
        if a0 <= 0.0:
            break

        # si estraggono due numeri casuali uniformi in (0,1)
        r1 = random.random()
        r2 = random.random()

        # tempo fino alla prossima reazione (distribuzione esponenziale)
        tau = (1.0 / a0) * math.log(1.0 / r1)
        t = t + tau
        if t > max_time:
            break

        # si seleziona quale reazione avviene
        threshold = r2 * a0
        cumulative = 0.0
        mu = None
        for i, a_mu in enumerate(a):
            cumulative += a_mu
            if cumulative >= threshold:
                mu = i
                break

        # si aggiorna lo stato con il vettore stechiometrico della reazione scelta
        for s in range(len(x)):
            x[s] += stoichiometry[mu][s]

        # si salva tempo e stato
        times.append(t)
        states.append(x[:])

    return times, states


def example_reversible():
    """
    Esempio del sistema A ⇄ B
      Reazione 1: A -> B con rate k1
      Reazione 2: B -> A con rate k2
    """
    # parametri
    k1 = 1.0   # rate A -> B
    k2 = 0.5   # rate B -> A

    # stato iniziale in [A, B]
    initial_state = [100, 0]

    # vettori stechiometrici per ciascuna reazione
    # reazione 0: A -> B  ->  A: -1, B: +1
    # reazione 1: B -> A  ->  A: +1, B: -1
    stoichiometry = [
        [-1, +1],   # reazione 0
        [+1, -1]    # reazione 1
    ]

    # funzione di propensità
    def propensity(state):
        A, B = state
        return [
            k1 * A,   # A -> B
            k2 * B    # B -> A
        ]

    max_time = 20.0

    times, states = gillespie_direct(initial_state,
                                     stoichiometry,
                                     propensity,
                                     max_time)

    # stampa di alcuni risultati
    print(f"Numero di eventi simulati: {len(times) - 1}")
    print(f"Tempo finale: {times[-1]:.3f}")
    print(f"Stato iniziale: A={initial_state[0]}, B={initial_state[1]}")
    print(f"Stato finale:   A={states[-1][0]}, B={states[-1][1]}")

    print("\nPrimi 10 punti della traiettoria (t, A, B):")
    for t, (A, B) in list(zip(times, states))[:10]:
        print(f"t = {t:8.3f}   A = {A:4d}   B = {B:4d}")


if __name__ == "__main__":
    example_reversible()
