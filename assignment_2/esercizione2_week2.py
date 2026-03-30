# =========================
# Golub Oleksandra 856706
#
# Codice che crea un simulatore stocastico della popolazione di cellule staminali
#
# - Si crea un modello a eventi discreti in tempo continuo basato su dadi stocastici
# - Simula l’evoluzione di alcuni tipi di cellule secondo i tassi di divisione, morte
#      e differenziazione (sono definiti in __init__ b dentro self.k, usati in _exp_time + nelle chiamate in _step, 
#      alla fine danno effetto nella scelta del prossimo evento (in base al minimo tempo) e nell’aggiornamento delle popolazioni) 
# - Restituisce la traccia temporale delle popolazioni
# =========================

import numpy as np

# Tipi di cellula:
# S è cellula staminale
# P è cellula progenitrice
# M è macrofago
# G è granulocita 
STEM, PROG, MACRO, GRAN = "S", "P", "M", "G"

class StemCellSim:     
    # i valori dei tassi e delle popolazioni iniziali sono scelti arbitrariamente per scopi di simulazione
    def __init__(self, Ns0=100, Np0=50,        # si crea un oggetto simulatore con le condizioni iniziali (Ns0, Np0) e e tassi di eventi (vari k...)
                 k_2S=0.3, k_AS=0.5, k_2P=0.15, k_die=0.05,
                 k_P2M=0.4, k_P2G=0.6):
        self.N = {STEM: Ns0, PROG: Np0, MACRO: 0, GRAN: 0}    # qua c'è il numero di cellule di ogni tipo
        self.k = {"2S": k_2S, "AS": k_AS, "2P": k_2P, "DIE": k_die,   # qua ci sono i tassi di reazione ovvero quanto rapidamente ogni evento può accadere
                  "P2M": k_P2M, "P2G": k_P2G}
        self.t = 0.0          # si inizializza il tempo e la prima riga dello stato iniziale
        self.trace = [(self.t, self.N[STEM], self.N[PROG], self.N[MACRO], self.N[GRAN])]

    # =========================
    #  parte WHEN
    # =========================

    def _exp_time(self, rate, pop):   # si calcola quando avverrà un certo evento (pop è popolazione corrente)
        return np.inf if pop <= 0 or rate <= 0 else np.random.exponential(1.0/(rate*pop))
     # se non ci sono cellule (pop ≤ 0), ritorna infinito ovvero ci risulta che l'evento è impossibile
     # altrimenti, estrae un tempo casuale esponenziale con media 1/(rate*pop)

     # A COSA SERVE FARE TUTTO CIò:
     # ogni cellula ha probabilità proporzionale al tasso, percio
     # il sistema “lancia i dadi” per vedere quando succederà il prossimo evento

    def _step(self):       # qua viene eseguito un singolo evento stocastico
        
        # tempi per eventi da S ovvero si calcola il tempo casuale per ciascun evento possibile delle cellule staminali
        t_2S = self._exp_time(self.k["2S"], self.N[STEM])
        t_AS = self._exp_time(self.k["AS"], self.N[STEM])
        t_2P = self._exp_time(self.k["2P"], self.N[STEM])
        t_die = self._exp_time(self.k["DIE"], self.N[STEM])
        
        # tempi per eventi da P ovvero che si fa stessa cosa er i progenitori
        t_P2M = self._exp_time(self.k["P2M"], self.N[PROG])
        t_P2G = self._exp_time(self.k["P2G"], self.N[PROG])

        times = {    # qua avviene il collegamento tra ogni evento con il suo tempo di accadimento
            ("S","2S"): t_2S, ("S","AS"): t_AS, ("S","2P"): t_2P, ("S","DIE"): t_die,
            ("P","P2M"): t_P2M, ("P","P2G"): t_P2G
        }
        evt, dt = min(times.items(), key=lambda kv: kv[1])    # qua si cerca l’evento con il tempo minimo cioè quello che avviene per primo
        
        if not np.isfinite(dt):    # se tutti i tempi erano infiniti, non ci sono eventi possibili, quindi la simulazione si ferma
            return False  # nessun evento abilitato

        # avanza il tempo
        self.t += dt       # si avanza il tempo globale di dt
        src, kind = evt    # src = tipo di cellula, kind = tipo di evento

        # =========================
        #  parte HOW
        # =========================
        
        # si aggiorna il numero di cellule in base all’evento scelto
        if src == "S":
            if kind == "2S":      # S -> 2S
                self.N[STEM] += 1
            elif kind == "AS":    # S -> S + P
                self.N[PROG] += 1
            elif kind == "2P":    # S -> 2P
                self.N[STEM] -= 1
                self.N[PROG] += 2
            elif kind == "DIE":   # S -> Ø
                self.N[STEM] -= 1
        else:  # stessa cosa per i progenitori
            if kind == "P2M":     # P -> M
                self.N[PROG] -= 1
                self.N[MACRO] += 1
            elif kind == "P2G":   # P -> G
                self.N[PROG] -= 1
                self.N[GRAN] += 1

        # pezzo per evitare conteggi negativi
        for k in self.N:
            if self.N[k] < 0: self.N[k] = 0

        self.trace.append((self.t, self.N[STEM], self.N[PROG], self.N[MACRO], self.N[GRAN]))
        return True  # si registra il nuovo stato (tempo e popolazioni) e viene subito indicato che lo step è riuscito


    # il programma cicla finche:
    # non raggiunge t_max (tempo massimo)
    # non supera max_events (numero massimo di eventi)
    # o finché ci sono cellule staminali vive
    def run(self, t_max=50.0, max_events=5000):   # 
        n = 0
        while self.t < t_max and n < max_events and self.N[STEM] > 0:
            if not self._step():    # chiama _step() a ogni iterazione e accumula la traccia dell’evoluzione
                break
            n += 1
        return self.trace

if __name__ == "__main__":      # si inizia la simulazione con seme fisso per riproducibilit
    np.random.seed(0)  # opzionale
    sim = StemCellSim()
    trace = sim.run()
    
    # stampa prime e ultime 5 righe della simulazione
    print("t, Ns, Np, Nm, Ng (prime 5):", trace[:5])
    print("t, Ns, Np, Nm, Ng (ultime 5):", trace[-5:])
