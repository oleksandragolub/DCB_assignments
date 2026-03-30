# =========================
# Golub Oleksandra 856706
#
# Codice per generare automaticamente le equazioni differenziali
# di un sistema di reazioni chimiche secondo la legge d’azione di massa
#
# - Ogni reazione è definita da reagenti, prodotti e costanti di velocità (diretta e inversa)
# - Il programma costruisce le ODE simboliche per ogni specie e la matrice stechiometrica del sistema
# =========================

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Union


# =========================
# Strutture dati di base
# =========================
@dataclass
class Reaction:
    reactants: Dict[str, int] # reagenti 
    products: Dict[str, int]  # prodotti
    k_forward: float          # costante di velocità diretta
    reversible: bool = False  # se vale true, la reazione è reversibile
    k_reverse: float = 0.0    # costante di velocità inversa
    kf_symbol: str = ""       # etichetta per kf nella stampa delle ODE
    kr_symbol: str = ""       # etichetta per kr, se la reazione è reversibile


class ReactionSystem:
    """
    Costruisce le specie, converte reazioni (legge d'azione di massa) in ODE testuali
    e calcola la matrice stechiometrica
    """
    def __init__(self):
        self.reactions: List[Reaction] = []      # lista delle reazioni aggiunte
        self.species_set = set()                 # insieme delle specie incontrate
        self.species_list: List[str] = []        # lista ordinata di specie, riempita da _finalize_species di sotto
        self.ode_system: Dict[str, str] = {}     # ODE convertite in formato testo

    def _finalize_species(self):   # ordina l’elenco specie
        self.species_list = sorted(self.species_set)

    def add_reaction(       # aggiunge una reazione al sistema
        self,
        reactants: Union[Dict[str, int], List[Tuple[str, int]]],     # reagenti con coefficienti
        rate: float,                                                 # valore numerico di kf
        products: Union[Dict[str, int], List[Tuple[str, int]]],      # prodotti con coefficienti
        reversible: bool = False,                                    # vale true se la reazione è reversibile
        reverse_rate: float = 0.0,                                   # valore numerico di kr, usata solo se la reazione è reversibile
        rxn_index: int = None,                                       # indice esplicito della reazione per generare etichette kf e kr (se omesso, l’indice viene assegnato in automatico)
    ):
        # Normalizza input (reactants/products) a dict
        if isinstance(reactants, list):
            reactants = dict(reactants)
        if isinstance(products, list):
            products = dict(products)

        # Indici/etichette parametri (kf che sia costante di velocità diretta, kr che sia costante di velocità inversa)
        if rxn_index is None:
            rxn_index = len(self.reactions) + 1
        kf_symbol = f"k{rxn_index}"
        kr_symbol = f"k_minus{rxn_index}" if reversible else ""

        # Registra reazione
        r = Reaction(
            reactants=reactants,
            products=products,
            k_forward=rate,
            reversible=reversible,
            k_reverse=reverse_rate if reversible else 0.0,
            kf_symbol=kf_symbol,
            kr_symbol=kr_symbol,
        )
        self.reactions.append(r)

        # Aggiorna specie
        self.species_set.update(reactants.keys())
        self.species_set.update(products.keys())

    def generate_odes(self) -> Dict[str, str]:  # costruisce le equazioni differenziali (testo) per ogni specie secondo la legge di azione di massa
        self._finalize_species()
        deriv_terms: Dict[str, List[str]] = {s: [] for s in self.species_list}

        for r in self.reactions:
            # termine diretto v_f = kf * Π [X]^stoich
            v_f = r.kf_symbol
            for sp, nu in r.reactants.items():
                v_f += f" * [{sp}]**{nu}" if nu != 1 else f" * [{sp}]"

            # consumo reagenti, produzione prodotti (diretta)
            for sp, nu in r.reactants.items():
                deriv_terms[sp].append(f"-{nu}*{v_f}" if nu != 1 else f"-{v_f}")
            for sp, nu in r.products.items():
                deriv_terms[sp].append(f"+{nu}*{v_f}" if nu != 1 else f"+{v_f}")

            # termine inverso (se reversibile) v_r = kr * Π [Prod]^stoich
            if r.reversible:
                v_r = r.kr_symbol
                for sp, nu in r.products.items():
                    v_r += f" * [{sp}]**{nu}" if nu != 1 else f" * [{sp}]"

                for sp, nu in r.reactants.items():
                    deriv_terms[sp].append(f"+{nu}*{v_r}" if nu != 1 else f"+{v_r}")
                for sp, nu in r.products.items():
                    deriv_terms[sp].append(f"-{nu}*{v_r}" if nu != 1 else f"-{v_r}")

        # si creano le stringhe finali
        self.ode_system = {
            s: ("d[" + s + "]/dt = " + " ".join(terms)).replace("+-", "-").replace("-+", "-")
            if terms else f"d[{s}]/dt = 0"
            for s, terms in deriv_terms.items()
        }
        return self.ode_system

    def print_odes(self):    # stampa su terminale le ODE già generate
        if not self.ode_system:
            self.generate_odes()
            
        print("\n" + "="*60)
        print("SISTEMA DI EQUAZIONI DIFFERENZIALI (Legge d'Azione di Massa)")
        print("="*60)
        
        for s in self.species_list:
            print(self.ode_system[s])
        print("="*60)

    def get_stoichiometric_matrix(self) -> np.ndarray:  # si calcola la matrice stechiometrica
        """
        Matrice stechiometrica S (numero di specie x numero di reazioni_effettive),
        contando una colonna per la diretta e una per l'inversa, se la reazione è reversibile
        """
        self._finalize_species()
        n_cols = len(self.reactions) + sum(1 for r in self.reactions if r.reversible)
        S = np.zeros((len(self.species_list), n_cols))
        c = 0
        
        for r in self.reactions:
            
            # diretta
            for i, s in enumerate(self.species_list):
                if s in r.reactants:
                    S[i, c] -= r.reactants[s]
                if s in r.products:
                    S[i, c] += r.products[s]
            c += 1
            
            # inversa
            if r.reversible:
                for i, s in enumerate(self.species_list):
                    if s in r.products:
                        S[i, c] -= r.products[s]
                    if s in r.reactants:
                        S[i, c] += r.reactants[s]
                c += 1
                
        return S


def reaction_to_odes(   # serve per costruire ODE da una lista di reazioni
    reactions_list: List[Tuple],
    reversible: bool = False
) -> Dict[str, str]:
    
    """
    Ogni reazione è (reactants, kf, products)  oppure (reactants, kf, products, kr)
    Se kr manca (si tratta di tuple a 3 elementi) e reversible=True, assume di default kr = 0.1*kf.
    """
    system = ReactionSystem()  # si crea un ReactionSystem, in cui si aggiunge ogni reazione, generando e poi ritornando le ode
    
    for idx, item in enumerate(reactions_list, start=1):
        if len(item) == 3:
            reactants, kf, products = item
            is_rev = reversible
            kr = 0.1 * kf if is_rev else 0.0
            
        elif len(item) == 4:
            reactants, kf, products, kr = item
            is_rev = (kr > 0.0)
            
        else:
            raise ValueError("Ogni reazione deve avere 3 o 4 elementi.")
            
        system.add_reaction(reactants, kf, products,
                           reversible=is_rev, reverse_rate=kr, rxn_index=idx)
                           
    return system.generate_odes()


# =========================
# Alcuni esempi
# =========================
if __name__ == "__main__":
    # Esempio 1: A + B -> C
    print("\nESEMPIO 1: A + B -> C")  # ODE con solo il termine diretto
    sys1 = ReactionSystem()
    sys1.add_reaction({"A":1,"B":1}, rate=0.5, products={"C":1}, rxn_index=1)
    sys1.generate_odes()
    sys1.print_odes()

    # Esempio 2: A + B <-> C
    print("\nESEMPIO 2: A + B <-> C")   # ODE con due termini: diretto e inverso
    sys2 = ReactionSystem()
    sys2.add_reaction({"A":1,"B":1}, rate=0.5, products={"C":1},
                      reversible=True, reverse_rate=0.1, rxn_index=1)
    sys2.generate_odes()
    sys2.print_odes()

    # Esempio 3: Michaelis–Menten  S + E <-> C -> P + E
    print("\nESEMPIO 3: S + E <-> C -> P + E")   # test di due reazioni: una reversibile (formazione del complesso C) e una irreversibile (formazione di P e rilascio E)
    sys3 = ReactionSystem()
    sys3.add_reaction({"S":1,"E":1}, rate=0.01, products={"C":1},
                      reversible=True, reverse_rate=0.001, rxn_index=1)
    sys3.add_reaction({"C":1}, rate=0.1, products={"P":1,"E":1}, rxn_index=2)
    sys3.generate_odes()
    sys3.print_odes()

    # Esempio 4: reaction_to_odes
    print("\nESEMPIO 4: reaction_to_odes")     # test di stampa dei tre reazioni semplici convertitit in odes
    reactions = [
        ({"A":2, "B":1}, 0.3, {"C":1}),   # 2A + B -> C
        ({"C":1}, 0.1, {"D":2}),         # C -> 2D
        ({"D":1, "A":1}, 0.05, {"E":1})  # D + A -> E
    ]
    odes = reaction_to_odes(reactions)
    for s in sorted(odes.keys()):
        print(odes[s])

    # Matrice stechiometrica per l'esempio 3
    print("\nMatrice stechiometrica dell'esempio 3:")  # stampa l’elenco delle specie e la matrice S (colonna per diretta e, se presente, per inversa)
    S = sys3.get_stoichiometric_matrix()
    print("Specie:", sys3.species_list)
    print(S)
