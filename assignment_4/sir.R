#######################################################
# Oleksandra Golub (856706) 
# 
# Modello SIR (Susceptible-Infected-Recovered) in R, 
# usando metodo di Eulero per risolvere ODE:
#    dS/dt = -beta * S * I / N
#    dI/dt =  beta * S * I / N - gamma * I
#    dR/dt =  gamma * I
# dove parametri sono:
#    beta: tasso di trasmissione/contatto (ovvero la probabilità di contagio per contatto tra un suscettibile e un infetto)
#    gamma: tasso di guarigione
#    R = beta/gamma: numero riproduttivo di base
#######################################################

# -----------
# Parametri
# -----------

N     <- 1000 # popolazione totale
beta  <- 0.5  # tasso di trasmissione
gamma <- 0.1  # tasso di guarigione

I0 <- 1  # numero iniziale di individui infetti all’inizio della simulazione
R0 <- 0  # numero iniziale di individui rimossi/guariti all’inizio della simulazione
S0 <- N - I0 - R0 # numero iniziale di individui suscettibili

t_max <- 160  # il tempo totale di simulazione
dt    <- 0.1  # il passo temporale del metodo di Eulero
times <- seq(0, t_max, by = dt)  # il vettore dei tempi della simulazione che va da 0 a t_max con incremento dt

R0_basic <- beta / gamma # indica il numero medio di nuovi contagi causati da un infetto in una pop suscettibile
cat("R0 =", R0_basic, "\n") # stampa a video il valore di R0 per verificare il regime epidemico (se maggiore di 1, c'e' un'epidemia)


# ------------------
# Funzione derivate
# ------------------

# funzione che definisce le equazioni differenziali del modello SIR
sir_deriv <- function(S, I, R, beta, gamma, N) {
  dS <- -beta * S * I / N  # rappresenta la variazione nel tempo dei suscettibili
  dI <-  beta * S * I / N - gamma * I  # rappresenta la variazione nel tempo degli infetti
  dR <-  gamma * I  # rappresenta la variazione nel tempo dei rimossi/guariti
  c(dS, dI, dR) # si restituiscono le tre derivate come vettore
}


# ---------------------------
# Integrazione numerica con il metodo di Eulero
# ---------------------------
S <- numeric(length(times)) # vettore che conterrà il numero di suscettibili per ogni istante di tempo
I <- numeric(length(times)) # vettore che conterrà il numero di infetti per ogni istante di tempo
R <- numeric(length(times)) # vettore che conterrà il numero di rimossi/guariti per ogni istante di tempo

# le condizioni iniziali di tre categorie al tempo t
S[1] <- S0
I[1] <- I0
R[1] <- R0

for (k in 2:length(times)) {
  deriv <- sir_deriv(S[k-1], I[k-1], R[k-1], beta, gamma, N) # si calcolano le derivate allo stato precedente

  # si aggiornano le tre categorie con il metodo di Eulero
  S[k] <- S[k-1] + deriv[1] * dt
  I[k] <- I[k-1] + deriv[2] * dt
  R[k] <- R[k-1] + deriv[3] * dt
}


# ----------
# Risultati
# ----------
peak_idx <- which.max(I) # si individua l’indice temporale in cui il numero di infetti è massimo
cat("Picco infetti =", round(I[peak_idx]), "a t =", times[peak_idx], "giorni\n")
cat("Recovered finali =", round(R[length(R)]), "\n")
cat("Suscettibili finali =", round(S[length(S)]), "\n")


# ----------------------
# Grafico + salvataggio
# ----------------------
dir.create("outputs", showWarnings = FALSE) # si crea la cartella outputs se non esiste già

png("outputs/sir_r_minimal.png", width = 1000, height = 500)

plot(times, S, type = "l", col = "blue", lwd = 2,
     xlab = "Tempo (giorni)", ylab = "Numero di individui",
     main = "Modello SIR (metodo di Eulero)")

lines(times, I, col = "red", lwd = 2)
lines(times, R, col = "green", lwd = 2)

legend("right",
       legend = c("Suscettibili (S)", "Infetti (I)", "Guariti (R)"),
       col = c("blue", "red", "green"),
       lwd = 2,
       bty = "n")

dev.off()

cat("Grafico salvato in: outputs/sir_r_minimal.png\n")
