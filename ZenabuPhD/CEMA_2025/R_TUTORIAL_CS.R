library(deSolve)

delta <- 0.01
timeHorizon <- 100
timesteps <- seq(0,timeHorizon, delta)

S=I=rep(0,length(timesteps))
N=500
I0=1
S[1]=N-I0
I[1]=I0
beta=0.3
gamma=0.1

#populate S and I using the loop
for (t in c(1:(length(timesteps)-1))) {
    S[t+1]= S[t] + delta*(-beta*I[t]/N*S[t] + gamma*I[t])
    I[t+1]= I[t] + delta*(beta*I[t]/N*S[t] - gamma*I[t])
    
}

dat <- data.frame(timesteps=timesteps ,S=S, I=I)

plot(x=timesteps, y=S, col="blue", type = "l", ylim = c(0,N), ylab = "Number of hosts")

lines(timesteps, I, col="red")
legend(x="topright",
       legend = c("Susceptible", "Infected"),
       col = c("blue", "red"), lwd = 1)


adansi_akrofuom,adansi_asokwa,adansi_south,adansi_north,amansie_central,banda,tain,
dormaa_west,sunyani_west,atebubu_amantin_municipal,nkoranza_north,techiman_north,
pru_west,kintampo_south,nkoranza_south_municipal,keta_municipal,akatsi_north,sefwi_akontombra,
suaman,ellembelle,tarkwa_nsuaem_municipal,