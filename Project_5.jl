using FFTW, Plots, StatGeochem
data = importdataset("BarHarborWL2013.csv", ',', importas=:Tuple)
# Signal: water Level
wl = data.Water_Level
# Number of sample points
N = length(data.Hour_Count)
# Sample period
Ts = 24/N
# Start and end time
t0 = 0
tmax = t0 + N * Ts
# Time coordinate
t = t0:Ts:tmax
# Fast Fourier Transform
ft = fft(wl)
F = fftshift(ft)
freqs = fftshift(fftfreq(N, 1.0/Ts))
# Plot results
time_domain = plot(t[1:end-1], wl, title="Signal", framestyle=:box)
display(time_domain)
freq_domain = plot(freqs, abs.(F), title="Spectrum", yscale=:log10, framestyle=:box)
display(freq_domain)


