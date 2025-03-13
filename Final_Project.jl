using Measurements
using Plots
using Statistics
using OptimalTransport, Distributions
## add @. in front of the function
## --- field values in column 1, magnetization values in column 2
## --- this is sample data from Maxunmix_example_data

goethite_dist = rand(Normal(3.30, 0.5),1000)
histogram(mostgoethite)
1   
hematite_dist = rand(Normal(2.3, 0.51961524227), 1000)
histogram(mosthematite)


using DelimitedFiles
    delimiter = ','
    data = readdlm("Maxunmix_example_data.csv", delimiter)

    using StatGeochem
    data = importdataset("Maxunmix_example_data.csv", delimiter, importas=:Tuple);

x = 98  ## -- this is the number of sample points that makes up each component (B)    
mu = data.Bh
sigma = data.DP
q = data.S
p = data.P
SGG(x,mu,sigma,q,p) = (1/(2^((1+1)/p)*sigma*roh*((1+1)/p)))*((abs((q*e^(qx))+q^(-1)*e^(e^(x/q)))))*e^(-0.5*((log(((e^(qx))+e^(x/q))/2)^p)))


x = 0:0.1:4
pdf1 = SSG.(Normal(1,1), x)
pdf2 = SSG.(Normal(2,1), x)

dist1 = DiscreteNonParametric(x, pdf1./sum(pdf1))
dist2 = DiscreteNonParametric(x, pdf2./sum(pdf2))

emd = wasserstein(dist1, dist2)

println("Earth Mover's Distance", emd)

## --- this is actual data from sample 16, 1A-309Y-1

# using DelimitedFiles
# delimiter = ','
# data = readdlm("16, 1A-309Y-1.csv", delimiter)

using StatGeochem
using SpecialFunctions
data = importdataset("16, 1A-309Y-1.csv", delimiter, importas=:Tuple);

x = 98  ## -- this is the number of sample points that makes up each component (B)    
mu = data.Bh
sigma = data.DP
q = data.S
p = data.P
function skewgg(x,mu,sigma,q,p)
    xs = (x - mu)/sigma
    result = 1/(2^(1+1/p) * sigma * gamma(1+1/p)) 
    result *= abs(q*exp(xs*q) + exp(xs/q)/q) / (exp(q*xs) + exp(xs/q))
    result *= exp(-0.5*abs(log(0.5*(exp(xs*q) + exp(xs/q))))^p)
    return result
end

x = 0:0.1:4
pdf1 = skewgg.(x, mu[1], sigma[1], q[1], p[1])
pdf2 = skewgg.(x, mu[2], sigma[2], q[2], p[2])


h = plot(xlabel="", ylabel="")
plot!(x, pdf1./sum(pdf1), label="1")
plot!(x, pdf2./sum(pdf2), label="2")
display(h)

# Turn pdfs into Distributions.jl distributions
dist1 = DiscreteNonParametric(x, pdf1./sum(pdf1))
dist2 = DiscreteNonParametric(x, pdf2./sum(pdf2))

goethite_dist = Normal(3.3, 0.25)
hematite_dist = Normal(2.3, 0.26)
discrete_goethite_dist = DiscreteNonParametric(x, pdf.(goethite_dist, x)./sum(pdf.(goethite_dist, x)))
discrete_hematite_dist = DiscreteNonParametric(x, pdf.(hematite_dist, x)./sum(pdf.(hematite_dist, x)))

# # If you wanted to draw from this distribution
# goethite_samples = rand(Normal(3.30, 0.5),10000)
# hematite_samples = rand(Normal(2.3, 0.51961524227), 10000)
# histogram(goethite_samples)
# histogram(hematite_samples)

emd_goethite = wasserstein(dist1, discrete_goethite_dist)
emd_hematite = wasserstein(dist1, discrete_hematite_dist)
println("Earth Mover's Distance: Component 1 Similarity to Goethite ", emd_goethite)
println("Earth Mover's Distance: Component 1 Similarity to Hematite ", emd_hematite)
emd_goethite = wasserstein(dist2, discrete_goethite_dist)
emd_hematite = wasserstein(dist2, discrete_hematite_dist)
println("Earth Mover's Distance: Component 2 Similarity to Goethite ", emd_goethite)
println("Earth Mover's Distance: Component 2 Similarity to Hematite ", emd_hematite)


#plot(x, pdf.(goethite_dist, x))
plot(x, pdf.(discrete_goethite_dist, x), xlabel="Applied field (log(T))",
ylabel = "Magnetization (Am^2/kg)", title = "16, 1A-309Y-1 coercivity spectra", label= "Idealized goethite coercivity", color=:yellow)
plot!(x, pdf.(discrete_hematite_dist, x), label= "Idealized hematite coercivity", color=:red)
plot!(x, pdf.(dist1, x), label= "Component 1 coercivity", color=:green)
plot!(x, pdf.(dist2, x), label= "Component 2 coercivity", color=:blue)

