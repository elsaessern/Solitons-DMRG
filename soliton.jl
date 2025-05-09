#-----------kitaev paper chiral soliton------------
using ITensors, ITensorMPS
using Plots
using LaTeXStrings
using CSV
using DataFrames
gr()


let
    N = 100 # lenght of chains
    sites = siteinds("S=1/2",N)
    os = OpSum()
    J = 1 # J <-> coupling
    h_xy = 0 # magnetic field strength in xy direction(s)
    h_z = 0 # im pretty sure this is same as h_xy, but since theta_z = 0 the magnetic field term for z direction is zero anyways
    
    phi_xy = pi/4
    theta_z = 0
        
# kitaev spin chain hamiltonian (alternating interactions: odd i have Sx at i and i+1, even i have Sy at i and i+1)
for i in 1:N-1
    if isodd(i)
        os += J, "Sx", i, "Sx", i+1
    else
        os += J, "Sy", i, "Sy", i+1
    end
end


# kitaev magnetic field term 
for i in 1:N
            os -= cos(phi_xy)*cos(theta_z)*h_xy, "Sx", i
            os -= sin(phi_xy)*cos(theta_z)*h_xy, "Sy", i
            os -= sin(theta_z)*h_z, "Sz", i            
end
    

H = MPO(os,sites)

    psi0 = random_mps(sites;linkdims=10)
        
    nsweeps = 150 # increase number of sweeps until it converges. For lower h term, LOTS more sweeps are needed (500+)
    maxdim = [10, 20, 100, 400, 800, 1200, 1200] # maxdim for convergence. Kitaev paper uses "more than 1000" for the bond dim.
    cutoff = [1E-12] # adjust this to tune accuracy of convergence. 1E-12 is "almost exact." kitaev paper uses less than 1E-10.
    noise = [1E-6]
    print("------Ground State------\n")
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    print("\nGround State Energy =", energy,"\n") 


   
function spin_expectation(psi, op::AbstractString, site::Int)
    return real(expect(psi, op; sites=site))
end
Sx_vals = zeros(N)
Sy_vals = zeros(N)
    for i in 1:N
        Sx_vals[i] = spin_expectation(psi, "Sx", i)
        Sy_vals[i] = spin_expectation(psi, "Sy", i)
    end

odd_sites = 1:2:N
    
    p1 = plot(
        odd_sites, Sx_vals[odd_sites],
        label = L"\langle S^x_i \rangle", 
        linewidth = 2, 
        color = :green,
        xlabel = "Site i",
        ylabel = "Spin Expectation Value"
    )
    plot!(
        odd_sites, Sy_vals[odd_sites],
        label = L"\langle S^y_i \rangle",
        linewidth = 2,
        color = :red,
        linestyle = :dash
    )



display(p1)

    
end
