# testing hamiltonian construction 
using ITensors, ITensorMPS
using Plots
using LaTeXStrings
using CSV
using DataFrames
gr()

function second_derivative(y::Vector{Float64}, x::Vector{Float64})
    n = length(y)
    if n != length(x) || n < 3
        error("input vectors must be same length & have 3 entries")
    end
    d2y = zeros(n)
    dx = x[2] - x[1] 

    # central difference for interior only
    for i in 2:n-1
        d2y[i] = (y[i+1] - 2*y[i] + y[i-1]) / (dx^2)
    end


    # could set endpoints equal to nearest interior point, but this is inaccurate. going to omit the endpoints for now
    # d2y[1] = d2y[2]
    # d2y[n] = d2y[n-1]

    # omits endpoints
    return d2y[2:n-1]
end



let
    N = 99 # lenght of chain
    sites = siteinds("S=1/2",N)
    J = 1 # J <-> coupling
    K = 1
    theta_z = 0

    h_values = 0.5:0.025:0.75
    theta_values_denomenator = 3.5:0.05:4.5


    #-----DATA STORAGE FOR LATER-----
    # tuple of (h,theta,soliton_mass)
    calculated_data = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    susceptibility_results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    
    # the big loop (menacing)
for h in h_values
        print("calculating for h = ", h)
        print("\n")
        local_d = Float64[] # theta values for this h 
        local_soliton_masses = Float64[] # soliton masses for this h
    for d in theta_values_denomenator
        print("  theta = ", pi/d)
        print("\n")
        
        halfway::Int = ceil(N/2)
        
        os = OpSum()
        bondenergy = OpSum()
            
        a = cos(pi/d)*cos(theta_z)
        b = sin(pi/d)*cos(theta_z)
        c = sin(theta_z)


    
        # kitaev spin chain hamiltonian (alternating interactions: odd i have Sx at i and i+1, even i have Sy at i and i+1)
        for i in 1:N-1
            if isodd(i)
                os += J, "Sx", i, "Sx", i+1
            else
                os += K, "Sy", i, "Sy", i+1
            end
        end
        
        # kitaev magnetic field term, first half
        for i in 1:N
                os -= a*h, "Sx", i
                os -= b*h, "Sy", i
                os -= c*h, "Sz", i            
        end
        
          
        
        H = MPO(os,sites)
        
        
            states  = [i <= halfway ? "Up" : "Dn" for i in 1:N]
        
            # initial MPS with randomness in bonds
            psi0    = random_mps(sites, states; linkdims = 30)
            orthogonalize!(psi0, halfway)
            
        
            
            # DMRG parameters
            nsweeps = 1000 # increase number of sweeps until it converges. For lower h term, LOTS more sweeps are needed (500+)
            maxdim = [600, 600, 800, 800, 1200, 1200, 1200,2000,2000,3000] # maxdim for convergence. Kitaev paper uses "more than 1000" for the bond dim.
            cutoff = [1E-10] # adjust this to tune accuracy of convergence. 1E-12 is "almost exact." kitaev paper uses less than 1E-10.
            noise = [1E-6]
            print("\n------Ground State------\n")
            energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff, noise, outputlevel=0)
            print("\nGround State Energy =", energy,"\n") 
        
        e0 = energy/N # average energy per site
        
        # computing bond energy e_i
        ei = zeros(Float64, N-1)
        
                    
                    
        # new opsum for local hamiltonian terms
        for i in 1:N-1
            if isodd(i)
                bondenergy += J, "Sx", i, "Sx", i+1
            else
                bondenergy += K, "Sy", i, "Sy", i+1
            end
               bondenergy -= a*h, "Sx", i
               bondenergy -= b*h, "Sy", i
        end
        
        
        # initalize endpoints
        #ei[1] = real(inner(psi' , MPO(0.5*bondenergy[1] + bondenergy[2] + bondenergy[3], sites) , psi))
        #ei[N] = real(inner(psi' , MPO(0.5*bondenergy[length(bondenergy)-2] + bondenergy[length(bondenergy)-1] + bondenergy[length(bondenergy)], sites) , psi))
        
        # magnetic field for endpoint of chain    
       # bondenergy -= a*h, "Sx", N
       # bondenergy -= b*h, "Sy", N
        
        # construct e_i's (with expectation value)
        for i in 4:3:length(bondenergy)-3
            k = Int(ceil(i/3)) 
            ei[k] = real(inner(psi' , MPO(0.5*bondenergy[i] + bondenergy[i+1] + bondenergy[i+2] + 0.5*bondenergy[i-3], sites) , psi))
        end
            
        
        
        
        # soliton mass
        soli_mass = sum(ei) 
        
        # --- store results for this theta ---
        push!(local_d, d)
        push!(local_soliton_masses, soli_mass)
        
        # clear opsums for next iteration
        os = OpSum()
        bondenergy = OpSum()
        
                end # end of theta loop
                
            calculated_data[h] = (local_d, local_soliton_masses) # store data in array for later
                
            end # end of h loop
    
# -------- END OF LOOP --------
    
print("DMRG FINISHED. CALCULATING SUSCEPTIBILITY...")
    
h_plot_values = sort(collect(keys(calculated_data))) # order h (idk if i need to do this but just in case)

for h in h_plot_values
    thetas_, masses = calculated_data[h]
    thetas = pi ./ thetas_
    if length(masses) < 3
        println("not enough data points for h = $h to calculate second derivative. skipping...")
        continue
    end
    # 2nd derivative
    d2_masses_dtheta2 = second_derivative(masses, thetas)

    # calculate susceptibility chi = - d^2(mass)/dtheta^2
    susceptibility = -d2_masses_dtheta2

    # store results
    thetas_for_derivative = thetas[2:end-1]
    susceptibility_results[h] = (thetas_for_derivative, susceptibility)
    print("calculated susceptibility for h = %.3f\n", h)
end
 h_plot_values_sorted = sort(collect(keys(calculated_data)))
    
    for h_val in h_plot_values_sorted
        denominators, masses = calculated_data[h_val]
        thetas_phi_xy = pi ./ denominators # these are the actual phi_xy values in radians
        
        if length(masses) < 3
            println("not enough data points for h = $h_val to calculate second derivative for chi_e. skipping...")
            continue
        end

        d2_masses_dtheta2 = second_derivative(masses, thetas_phi_xy)
        susceptibility_chi_e = -d2_masses_dtheta2
        
        thetas_for_derivative = thetas_phi_xy[2:end-1]
        susceptibility_results[h_val] = (thetas_for_derivative, susceptibility_chi_e)
        print("calculated energetic susceptibility chi_e for h = ", round(h_val, digits=3), "\n")
    end
        
    print("\nPLOTTING ENERGETIC SUSCEPTIBILITY...\n")
    
    # define colors for the plot
    colors_chi_e = cgrad(:viridis, length(h_plot_values_sorted), categorical=true)

    final_plot_chi_e = plot(
        xlabel = L"\phi_{xy} \text{ (degrees)}", 
        ylabel = L"\chi^e_{\phi_{xy}}",
        title = "Energetic Susceptibility vs. Field Angle",
        legend = :topright 
    )
    
    plot_indx = 1
    for h_val in h_plot_values_sorted
        if haskey(susceptibility_results, h_val)
            thetas_rad, susceptibility = susceptibility_results[h_val]
            thetas_deg = thetas_rad .* (180 / pi) # convert radians to degrees

            plot!(
                final_plot_chi_e,
                thetas_deg,
                susceptibility,
                label = "h = $(round(h_val, digits=3))", 
                linewidth = 2.0, 
                color = colors_chi_e[plot_idx]
            )
            plot_indx += 1
        end
    end
    display(final_plot_chi_e) 
    savefig(final_plot_chi_e, "chi_e_vs_phi_xy.png")


end
