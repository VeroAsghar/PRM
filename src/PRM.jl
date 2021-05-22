module PRM

    using Plots

    export run_sim, SimParam, SimState, SimOutput

    Base.@kwdef struct SimParam
        Tm::Float32 # Time constant of the PT1 Block [s]
        epsilon::Float32 # Local tolerance of update step size function
        t0::Float32 # Simulation start time [s]
        tf::Float32 # Simulation end time [s]
        ts::Float32 # Time at which input goes high
        x0::Float32 # Initial state of PT1-Block 
        h0::Float32 # Intial step size [s]
        u0::Float32 # Rising edge size
    end

    mutable struct SimState
        x::Float32 # Current state of PT1 Block
        u1::Float32 # Current input
        u2::Float32
        u3::Float32
        y::Array{Float32} # Current output
        d::Float32 # Current local discretization error
        h::Float32 # Current step size
        t::Float32 # Current time
        i::Int32 # Current loop iteration
        hys_last_run::Int32 # Hysteresis memory
        hys_last_loop::Int32
    end
    SimState() = SimState(0.0, 0.0, 0.0, 0.0, [0.0, 0.0, 0.0], 0.0, 0.0, 0.0, 0, 0, 0)


    
    struct SimOutput{T}
        x_values::Vector{T} # State of PT1-Block
        u_values::Vector{T} # Input
        y_values::Vector{T} # Outputs of all blocks
        d_values::Vector{T} # Local Discretization Error
        t_values::Vector{T} # Time
        h_values::Vector{T} # Step size
    end
    SimOutput(T) = SimOutput{T}(T[], T[], T[], T[], T[], T[])
    SimOutput() = SimOutput(Float32)


    # PT1 Block
    function pt1(x, u)
        Tm = 10 # Time constant
        xdot = -1/Tm * x + 1/Tm * u # State equation
        return (; xdot, x)
    end

    # Hysteresis block
    function hysteresis(u, hys_last)
        h_e = 0.085 # Rising edge limit
        h_a = 0.065 # Falling edge limit

        if u >= h_e
            y = 1
        elseif u <= -h_e
            y = -1
        elseif -h_a <= u <= h_a
            y = 0
        end
    
        if h_a < u < h_e
            if hys_last == 1
                y = 1
            else
                y = 0
            end
        elseif -h_e < u < -h_a
            if hys_last == -1
                y = 1
            else
                y = 0
            end
        end
        return y
    end

    function step(u, t, t_step)
        if t < t_step
            return 0
        else
            return u
        end
    end

    function sim_prm!(x, u, t, state::SimState, output::Bool)
        # Block 3, PT1, only output
        _, y3 = pt1(x, 0)
        # Block 1, Subtractor
        y1 = -(u, y3)
        # Block 2, Hysteresis
        y2 = hysteresis(y1, state.hys_last_run)
        state.hys_last_run = y2
        # Block 3, PT1, only state
        xdot, _ = pt1(x, y2)
        
        if output
            state.y[1] = y1
            state.y[2] = y2
            state.y[3] = y3
        end

        return xdot
    end



      # Integration method
      function midpoint_method!(sim_topology::Function, state::SimState)
        x = state.x
        t = state.t
        h = state.h

        k1 = sim_topology(x, state.u1, t, state, true)
        k2 = sim_topology(x + h/2*k1, state.u2, t + h/2, state, false)
        k3 = sim_topology(x - h*k1 + 2*h*k2, state.u3, t + h, state, false)
        x = x + h*k2
        d = h/6*(k1 - 2*k2 + k3)

        state.x = x
        state.d = d
    end



    function log_output!(state::SimState, output::SimOutput)
        append!(output.u_values, state.u1)
        append!(output.x_values, state.x)
        append!(output.y_values, state.y)
        append!(output.t_values, state.t)
        append!(output.d_values, state.d)
        append!(output.h_values, state.h)
    end



    function render_output(output::SimOutput)
        y_values = reshape(output.y_values, (3, :))
        
        p1 = plot(output.t_values, output.u_values, label="Input")
        p2 = plot(output.t_values, y_values[1, :], label="Output Subtractor")
        p3 = plot(output.t_values, y_values[2, :], label="Output Hysteresis")
        p4 = plot(output.t_values, y_values[3, :], label="Output PT1")
        p5 = plot(output.t_values, output.h_values, label="Step Size")
        p6 = plot(output.t_values, output.d_values, label="LDE")

        plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), xlabel="Time [s]")
    end


    function run_sim(simp::SimParam)
        output = SimOutput();
        state = SimState();

        state.x = simp.x0
        state.h = simp.h0
        state.i = 1

        state.t = simp.t0

        h_min = 12*simp.epsilon
        h_max = simp.Tm*2
        
        while state.t < simp.tf

            state.u1 = step(simp.u0, state.t, simp.ts)
            state.u2 = step(simp.u0, state.t + state.h/2, simp.ts)
            state.u3 = step(simp.u0, state.t + state.h, simp.ts)

            state.hys_last_run = state.hys_last_loop
            midpoint_method!(sim_prm!, state)

            

            if abs(state.d) > 0
                h_new = state.h*(simp.epsilon/abs(state.d))^(1/3)
                
                h_new = min(max(h_new, h_min), 0.99*h_max)

                if h_new > 2*state.h
                    state.h = h_new
                elseif h_new <= state.h
                    state.h = 0.75*h_new
                    continue
                end
            end

            state.hys_last_loop = state.hys_last_run
            log_output!(state, output)

            state.t = state.t + state.h
            state.i = state.i + 1
        end

        render_output(output)
        
    end



end # module
simp = PRM.SimParam(Tm=10, epsilon=1e-10, x0=0, h0=0.01, t0=0, tf=20, ts=1, u0=0.17);

PRM.run_sim(simp)