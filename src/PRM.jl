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
    end

    mutable struct SimState
        x::Float32 # Current state of PT1 Block
        u::Float32 # Current input
        y::Vector{Float32} # Current output
        d::Float32 # Current local discretization error
        h::Float32 # Current step size
        t::Float32 # Current time
        i::Int32 # Current loop iteration
        hys_state::Vector{Float32} # Hysteresis memory
    end
    SimState() = SimState(0.0, 0.0, [], 0.0, 0.0, 0.0, 0, [])


    
    mutable struct SimOutput
        x_values::Vector{Float32} # State of PT1-Block
        u_values::Vector{Float32} # Input
        y_values::Vector{Float32} # Outputs of all blocks
        d_values::Vector{Float32} # Local Discretization Error
        t_values::Vector{Float32} # Time
        h_values::Vector{Float32} # Step size
    end
    SimOutput() = SimOutput([], [], [], [], [], [])


    # PT1 Block
    function pt1(x, u)
        Tm = 10 # Time constant
        xdot = -1/Tm * x + 1/Tm * u # State equation
        return (; xdot, x)
    end

    # Hysteresis block
    function hysteresis(u, hys_state, i)
        h_e = 0.085 # Rising edge limit
        h_a = 0.065 # Falling edge limit

        if u >= h_e
            y = 1;
        elseif u <= -h_e
            y = -1;
        else
            y = 0;
        end
    
        if i > 1
            if hys_state[i - 1] == 1 && u <= h_a
                y = 0;
            elseif hys_state[i - 1] == -1 && u >= -h_a
                y = 0;
            end
        end

        append!(hys_state, y)
        return y
    end



    function sim_prm(x, u, t, hys_state, i)
        # Block 3, PT1, only output
        _, y3 = pt1(x, 0)
        # Block 1, Subtractor
        y1 = -(u, y3)
        # Block 2, Hysteresis
        y2 = hysteresis(y1, hys_state, i)
        # Block 3, PT1, only state
        xdot, _ = pt1(x, y2)
        return xdot, [y1, y2, y3]
    end



      # Integration method
      function midpoint_method!(sim_topology::Function, state::SimState)
        x = state.x
        u = state.u
        t = state.t
        h = state.h

        k1, y = sim_topology(x, u, t, state.hys_state, state.i)
        k2, _ = sim_topology(x + h/2*k1, u, t + h/2, state.hys_state, state.i)
        k3, _ = sim_topology(x - h*k1 + 2*h*k2, u, t + h, state.hys_state, state.i)
        x = x + h*k2
        d = h/6*(k1 - 2*k2 + k3)

        state.x = x
        state.y = y
        state.d = d
    end



    function update_step_size!(state::SimState, simp::SimParam)
        h = state.h
        h_min = 12*simp.epsilon
        h_max = simp.Tm*2
        h = 0.9*h*min(max((simp.epsilon/abs(state.d))^(1/3), 0.3), 2)
        h = min(max(h, h_min), h_max)
        state.h = h
    end



    function log_output!(state::SimState, output::SimOutput)
        append!(output.u_values, state.u)
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

        plot(p1, p2, p3, p4, p5, p6, layout=(3, 2))
    end


    function run_sim(simp::SimParam, state::SimState, output::SimOutput)
        state.x = simp.x0
        state.h = simp.h0
        state.i = 1

        state.t = simp.t0

        while state.t < simp.tf

            if state.t < simp.ts
                state.u = 0
            else
                state.u = 0.49
            end

            midpoint_method!(sim_prm, state)
            log_output!(state, output)
            update_step_size!(state, simp)

            state.t = state.t + state.h
            state.i = state.i + 1
        end

        render_output(output)
        
    end



    function main()
        
        

    end




end # module
