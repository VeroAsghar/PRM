using Plots

function midpoint_method(f::Function, x, u, t, h)
    k1, y = f(x, u, t)
    k2, _ = f(x + h/2*k1, u, t + h/2)
    k3, _ = f(x - h*k1 + 2*h*k2, u, t + h)
    x_VPG = x + h*k2
    d_VPG = h/6*(k1 - 2*k2 + k3)
    return x_VPG, y, d_VPG
end

function pt1(x, u)
    Tm = 10
    xdot = -1/Tm * x + 1/Tm * u
    y = x
    return xdot, y
end

function subtractor(u1, u2)
    return u1 - u2
end

function hysteresis(u)
    h_e = 0.085
    h_a = 0.065
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

function sys_PRM(x, u, t)
    _, y1::Float64 = pt1(x, 0)
    y2::Float64 = subtractor(u, y1)
    y3::Float64 = hysteresis(y2)
    xdot1, _ = pt1(x, y3)
    return xdot1, [y1 y2 y3]
end

t_step = 1;
Tm = 10
t0 = 0
tf = 20
epsilon = 1e-3
h_max = 2*Tm
h_min = 12*epsilon
u_values = []
x_values = []
y_values = Vector{Float64}()
d_values = []
t_values = []
h_values = []
hys_state = []
x = 0
h = 0.2
t = 0
i = 1
while t <= tf
    if t < t_step
        global u = 0
    else
        global u = 0.17
    end

    global x, y, d = midpoint_method(sys_PRM, x, u, t, h);
    global h = 0.9*h*min(max((epsilon/max(abs(d)))^(1/3), 0.3), 2)
    global h = min(max(h, h_min), h_max)

    append!(x_values, x)
    append!(y_values, y)
    append!(d_values, d)
    append!(h_values, h)
    append!(t_values, t)

    global t = t + h
    global i = i + 1
end

y_values = reshape(y_values, (3, :))

p1 = plot(t_values, y_values[1, :])
p2 = plot(t_values, y_values[2, :])
p3 = plot(t_values, y_values[3, :])
p4 = plot(t_values, h_values)
p5 = plot(t_values, d_values)

plot(p1, p2, p3, p4, p5, layout=(5, 1))
