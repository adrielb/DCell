

dt = 1e-6 
dc = 2e-7  
iterations = 600


sf = 1e8
viz = 3
timestep = dt


R_int = 0.003

side_length = 0.1357 * 1.5 
half_box = side_length*0.5
volume = side_length*side_length*side_length 
liters = volume*1e-15
Na = 6.022e23

main_scale = 8.0

A = 0.01 
B = 1.0 
in_X = 0.003 / main_scale
in_Y = 0.02 / main_scale
in_Z = 0.00025 / main_scale

k1 = 0.025 * sf 
k2 = 1 * sf
k3 = 1 * sf
k4 = 1 * sf
k5 = 0.0004 * sf

LNa = liters*Na

x_count = in_X * LNa 
y_count = in_Y * LNa
z_count = in_Z * LNa




