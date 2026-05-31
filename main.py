from mesh import generate_mesh
from initialization import initialize_flow
from solver import maccormack_step
from postprocess import plot_results

# constants
GAMMA = 1.4
N = 61
LENGTH = 3.0
CFL = 0.5
CX = 0.2
ITERATIONS = 2000
EXIT_PRESSURE = 0.674


# mesh
x, dx, area = generate_mesh(N,LENGTH)

# initialize
U1, U2, U3 = initialize_flow(x,area,GAMMA)

# solve
for i in range(ITERATIONS):
    U1, U2, U3 = maccormack_step(U1,U2,U3,area,dx,GAMMA,CFL,CX,EXIT_PRESSURE)

# plots
plot_results(x,area,U1,U2,U3,GAMMA)