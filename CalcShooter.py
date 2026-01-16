from matplotlib.ticker import MultipleLocator
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed)


# Inputs
P_x = 3
P_y = 4
O_r = 0
w = 2333
O_b = 73/180*m.pi
Z_sh = 0.25
C_D = 0.5
C_L = 0
C_roll = 0.5
R_w = 2*2.54/100


# Constants
R_b = 0.15/2
g = 9.81
m_ball = 0.225
A = m.pi*R_b**2
p = 1.2
enable2d = False

# Simulation parameters
duration = 2
n_steps = 500
dt = duration/n_steps

# Field Dimensions
L_x = 158.611*2.54/100
P_x = L_x - 1.15 + 0.1
L_y = 318.188*2.54/100
W_H = 60*2.54/100
H_H = 60*2.54/100
H_H2 = 13.5*2.54/100
H_W = 20*2.54/100
R_1 = 25/2*2.54/100
R_2 = 45/2*2.54/100
H_Shoot = 0.1

# Sliders
init_angle = O_b
init_w = w
init_C_roll = C_roll
init_C_D = C_D
init_C_L = C_L 
init_P_x = P_x
init_P_y = P_y
init_O_r = O_r

# 2D Plot

ax2d = None

def setup_2d_plot():
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.3)

    ax_angle = plt.axes([0.2, 0.24, 0.65, 0.03])
    ax_wheel = plt.axes([0.2, 0.18, 0.65, 0.03])
    ax_C_roll = plt.axes([0.2, 0.12, 0.65, 0.03])
    ax_C_D = plt.axes([0.2, 0.06, 0.65, 0.03])
    ax_C_L = plt.axes([0.2, 0.0, 0.65, 0.03])

    angle_slider = Slider(
        ax_angle, 'Angle (deg)', 10, 80,
        valinit=init_angle * 180 / m.pi
    )

    wheel_slider = Slider(
        ax_wheel, 'Wheel RPM', 1000, 5000,
        valinit=init_w
    )

    C_roll_slider = Slider(
        ax_C_roll, 'C_roll', 0, 1,
        valinit=init_C_roll
    )   

    C_D_slider = Slider(
        ax_C_D, 'C_D', 0, 1,
        valinit=init_C_D
    )
    C_L_slider = Slider(
        ax_C_L, 'C_L', 0, 1,
        valinit=init_C_L
    )

    
        
    def on_change(val):
        angle = angle_slider.val * m.pi / 180
        rpm = wheel_slider.val
        C_r = C_roll_slider.val
        C_d = C_D_slider.val
        C_l = C_L_slider.val
        update_plot(angle, rpm, C_r, C_d, C_l)
        fig.canvas.draw_idle()

    angle_slider.on_changed(on_change)
    wheel_slider.on_changed(on_change)
    C_roll_slider.on_changed(on_change)
    C_D_slider.on_changed(on_change)
    C_L_slider.on_changed(on_change)
    
    return fig, ax

# Fig3D

fig3d = plt.figure()

gs = GridSpec(
    nrows=12,
    ncols=2,
    width_ratios=[4, 1],   # LEFT big plot, RIGHT sliders
    figure=fig3d
)

ax3d = fig3d.add_subplot(gs[:, 0], projection='3d')

ax3d_P_x = fig3d.add_subplot(gs[0, 1])
ax3d_P_y = fig3d.add_subplot(gs[1, 1])
ax3d_O_r = fig3d.add_subplot(gs[2, 1])
ax3d_angle = fig3d.add_subplot(gs[3, 1])
ax3d_wheel = fig3d.add_subplot(gs[4, 1])
ax3d_wradius = fig3d.add_subplot(gs[5, 1])
ax3d_zShoot = fig3d.add_subplot(gs[6, 1])
ax3d_C_roll = fig3d.add_subplot(gs[7, 1])
ax3d_C_D = fig3d.add_subplot(gs[8, 1])
ax3d_C_L = fig3d.add_subplot(gs[9, 1])
ax3d_duration = fig3d.add_subplot(gs[10, 1])
ax3d_nsteps = fig3d.add_subplot(gs[11, 1])

P_x_slider = Slider(
    ax3d_P_x, 'P_x (m)', 0, L_x,
    valinit=init_P_x
)

P_y_slider = Slider(
    ax3d_P_y, 'P_y (m)', 0, L_y,
    valinit=init_P_y
)

O_r_slider = Slider(
    ax3d_O_r, 'O (deg)', -90, 90,
    valinit=init_O_r * 180 / m.pi
)

angle_3d_slider = Slider(
    ax3d_angle, 'O_b (deg)', 10, 80,
    valinit=init_angle * 180 / m.pi
)

wheel_3d_slider = Slider(
    ax3d_wheel, 'w (RPM)', 1000, 5000,
    valinit=init_w
)

wheel_3d_radius_slider = Slider(
    ax3d_wradius, 'r_w (m)', 0.01, 0.1,
    valinit=R_w
)

zShoot_3d_slider = Slider(
    ax3d_zShoot, 'Z_sh (m)', 0, 1,
    valinit=Z_sh
)

C_roll_3d_slider = Slider(
    ax3d_C_roll, 'C_roll', 0, 1,
    valinit=init_C_roll
)

C_D_3d_slider = Slider(
    ax3d_C_D, 'C_D', 0, 1,
    valinit=init_C_D
)

C_L_3d_slider = Slider(
    ax3d_C_L, 'C_L', 0, 1,
    valinit=init_C_L
)

duration_3d_slider = Slider(
    ax3d_duration, 'Duration (s)', 0.1, 5,
    valinit=duration
)

nsteps_3d_slider = Slider(
    ax3d_nsteps, 'N Steps', 10, 1000,
    valinit=n_steps
)

def on_3d_change(val):
    global P_x, P_y, O_r, angle, wheel, R_w, Z_sh, C_roll, C_D, C_L, time, duration, n_steps, dt
    
    P_x = P_x_slider.val
    P_y = P_y_slider.val
    O_r = O_r_slider.val * m.pi / 180
    angle = angle_3d_slider.val * m.pi / 180
    wheel = wheel_3d_slider.val
    R_w = wheel_3d_radius_slider.val
    Z_sh = zShoot_3d_slider.val
    C_roll = C_roll_3d_slider.val
    C_D = C_D_3d_slider.val
    C_L = C_L_3d_slider.val
    duration = duration_3d_slider.val
    n_steps = int(nsteps_3d_slider.val)

    dt = duration / n_steps
    time = np.linspace(0, duration, n_steps)

    update_plot()
    fig3d.canvas.draw_idle()

P_x_slider.on_changed(on_3d_change)
P_y_slider.on_changed(on_3d_change)
O_r_slider.on_changed(on_3d_change)
angle_3d_slider.on_changed(on_3d_change)
wheel_3d_slider.on_changed(on_3d_change)
wheel_3d_radius_slider.on_changed(on_3d_change)
zShoot_3d_slider.on_changed(on_3d_change)
C_roll_3d_slider.on_changed(on_3d_change)
C_D_3d_slider.on_changed(on_3d_change)
C_L_3d_slider.on_changed(on_3d_change)
duration_3d_slider.on_changed(on_3d_change)
nsteps_3d_slider.on_changed(on_3d_change)

# Calculations Ideal
V_s = (w/60)*2*m.pi*R_w
V_b0 = V_s*C_roll
w_b = V_s/R_b*(1 - C_roll)
V_u0 = V_b0*m.cos(O_b)
V_n0 = V_b0*m.sin(O_b)

def calculate_iV_n(t):
    v = V_n0 - g*t
    return v 

def calculate_iN(t):
    iV_n = calculate_iV_n(t)
    n = Z_sh + iV_n*t + (V_n0-iV_n)*t/2
    if n < 0:
        n = 0
    return n

def calculate_iZ(t):
    z = calculate_iN(t)
    return z

def calculate_iU(t):
    iV_u = V_u0
    u = iV_u*t
    return u

def calculate_iY(t):
    y = P_y + calculate_iU(t)*m.sin(O_r)
    if y < 0:
        y = 0
    return y

def calculate_iX(t):
    x = P_x + calculate_iU(t)*m.cos(O_r)
    return x

def calculate_iPoints():
    iN_vals = []
    iU_vals = []
    ideal_point = []
    for t in time:
        iN_val = calculate_iN(t)
        iU_val = calculate_iU(t)
        iN_vals.append(iN_val)
        iU_vals.append(iU_val)
        if abs(iN_val-(H_H+H_H2+H_Shoot)) < 0.01:
            ideal_point = [iU_val, iN_val, t] 

    iPoints = [iU_vals, iN_vals, time]

    return iPoints, ideal_point

# Calculations Friction
def f_drag(v):
    F_D = 0.5*C_D*p*A*v**2
    return F_D/m_ball

def f_lift(v):
    F_L = 0.5*p*A*R_b*w_b*v*C_L
    return F_L/m_ball

def calculate_Re(V):
    Re = V*2*R_b/(1.5e-5)
    return Re

def calculate_frPoints():
    N_vals = []
    U_vals = []
    V_u = V_u0
    V_n = V_n0
    N = Z_sh
    U = 0
    frIdealPoint = []
    for t in time:
        V = m.sqrt(V_u**2 + V_n**2)
        a_drag = -1*f_drag(V)
        a_lift = f_lift(V)
        theta = m.atan2(V_n, V_u)


        V_u = V_u + dt*(m.cos(theta)*a_drag + m.cos(theta+ m.pi/2)*a_lift)
        V_n = V_n + dt*(m.sin(theta)*a_drag + m.sin(theta+ m.pi/2)*a_lift - g)

        U = U + V_u*dt
        N = N + V_n*dt

        if N < 0:
            N = 0

        N_vals.append(N)
        U_vals.append(U)

        if abs(N-(H_H+H_H2+H_Shoot)) < 0.01:
            frIdealPoint = [U, N, t] 
    
    frPoints = [U_vals, N_vals, time]
    return frPoints, frIdealPoint

# Coordinate Calculations
def calculate_XYZ_fromPoints(points):
    X_vals = []
    Y_vals = []
    Z_vals = []
    for i in range(len(points[0])):
        U = points[0][i]
        N = points[1][i]
        X = P_x + U*m.cos(O_r)
        Y = P_y + U*m.sin(O_r)
        Z = N
        X_vals.append(X)
        Y_vals.append(Y)
        Z_vals.append(Z)

    points_xyz = [X_vals, Y_vals, Z_vals, points[2]]
    return points_xyz
    

# Plotting 2D
time = np.linspace(0, duration, n_steps)

def plotIdeal(iPoints, ideal_point):

    ax2d.plot(iPoints[0], iPoints[1])
    ax2d.axline((P_x, H_H+H_H2), slope=0, color='g', linestyle='--')
    ax2d.axline((P_x, H_H+H_H2+H_Shoot), slope=0, color='g', linestyle='--')
    if ideal_point:
        ax2d.axvline(ideal_point[0], color='r', linestyle='--', label='Ideal Shot Point') 
        ax2d.text(
        ideal_point[0], 
        ax2d.get_ylim()[0],          # bottom of plot
        f'x = {ideal_point[0]:.2f}', 
        rotation=90,
        verticalalignment='bottom',
        horizontalalignment='right',
        color='r'
        )   
        ax2d.axvline(ideal_point[0]+R_2*2, color='r', linestyle='--', label='End Of Hub')
def plotWithFriction(frPoints, frIdealPoint):

    ax2d.plot(frPoints[0], frPoints[1], color='orange')

    if frIdealPoint:
        ax2d.axvline(frIdealPoint[0], color='m', linestyle='--', label='Fr Ideal Shot Point') 
        ax2d.text(
        frIdealPoint[0], 
        ax2d.get_ylim()[0],          # bottom of plot
        f'x = {frIdealPoint[0]:.2f}', 
        rotation=90,
        verticalalignment='bottom',
        horizontalalignment='right',
        color='m'
        )   
        ax2d.axvline(frIdealPoint[0]+R_2*2, color='m', linestyle='--', label='Fr End Of Hub')
# Plotting 3D

def plotIdeal3d(iPoints_xyz):
    ax3d.plot(iPoints_xyz[0], iPoints_xyz[1], iPoints_xyz[2], 
              color='black', alpha=1)

def plotWithFriction3d(frPoints_xyz):
    ax3d.plot(frPoints_xyz[0], frPoints_xyz[1], frPoints_xyz[2],
            color='green', alpha=1)

def draw_ground():
    X = np.array([[0, L_x], [0, L_x]])
    Y = np.array([[0, 0], [L_y, L_y]])
    Z = np.zeros_like(X)

    ax3d.plot_surface(
        X, Y, Z,
        color='gray',
        alpha=0.4,
        shade=True
    )

def draw_field_walls():
    X_1 = np.array([[0, L_x], [0, L_x]])
    Z_1 = np.array([[0, 0], [H_W, H_W]])
    Y_1 = np.full_like(X_1, 0)

    X_2 = np.array([[0, L_x], [0, L_x]])
    Z_2 = np.array([[0, 0], [H_W, H_W]])
    Y_2 = np.full_like(X_1, L_y)

    X_3 = np.full_like(X_1, 0)
    Z_3 = np.array([[0, H_W], [0, H_W]])
    Y_3 = np.array([[0, 0], [L_y, L_y]])

    ax3d.plot_surface(
        X_1, Y_1, Z_1,
        color='gray',
        alpha=0.4,
        shade=True
    )

    ax3d.plot_surface(
        X_2, Y_2, Z_2,
        color='gray',
        alpha=0.4,
        shade=True
    )

    ax3d.plot_surface(
        X_3, Y_3, Z_3,
        color='gray',
        alpha=0.4,
        shade=True
    )

def draw_hub_post():
    X_1 = np.array([[L_x, L_x + W_H], [L_x, L_x + W_H]])
    Y_1 = np.array([[(L_y-W_H)/2, (L_y-W_H)/2], [(L_y+W_H)/2, (L_y+W_H)/2]])
    Z_1 = np.full_like(X_1, 0)

    X_2 = np.array([[L_x, L_x + W_H], [L_x, L_x + W_H]])
    Y_2 = np.array([[(L_y-W_H)/2, (L_y-W_H)/2], [(L_y+W_H)/2, (L_y+W_H)/2]])
    Z_2 = np.full_like(X_1, H_H)

    X_3 = np.full_like(X_1, L_x)
    Y_3 = np.array([[(L_y-W_H)/2, (L_y-W_H)/2], [(L_y+W_H)/2, (L_y+W_H)/2]])
    Z_3 = np.array([[0, H_H], [0, H_H]])

    X_4 = np.full_like(X_1, L_x + W_H)
    Y_4 = np.array([[(L_y-W_H)/2, (L_y-W_H)/2], [(L_y+W_H)/2, (L_y+W_H)/2]])
    Z_4 = np.array([[0, H_H], [0, H_H]])

    X_5 = np.array([[L_x, L_x + W_H], [L_x, L_x + W_H]])
    Y_5 = np.full_like(X_1, (L_y - W_H)/2)
    Z_5 = np.array([[0, 0], [H_H, H_H]])

    X_6 = np.array([[L_x, L_x + W_H], [L_x, L_x + W_H]])
    Y_6 = np.full_like(X_1, (L_y + W_H)/2)
    Z_6 = np.array([[0, 0], [H_H, H_H]])

    ax3d.plot(
        [L_x, L_x], [(L_y-W_H)/2, (L_y-W_H)/2], [0, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x], [(L_y+W_H)/2, (L_y+W_H)/2], [0, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x + W_H, L_x + W_H], [(L_y-W_H)/2, (L_y-W_H)/2], [0, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x + W_H, L_x + W_H], [(L_y+W_H)/2, (L_y+W_H)/2], [0, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x + W_H], [(L_y-W_H)/2, (L_y-W_H)/2], [H_H, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x + W_H], [(L_y+W_H)/2, (L_y+W_H)/2], [H_H, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x + W_H], [(L_y-W_H)/2, (L_y-W_H)/2], [0, 0],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x + W_H], [(L_y+W_H)/2, (L_y+W_H)/2], [0, 0],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x], [(L_y-W_H)/2, (L_y+W_H)/2], [0, 0],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x + W_H, L_x + W_H], [(L_y-W_H)/2, (L_y+W_H)/2], [0, 0],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x, L_x], [(L_y-W_H)/2, (L_y+W_H)/2], [H_H, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot(
        [L_x + W_H, L_x + W_H], [(L_y-W_H)/2, (L_y+W_H)/2], [H_H, H_H],
        color='black',
        alpha=1,
    )

    ax3d.plot_surface(
        X_1, Y_1, Z_1,
        color='blue',
        alpha=0.2,
        shade=True
    )

    ax3d.plot_surface(
        X_2, Y_2, Z_2,
        color='blue',
        alpha=0.2,
        shade=True
    )

    ax3d.plot_surface(
        X_3, Y_3, Z_3,
        color='blue',
        alpha=0.2,
        shade=True
    )

    ax3d.plot_surface(
        X_4, Y_4, Z_4,
        color='blue',
        alpha=0.2,
        shade=True
    )

    ax3d.plot_surface(
        X_5, Y_5, Z_5,
        color='blue',
        alpha=0.2,
        shade=True
    )

    ax3d.plot_surface(
        X_6, Y_6, Z_6,
        color='blue',
        alpha=0.2,
        shade=True
    )

def draw_hub():
    # Center
    x_c = L_x + W_H / 2
    y_c = L_y / 2

    # Parameters
    theta = np.linspace(0, 2*np.pi, 200)
    z = np.linspace(H_H, H_H + H_H2, 100)

    Theta, Z = np.meshgrid(theta, z)

    # Radius as function of z
    R = (R_2 - R_1) / H_H2 * (Z - H_H) + R_1

    # Parametric surface
    X = x_c + R * np.cos(Theta)
    Y = y_c + R * np.sin(Theta)

    ax3d.plot_surface(
        X, Y, Z,
        color='purple',
        alpha=0.3,
        linewidth=0,
        shade=True
    )

    x_1 = x_c + R_1 * np.cos(theta)
    y_1 = y_c + R_1 * np.sin(theta)
    z_1 = np.full_like(x_1, H_H)
    ax3d.plot(x_1, y_1, z_1, color='black', alpha=0.6)

    x_2 = x_c + R_2 * np.cos(theta)
    y_2 = y_c + R_2 * np.sin(theta)
    z_2 = np.full_like(x_2, H_H + H_H2)
    ax3d.plot(x_2, y_2, z_2, color='black', alpha=0.6)

def draw_field():
    draw_ground()
    draw_field_walls()
    draw_hub_post()
    draw_hub()

# Main Stuff

def plot2D():

    ax2d.cla()  # clear axes
    ax2d.set_xlim(0, L_x*1.5)
    ax2d.set_ylim(0, H_H + H_H2 + 1)
    ax2d.set_aspect('equal', adjustable='box') 

    ax2d.xaxis.set_major_locator(MultipleLocator(0.5))
    ax2d.yaxis.set_major_locator(MultipleLocator(0.5))
    ax2d.grid(True)

    idealPoints, ideal_point = calculate_iPoints()
    frPoints, frIdealPoint = calculate_frPoints()

    plotIdeal(idealPoints, ideal_point)
    plotWithFriction(frPoints, frIdealPoint)

def plot3D():

    ax3d.cla()

    ax3d.set_xlim(0, L_x * 1.5)
    ax3d.set_ylim(0, L_y)
    ax3d.set_zlim(0, H_H + H_H2 + 1)

    ax3d.set_xlabel('X (Forward)')
    ax3d.set_ylabel('Y (Side)')
    ax3d.set_zlabel('Z (Height)')

    draw_field()

    idealPoints, ideal_point = calculate_iPoints()
    frPoints, frIdealPoint = calculate_frPoints()

    idealPoints_xyz = calculate_XYZ_fromPoints(idealPoints)
    frPoints_xyz = calculate_XYZ_fromPoints(frPoints)

    plotIdeal3d(idealPoints_xyz)
    plotWithFriction3d(frPoints_xyz)

def update_plot():
    global V_s, V_b0, w_b, V_u0, V_n0

    # --- recompute dependent values ---
    V_s = (w / 60) * 2 * m.pi * R_w
    V_b0 = V_s * C_roll
    w_b = V_s / R_b * (1 - C_roll)
    V_u0 = V_b0 * m.cos(O_b)
    V_n0 = V_b0 * m.sin(O_b)    

    if enable2d:
        plot2D()
    plot3D()
    


ax3d.view_init(elev=40, azim=235)

if enable2d:
    _, ax2d = setup_2d_plot()
    ax2d.set_title('Projectile Motion of the Ball')
    ax2d.set_xlabel('Horizontal Distance (m)')
    ax2d.set_ylabel('Vertical Distance (m)')

update_plot()


plt.show()