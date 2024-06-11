import sys, numpy as np
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from particles import *
from kepler import *
from load import *
from visualization import *

output_path = './results/'

nbody = Particles()

window_width  = 640
window_height = 480
w_margin = 15
h_margin = 15
column_gap = 20

# Create application and window
app = QApplication(sys.argv)
w = QWidget()
w.resize(window_width, window_height)
w.setWindowTitle('Settings')

#------------
# New object
#------------

vec_box_size = 30
vec_box_offset = 6
vec_box_w_gap = 5
vec_box_h_gap = 10

params_title = QLabel(w); params_title.setText('NEW OBJECT'); params_title.setFont(QFont("Arial", 13, QFont.Bold))
params_title.move(w_margin, h_margin)

name_label = QLabel(w); name_label.setText('Name')
name_label.move(params_title.pos().x(), params_title.pos().y() + params_title.height())

name_box = QLineEdit(w); name_box.setText('Particle 1')
name_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
name_box.move(name_label.pos().x() + name_label.width(), name_label.pos().y() - vec_box_offset)

mass_label = QLabel(w); mass_label.setText('Mass [kg]')
mass_label.move(name_label.pos().x(), name_box.pos().y() + name_box.height() + vec_box_h_gap)

mass_box = QLineEdit(w); mass_box.setText('1e23')
mass_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
mass_box.move(mass_label.pos().x() + mass_label.width(), mass_label.pos().y() - vec_box_offset)

divider1 = QLabel(w); divider1.setText('----------------------------------')
divider1.move(mass_label.pos().x(), mass_box.pos().y() + mass_box.height())

# Orbital elements

ecc_label = QLabel(w); ecc_label.setText('Eccentricity')
ecc_label.move(divider1.pos().x(), divider1.pos().y() + vec_box_h_gap * 2)

ecc_box = QLineEdit(w); ecc_box.setText('0.2')
ecc_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
ecc_box.move(ecc_label.pos().x() + ecc_label.width(), ecc_label.pos().y() - vec_box_offset)

a_label = QLabel(w); a_label.setText('a [m]') # semi-major axis
a_label.move(ecc_label.pos().x(), ecc_box.pos().y() + ecc_box.height() + vec_box_h_gap)

a_box = QLineEdit(w); a_box.setText('1e10')
a_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
a_box.move(a_label.pos().x() + a_label.width(), a_label.pos().y() - vec_box_offset)

inc_label = QLabel(w); inc_label.setText('Inclination [deg]')
inc_label.move(a_label.pos().x(), a_box.pos().y() + a_box.height() + vec_box_h_gap)

inc_box = QLineEdit(w); inc_box.setText('10')
inc_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
inc_box.move(inc_label.pos().x() + inc_label.width(), inc_label.pos().y() - vec_box_offset)

lan_label = QLabel(w); lan_label.setText('LAN [deg]') # longitude of the ascending node
lan_label.move(inc_label.pos().x(), inc_box.pos().y() + inc_box.height() + vec_box_h_gap)

lan_box = QLineEdit(w); lan_box.setText('48')
lan_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
lan_box.move(lan_label.pos().x() + lan_label.width(), lan_label.pos().y() - vec_box_offset)

ap_label = QLabel(w); ap_label.setText('AP [deg]') # argument of periapsis
ap_label.move(lan_label.pos().x(), lan_box.pos().y() + lan_box.height() + vec_box_h_gap)

ap_box = QLineEdit(w); ap_box.setText('29')
ap_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
ap_box.move(ap_label.pos().x() + ap_label.width(), ap_label.pos().y() - vec_box_offset)

ma_label = QLabel(w); ma_label.setText('MA [deg]') # mean anomaly
ma_label.move(ap_label.pos().x(), ap_box.pos().y() + ap_box.height() + vec_box_h_gap)

ma_box = QLineEdit(w); ma_box.setText('174')
ma_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
ma_box.move(ma_label.pos().x() + ma_label.width(), ma_label.pos().y() - vec_box_offset)

larger_mass_label = QLabel(w); larger_mass_label.setText('Larger mass [kg]') # mass of the central body
larger_mass_label.move(ma_label.pos().x(), ma_box.pos().y() + ma_box.height() + vec_box_h_gap)

larger_mass_box = QLineEdit(w); larger_mass_box.setText('1.9884e30')
larger_mass_box.resize(vec_box_size * 3 + vec_box_w_gap * 2, vec_box_size)
larger_mass_box.move(larger_mass_label.pos().x() + larger_mass_label.width(), larger_mass_label.pos().y() - vec_box_offset)

divider2 = QLabel(w); divider2.setText('----------------------------------')
divider2.move(larger_mass_label.pos().x(), larger_mass_box.pos().y() + larger_mass_box.height())

r_vec_label = QLabel(w); r_vec_label.setText('Position [m]')
r_vec_label.move(divider2.pos().x(), divider2.pos().y() + vec_box_h_gap * 2)

r_vec_x_box = QLineEdit(w); r_vec_x_box.setText('-1e10')
r_vec_x_box.resize(vec_box_size, vec_box_size)
r_vec_x_box.move(r_vec_label.pos().x() + r_vec_label.width(), r_vec_label.pos().y() - vec_box_offset)

r_vec_y_box = QLineEdit(w); r_vec_y_box.setText('-1e9')
r_vec_y_box.resize(vec_box_size, vec_box_size)
r_vec_y_box.move(r_vec_x_box.pos().x() + r_vec_x_box.width() + vec_box_w_gap, r_vec_x_box.pos().y())

r_vec_z_box = QLineEdit(w); r_vec_z_box.setText('1e9')
r_vec_z_box.resize(vec_box_size, vec_box_size)
r_vec_z_box.move(r_vec_y_box.pos().x() + r_vec_y_box.width() + vec_box_w_gap, r_vec_y_box.pos().y())

v_vec_label = QLabel(w); v_vec_label.setText('Velocity [m/s]')
v_vec_label.move(r_vec_label.pos().x(), r_vec_z_box.pos().y() + r_vec_z_box.height() + vec_box_h_gap)

v_vec_x_box = QLineEdit(w); v_vec_x_box.setText('-1e3')
v_vec_x_box.resize(vec_box_size, vec_box_size)
v_vec_x_box.move(r_vec_x_box.pos().x(), v_vec_label.pos().y() - vec_box_offset)

v_vec_y_box = QLineEdit(w); v_vec_y_box.setText('1e5')
v_vec_y_box.resize(vec_box_size, vec_box_size)
v_vec_y_box.move(r_vec_y_box.pos().x(), v_vec_x_box.pos().y())

v_vec_z_box = QLineEdit(w); v_vec_z_box.setText('-1e4')
v_vec_z_box.resize(vec_box_size, vec_box_size)
v_vec_z_box.move(r_vec_z_box.pos().x(), v_vec_y_box.pos().y())

convert_button = QPushButton('Convert', w); convert_button.setToolTip('Convert orbital elements to state vectors')
convert_button.resize(convert_button.sizeHint())
convert_button.move(v_vec_x_box.pos().x() - convert_button.width(), v_vec_x_box.pos().y() + v_vec_x_box.height())

def on_click_convert():
    try:
        e   = float(ecc_box.text())
        i   = float(inc_box.text())
        a   = float(a_box.text())
        LAN = float(lan_box.text())
        AP  = float(ap_box.text())
        MA  = float(ma_box.text())
        ms  = float(mass_box.text())
        ml  = float(larger_mass_box.text())
        
        X, V = kep2cart(e, a, i, LAN, ms + ml, ω=AP, M=MA)
        r_vec_x_box.setText(str(X[0]))
        r_vec_y_box.setText(str(X[1]))
        r_vec_z_box.setText(str(X[2]))
        v_vec_x_box.setText(str(V[0]))
        v_vec_y_box.setText(str(V[1]))
        v_vec_z_box.setText(str(V[2]))

        return 0

    except ValueError:
        print("Error: invalid parameters.")
        return 1

convert_button.clicked.connect(on_click_convert) # Connect the signals to the slots

add_button = QPushButton('Add', w); add_button.setToolTip('Add the object to list')
add_button.resize(add_button.sizeHint())
add_button.move(v_vec_x_box.pos().x() + vec_box_w_gap, v_vec_x_box.pos().y() + v_vec_x_box.height())

def on_click_add():
    try:
        name = str(name_box.text())
        mass = float(mass_box.text())
        x  = float(r_vec_x_box.text())
        y  = float(r_vec_y_box.text())
        z  = float(r_vec_z_box.text())
        vx = float(v_vec_x_box.text())
        vy = float(v_vec_y_box.text())
        vz = float(v_vec_z_box.text())
        
    except ValueError:
        print("Error: invalid parameters.")
        return 1

    nbody.add(mass, x, y, z, vx, vy, vz, name)
    update_table(nbody.obj) # Updates table

    return 0

add_button.clicked.connect(on_click_add)

#--------------
# Load presets
#--------------

column2_x = name_box.pos().x() + name_box.width() + column_gap

load_title = QLabel(w); load_title.setText('LOAD PRESETS'); load_title.setFont(QFont("Arial", 13, QFont.Bold))
load_title.move(column2_x, h_margin)

# Drop-down list for selecting preset planetary systems
presets_box1a = QComboBox(w); presets_box1a.move(column2_x, load_title.pos().y() + load_title.height() - 2 * vec_box_offset)
presets = preset_list(preset_path)
for x in presets:
    presets_box1a.addItem(x)

def on_activated_presets_box1a(text):
    presets_box1b.clear()
    global current_preset
    current_preset = load_preset(preset_path + str(text) + '.txt')
    for x in current_preset.keys():
        presets_box1b.addItem(x)

presets_box1a.activated[str].connect(on_activated_presets_box1a)

# Drop-down list for selecting individual objects
presets_box1b = QComboBox(w); presets_box1b.move(presets_box1a.pos().x(), presets_box1a.pos().y() + presets_box1a.height())
current_preset = load_preset(preset_path + presets[0] + '.txt') # Default preset
for x in current_preset.keys():
    presets_box1b.addItem(x)

def on_activated_presets_box1b(text):
    text_str = str(text)
    if text_str[0] == '*': # '*' indicates central body of system
        x = y = z = vx = vy = vz = 0
        m = current_preset[text_str][6]
        ecc_box.setText('N/A')
        a_box.setText('N/A')
        inc_box.setText('N/A')
        lan_box.setText('N/A')
        ap_box.setText('N/A')
        ma_box.setText('N/A')
        larger_mass_box.setText('N/A')
    else:
        e, a, i, Om, w, M, m = current_preset[text_str]
        m0 = current_preset[list(current_preset.keys())[0]][6]
        ecc_box.setText(str(e))
        a_box.setText(str(a))
        inc_box.setText(str(i))
        lan_box.setText(str(Om))
        ap_box.setText(str(w))
        ma_box.setText(str(M))
        larger_mass_box.setText(str(m0))
        (x, y, z), (vx, vy, vz) = kep2cart(e, a, i, Om, m + m0, ω=w, M=M)

    name_box.setText(text)
    mass_box.setText(str(m))

    r_vec_x_box.setText(str(x))
    r_vec_y_box.setText(str(y))
    r_vec_z_box.setText(str(z))
    v_vec_x_box.setText(str(vx))
    v_vec_y_box.setText(str(vy))
    v_vec_z_box.setText(str(vz))

presets_box1b.activated[str].connect(on_activated_presets_box1b)

# Drop-down list for selecting random or test particles
presets_box2 = QComboBox(w); presets_box2.move(presets_box1b.pos().x(), presets_box1b.pos().y() + presets_box1b.height())
test_cases = load_test_cases(test_cases)
for x in test_cases.keys():
    presets_box2.addItem(x)
presets_box2.addItem('Random object')

def on_activated_presets_box2(text):
    text_str = str(text)
    if text_str == 'Random object':
        name = 'Random object'
        m, x, y, z, vx, vy, vz = randomize(1e11, 1e27, 1e3)
    else:
        name = text_str
        m, x, y, z, vx, vy, vz = test_cases[text_str]
    
    name_box.setText(name)
    mass_box.setText(str(m))

    ecc_box.setText('N/A')
    a_box.setText('N/A')
    inc_box.setText('N/A')
    lan_box.setText('N/A')
    ap_box.setText('N/A')
    ma_box.setText('N/A')
    larger_mass_box.setText('N/A')

    r_vec_x_box.setText(str(x))
    r_vec_y_box.setText(str(y))
    r_vec_z_box.setText(str(z))
    v_vec_x_box.setText(str(vx))
    v_vec_y_box.setText(str(vy))
    v_vec_z_box.setText(str(vz))
    
presets_box2.activated[str].connect(on_activated_presets_box2)

quick_add_button = QPushButton('Quick Add', w); quick_add_button.setToolTip('Add all objects in the preset planetary system')
quick_add_button.resize(quick_add_button.sizeHint())
quick_add_button.move(presets_box2.pos().x(), presets_box2.pos().y() + presets_box2.height())

def on_click_quick_add():
    nbody.reset(N=len(current_preset))
    nbody.label = presets_box1a.currentText().replace(' ', '')
    for x in current_preset.keys():
        on_activated_presets_box1b(x)
        on_click_add()

quick_add_button.clicked.connect(on_click_quick_add)

divider3 = QLabel(w); divider3.setText('--------------------------------')
divider3.move(load_title.pos().x(), quick_add_button.pos().y() + quick_add_button.height() - vec_box_offset)

#-----------------
# Random ensemble
#-----------------

rand_title = QLabel(w); rand_title.setText('RANDOM ENSEMBLE'); rand_title.setFont(QFont("Arial", 13, QFont.Bold))
rand_title.move(divider3.pos().x(), divider3.pos().y() + vec_box_h_gap * 2)

nobj_label = QLabel(w); nobj_label.setText('Number')
nobj_label.move(rand_title.pos().x(), rand_title.pos().y() + rand_title.height())

nobj_box = QLineEdit(w); nobj_box.setText('12')
nobj_box.resize(vec_box_size * 3, vec_box_size)
nobj_box.move(nobj_label.pos().x() + nobj_label.width(), nobj_label.pos().y() - vec_box_offset)

radius_label = QLabel(w); radius_label.setText('Radius [m]')
radius_label.move(nobj_label.pos().x(), nobj_box.pos().y() + nobj_box.height() + vec_box_h_gap)

radius_box = QLineEdit(w); radius_box.setText('1e11')
radius_box.resize(vec_box_size * 3, vec_box_size)
radius_box.move(radius_label.pos().x() + radius_label.width(), radius_label.pos().y() - vec_box_offset)

vmax_label = QLabel(w); vmax_label.setText('Max speed [m/s]')
vmax_label.move(radius_label.pos().x(), radius_box.pos().y() + radius_box.height() + vec_box_h_gap)

vmax_box = QLineEdit(w); vmax_box.setText('1e3')
vmax_box.resize(vec_box_size * 3, vec_box_size)
vmax_box.move(vmax_label.pos().x() + vmax_label.width(), vmax_label.pos().y() - vec_box_offset)

mmax_label = QLabel(w); mmax_label.setText('Max mass [kg]')
mmax_label.move(vmax_label.pos().x(), vmax_box.pos().y() + vmax_box.height() + vec_box_h_gap)

mmax_box = QLineEdit(w); mmax_box.setText('1e27')
mmax_box.resize(vec_box_size * 3, vec_box_size)
mmax_box.move(mmax_label.pos().x() + mmax_label.width(), mmax_label.pos().y() - vec_box_offset)

add_random_button = QPushButton('Add Random', w); add_random_button.setToolTip('Add an ensemble of random objects')
add_random_button.resize(add_random_button.sizeHint())
add_random_button.move(mmax_box.pos().x() - int(add_random_button.width() / 2), mmax_box.pos().y() + mmax_box.height())
def on_click_add_random():
    try:
        N    = int(nobj_box.text())
        r0   = float(radius_box.text())
        vmax = float(vmax_box.text())
        mmax = float(mmax_box.text())
    except ValueError:
        N    = 12
        r0   = 1e11
        vmax = 1e3
        mmax = 1e27
        print("Warning: initial parameters not entered; set to default values.")

    nbody.gen_random(N, r0, mmax, vmax)
    update_table(nbody.obj)

add_random_button.clicked.connect(on_click_add_random)

divider4 = QLabel(w); divider4.setText('--------------------------------')
divider4.move(divider3.pos().x(), add_random_button.pos().y() + add_random_button.height() - vec_box_offset)

#------------------------
# Visualization settings
#------------------------

visual_title = QLabel(w); visual_title.setText('VISUALIZATION'); visual_title.setFont(QFont("Arial", 13, QFont.Bold))
visual_title.move(divider4.pos().x(), divider4.pos().y() + vec_box_h_gap * 2)

# Checkbox for setting range limits to solar system
range_check = QCheckBox(w); range_check.setText('Set radius at          [m]'); range_check.setChecked(True)
range_check.move(visual_title.pos().x(), visual_title.pos().y() + visual_title.height() - vec_box_offset)

range_box = QLineEdit(w); range_box.setText('6e12')
range_box.resize(vec_box_size, vec_box_size)
range_box.move(range_check.pos().x() + range_check.width(), range_check.pos().y() - vec_box_offset + 1)

# Checkbox for COM frame
com_check = QCheckBox(w); com_check.setText('Use COM frame'); com_check.setChecked(True)
com_check.move(range_check.pos().x(), range_check.pos().y() + range_check.height())

#-----------------
# List of objects
#-----------------

column3_x = nobj_box.pos().x() + nobj_box.width() + column_gap

table_width = window_width - w_margin - column3_x
table_height = 200

# Table of all particles created for the simulation
table_title = QLabel(w); table_title.setText('LIST OF OBJECTS'); table_title.setFont(QFont("Arial", 13, QFont.Bold))
table_title.move(column3_x, h_margin)

table = QTableWidget(w)
table.resize(table_width, table_height)
table.move(table_title.pos().x(), table_title.pos().y() + table_title.height() - vec_box_offset)

def update_table(objects): # update table after adding object
    num_items = objects['xi'].shape[0]
    table.setColumnCount(1)
    table.setRowCount(num_items)
    count = 0
    for i in range(num_items):
        Name = objects['id'][i]
        table.setItem(count, 0, QTableWidgetItem(Name))
        count += 1

divider5 = QLabel(w); divider5.setText('-----------------------------')
divider5.move(table.pos().x(), table.pos().y() + table.height())

#------------
# Simulation
#------------

sim_title = QLabel(w); sim_title.setText('SIMULATION'); sim_title.setFont(QFont("Arial", 13, QFont.Bold))
sim_title.move(divider5.pos().x(), divider5.pos().y() + vec_box_h_gap * 2)

# Entering epoch at start in Julian date; default: J2000
epoch_label = QLabel(w); epoch_label.setText('Start epoch (JD)')
epoch_label.move(sim_title.pos().x(), sim_title.pos().y() + sim_title.height())

epoch_box = QLineEdit(w); epoch_box.setText('2451545.0')
epoch_box.resize(vec_box_size * 3, vec_box_size)
epoch_box.move(epoch_label.pos().x() + epoch_label.width(), epoch_label.pos().y() - vec_box_offset)

# Entering number of steps per day
res_label = QLabel(w); res_label.setText('Resolution [1/d]')
res_label.move(epoch_label.pos().x(), epoch_box.pos().y() + epoch_box.height() + vec_box_h_gap)

res_box = QLineEdit(w); res_box.setText('10') # 10 steps per day gives reasonable accuracy
res_box.resize(vec_box_size * 3, vec_box_size)
res_box.move(res_label.pos().x() + res_label.width(), res_label.pos().y() - vec_box_offset)

# Entering number of years
duration_label = QLabel(w); duration_label.setText('Duration [d]')
duration_label.move(res_label.pos().x(), res_box.pos().y() + res_box.height() + vec_box_h_gap)

duration_box = QLineEdit(w); duration_box.setText('3652.5')
duration_box.resize(vec_box_size * 3, vec_box_size)
duration_box.move(duration_label.pos().x() + duration_label.width(), duration_label.pos().y() - vec_box_offset)

# Entering force law
force_label = QLabel(w); force_label.setText('Force index')
force_label.move(duration_label.pos().x(), duration_box.pos().y() + duration_box.height() + vec_box_h_gap)

force_box = QLineEdit(w); force_box.setText('2')
force_box.resize(vec_box_size * 3, vec_box_size)
force_box.move(force_label.pos().x() + force_label.width(), force_label.pos().y() - vec_box_offset)

start_button = QPushButton('Start', w); start_button.setToolTip('Start simulation')
start_button.resize(start_button.sizeHint())
start_button.move(force_box.pos().x() - start_button.width(), add_button.pos().y())
def on_click_start():
    try:
        epoch       = float(epoch_box.text())
        sample_rate = float(res_box.text())
        duration    = float(duration_box.text())
        force       = float(force_box.text())
    except ValueError: # If these not entered, set default values
        sample_rate = 1. # [1/d]
        duration    = 3652.5 # [d]
        epoch       = 2451545.0
        force       = 2.
        print("Warning: number of steps not entered; set to default values.")

    Nsteps = int(sample_rate * duration)

    if nbody.n_obj == 0:
        print("Error: no objects entered.")
        return 1
    
    nbody.epoch = epoch
    results = nbody.run(Nsteps, sample_rate)
    results['kepler'] = current_preset

    if '*Sun' in results['id']:
        cm = 'solar'
    elif '*Kerbol' in results['id']:
        cm = 'ksp'
    elif '*Jool' in results['id']:
        cm = 'jool'
    else:
        cm = 'default'

    # show preview
    animate(results, set_range=range_check.isChecked(), plot_range=float(range_box.text()), 
            COM=com_check.isChecked(), cm=cm, show=True)

    output_name = '%s_%dsteps_%.1fd_n_%.1f'%(nbody.label, sample_rate, duration, force)

    # plot final state
    fig1, fig2 = plot(results, set_range=range_check.isChecked(), plot_range=float(range_box.text()), 
                      COM=com_check.isChecked(), cm=cm)
    fig1.savefig(plot_dir+output_name+'_3d.pdf', bbox_inches='tight')
    fig2.savefig(plot_dir+output_name+'_2d.pdf', bbox_inches='tight')
    
    # save results
    np.save(output_path+output_name+'.npy', results)

start_button.clicked.connect(on_click_start)

# Button for resetting the simulation
reset_button = QPushButton('Reset', w); reset_button.setToolTip('Reset simulation')
reset_button.resize(reset_button.sizeHint())
reset_button.move(start_button.pos().x() + start_button.width(), start_button.pos().y())
def on_click_reset():
    nbody.reset()
    name_box.setText('Particle 1')
    mass_box.setText('1e23')
    ecc_box.setText('0.2')
    a_box.setText('1e10')
    inc_box.setText('10')
    lan_box.setText('48')
    ap_box.setText('29')
    ma_box.setText('174')
    larger_mass_box.setText('1.9884e30')
    r_vec_x_box.setText('-1e10')
    r_vec_y_box.setText('-1e9')
    r_vec_z_box.setText('1e9')
    v_vec_x_box.setText('-1e3')
    v_vec_y_box.setText('-1e5')
    v_vec_z_box.setText('-1e4')
    nobj_box.setText('12')
    radius_box.setText('1e11')
    vmax_box.setText('1e3')
    mmax_box.setText('1e27')
    epoch_box.setText('2451545.0')
    res_box.setText('10')
    duration_box.setText('3652.5')
    force_box.setText('2')
    table.setColumnCount(0)
    table.setRowCount(0)
    plt.close(); plt.close()

reset_button.clicked.connect(on_click_reset)

w.show()

sys.exit(app.exec_())
