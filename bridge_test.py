import numpy as np
import math
# This Python Code Does the Necessary Calculations Needed after obtaining the SFE and BME
# After the maximum force and moment and their respective locations have been determined.

file_path = "data.txt"  

# Open and read the file
try:
    with open(file_path, 'r') as file:
        # Read the file content
        content = file.read()
except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")

# Loop Through the TextFile to extract the values
text = content.splitlines()
# Store the shapes
shapes = {}
for x in range (len(text)):
    if (text[x][0] == "/"):
        pass
    elif (text[x][0] == "S"):
        shape = text[x]
        x += 1
        b = float(text[x][(text[x].index(" ") + 1):])
        x += 1
        h = float(text[x][(text[x].index(" ") + 1):])
        x += 1
        bot = float(text[x][(text[x].index(" ") + 1):])
        shapes.update({shape: [b, h, bot]})

# add the centroid location to the shapes
for x in shapes.values():
    x.append(x[2] + x[1]/2)

# calculate y bar
sumtop = 0
sumarea = 0
for x in shapes.values():
    sumtop += (x[0]*x[1]*x[3])
    sumarea += (x[0]*x[1])
y_bar = sumtop / sumarea

# calculate I
I = 0
for x in shapes.values():
    # check y bar of the centroid of each shape
    if (x[3] > y_bar):
        I += x[0]*x[1]*((x[3]-y_bar)**2) + (x[0]*x[1]**3) / 12
    else:
        I += x[0]*x[1]*((y_bar-x[3])**2) + (x[0]*x[1]**3) / 12

# read sfe & bme
file_path = "SFE_BME.txt"  

# Open and read the file
try:
    with open(file_path, 'r') as file:
        # Read the file content
        content = file.read()
except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")

text = content.splitlines()
load1_SFE = text[0]
load1_SFE = load1_SFE.split()
load1_BME = text[1]
load1_BME = load1_BME.split()
load2_SFE = text[2]
load2_SFE = load2_SFE.split()
load2_BME = text[3]
load2_BME = load2_BME.split()

# Q calculations
Q_centroid = 0
# check for all shapes above and below the y bar to do the respective calculations 
for x in shapes.values():
    # check bottom above
    if x[2] > y_bar:
        Q_centroid += x[0]*x[1]*(x[3]-y_bar)
    else:
        # check if top is above the y_bar
        if (x[2]+x[1] > y_bar):
            # create the cross section within this difference
            b = x[0]
            bot = y_bar
            h = x[2] + x[1] - y_bar
            centroid = bot + h / 2
            Q_centroid += b*h*(centroid-y_bar)

Q_glue = 0
glueLoc = 75
# check for all shapes above the glue location
for x in shapes.values():
    # check bottom
    if x[2] >= glueLoc:
        if (x[3] > y_bar):
            Q_glue += x[0]*x[1]*(x[3]-y_bar)
        else:
            Q_glue += x[0]*x[1]*(y_bar-x[3])

# total height adjust if needed
height = 76.27
# array of tensile stresses
load1_tensile_stress = []
for x in load1_BME:
    load1_tensile_stress.append((float(x)*(y_bar))/I)
load2_tensile_stress = []
for x in load2_BME:
    load2_tensile_stress.append((float(x)*(y_bar))/I)

# array of compressive stresses
load1_compressive_stress = []
for x in load1_BME:
    load1_compressive_stress.append((float(x)*(height-y_bar))/I)
load2_compressive_stress = []
for x in load2_BME:
    load2_compressive_stress.append((float(x)*(height-y_bar))/I)

# array of centroid shear stresses
# make sure to specify b
b = 2.57
load1_centroid_shear_stress = []
for x in load1_SFE:
    load1_centroid_shear_stress.append((float(x)*Q_centroid)/(I*b))
load2_centroid_shear_stress = []
for x in load2_SFE:
    load2_centroid_shear_stress.append((float(x)*Q_centroid)/(I*b))

# array of glue shear stresses
# make sure to specify b
b = 12.54
load1_glue_shear_stress = []
for x in load1_SFE:
    load1_glue_shear_stress.append((float(x)*Q_glue)/(I*b))
load2_glue_shear_stress = []
for x in load2_SFE:
    load2_glue_shear_stress.append((float(x)*Q_glue)/(I*b))

# Plate Buckling
# Case 1
# specify b depending on the scross section cut
b = 67.46
# specify t
t = 1.27
# specify k
k = 4
sigma_buck_1 = ((k * math.pi**2 * 4000) / (12 * (1 - 0.2**2))) * ((t / b) ** 2)

# Case 2
# specify b depending on the scross section cut
b = 10
# specify t
t = 1.27
# specify k
k = 0.425
sigma_buck_2 = ((k * math.pi**2 * 4000) / (12 * (1 - 0.2**2))) * ((t / b) ** 2)

# Case 3
# specify b depending on the scross section cut
b = 33.6
# specify t
t = 1.27
# specify k
k = 6
sigma_buck_3 = ((k * math.pi**2 * 4000) / (12 * (1 - 0.2**2))) * ((t / b) ** 2)

# find the smallest of the plate bucklings to calculate FOS
sigma_plate_buck = min(sigma_buck_1, sigma_buck_2, sigma_buck_3)

# Shear Buckling
# specify b depending on the scross section cut
a = 400
# specify t
t = 1.27
# specify h
h = 76.27
sigma_shear_buck = ((5*math.pi**2*4000)/(12*(1-0.2**2))) * ((t/h)**2 + (t/a)**2)
print (sigma_shear_buck)

# Loading Case 1 FOS
FOS1_tension = []
for x in load1_tensile_stress:
    if (x == 0) or (abs(30/x) > 100):
        FOS1_tension.append(100)
    else:
        FOS1_tension.append(abs(30/x))
FOS1_compression = []
for x in load1_compressive_stress:
    if (x == 0) or (abs(6/x) > 100):
        FOS1_compression.append(100)
    else:
        FOS1_compression.append(abs(6/x))
FOS1_buckling = []
for x in load1_compressive_stress:
    if (x == 0) or (abs(sigma_plate_buck/x) > 100):
        FOS1_buckling.append(100)
    else:
        FOS1_buckling.append(abs(sigma_plate_buck/x))
FOS1_centroid_shear = []
for x in load1_centroid_shear_stress:
    if (x == 0) or (abs(4/x) > 100):
        FOS1_centroid_shear.append(100)
    else:
        FOS1_centroid_shear.append(abs(4/x))
FOS1_glue_shear = []
for x in load1_glue_shear_stress:
    if ((x == 0) or abs(2/x) > 100):
        FOS1_glue_shear.append(100)
    else:
        FOS1_glue_shear.append(abs(2/x))
FOS1_shear_buck = []
for x in load1_centroid_shear_stress:
    if (x == 0) or (abs(sigma_shear_buck/x) > 100):
        FOS1_shear_buck.append(100)
    else:
        FOS1_shear_buck.append(abs(sigma_shear_buck/x))

# put into txt
with open("FOS_loadcase1.txt", "w") as file:
    for x in FOS1_tension:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS1_compression:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS1_buckling:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS1_centroid_shear:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS1_glue_shear:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS1_shear_buck:
        file.write(str(x) + " ")

# Loading Case 2 FOS
FOS2_tension = []
for x in load2_tensile_stress:
    if (x == 0) or (abs(30/x) > 100):
        FOS2_tension.append(100)
    else:
        FOS2_tension.append(abs(30/x))
FOS2_compression = []
for x in load2_compressive_stress:
    if (x == 0) or (abs(6/x) > 100):
        FOS2_compression.append(100)
    else:
        FOS2_compression.append(abs(6/x))
FOS2_buckling = []
for x in load2_compressive_stress:
    if (x == 0) or (abs(sigma_plate_buck > 100)):
        FOS2_buckling.append(100)
    else:
        FOS2_buckling.append(abs(sigma_plate_buck/x))
FOS2_centroid_shear = []
for x in load2_centroid_shear_stress:
    if (x == 0) or (abs(4/x) > 100):
        FOS2_centroid_shear.append(100)
    else:
        FOS2_centroid_shear.append(abs(4/x))
FOS2_glue_shear = []
for x in load2_glue_shear_stress:
    if (x == 0) or (abs(2/x) > 100):
        FOS2_glue_shear.append(100)
    else:
        FOS2_glue_shear.append(abs(2/x))
FOS2_shear_buck = []
for x in load2_centroid_shear_stress:
    if (x == 0) or (abs(sigma_shear_buck/x) > 100):
        FOS2_shear_buck.append(100)
    else:
        FOS2_shear_buck.append(abs(sigma_shear_buck/x))

# put into txt
with open("FOS_loadcase2.txt", "w") as file:
    for x in FOS2_tension:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS2_compression:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS2_buckling:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS2_centroid_shear:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS2_glue_shear:
        file.write(str(x) + " ")
    file.write("\n")
    for x in FOS2_shear_buck:
        file.write(str(x) + " ")