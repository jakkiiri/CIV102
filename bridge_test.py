# This Python Code Does the Necessary Calculations Needed after obtaining the SFE and BME
# After the maximum force and moment and their respective locations have been determined.

file_path = "data.txt"  # Replace 'your_file.txt' with the path to your file

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

print (shapes)
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
print(I)