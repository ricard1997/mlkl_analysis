


file = open("rotated_ref_1.pdb")


lines = file.readlines()

new_file = []

for line in lines:
    items = line.split()
    if items[0] == "ATOM":
        if items[2] == "CA" and int(items[5]) > 180:
            new_file.append(line)
    else:
        new_file.append(line)

file.close()

new = open("new_rotated.pdb","w")
for line in new_file:
    new.write(line)
new.close()

        

