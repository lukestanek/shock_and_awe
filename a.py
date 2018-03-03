import numpy as np

fp  = open("0_case_movie_3_pos_vel_KE.txt")
fpOut = open("0_case_movie_3_pos_vel_KE_corrected.txt", "w")

loop = True
j = 0

while loop:
    print "{0}/{1}".format(j,5000)
    line = fp.readline()
    numParts = int(line)
    comm = fp.readline()

    fpOut.write(str(numParts)+"\n")
    fpOut.write(comm[0:-1] + " |v|\n")
    
    for i in range(numParts):
        line = fp.readline()
        vals = line.split()
        mag = str(np.sqrt(float(vals[4])**2 + float(vals[5])**2 + float(vals[6])**2))

        vals.append(mag)
        fpOut.write(" ".join(vals)+"\n")
    j += 1

print "Done!"
fp.close()
fpOut.close()



