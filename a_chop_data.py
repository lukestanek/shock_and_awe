import numpy as np

startFrame = 900
endFrame = 1400

fp  = open("0_case_movie_3_pos_vel_KE_corrected.txt")
fpOut = open("0_case_movie_3_pos_vel_KE_chopped.txt", "w")


frame = 0

while frame <= endFrame+1:
    print "{0}/{1}".format(frame,5000)
    line = fp.readline()
    numParts = int(line)
    comm = fp.readline()

    if frame >= startFrame and frame <= endFrame:
        # save to output file
        fpOut.write(str(numParts)+"\n")
        fpOut.write(comm)
    
    for i in range(numParts):
        line = fp.readline()

        if frame >= startFrame and frame <= endFrame:
            # save to output file
            fpOut.write(line)

    frame += 1

print "Done!"
fp.close()
fpOut.close()



