# Note that the code will pass through the first 1746 lines which include all header lines and the "INFO" column header line. Since we add the column header, we get a difference of 1745 lines.

with open('info.txt', 'r') as f1, open("IC_values.txt", 'w') as f2, open("QD.txt", 'w') as f3, open("FS.txt", 'w') as f4, open("MQ.txt", 'w') as f5:
    f2.write('InbreedingCoeff\n')
    f3.write('QualityByDepth\n')
    f4.write('FisherStrandBias\n')
    f5.write("RMSMappingQuality\n")
    for line in f1:
        if line.startswith('#') or line.startswith("INFO"):
            continue
        info = line.strip().split(';')
        for i in info:
            if i.startswith("InbreedingCoeff"):
                ic = i.split('=')
                f2.write(f"{ic[1]}\n")
            elif i.startswith("QD"):
                qd = i.split('=')
                f3.write(f"{qd[1]}\n")
            elif i.startswith("FS"):
                fs = i.split('=')
                f4.write(f"{fs[1]}\n")
            elif i.startswith("MQ="):
                mq= i.split('=')
                f5.write(f"{mq[1]}\n")
