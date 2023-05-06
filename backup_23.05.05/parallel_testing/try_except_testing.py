a = [2,4,6,'apple','pear',10,12,'orange',14]


for i in a:
    try:
        #Try the below code which may fail to run
        num = int(i)
        print(num)
    except:
        #If an error occurs, record that iteration into a file and continue with the rest of the script using the pass statement
        with open('fail_log.txt', 'a+') as fh:
            fh.write(f"{i} not an interger\n")
        pass
