def Timer(Outfile, Time, Step, i):
    
    if (100*i / (Time/Step)) % 5 == 0:
        Per = int(100*i / (Time/Step))
        with open(Outfile, "a") as f:
            f.write(
                "Currently performed %.3f%% of the calculation.\n" % (Per))
            f.write('['+int(Per/5)*'##'+(20-int(Per/5))*'  '+']\n')