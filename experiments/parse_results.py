import sys
import os
import re

def main():
    str_len = len(sys.argv[1])
    apps = sys.argv[1][1:str_len-1]
    apps = apps.split(',')

    f_options = ['full', 'layer', 'google', 'opt', 'uniform', 'naive']
    out_str = " , full, layer, google, opt, uniform, naive\n" 
    for o in apps:
        print("Collecting results for " + o)
        avg_success = {}
        avg_compile = {}
        avg_decohere = {}
        avg_depth = {}
        path = o+'_res/'
        #for q in [4,9,16,25]:
        for q in [4]:
            filepath = path + "q" + str(q) + ".out"
            if not os.path.isfile(filepath):
                print("File path {} does not exist. Exiting...".format(filepath))
                sys.exit()
            with open(filepath) as fp:
                out_str += o + "(" + str(q) + "), "
                success_rate, compile_time, decohere_err, depth = record_data(fp, f_options)
                avg_success = {k:mean(v) for (k,v) in success_rate.items()}
                avg_compile = {k:mean(v) for (k,v) in compile_time.items()}
                avg_decohere = {k:mean(v) for (k,v) in decohere_err.items()}
                avg_depth = {k:int(mean(v)) for (k,v) in depth.items()}
                #print(o + "(" + str(q) + ")")
                #print("Success_rate:")
                #print(avg_success)
                #print("Compile_time:")
                #print(avg_compile)
                #print("Decohere_error:")
                #print(avg_decohere)
                out_str += str(avg_success['full']) + ", "
                out_str += str(avg_success['layer']) + ", "
                out_str += str(avg_success['google']) + ", "
                out_str += str(avg_success['opt']) + ", "
                out_str += str(avg_success['uniform']) + ", "
                out_str += str(avg_success['naive']) + "\n"

                #out_str += str(avg_depth['full']) + ", "
                #out_str += str(avg_depth['layer']) + ", "
                #out_str += str(avg_depth['google']) + ", "
                #out_str += str(avg_depth['opt']) + ", "
                #out_str += str(avg_depth['uniform']) + ", "
                #out_str += str(avg_depth['naive']) + "\n"
    print(out_str)

def mean(s):
    return sum(s)/len(s)


def record_data(fp, f_options):
    lines = []
    for line in fp:
        lines.append(line)
    success_rate = {f:[] for f in f_options}
    header_success = "Final success rate:"
    compile_time = {f:[] for f in f_options}
    header_compile = "Compilation time:"
    decohere_err = {f:[] for f in f_options}
    header_decohere = "Decoherence error:"
    depth = {f:[] for f in f_options}
    header_depth = "Circuit depth:"
    for (i, line) in enumerate(lines):
        #print("line {} contents {}".format(cnt, line))
        line = lines[i].strip().split(' ')
        if "---" in line: # filename line
            f = line[2]
            assert(f in f_options)
            noisy = []
            if header_success in lines[i+1]:
                success_rate[f].append(float(lines[i+1][len(header_success)+1:]))
            if header_compile in lines[i+7]:
                compile_time[f].append(float(lines[i+7][len(header_compile)+1:-2]))
            if header_decohere in lines[i+6]:
                decohere_err[f].append(float(lines[i+6][len(header_decohere)+1:]))
            if header_depth in lines[i+4]:
                tmp = lines[i+4][len(header_depth)+1:].strip().split(',')
                depth[f].append(int(tmp[1]))
    #print(success_rate)
    #print(compile_time)
    return success_rate, compile_time, decohere_err, depth

if __name__ == '__main__':
    main()
