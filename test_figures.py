import argparse
import subprocess
argparser = argparse.ArgumentParser()
argparser.add_argument('figures',help='List of figures to generate',nargs='+')
argparser.add_argument('-a',"--all",help='Generate all figures',action='store_true')

args = argparser.parse_args()
spec_files = {"1":[
                    ["fig1_1D.py","s"],
                    ["fig1_size_lg.py","s"],
                    ["fig1_size_sm.py","s"],
                    ["fig1_stdn_pir.py","m"],
                    ["fig1_stdn_res.py","s"],
                    ["fig1_stdn_sag.py","m"],
                    ["fig1_stdn_subthresh.py","s"]],
                "1_lg":[["fig1_2D.py","s"]],
                "2":[["fig2_c.py","m"]],
                "3":[["fig3_abc.py","s"],
                     ["fig3_de.py","s"]]}

for i,fig in enumerate(args.figures):
    print(f"Generating figure {fig}")
    for fig_num in spec_files[fig]:
        print(f"Running sim {fig_num[0]}")
        proc=subprocess.run(["python",f"{fig_num[1]}_sim_setup.py",f"specs/fig{fig}/{fig_num[0]}","-v","-o"],capture_output=True)
        if proc.returncode != 0:
            print(proc.stderr)
            raise Exception(f"Error in sim {fig_num[0]}")
    if args.all:
        for fig_num in spec_files[f"{fig}_lg"]:
            print(f"Running sim {fig_num[0]}")
            proc=subprocess.run(["python",f"{fig_num[1]}_sim_setup.py",f"specs/fig{fig}/{fig_num[0]}","-v","-o"],capture_output=True)
            if proc.returncode != 0:
                print(proc.stderr)
                raise Exception(f"Error in sim {fig_num[0]}")