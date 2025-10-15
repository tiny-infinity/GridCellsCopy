import subprocess

sim_names = ["THETA_CONST","THETA_OSC","THETA_MOD"]

for name in sim_names:
    print(f"Running {name}")
    result = subprocess.run(["python","m_sim_setup.py",f"specs/theta_mod/{name}.py", "-o"],capture_output=True,text=True,check=True)
    print(f"Success",result.stdout)
