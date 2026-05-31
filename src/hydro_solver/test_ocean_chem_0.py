from phreeqpy.iphreeqc.phreeqc_dll import IPhreeqc

iphreeqc = IPhreeqc()
print("phreeqpy import OK")

# Load database (important!)
iphreeqc.load_database("./external/phreeqc/database/phreeqc.dat")

# Define solution
input_str = """
SOLUTION 1
    temp 25
    pH 7
    Na 1
    Cl 1

SELECTED_OUTPUT
    -reset false
    -pH true
    -temperature true
    -ionic_strength true

USER_PUNCH
    -headings Na Cl
    PUNCH TOT("Na"), TOT("Cl")

END
"""

iphreeqc.run_string(input_str)

# Get results
output = iphreeqc.get_selected_output_array()
print(output)