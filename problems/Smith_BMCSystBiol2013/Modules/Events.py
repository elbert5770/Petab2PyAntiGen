# Events.py: none.
#
# The SBML has three events, an insulin pulse that ends (Ins = 0 at t_ins) and a second one (Ins = 499999 at 2880, 0 at 2895).
# Every condition's last measurement is at or before its t_ins (15, 30, 240 or 960 min, the pulse ends when the data do), and
# the second pulse is at 2880: no event can change a simulated value that is measured, so there is nothing to simulate.
# (tests/test_smith.py checks this against PEtab's simulation, which does include the events.)

def generate_no_events(replicate, df_dict):
    return ''
